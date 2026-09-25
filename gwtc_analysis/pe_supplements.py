"""Fill products missing from official PE releases with public supplementary data.

A few catalog PE files ship without noise PSDs (empty ``psds`` groups), which
leaves the strain overlay without its whitening PSD, or without skymaps. For
those events the products are taken from another public release of the same
event, usually the data release of its discovery paper on the DCC or Zenodo,
and cached locally.

Policy: if a PE label has no PSD (or no skymap), look the event up in
``PSD_SUPPLEMENTS``; if registered, fetch the product there (for PSDs only the
``psds`` group is read, over HTTP range requests when the server allows it)
and attach it to every label that lacks one. The source and its caveat are
written to the log, since these products come from a different analysis than
the catalog run.
"""
from __future__ import annotations

import json
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Optional

from .gwpe_utils import label_has_psd

USER_AGENT = "gwtc_analysis (https://github.com/danielsentenac/gwtc_analysis)"


@dataclass(frozen=True)
class PSDSupplement:
    event: str                   # full event name, e.g. GW230529_181500
    url: str                     # public PESummary file holding the PSDs
    reference: str               # release the PSDs come from, for the log
    run: Optional[str] = None    # label to read in that file (None: first with PSDs)
    note: str = ""               # caveat logged when the PSDs are used
    range_read: bool = True      # False: download the whole file (Zenodo rejects range reads)
    skymap_url: Optional[str] = None  # FITS skymap for labels without one


PSD_SUPPLEMENTS: dict[str, PSDSupplement] = {
    s.event: s
    for s in (
        PSDSupplement(
            event="GW230529_181500",
            url="https://zenodo.org/records/10845779/files/posterior_samples.h5?download=1",
            reference="GW230529 discovery-paper data release (Zenodo 10845779, LIGO-P2300352)",
            note="The GWTC-4 PE release has empty PSD groups for this event; the L1 PSD is identical "
                 "in all 15 discovery runs.",
            range_read=False,
            skymap_url="https://zenodo.org/records/10845779/files/skymap_combined_PHM_high_spin.fits?download=1",
        ),
        PSDSupplement(
            event="GW190425_081805",
            url="https://dcc.ligo.org/public/0165/P2000026/002/posterior_samples.h5",
            reference="GW190425 discovery-paper data release (LIGO-P2000026)",
            run="PhenomPNRT-HS",
            note="These are the PSDs of the earlier LALInference discovery analysis, not the PE_v3 PSDs "
                 "of the GWTC-2.1 bilby runs (those are not public).",
        ),
        PSDSupplement(
            event="GW200105_162426",
            url="https://dcc.ligo.org/public/0175/P2100143/002/GW200105_162426_posterior_samples_v2.h5",
            reference="GW200105/GW200115 discovery-paper data release (LIGO-P2100143)",
            run="C01:PhenomXPHM_high_spin",
            note="Deglitched BayesWave PSDs of the discovery analysis; GWTC-3 PE v1/v2 lack PSDs, v3 "
                 "includes them (they agree with these to ~2% median, differing mostly on lines).",
        ),
    )
}


def _base_name(name: str) -> str:
    return re.sub(r"-v\d+$", "", str(name).strip())


def get_psd_supplement(src_name: str) -> Optional[PSDSupplement]:
    """Registry lookup by full name (GW230529_181500) or short name (GW230529)."""
    name = _base_name(src_name)
    if name in PSD_SUPPLEMENTS:
        return PSD_SUPPLEMENTS[name]
    matches = [s for key, s in PSD_SUPPLEMENTS.items() if key.split("_")[0] == name]
    return matches[0] if len(matches) == 1 else None


def _cache_dir() -> Path:
    return Path.home() / ".gwcache" / "psd_supplements"


def _read_psds(h5, run: Optional[str]) -> dict:
    import numpy as np

    labels = [run] if run else [k for k in h5 if hasattr(h5[k], "keys") and "psds" in h5[k]]
    for lab in labels:
        grp = h5[lab]["psds"]
        psds = {det: np.asarray(grp[det][()], dtype=float) for det in grp}
        if psds:
            return psds
    raise KeyError(f"no PSDs found in run(s) {labels}")


def load_supplementary_psds(
    spec: PSDSupplement,
    *,
    log_cb: Callable[[str], None] | None = None,
    cache_dir: str | Path | None = None,
) -> Optional[dict]:
    """PSDs {det: (N, 2) array} of a registered supplement, cached as .npz."""
    import numpy as np

    cache = Path(cache_dir or _cache_dir())
    npz = cache / f"{spec.event}_psds.npz"
    meta = cache / f"{spec.event}_psds.json"
    if npz.exists():
        try:
            cached_url = json.loads(meta.read_text(encoding="utf-8")).get("url")
        except (OSError, ValueError):
            cached_url = None
        if cached_url == spec.url:
            with np.load(npz) as z:
                return {det: z[det] for det in z.files}

    import h5py

    if log_cb:
        log_cb(f"ℹ️ [DOWNLOAD] PSDs for {spec.event} from {spec.reference}: {spec.url}")
    try:
        if spec.range_read:
            import fsspec

            fs = fsspec.filesystem("http", block_size=2**20, headers={"User-Agent": USER_AGENT})
            with fs.open(spec.url, "rb") as fo, h5py.File(fo, "r") as h5:
                psds = _read_psds(h5, spec.run)
        else:
            import requests

            cache.mkdir(parents=True, exist_ok=True)
            tmp = cache / f"{spec.event}_source.h5.part"
            try:
                with requests.get(spec.url, stream=True, timeout=(10, 600),
                                  headers={"User-Agent": USER_AGENT}) as r:
                    r.raise_for_status()
                    with tmp.open("wb") as out:
                        for chunk in r.iter_content(chunk_size=2**20):
                            out.write(chunk)
                with h5py.File(tmp, "r") as h5:
                    psds = _read_psds(h5, spec.run)
            finally:
                tmp.unlink(missing_ok=True)  # keep only the extracted PSDs
    except Exception as e:
        if log_cb:
            log_cb(f"⚠️ [WARN] Could not fetch supplementary PSDs for {spec.event}: {type(e).__name__}: {e}")
        return None

    cache.mkdir(parents=True, exist_ok=True)
    np.savez(npz, **psds)
    meta.write_text(json.dumps({"url": spec.url, "run": spec.run, "reference": spec.reference}), encoding="utf-8")
    return psds


def fill_missing_psds(
    pedata,
    src_name: str,
    *,
    log_cb: Callable[[str], None] | None = None,
    cache_dir: str | Path | None = None,
) -> list[str]:
    """Attach supplementary PSDs to the labels of `pedata` that have none.

    Returns the labels that were filled (empty if nothing was missing, the event
    has no registered supplement, or the supplement could not be fetched).
    """
    labels = list(getattr(pedata, "labels", None) or [])
    missing = [lab for lab in labels if not label_has_psd(pedata, lab)]
    if not missing:
        return []

    spec = get_psd_supplement(src_name)
    if spec is None:
        if log_cb:
            log_cb(
                f"⚠️ [WARN] {len(missing)} PE label(s) of {src_name} have no PSD and no supplementary "
                "public release is registered for this event; the overlay will whiten with a PSD "
                "estimated from the strain."
            )
        return []

    psds = load_supplementary_psds(spec, log_cb=log_cb, cache_dir=cache_dir)
    if not psds:
        return []

    from pesummary.gw.file.psd import PSDDict

    if getattr(pedata, "psd", None) is None:
        pedata.psd = {}
    for lab in missing:
        pedata.psd[lab] = PSDDict({det: arr for det, arr in psds.items()})

    if log_cb:
        log_cb(
            f"ℹ️ [INFO] {src_name}: the PE file has no PSD for {len(missing)} label(s); using "
            f"{'/'.join(sorted(psds))} PSDs from the {spec.reference}. {spec.note}"
        )
    return missing


def fill_missing_skymaps(
    pedata,
    src_name: str,
    *,
    log_cb: Callable[[str], None] | None = None,
    cache_dir: str | Path | None = None,
) -> list[str]:
    """Attach the registered supplementary skymap to the labels without one."""
    labels = list(getattr(pedata, "labels", None) or [])
    container = getattr(pedata, "skymap", None)
    missing = [lab for lab in labels if container is None or lab not in container]
    spec = get_psd_supplement(src_name)
    if not missing or spec is None or not spec.skymap_url or container is None:
        return []

    cache = Path(cache_dir or _cache_dir())
    path = cache / f"{spec.event}_skymap.fits"
    try:
        if not (path.exists() and path.stat().st_size > 0):
            import requests

            if log_cb:
                log_cb(f"ℹ️ [DOWNLOAD] Skymap for {spec.event} from {spec.reference}: {spec.skymap_url}")
            r = requests.get(spec.skymap_url, timeout=(10, 300), headers={"User-Agent": USER_AGENT})
            r.raise_for_status()
            cache.mkdir(parents=True, exist_ok=True)
            path.write_bytes(r.content)

        from ligo.skymap.io import read_sky_map
        from pesummary.gw.file.skymap import SkyMap

        prob, meta = read_sky_map(str(path), nest=True)  # also flattens multi-order maps
        skymap = SkyMap(prob, meta_data={
            "nest": True,
            "objid": spec.event,
            "gps_time": float(meta.get("gps_time", 0.0)),
            "creator": meta.get("creator", ""),
            "origin": meta.get("origin", ""),
            "distmean": float(meta.get("distmean", 0.0)),
            "diststd": float(meta.get("diststd", 0.0)),
        })
    except Exception as e:
        if log_cb:
            log_cb(f"⚠️ [WARN] Could not load the supplementary skymap for {spec.event}: {type(e).__name__}: {e}")
        return []

    for lab in missing:
        container[lab] = skymap
    if log_cb:
        log_cb(
            f"ℹ️ [INFO] {src_name}: the PE file has no skymap for {len(missing)} label(s); using the "
            f"skymap of the {spec.reference}."
        )
    return missing
