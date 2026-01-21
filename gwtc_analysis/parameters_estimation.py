"""GW event parameter-estimation (PE) plotting pipeline.

This module is a refactor of the original notebook into a Galaxy/MMODA-friendly
entrypoint:

    run_parameters_estimation(args) -> int

`args` must provide exactly these fields:
  - src_name
  - sample_method
  - strain_approximant
  - start
  - stop
  - fs_low
  - fs_high

The function produces plot files in the working directory and, when available,
wraps them as ODA PictureProduct objects.

Notes
-----
- This code intentionally mirrors the notebook control-flow (with `go_next_cell`)
  to preserve behaviour and logs.
- The module stores the last run outputs in the module global `LAST_RESULT` for
  callers that want access despite the `-> int` signature.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple, Union
from pathlib import Path


@dataclass
class PEResult:
    """Collected outputs from a PE run."""

    fig_distribution: List[Any] = field(default_factory=list)
    fig_strain: List[Any] = field(default_factory=list)
    fig_psd: List[Any] = field(default_factory=list)
    fig_skymap: List[Any] = field(default_factory=list)
    tool_log: List[str] = field(default_factory=list)

    # Convenience: raw filenames (always populated, even if ODA is unavailable)
    files_distribution: List[str] = field(default_factory=list)
    files_strain: List[str] = field(default_factory=list)
    files_psd: List[str] = field(default_factory=list)
    files_skymap: List[str] = field(default_factory=list)


LAST_RESULT: Optional[PEResult] = None

import html
import json
import re
import requests


ZENODO_PE_RECORD_IDS = [6513631, 8177023, 17014085]

def _extract_event_id_from_filename(name: str) -> str | None:
    # Strong match used across GWTC releases
    m = re.search(r"(GW\d{6}_\d{6})", name)
    return m.group(1) if m else None


def _zenodo_record_files(record_id: int) -> list[dict[str, Any]]:
    url = f"https://zenodo.org/api/records/{record_id}"
    r = requests.get(url, timeout=(10, 120))
    r.raise_for_status()
    j = r.json()
    return j.get("files", []) or []


def build_zenodo_pe_index(
    *,
    cache_dir: str | Path = ".cache_gwosc",
    record_ids: list[int] | None = None,
    force_refresh: bool = False,
) -> dict[str, list[dict[str, Any]]]:
    """
    Returns:
      index[event_id] = [{record_id, filename, url}, ...]
    Cached in: <cache_dir>/zenodo_pe_index.json
    """
    record_ids = record_ids or ZENODO_PE_RECORD_IDS
    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)
    cache_path = cache_dir / "zenodo_pe_index.json"

    if (not force_refresh) and cache_path.exists() and cache_path.stat().st_size > 0:
        return json.loads(cache_path.read_text(encoding="utf-8"))

    index: dict[str, list[dict[str, Any]]] = {}

    for rid in record_ids:
        for f in _zenodo_record_files(rid):
            # Zenodo API uses 'key' for filename
            fname = f.get("key") or f.get("filename") or ""
            ev = _extract_event_id_from_filename(fname)
            if not ev:
                continue

            links = f.get("links", {}) or {}
            # Usually links.self is an API URL like .../api/records/<rid>/files/<fname>
            url = links.get("self") or links.get("download") or ""
            if url and "/api/" in url and "/files/" in url and not url.endswith("/content"):
                url = url.rstrip("/") + "/content"

            index.setdefault(ev, []).append(
                {"record_id": rid, "filename": fname, "url": url}
            )

    cache_path.write_text(json.dumps(index, indent=2), encoding="utf-8")
    return index


def choose_best_pe_file(candidates: list[dict[str, Any]]) -> dict[str, Any] | None:
    """
    Rule:
      - prefer 'mixed_cosmo' when exists
      - otherwise 'recombined'
      - fallback: prefer .hdf5/.h5, else first
    """
    if not candidates:
        return None

    def n(s: str) -> str:
        return (s or "").lower()

    for c in candidates:
        if "mixed_cosmo" in n(c.get("filename", "")):
            return c

    for c in candidates:
        if "recombined" in n(c.get("filename", "")):
            return c

    for ext in (".hdf5", ".h5"):
        for c in candidates:
            if n(c.get("filename", "")).endswith(ext):
                return c

    return candidates[0]


def _download_url_to_path(url: str, out_path: Path, *, chunk_size: int = 1024 * 1024) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with requests.get(url, stream=True, timeout=(10, 300)) as r:
        r.raise_for_status()
        with out_path.open("wb") as f:
            for chunk in r.iter_content(chunk_size=chunk_size):
                if chunk:
                    f.write(chunk)


def download_zenodo_pe_file(
    chosen: dict[str, Any],
    *,
    outdir: Path,
    progress_cb=None,
    log_cb=None,
) -> Path:
    """Download a single PE file described by `chosen` into outdir (cached).

    `chosen` must have keys: filename, url.
    """
    fname = str(chosen.get("filename") or "")
    url = str(chosen.get("url") or "")
    if not fname or not url:
        raise ValueError("Invalid Zenodo PE descriptor: missing filename/url")

    out_path = outdir / Path(fname).name
    out_path.parent.mkdir(parents=True, exist_ok=True)

    if out_path.exists() and out_path.stat().st_size > 0:
        if log_cb:
            log_cb(f"ℹ️ [CACHE] Using existing Zenodo PE file: {out_path}")
        return out_path

    if log_cb:
        log_cb(f"ℹ️ [DOWNLOAD] Zenodo PE: {fname}")
    if progress_cb:
        try:
            progress_cb("Download data", 15, "Zenodo PE download")
        except Exception:
            pass

    try:
        _download_url_to_path(url, out_path)
    except Exception as e:
        try:
            out_path.unlink(missing_ok=True)
        except Exception:
            pass
        raise RuntimeError(f"Failed downloading Zenodo PE file {fname}: {e}") from e

    return out_path

def _is_skymap_plot(p: Path) -> bool:
    name = p.name.lower()
    return any(k in name for k in (
        "skymap", "localization", "allsky", "radec", "ra_dec", "sky"
    ))

def _pe_plot_sort_key(p: Path) -> tuple[int, int, str]:
    """
    PE report plot ordering:

      0) Posteriors (corner/posterior, etc.)
      1) PSD
      2) Strain + waveform + Q-transform (grouped by detector)
      3) Skymaps / localization / RA-Dec (last)
      9) Everything else
    """
    name = p.name.lower()

    # 0) Posteriors
    if any(k in name for k in (
        "corner", "posterior", "posteriors",
        "credible", "ci", "marginal"
    )):
        group = 0

    # 1) PSD
    elif any(k in name for k in ("psd", "noise", "spectrum", "asd")):
        group = 1

    # 2) Strain / waveform / time-frequency
    elif any(k in name for k in (
        "strain", "waveform",
        "qtransform", "q-transform", "q_transform",
        "spectrogram", "timefreq", "time-freq", "tf"
    )):
        group = 2

    # 3) Skymaps / localization (last)
    elif any(k in name for k in (
        "skymap", "localization", "allsky", "sky", "radec", "ra_dec"
    )):
        group = 3

    else:
        group = 9

    # Secondary ordering: keep per-detector plots grouped (H1, L1, V1, K1)
    det_order = 9
    for i, det in enumerate(("h1", "l1", "v1", "k1")):
        if det in name:
            det_order = i
            break

    return (group, det_order, name)


def _write_parameters_estimation_report(
    *,
    out_path: Path,
    title: str,
    plots: list[Path],
    files: list[Path],
    params: dict[str, Any],
) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)

    def rel(p: Path) -> str:
        try:
            return str(p.relative_to(out_path.parent))
        except Exception:
            return str(p)

    def _norm_key_for_pairing(stem: str) -> str:
        """
        Build a pairing key by removing common RA/DEC tokens.
        Example:
          "event_ra"  -> "event"
          "event_dec" -> "event"
          "ra_event"  -> "event"
        """
        s = stem.lower()
        for tok in ("_ra", "_dec", "ra_", "dec_", "-ra", "-dec", "ra-", "dec-"):
            s = s.replace(tok, "_")
        # remove standalone ra/dec tokens
        parts = [p for p in s.replace("-", "_").split("_") if p and p not in ("ra", "dec")]
        return "_".join(parts) if parts else s

    def _is_ra_plot(p: Path) -> bool:
        n = p.name.lower()
        return ("_ra" in n) or n.startswith("ra_") or n.endswith("_ra.png") or (n == "ra.png")

    def _is_dec_plot(p: Path) -> bool:
        n = p.name.lower()
        return ("_dec" in n) or n.startswith("dec_") or n.endswith("_dec.png") or (n == "dec.png")

    def _is_skymap_plot(p: Path) -> bool:
        name = p.name.lower()
        return any(k in name for k in ("skymap", "localization", "allsky", "sky"))

    lines: list[str] = []
    lines.append("<html><head><meta charset='utf-8'>")
    lines.append(f"<title>{html.escape(title)}</title>")
    lines.append("</head><body>")
    lines.append(f"<h1>{html.escape(title)}</h1>")

    # Parameters block
    lines.append("<h2>Run parameters</h2>")
    lines.append("<table border='1' cellspacing='0' cellpadding='6'>")
    lines.append("<tr><th>Parameter</th><th>Value</th></tr>")
    for k, v in params.items():
        lines.append(
            f"<tr><td><code>{html.escape(str(k))}</code></td><td>{html.escape(str(v))}</td></tr>"
        )
    lines.append("</table>")

    # Output summary
    lines.append("<h2>Outputs</h2>")
    lines.append("<ul>")
    lines.append(f"<li>Plots: {len(plots)}</li>")
    lines.append(f"<li>Other files: {len(files)}</li>")
    lines.append("</ul>")

    if files:
        lines.append("<h3>Files</h3>")
        lines.append("<ul>")
        for p in files:
            rp = html.escape(rel(p))
            lines.append(f"<li><a href='{rp}'><code>{rp}</code></a></li>")
        lines.append("</ul>")

    if not plots:
        lines.append("<p><b>No plots found</b> in the plots directory.</p>")
        lines.append("</body></html>\n")
        out_path.write_text("\n".join(lines), encoding="utf-8")
        return

    # -----------------------
    # Plots
    # -----------------------
    lines.append("<h2>Plots</h2>")

    # 1) Pair RA/DEC plots explicitly and show them on the same row
    ra_plots = [p for p in plots if _is_ra_plot(p)]
    dec_plots = [p for p in plots if _is_dec_plot(p)]

    ra_by_key: dict[str, Path] = {_norm_key_for_pairing(p.stem): p for p in ra_plots}
    dec_by_key: dict[str, Path] = {_norm_key_for_pairing(p.stem): p for p in dec_plots}

    paired_keys = sorted(set(ra_by_key) | set(dec_by_key))

    used: set[Path] = set()

    if paired_keys:
        lines.append("<h3>RA / DEC</h3>")
        lines.append("<table><tr><th>RA</th><th>DEC</th></tr>")

        for k in paired_keys:
            ra_p = ra_by_key.get(k)
            dec_p = dec_by_key.get(k)

            lines.append("<tr>")

            # RA cell
            if ra_p is not None:
                used.add(ra_p)
                rp = html.escape(rel(ra_p))
                lines.append(
                    "<td style='padding:10px; vertical-align:top;'>"
                    f"<a href='{rp}'><img src='{rp}' style='max-width:450px; height:auto;'/></a>"
                    f"<br/><code>{rp}</code>"
                    "</td>"
                )
            else:
                lines.append("<td style='padding:10px;'><i>missing RA plot</i></td>")

            # DEC cell
            if dec_p is not None:
                used.add(dec_p)
                dp = html.escape(rel(dec_p))
                lines.append(
                    "<td style='padding:10px; vertical-align:top;'>"
                    f"<a href='{dp}'><img src='{dp}' style='max-width:450px; height:auto;'/></a>"
                    f"<br/><code>{dp}</code>"
                    "</td>"
                )
            else:
                lines.append("<td style='padding:10px;'><i>missing DEC plot</i></td>")

            lines.append("</tr>")

        lines.append("</table>")

    # Remaining plots excluding the RA/DEC ones already displayed
    remaining = [p for p in plots if p not in used]

    # Skymaps (other than RA/DEC) — render full-width like other plots
    skymaps = [p for p in remaining if _is_skymap_plot(p)]
    others = [p for p in remaining if p not in skymaps]

    # Normal plots: one per row
    for p in others:
        rp = html.escape(rel(p))
        lines.append(f"<p><a href='{rp}'><code>{rp}</code></a></p>")
        lines.append(
            f"<p><a href='{rp}'><img src='{rp}' style='max-width:900px; height:auto;'/></a></p>"
        )

    # Skymaps: one per row, full width
    if skymaps:
        lines.append("<h3>Sky localization</h3>")
        for p in skymaps:
            rp = html.escape(rel(p))
            lines.append(f"<p><a href='{rp}'><code>{rp}</code></a></p>")
            lines.append(
                f"<p><a href='{rp}'><img src='{rp}' style='max-width:900px; height:auto;'/></a></p>"
            )

    lines.append("</body></html>\n")
    out_path.write_text("\n".join(lines), encoding="utf-8")

def _get_arg(args: Any, name: str, default: Any = None, required: bool = True) -> Any:
    """Read an argument from dict/Namespace/object."""
    if isinstance(args, dict):
        if name in args:
            return args[name]
        if required:
            raise KeyError(f"Missing required argument '{name}'")
        return default

    # argparse.Namespace or similar
    if hasattr(args, name):
        val = getattr(args, name)
        if val is None and required:
            raise ValueError(f"Argument '{name}' is required but None")
        return default if val is None else val

    if required:
        raise AttributeError(f"Missing required argument '{name}'")
    return default


def _set_output(args: Any, name: str, value: Any) -> None:
    """Best-effort: attach outputs back to args (Galaxy/MMODA style)."""
    try:
        if isinstance(args, dict):
            args[name] = value
        else:
            setattr(args, name, value)
    except Exception:
        # Non-fatal: outputs still available via LAST_RESULT
        pass
def run_parameters_estimation(
    *,
    src_name: str,
    plots_dir: str | None = None,
    start: float = 0.5,
    stop: float = 0.1,
    fs_low: float = 20.0,
    fs_high: float = 300.0,
    sample_method: str = "Mixed",
    strain_approximant: str = "IMRPhenomXPHM",
    out_report_html: str | None = None,
    data_repo: str = "local",
    catalog: str | None = None,
) -> int:
    """Entry point: run the PE plotting pipeline.

    Parameters
    ----------
    src_name:
        Event name (e.g. GW200220_061928).
    plots_dir:
        Directory where plots and downloaded inputs are written.
    out_report_html:
        Optional HTML report path.
    data_repo:
        local | s3 | zenodo (others can be added later).
    catalog:
        Optional catalog (not required for zenodo PE selection).

    Returns
    -------
    int
        0 on success, non-zero on failure.
    """
    global LAST_RESULT

    from pathlib import Path
    from typing import Any, Dict, List, Optional
    import os
    import json
    import warnings

    # ---------------------------------------------------------------------
    # Normalize output dir
    # ---------------------------------------------------------------------
    if plots_dir is None:
        if out_report_html:
            plots_dir = str(Path(out_report_html).parent / "pe_plots")
        else:
            plots_dir = "pe_plots"

    outdir = Path(plots_dir) if plots_dir else Path(".")
    outdir.mkdir(parents=True, exist_ok=True)

    # --- imports kept inside to make module import cheap and Galaxy-friendly ---
    warnings.filterwarnings("ignore", category=UserWarning, append=True)
    try:
        from astropy.wcs import FITSFixedWarning  # type: ignore

        warnings.simplefilter("ignore", category=FITSFixedWarning)
    except Exception:
        pass

    warnings.filterwarnings("ignore", "Wswiglal-redir-stdio")

    import lal  # type: ignore

    lal.swig_redirect_standard_output_error(False)

    import matplotlib  # type: ignore

    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt  # type: ignore
    from matplotlib.colorbar import Colorbar  # type: ignore

    from pesummary.io import read  # type: ignore

    from .gwpe_utils import (
        load_strain,
        generate_projected_waveform,
        plot_whitened_overlay,
        plot_time_frequency,
        plot_basic_posteriors,
        select_label,
        label_report,
    )

    # ODA/MMODA objects are optional
    try:
        from oda_api.data_products import PictureProduct  # type: ignore
        from oda_api.api import ProgressReporter  # type: ignore

        oda_available = True
    except Exception:
        PictureProduct = None  # type: ignore
        ProgressReporter = None  # type: ignore
        oda_available = False

    pr = ProgressReporter() if oda_available else None

    # --- outputs ---
    result = PEResult()
    event_logs: List[str] = [""]

    fig_distributionList: List[Any] = []
    fig_strainList: List[Any] = []
    fig_psdList: List[Any] = []
    fig_skymapList: List[Any] = []

    go_next_cell = True

    def _log(msg: str) -> None:
        event_logs[-1] += f"\n{msg}\n"

    def _progress(stage: str, progress: int, substage: str) -> None:
        if pr is not None:
            try:
                pr.report_progress(stage=stage, progress=progress, substage=substage)
            except Exception:
                pass

    # ---------------------------------------------------------------------
    # Resolve PE file (local / S3 / Zenodo)
    # ---------------------------------------------------------------------
    _progress("Download data", 10, "step 1")

    data = None
    label_waveform: Optional[str] = None
    local_pe_path: Optional[str] = None

    try:
        if data_repo == "zenodo":
            from difflib import get_close_matches

            # 1) load cached index
            index = build_zenodo_pe_index(cache_dir=".cache_gwosc", force_refresh=False)
            index = build_zenodo_pe_index(cache_dir=".cache_gwosc", force_refresh=False)

            # 2) if missing, refresh once (cache may be stale)
            if src_name not in index:
                _log(f"⚠️ [WARN] Event {src_name} not found in cached Zenodo index. Refreshing index…")
                index = build_zenodo_pe_index(cache_dir=".cache_gwosc", force_refresh=True)

            cands = index.get(src_name, [])
            chosen = choose_best_pe_file(cands)

            if chosen is None:
                # Suggestions to help user (and to diagnose naming mismatch)
                keys = list(index.keys())
                sugg = get_close_matches(src_name, keys, n=5, cutoff=0.6)
                msg = (
                    f"No Zenodo PE file found for event {src_name}. "
                    f"Closest matches: {', '.join(sugg) if sugg else '(none)'}"
                )
                _log(f"❌ [ERROR] {msg}")
                raise ValueError(msg)

            local_path = download_zenodo_pe_file(
                chosen,
                outdir=outdir,
                progress_cb=_progress if pr is not None else None,
                log_cb=_log,
            )
            local_pe_path = str(local_path)


        elif data_repo == "s3":
            from minio import Minio  # type: ignore

            credentials_env = os.environ.get("S3_CREDENTIALS")
            if credentials_env:
                try:
                    credentials = json.loads(credentials_env)
                except Exception as e:
                    credentials = {"endpoint": "minio-dev.odahub.fr", "secure": True}
                    _log(f"⚠️ [WARN] Could not parse S3_CREDENTIALS JSON: {e}")
            else:
                credentials = {"endpoint": "minio-dev.odahub.fr", "secure": True}

            client = Minio(
                endpoint=credentials["endpoint"],
                secure=credentials.get("secure", True),
                access_key=credentials.get("access_key"),
                secret_key=credentials.get("secret_key"),
            )

            bucket = "gwtc"
            pe_prefixes = ["GWTC-2.1/PE/", "GWTC-3/PE/", "GWTC-4/PE/"]
            candidates: List[str] = []

            for prefix in pe_prefixes:
                for obj in client.list_objects(bucket, prefix=prefix, recursive=True):
                    key = obj.object_name
                    low = key.lower()
                    if src_name not in key:
                        continue
                    if not (low.endswith(".h5") or low.endswith(".hdf5")):
                        continue
                    candidates.append(key)

            if not candidates:
                msg = (
                    f"❌ [ERROR] No S3 PE file found for event {src_name} "
                    f"under {', '.join(pe_prefixes)}"
                )
                _log(msg)
                go_next_cell = False
            else:
                def _rank(k: str) -> Tuple[int, str]:
                    kl = k.lower()
                    if "mixed_cosmo" in kl:
                        return (0, kl)
                    if "recombined" in kl:
                        return (1, kl)
                    return (9, kl)

                remote_file = sorted(candidates, key=_rank)[0]
                _log(f"ℹ️ [INFO] Selected remote PE file: {remote_file}")

                local_pe_path = str(outdir / Path(remote_file).name)

                if os.path.isfile(local_pe_path):
                    _log(f"ℹ️ [CACHE] Using existing local PE file: {local_pe_path}")
                else:
                    _log(f"ℹ️ [DOWNLOAD] Fetching from S3: {remote_file}")
                    client.fget_object(bucket, remote_file, local_pe_path)
        else:
            # local (or galaxy): expect user already has the PE file locally or
            # implement your local discovery logic here.
            # For now, keep existing behavior: look for a PE file in outdir.
            candidates = sorted(outdir.glob(f"*{src_name}*PEDataRelease*.h5")) + sorted(
                outdir.glob(f"*{src_name}*PEDataRelease*.hdf5")
            )
            if not candidates:
                _log(
                    "❌ [ERROR] local mode: no PE file found in output directory. "
                    "Provide a PEDataRelease .h5/.hdf5 file or use --data-repo s3/zenodo."
                )
                go_next_cell = False
            else:
                local_pe_path = str(candidates[0])
                _log(f"ℹ️ [INFO] local mode: using PE file {local_pe_path}")

        if go_next_cell and local_pe_path is not None:
            _log(f"ℹ️ [INFO] Reading PE data from: {local_pe_path}")
            _progress("Read data", 25, "step 2")
            try:
                data = read(local_pe_path)
                label_report(data, sample_method=sample_method, strain_approximant=strain_approximant)
            except Exception as e:
                _log(f"❌ [ERROR] Failed reading PE data: {e}")
                go_next_cell = False

    except Exception as e:
        _log(f"❌ [ERROR] Failed resolving PE input ({data_repo}): {e}")
        go_next_cell = False

    # ---------------------------------------------------------------------
    # Posterior distributions
    # ---------------------------------------------------------------------
    if go_next_cell and data is not None:
        _progress("Posterior distributions", 50, "step 3")
        try:
            samples_dict = data.samples_dict
            labels = list(data.labels)

            label = select_label(
                data,
                sample_method,
                event_logs=event_logs,
                require_psd=False,
                show_labels=True,
            )

            if label is None:
                _log(f"❌ [ERROR] PE samples: no suitable label found; labels present: {labels}")
                go_next_cell = False
            else:
                label_waveform = label
                _log(f"ℹ️ [INFO] PE samples: final label = {label}")

                posterior_samples = samples_dict[label]
                approximant_for_title = label.split(":", 1)[1] if ":" in label else label

                posterior_files = plot_basic_posteriors(
                    posterior_samples,
                    src_name,
                    outdir,
                    approximant=approximant_for_title,
                )

                for _par_name, fname in posterior_files.items():
                    result.files_distribution.append(fname)
                    if oda_available:
                        fig_distributionList.append(PictureProduct.from_file(fname))

        except Exception as e:
            _log(f"❌ [ERROR] Failed creating posterior samples: {e}")
            go_next_cell = False

    # ---------------------------------------------------------------------
    # Detector strains + waveform overlay
    # ---------------------------------------------------------------------
    strain_data: Dict[str, Any] = {}
    merger_times: Dict[str, float] = {}
    if go_next_cell and data is not None:
        _progress("Download detector strains", 55, "step 4")
        try:
            samples_dict = data.samples_dict

            label_waveform = select_label(
                data,
                strain_approximant,
                event_logs=event_logs,
                require_psd=True,
                show_labels=False,
            )

            if label_waveform is None:
                _log(
                    "❌ [ERROR] Waveform model: no suitable label found in the PE file; "
                    "unable to load strain / build waveform."
                )
                go_next_cell = False
            else:
                _log(f"ℹ️ [INFO] Waveform model: using label_waveform = {label_waveform}")

            if go_next_cell:
                posterior_samples_wave = samples_dict[label_waveform]

                try:
                    dets_available = list(data.psd[label_waveform].keys())
                except Exception:
                    dets_available = ["H1", "L1", "V1"]

                _log(f"ℹ️ [INFO] Detectors in PE file (label {label_waveform}): {dets_available}")

                for det in dets_available:
                    det_time_key = f"{det}_time"
                    t0_det: Optional[float] = None
                    try:
                        t0_det = float(posterior_samples_wave.maxL[det_time_key][0])
                        _log(f"ℹ️ [INFO] Using maxL {det_time_key} = {t0_det} for {det}")
                    except Exception:
                        try:
                            t0_det = float(posterior_samples_wave.maxL["geocent_time"][0])
                            _log(f"ℹ️ [INFO] Using geocent_time = {t0_det} for {det}")
                        except Exception:
                            _log(f"⚠️ [WARN] No time information for {det} in maxL table")
                            continue

                    try:
                        strain_data[det] = load_strain(src_name, t0_det, det)
                        merger_times[det] = t0_det
                    except Exception as e:
                        _log(f"⚠️ [WARN] Could not load strain for {det}: {e}")

                if not strain_data:
                    _log("⚠️ [WARN] No strain data could be loaded for any detector.")
                    go_next_cell = False

        except Exception as e:
            _log(f"❌ [ERROR] Failed loading strain data: {e}")
            go_next_cell = False

    # Generate projected waveforms / q-transforms
    if go_next_cell and data is not None and label_waveform is not None:
        _progress("Generate projected waveforms", 90, "step 5")
        try:
            start_before = float(start)
            stop_after = float(stop)

            for det, strain in strain_data.items():
                t0 = merger_times[det]

                bp_cropped, crop_temp, used_aprx = generate_projected_waveform(
                    strain=strain,
                    event=src_name,
                    det=det,
                    t0=t0,
                    pedata=data,
                    label=label_waveform,
                    freqrange=(fs_low, fs_high),
                    time_window=(start_before, stop_after),
                )

                if bp_cropped is not None and crop_temp is not None:
                    fname_overlay = plot_whitened_overlay(
                        bp_cropped,
                        crop_temp,
                        src_name,
                        det,
                        outdir,
                        approximant=used_aprx,
                        t0=t0,
                    )
                    result.files_strain.append(fname_overlay)
                    if oda_available:
                        fig_strainList.append(PictureProduct.from_file(fname_overlay))
                else:
                    _log(f"⚠️ [WARN] Skipping overlay for {det} (no waveform).")

                try:
                    fname_q = plot_time_frequency(
                        strain=strain,
                        t0=t0,
                        event=src_name,
                        det=det,
                        outdir=outdir,
                        outseg=(-start_before, stop_after),
                        frange=(fs_low, fs_high),
                        approximant=used_aprx,
                    )
                    result.files_strain.append(fname_q)
                    if oda_available:
                        fig_strainList.append(PictureProduct.from_file(fname_q))
                except Exception as e_q:
                    _log(f"⚠️ [WARN] Could not create q-transform for {det}: {e_q}")

        except Exception as e:
            _log(f"❌ [ERROR] Failed creating projected waveforms / q-transforms: {e}")
            go_next_cell = False

    # ---------------------------------------------------------------------
    # PSD plot
    # ---------------------------------------------------------------------
    if go_next_cell and data is not None and label_waveform is not None:
        try:
            psd = data.psd[label_waveform]
            fig = psd.plot(fmin=20)
            ax = fig.gca()
            ax.set_ylim(1e-48, 1e-40)
            ax.set_title(f"PSD model for waveform: {label_waveform}", fontsize=12)
            plt.tight_layout()

            waveform = label_waveform.split(":", 1)[1] if ":" in label_waveform else label_waveform
            fname_psd = os.path.join(outdir, f"{src_name}_{waveform}_psd.png")
            fig.savefig(fname_psd, dpi=150)
            plt.close(fig)

            result.files_psd.append(fname_psd)
            if oda_available:
                fig_psdList.append(PictureProduct.from_file(fname_psd))

        except Exception as e:
            _log(f"❌ [ERROR] Failed creating PSD plot: {e}")
            go_next_cell = False

    # ---------------------------------------------------------------------
    # Skymap plot
    # ---------------------------------------------------------------------
    if go_next_cell and data is not None and label_waveform is not None:
        try:
            fig, ax = data.skymap[label_waveform].plot(contour=[50, 90])
            ax.set_title(f"Skymap for waveform: {label_waveform}", fontsize=9)
            ax.tick_params(axis="both", labelsize=7)
            ax.set_xlabel(ax.get_xlabel(), fontsize=7)
            ax.set_ylabel(ax.get_ylabel(), fontsize=7)

            try:
                for txt in ax.texts:
                    txt.set_fontsize(7)
            except Exception:
                pass

            for obj in fig.get_children():
                if isinstance(obj, Colorbar):
                    obj.ax.tick_params(labelsize=6)
                    if obj.ax.yaxis.label is not None:
                        obj.ax.yaxis.label.set_fontsize(7)

            waveform = label_waveform.split(":", 1)[1] if ":" in label_waveform else label_waveform
            fname_sky = os.path.join(outdir, f"{src_name}_{waveform}_skymap.png")
            fig.savefig(fname_sky, dpi=150)
            plt.close(fig)

            result.files_skymap.append(fname_sky)
            if oda_available:
                fig_skymapList.append(PictureProduct.from_file(fname_sky))

        except Exception as e:
            _log(f"❌ [ERROR] Failed creating skymap plot: {e}")
            go_next_cell = False

    _progress("Finish", 100, "step 6")


    # ---------------------------------------------------------------------
    # HTML report
    # ---------------------------------------------------------------------
    if out_report_html:
    
        out_path = Path(out_report_html)
        out_path.parent.mkdir(parents=True, exist_ok=True)

        plot_files = sorted(outdir.glob("*.png"), key=_pe_plot_sort_key)
        other_files = sorted(
            [p for p in outdir.iterdir() if p.is_file() and p.suffix.lower() not in {".png"}]
        )

        _write_parameters_estimation_report(
            out_path=out_path,
            title=f"Parameter estimation: {src_name}",
            plots=plot_files,
            files=other_files,
            params=dict(
                src_name=src_name,
                start=start,
                stop=stop,
                fs_low=fs_low,
                fs_high=fs_high,
                sample_method=sample_method,
                strain_approximant=strain_approximant,
                data_repo=data_repo,
            ),
        )

    # ---------------------------------------------------------------------
    # Finalize outputs
    # ---------------------------------------------------------------------
    result.tool_log = event_logs
    result.fig_distribution = fig_distributionList
    result.fig_strain = fig_strainList
    result.fig_psd = fig_psdList
    result.fig_skymap = fig_skymapList

    LAST_RESULT = result
 
    if not go_next_cell or data is None:
        # Stop pretending success: show the user the real reason
        last_msg = event_logs[-1].strip()
        raise ValueError(last_msg or "Parameter estimation failed before producing outputs.")
    # Return a plain dict so callers can see exactly what was produced
    return {
        "fig_distribution": result.files_distribution,
        "fig_strain": result.files_strain,
        "fig_psd": result.files_psd,
        "fig_skymap": result.files_skymap,
        "tool_log": event_logs,
    }

   
