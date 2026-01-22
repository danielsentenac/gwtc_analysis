#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
gw_stat_zenodo.py
=================

Extension of gw_stat to compute credible area for *all* confident-catalog events by using the
official Zenodo tarballs that contain skymap FITS files.

Why this exists
---------------
The GWOSC v2 parameters endpoint often provides only a scalar `sky_area` for a subset
of events and may not expose skymap links uniformly. However, the LVK "confident"
catalog PE releases on Zenodo include comprehensive skymap tarballs, e.g.:
- GWTC-4 archived skymaps tarball (Zenodo record 17014085) citeturn26search0
- GWTC-3.0 PE skylocalizations tarball (Zenodo record 8177023) citeturn26search1
- GWTC-2.1 PE skymaps tarball (Zenodo record 6513631) citeturn26search2

This module provides:
- download + cache of the tarball once per catalog
- per-event credible area extraction by reading the event's FITS from the tarball

Notes
-----
- Requires: `pip install ligo.skymap astropy healpy`
- The tarballs are large; caching avoids repeated downloads.
- credible area is computed by integrating the % credible region on the HEALPix map.
"""

from __future__ import annotations

import io
import os
import re
import time
import tarfile
from dataclasses import dataclass
from typing import Any, Dict, Iterable, Optional, Tuple, Union

import numpy as np
import pandas as pd
import requests

from .data_repo import zenodo_skymap_url, zenodo_cache_dir

FIGSIZE = (6, 4)

ALLOWED_CATALOGS = (
    "GWTC-1",
    "GWTC-2.1",
    "GWTC-3",
    "GWTC-4",
    "ALL",
)
def _download_with_byte_progress(
    url: str,
    dest: Path,
    *,
    progress: bool = True,
    verbose: bool = False,
    timeout: tuple[int, int] = (10, 600),
    chunk_size: int = 1024 * 1024,  # 1 MiB
    allow_resume: bool = True,
) -> Path:
    """
    Download `url` to `dest` with streaming + optional resume + tqdm byte progress.
    """
    dest = Path(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)

    try:
        from tqdm.auto import tqdm  # type: ignore
    except Exception:
        tqdm = None  # type: ignore

    # Resume if partial file exists
    existing = dest.stat().st_size if dest.exists() else 0
    headers: dict[str, str] = {}
    mode = "wb"

    if allow_resume and existing > 0:
        headers["Range"] = f"bytes={existing}-"
        mode = "ab"
        if verbose:
            print(f"[gw_stat] resuming download at byte {existing}: {dest}")

    with requests.get(url, stream=True, headers=headers, timeout=timeout) as r:
        # If server doesn't support Range, it may return 200; then restart.
        if existing > 0 and r.status_code == 200 and "Range" in headers:
            if verbose:
                print("[gw_stat] server ignored Range; restarting full download")
            existing = 0
            mode = "wb"

        r.raise_for_status()

        # Determine total size (best effort)
        total = None
        cl = r.headers.get("Content-Length")
        if cl is not None:
            try:
                total = int(cl)
            except ValueError:
                total = None

        # If resuming and server replied 206, total is remaining bytes; add existing for full bar
        if existing > 0 and r.status_code == 206 and total is not None:
            total = existing + total

        pbar = None
        if progress and tqdm is not None:
            pbar = tqdm(
                total=total,
                initial=existing if total is not None else 0,
                unit="B",
                unit_scale=True,
                unit_divisor=1024,
                desc=dest.name,
                leave=False,
            )

        tmp = dest.with_suffix(dest.suffix + ".part")

        # Write to .part then move into place at end (safer)
        # If resuming, we append to .part if it exists, else start from dest itself
        out_path = tmp
        if mode == "ab" and tmp.exists():
            # keep appending to existing .part
            pass
        elif mode == "ab" and dest.exists() and not tmp.exists():
            # move existing dest to .part so we append to the temp file
            dest.replace(tmp)

        with open(out_path, mode) as f:
            for chunk in r.iter_content(chunk_size=chunk_size):
                if not chunk:
                    continue
                f.write(chunk)
                if pbar is not None:
                    pbar.update(len(chunk))

        if pbar is not None:
            pbar.close()

        os.replace(out_path, dest)
        return dest
        
def _require_ligo_skymap():
    try:
        import ligo.skymap  # noqa: F401
    except ModuleNotFoundError as e:
        raise ModuleNotFoundError(
            "ligo.skymap is required for skymap/credible area functionality. "
            "Install the 'ligo.skymap' dependency (and its stack: astropy/healpy/scipy)."
        ) from e

import logging
logger = logging.getLogger(__name__)

USER_AGENT = "gw-stat/1.1 (+https://gwosc.org/)"


# ---------------------------------------------------------------------
# Basic GWOSC jsonfull download helpers (same as before)
# ---------------------------------------------------------------------
GWTC_URL_DEFAULT = "https://gwosc.org/eventapi/jsonfull/GWTC/"

def _safe_get(
    session: requests.Session,
    url: str,
    *,
    timeout=(10.0, 120.0),
    retries: int = 3,
    backoff_s: float = 0.6,
    verbose: bool = False,
    headers: Optional[dict[str, str]] = None,
    stream: bool = False,
) -> requests.Response:
    last: Optional[Exception] = None

    # Merge caller headers with our UA (caller headers win except UA is forced unless caller sets it)
    merged_headers = {"User-Agent": USER_AGENT}
    if headers:
        merged_headers.update(headers)

    for k in range(retries):
        try:
            r = session.get(url, timeout=timeout, headers=merged_headers, stream=stream)
            try:
                r.raise_for_status()
            except requests.HTTPError as e:
                raise requests.HTTPError(f"{e} (URL={r.url})") from e
            return r
        except Exception as e:  # noqa: BLE001
            last = e
            if verbose:
                print(f"[gw_stat] GET failed ({k+1}/{retries}) {type(e).__name__}: {e} :: {url}")
            time.sleep(backoff_s * (k + 1))

    raise last  # type: ignore[misc]


def catalog_key_to_jsonfull_url(catalog: str) -> str:
    cat = str(catalog).strip()
    if cat.startswith("http://") or cat.startswith("https://"):
        return cat if cat.endswith("/") else cat + "/"
    return f"https://gwosc.org/eventapi/jsonfull/{cat}/"


def fetch_gwtc_events(catalog: str):
    # Expand ALL into a merged dict of events from all catalogs (excluding ALL itself)
    if catalog == "ALL":
        merged: dict = {"events": {}}
        for c in ALLOWED_CATALOGS:
            if c == "ALL":
                continue
            raw = fetch_gwtc_events(c)
            # Defensive: accept either dict events or missing
            evs = raw.get("events", {}) if isinstance(raw, dict) else {}
            merged["events"].update(evs)
        return merged

    url = f"https://gwosc.org/eventapi/jsonfull/{catalog}/"
    s = requests.Session()
    try:
        return _safe_get(s, url, timeout=(10, 120), retries=3).json()
    except requests.HTTPError as e:
        status = getattr(e.response, "status_code", None)
        if status == 404:
            raise ValueError(
                f"Unknown catalog '{catalog}'. Allowed catalogs are: {', '.join(ALLOWED_CATALOGS)}"
            ) from e
        raise


def events_to_dataframe(raw_events: Dict[str, Any]) -> pd.DataFrame:
    rows = []
    for ev_id, ev in raw_events.items():
        rows.append(
            dict(
                event_id=ev_id,
                common_name=ev.get("commonName"),
                gps=ev.get("GPS"),
                catalog=ev.get("catalog.shortName"),
                version=ev.get("version"),
                mass_1_source=ev.get("mass_1_source"),
                mass_2_source=ev.get("mass_2_source"),
                total_mass_source=ev.get("total_mass_source"),
                final_mass_source=ev.get("final_mass_source"),
                chirp_mass_source=ev.get("chirp_mass_source"),
                redshift=ev.get("redshift"),
                luminosity_distance=ev.get("luminosity_distance"),
                chi_eff=ev.get("chi_eff"),
                chi_p=ev.get("chi_p"),
                snr=ev.get("network_matched_filter_snr"),
                far=ev.get("far"),
                p_astro=ev.get("p_astro"),
            )
        )
    return pd.DataFrame(rows)

def chirp_mass_from_masses(m1, m2):
    m1a = np.asarray(m1, dtype=float)
    m2a = np.asarray(m2, dtype=float)
    return (m1a * m2a) ** (3.0 / 5.0) / (m1a + m2a) ** (1.0 / 5.0)

def classify_binary(m1: float, m2: float, ns_threshold: float = 3.0) -> str:
    if np.isnan(m1) or np.isnan(m2):
        return "unknown"
    is_ns1 = m1 < ns_threshold
    is_ns2 = m2 < ns_threshold
    if is_ns1 and is_ns2:
        return "NS-NS"
    if is_ns1 ^ is_ns2:
        return "BH-NS"
    return "BBH"

def prepare_catalog_df(df: pd.DataFrame, ns_threshold: float = 3.0) -> pd.DataFrame:
    out = df.copy()
    out = out.dropna(subset=["mass_1_source", "mass_2_source"]).copy()
    swap = out["mass_2_source"] > out["mass_1_source"]
    if swap.any():
        out.loc[swap, ["mass_1_source", "mass_2_source"]] = out.loc[
            swap, ["mass_2_source", "mass_1_source"]
        ].to_numpy()
    out["total_mass_source"] = out["mass_1_source"] + out["mass_2_source"]
    out["q"] = out["mass_2_source"] / out["mass_1_source"]
    out["eta"] = (out["mass_1_source"] * out["mass_2_source"]) / (out["total_mass_source"] ** 2)
    missing = out["chirp_mass_source"].isna()
    if missing.any():
        out.loc[missing, "chirp_mass_source"] = chirp_mass_from_masses(
            out.loc[missing, "mass_1_source"].to_numpy(),
            out.loc[missing, "mass_2_source"].to_numpy(),
        )
    out["binary_type"] = [
        classify_binary(m1, m2, ns_threshold=ns_threshold)
        for m1, m2 in out[["mass_1_source", "mass_2_source"]].to_numpy()
    ]
    return out

def _print_skymap_members(members, max_lines=50):
    """
    Print FITS skymap filenames for inspection.
    """
    print(f"[gw_stat] Found {len(members)} FITS skymaps in tarball")
    for i, name in enumerate(members[:max_lines]):
        print("  ", name)
    if len(members) > max_lines:
        print(f"  ... ({len(members) - max_lines} more)")


# ---------------------------------------------------------------------
# Detector network (H1 / L1 / V1)
# ---------------------------------------------------------------------

GWOSC_V2_EVENTVER_URL = "https://gwosc.org/api/v2/event-versions/{name}-v{ver}"

def _parse_name_version(event_id: str, version: Optional[int]):
    """
    Extract (name, version) for GWOSC v2 endpoints.
    """
    if event_id is None:
        return None, None
    m = re.match(r"^(?P<name>.+)-v(?P<ver>\d+)$", str(event_id))
    if m:
        return m.group("name"), int(m.group("ver"))
    if version is not None:
        return str(event_id), int(version)
    return None, None

def add_detectors_and_virgo_flag(
    df: pd.DataFrame,
    *,
    column_detectors: str = "detectors",
    column_has_v1: str = "has_V1",
    progress: bool = True,
    verbose: bool = False,
    plot_network_pie: bool = True,
):
    """
    Query GWOSC v2 API to retrieve detector network per event.

    Adds columns:
      - detectors : list[str]
      - has_V1    : True if Virgo (V1) present

    Optionally plots a donut pie chart of detector-network composition
    using a legend (no labels on the pie itself).
    """

    import requests
    import pandas as pd
    from typing import Dict, Tuple, Optional

    out = df.copy()
    detectors_col = [None] * len(out)

    # ──────────────────────────────────────────────────────────────
    # Progress bar (optional)
    # ──────────────────────────────────────────────────────────────
    try:
        from tqdm.auto import tqdm  # type: ignore
    except Exception:
        tqdm = None  # type: ignore

    it = out.iterrows()
    if progress and tqdm is not None:
        it = tqdm(list(it), desc="Detector networks", total=len(out))

    # ──────────────────────────────────────────────────────────────
    # Query GWOSC (with caching)
    # ──────────────────────────────────────────────────────────────
    cache: Dict[Tuple[str, int], Optional[list]] = {}

    with requests.Session() as s:
        for idx, row in it:
            name, ver = _parse_name_version(row["event_id"], row.get("version"))
            if not name or ver is None:
                continue

            key = (name, ver)
            if key in cache:
                detectors_col[out.index.get_loc(idx)] = cache[key]
                continue

            url = GWOSC_V2_EVENTVER_URL.format(name=name, ver=ver)
            try:
                r = _safe_get(s, url, timeout=(10, 30), retries=3, verbose=verbose)
                dets = r.json().get("detectors")
                cache[key] = dets
                detectors_col[out.index.get_loc(idx)] = dets
            except Exception as e:  # noqa: BLE001
                cache[key] = None
                detectors_col[out.index.get_loc(idx)] = None
                if verbose:
                    print(f"[gw_stat] detector error for {name}-v{ver}: {e}")

    out[column_detectors] = detectors_col
    out[column_has_v1] = out[column_detectors].apply(
        lambda x: isinstance(x, list) and ("V1" in x)
    )

    # ──────────────────────────────────────────────────────────────
    # Plot detector-network pie (legend only)
    # ──────────────────────────────────────────────────────────────
    fig = None
    if plot_network_pie:
        import matplotlib.pyplot as plt

        allowed = {"H1", "L1", "V1", "G1"}

        def classify_network(dets):
            if not isinstance(dets, list):
                return None

            s = set(dets) & allowed
            if not s:
                return None

            if len(s) == 1:
                return "1 detector"
            if len(s) == 4:
                return "H1–L1–V1-G1"
            if {"H1", "L1", "V1"}.issubset(s):
                return "H1–L1–V1"
            if {"H1", "L1"}.issubset(s) and "V1" not in s:
                return "H1–L1"
            if "V1" in s and (("H1" in s) ^ ("L1" in s)) and len(s) == 2:
                return "H1/L1–V1"

            return None  # omit rare combinations

        order = [
            "1 detector",
            "H1–L1",
            "H1/L1–V1",
            "H1–L1–V1",
            "H1–L1–V1-G1",
        ]

        counts = (
            out[column_detectors]
            .apply(classify_network)
            .dropna()
            .value_counts()
            .reindex(order, fill_value=0)
        )

        counts = counts[counts > 0]
        if counts.sum() > 0:
            colors = {
                "1 detector": "#9467bd",
                "H1–L1": "#ff7f0e",
                "H1/L1–V1": "#2ca02c",
                "H1–L1–V1": "#1f77b4",
                "H1–L1–V1-G1": "#d62728",
            }

            fig, ax = plt.subplots(figsize=FIGSIZE)
            fig.subplots_adjust(right=0.78)
            wedges, _, autotexts = ax.pie(
                counts.values,
                labels=None,
                colors=[colors[k] for k in counts.index],
                startangle=90,
                counterclock=False,
                wedgeprops=dict(width=0.42),
                autopct=lambda p: f"{p:.1f}%" if p >= 3 else "",
                pctdistance=0.78,
            )

            # Center label
            total = int(counts.sum())
            ax.text(
                0, 0,
                f"N = {total}",
                ha="center", va="center",
                fontsize=13,
                fontweight="bold",
            )

            # Legend
            ax.legend(
                wedges,
                [f"{k} ({counts[k]})" for k in counts.index],
                title="Detector network",
                loc="center left",
                bbox_to_anchor=(1.02, 0.5),
                frameon=False,
            )

            ax.set_title("Events network detection", fontsize=14, pad=18)
            ax.axis("equal")

        elif verbose:
            print("[gw_stat] No detector-network data to plot.")

    return out, fig




# ---------------------------------------------------------------------
# Credible area from Zenodo tarballs
# ---------------------------------------------------------------------

# ---------------------------------------------------------------------
# Galaxy collection skymaps helpers (directory of symlinks -> staged blobs)
# ---------------------------------------------------------------------

from pathlib import Path

def _galaxy_inputs_blob_dir() -> Path:
    """Return the job-local staged inputs directory if present."""
    wd = Path(os.environ.get("GALAXY_WORKING_DIR", "."))
    job_dir = wd.parent  # .../jobs/<jobid>
    return job_dir / "inputs"

def build_skymap_index_from_directory(skymap_dir: str | Path, verbose: bool = False):
    """
    Build an index from a Galaxy collection directory containing symlinks to FITS skymaps.
    The symlink targets may point to an object store path not mounted in the runtime.
    We redirect those to the staged blobs in <job>/inputs/dataset_<uuid>.dat.
    Returns dict[(GWYYMMDD_HHMMSS, approximant_key)] -> Path(blob_file)
    """
    skymap_dir = Path(skymap_dir)
    if not skymap_dir.exists():
        raise FileNotFoundError(f"Skymap directory not found: {skymap_dir}")

    blob_dir = _galaxy_inputs_blob_dir()
    blob_map = {}
    if blob_dir.exists():
        for p in blob_dir.glob("dataset_*.dat"):
            blob_map[p.name] = p

    files = list(skymap_dir.rglob("*.fit*"))
    index = {}
    redirected = broken = 0

    for link_path in files:
        link_name = link_path.name  # contains GW... and approximant
        p = link_path

        if p.is_symlink():
            # First try to use the symlink target directly (local runs / staged symlinks)
            try:
                t = p.readlink()
                target = (p.parent / t) if not Path(t).is_absolute() else Path(t)
                if target.exists() and target.is_file():
                    p = target
                else:
                    # Fall back to Galaxy redirection: target name is dataset_<uuid>.dat
                    tname = Path(t).name
                    real = blob_map.get(tname)
                    if real is None or not real.exists():
                        broken += 1
                        continue
                    p = real
                    redirected += 1
            except Exception:
                broken += 1
                continue
        else:
            if not p.exists():
                continue

        m = re.search(r"(GW\d{6}_\d{6})", link_name)
        if not m:
            continue
        ev_key = m.group(1)

        low = link_name.lower()
        if "mixed" in low:
            approx = "Mixed"
        else:
            post = link_name.split(ev_key, 1)[-1]
            post = post.lstrip("_- .")
            token = re.split(r"[_.\- ]+", post)[0] if post else ""
            approx = token if token else "Unknown"

        index.setdefault((ev_key, approx), p)

    if verbose:
        print(f"[gw_stat] directory index entries: {len(index)} redirected={redirected} broken_links={broken}")

    if not index:
        raise FileNotFoundError(
            f"No readable skymap datasets found under {skymap_dir}. "
            f"(redirected={redirected}, broken_links={broken}, blob_dir={blob_dir})"
        )

    return index


def select_skymap_path(index: dict, ev_key: str, prefer: str = "Mixed") -> Path | None:
    if (ev_key, prefer) in index:
        return index[(ev_key, prefer)]
    for (e, _a), p in index.items():
        if e == ev_key:
            return p
    return None

def add_localization_area_from_directory(
    df: pd.DataFrame,
    *,
    skymap_dir: str | Path,
    cred: float = 0.9,
    column: Optional[str] = None,
    progress: bool = True,
    verbose: bool = False,
) -> pd.DataFrame:
    """
    Compute A{cred*100} sky area for events by reading skymaps from a Galaxy collection directory.
    The directory may contain symlinks; we redirect to staged blobs automatically.

    Returns a copy of df with `column` filled where possible.
    """
    index = build_skymap_index_from_directory(skymap_dir, verbose=verbose)

    out = df.copy()
    areas = [np.nan] * len(out)

    try:
        from tqdm.auto import tqdm  # type: ignore
    except Exception:
        tqdm = None  # type: ignore

    it = out.iterrows()
    if progress and tqdm is not None:
        it = tqdm(list(it), desc=f"A{int(cred*100)} from local skymaps", total=len(out))

    found = ok = miss = errors = 0

    for idx_row, row in it:
        ev_id = row.get("event_id")
        ev = _normalize_event_name(ev_id)
        if ev is None:
            miss += 1
            continue

        m = re.search(r"(GW\d{6}_\d{6})", str(ev))
        if not m:
            miss += 1
            continue
        ev_key = m.group(1)

        p = select_skymap_path(index, ev_key, prefer="Mixed")
        if p is None:
            miss += 1
            continue

        found += 1
        try:
            raw = Path(p).read_bytes()
            prob = _read_sky_map_from_bytes(raw)
            area = _credible_area_from_probmap(prob, cred=cred)
            areas[out.index.get_loc(idx_row)] = float(area)
            ok += 1
        except Exception as e:  # noqa: BLE001
            errors += 1
            print(f"[gw_stat] A{int(cred*100)} error for {ev_id} using {p}: {type(e).__name__}: {e}")

        if tqdm is not None and progress:
            try:
                it.set_postfix(found=found, ok=ok, miss=miss, errors=errors)  # type: ignore[attr-defined]
            except Exception:
                pass

    print(f"[gw_stat] A{int(cred*100)} summary (directory): found={found} ok={ok} miss={miss} errors={errors}")
    if column is None:
        level = int(round(100 * cred))
        column = f"A{level}_deg2"
    out[column] = areas
    return out

def _credible_area_from_probmap(prob: np.ndarray, cred: float = 0.9) -> float:
    _require_ligo_skymap()
    from ligo.skymap.postprocess import find_greedy_credible_levels  # type: ignore
    cls = find_greedy_credible_levels(prob)
    mask = cls <= cred
    pixel_area = 41252.96 / len(prob)  # deg^2
    return float(mask.sum() * pixel_area)

def _read_sky_map_from_bytes(data: bytes) -> np.ndarray:
    import gzip
    import io
    import numpy as np
    _require_ligo_skymap()
    from ligo.skymap.io import read_sky_map  # type: ignore

    # If gzipped FITS, decompress first (gzip magic bytes: 1f 8b)
    if len(data) >= 2 and data[0] == 0x1F and data[1] == 0x8B:
        data = gzip.decompress(data)

    bio = io.BytesIO(data)

    # Newer ligo.skymap supports ignore_missing_simple; older versions don't.
    try:
        prob, _meta = read_sky_map(bio, nest=None, ignore_missing_simple=True)
    except TypeError as e:
        if "ignore_missing_simple" not in str(e):
            raise
        bio.seek(0)
        prob, _meta = read_sky_map(bio, nest=None)

    return np.asarray(prob, dtype=float)

def _normalize_event_name(ev: str) -> str:
    # event_id in jsonfull looks like "GW200129_065458-v2"
    # skymap files usually include the base event name without "-vN"
    ev = str(ev)
    return re.sub(r"-v\d+$", "", ev)
    
def download_zenodo_skymaps_tarball(
    catalog_key: str,
    *,
    dest_path: Optional[str] = None,
    cache_dir: str = ".cache_gwosc",
    progress: bool = True,
    verbose: bool = False,
) -> str:
    """
    Download the official Zenodo skymap tarball for a confident catalog key.
    Returns local filepath.

    The file is cached on disk; if already present, it is not downloaded again.

    Adds streaming byte-progress via tqdm and supports resume (HTTP Range) when possible.
    """
    url = zenodo_skymap_url(catalog_key)

    # Respect explicit cache_dir argument (your previous code overwrote it)
    os.makedirs(cache_dir, exist_ok=True)

    if dest_path is None:
        safe = catalog_key.replace("/", "_")
        dest_path = os.path.join(cache_dir, f"{safe}_skymaps.tar.gz")

    # If already downloaded, use it
    if os.path.exists(dest_path) and os.path.getsize(dest_path) > 0:
        if verbose:
            print(f"[gw_stat] using cached tarball: {dest_path}")
        return dest_path

    try:
        from tqdm.auto import tqdm  # type: ignore
    except Exception:
        tqdm = None  # type: ignore

    part_path = dest_path + ".part"
    existing = os.path.getsize(part_path) if os.path.exists(part_path) else 0

    headers = {}
    if existing > 0:
        headers["Range"] = f"bytes={existing}-"
        if verbose:
            print(f"[gw_stat] resuming tarball download at byte {existing}: {part_path}")

    with requests.Session() as s:
        # NOTE: _safe_get must accept headers=... and stream=True, or you can pass stream=True via _safe_get.
        # If your _safe_get doesn't support these kwargs yet, add them there.
        r = _safe_get(
            s,
            url,
            timeout=(10, 1200),
            retries=3,
            verbose=verbose,
            headers=headers,   # e.g. {"Range": "bytes=123-"} or {}
            stream=True,
        )


        # If server ignored Range and returned 200, restart from scratch
        if existing > 0 and r.status_code == 200:
            if verbose:
                print("[gw_stat] server ignored Range; restarting full download")
            existing = 0
            headers = {}
            # restart request cleanly
            r.close()
            r = _safe_get(
                s,
                url,
                timeout=(10, 1200),
                retries=3,
                verbose=verbose,
                headers=headers,
                stream=True,
            )

        # Determine total size (best effort)
        total = None
        cl = r.headers.get("Content-Length")
        if cl:
            try:
                total = int(cl)
            except ValueError:
                total = None

        # If resuming and status 206, Content-Length is remaining; convert to full size for tqdm
        if existing > 0 and r.status_code == 206 and total is not None:
            total = existing + total

        pbar = None
        if progress and tqdm is not None:
            pbar = tqdm(
                total=total,
                initial=existing if (total is not None and existing > 0) else 0,
                unit="B",
                unit_scale=True,
                unit_divisor=1024,
                desc=f"Download {catalog_key} skymaps",
                leave=False,
            )

        # Write to .part then atomically move into place at the end
        mode = "ab" if existing > 0 else "wb"
        with open(part_path, mode) as f:
            for chunk in r.iter_content(chunk_size=1024 * 1024):
                if not chunk:
                    continue
                f.write(chunk)
                if pbar is not None:
                    pbar.update(len(chunk))

        if pbar is not None:
            pbar.close()

    # Finalize
    os.replace(part_path, dest_path)
    if verbose:
        print(f"[gw_stat] downloaded tarball: {dest_path}")

    return dest_path

def build_skymap_index_from_tar(
    tar_path: str,
    *,
    verbose: bool = False,
):
    import os
    import re
    import tarfile

    index = {}
    ev_re = re.compile(r"(GW\d{6}_\d{6}|GW\d{6}|S\d{6}[a-z]?)", re.IGNORECASE)

    with tarfile.open(tar_path, mode="r:gz") as tf:
        for mem in tf.getmembers():
            name = mem.name
            lname = name.lower()

            # Accept both FITS and gzipped FITS
            if not (lname.endswith(".fits") or lname.endswith(".fits.gz")):
                continue

            matches = ev_re.findall(name)
            if not matches:
                continue
            event_id = max(matches, key=len)  # prefer GWYYYYMM_DDHHMM over GWYYYYMM

            base = os.path.basename(name)

            # Strip extension for parsing
            lbase = base.lower()
            if lbase.endswith(".fits.gz"):
                stem = base[:-8]   # remove ".fits.gz"
            elif lbase.endswith(".fits"):
                stem = base[:-5]   # remove ".fits"
            else:
                stem = base

            parts = stem.split(":")

            approximant = "Unknown"
            if len(parts) >= 3:
                approximant = f"{parts[-2]}:{parts[-1]}"
            elif len(parts) == 2:
                approximant = parts[-1]

            key = (event_id, approximant)

            if key in index:
                if verbose:
                    kept = index[key]
                    print(
                        f"[gw_stat] duplicate skymap for {event_id} approximant {approximant}: "
                        f"kept='{kept}' skipped='{name}'"
                    )
                continue

            index[key] = name

            if verbose:
                print(f"[gw_stat] indexed: event={event_id} approx={approximant} member='{name}'")

    if verbose:
        n_events = len({ev for ev, _ in index.keys()})
        print(
            f"[gw_stat] Indexed {len(index)} skymaps for {n_events} events from {os.path.basename(tar_path)}"
        )

    return index



def select_skymap_member(
    index: dict,
    event_id: str,
    *,
    prefer: str = "Mixed",
):
    """
    Select one skymap tar member for an event from an index keyed by
    (event_id, approximant).

    Preference order:
    1) Exact match on preferred approximant
    2) Any available approximant for the event

    Returns
    -------
    member : str or None
    """
    # Preferred approximant
    key = (event_id, prefer)
    if key in index:
        return index[key]

    # Fallback: first available approximant
    for (ev, _approx), member in index.items():
        if ev == event_id:
            return member

    return None

    
import re
import json
import numpy as np
import pandas as pd
from minio import Minio
import os
import json

def _catalog_to_s3_prefix(catalog_key: str) -> str:
    # catalogs are already canonical: GWTC-2.1, GWTC-3, GWTC-4
    return f"{catalog_key.rstrip('/')}/"

def _s3_client_from_env() -> Minio:
    """
    Connect to S3/MinIO.

    - If S3_CREDENTIALS env var is set (JSON), use it.
    - Otherwise use anonymous development endpoint (public).
    """
    credentials_env = os.environ.get("S3_CREDENTIALS")
    if credentials_env:
        credentials = json.loads(credentials_env)
    else:
        # anonymous development server
        credentials = {
            "endpoint": "minio-dev.odahub.fr",
            "secure": True,
        }

    return Minio(
        endpoint=credentials["endpoint"],
        secure=credentials.get("secure", True),
        access_key=credentials.get("access_key"),
        secret_key=credentials.get("secret_key"),
    )



def _catalog_to_s3_prefix(catalog_key: str) -> str:
    """
    Map your CLI catalog keys to S3 folder names.
    Adjust if you want GWTC-3 to point to GWTC-3/, etc.
    """
    if catalog_key.startswith("GWTC-2.1"):
        return "GWTC-2.1/"
    if catalog_key.startswith("GWTC-3"):
        return "GWTC-3/"
    if catalog_key.startswith("GWTC-4"):
        return "GWTC-4/"
    # fallback: try exact
    return f"{catalog_key}/"


def _build_s3_skymap_index(mc, bucket: str, prefix: str, verbose: bool = False) -> dict[str, str]:
    """
    Return dict: ev_key ('GWYYYYMM_DDHHMM') -> object_name.
    If multiple matches exist for same event, keep the first (or refine later).
    """
    index: dict[str, str] = {}
    rx = re.compile(r"(GW\d{6}_\d{6})")

    for obj in mc.list_objects(bucket, prefix=prefix, recursive=True):
        name = obj.object_name
        if not (name.endswith(".fits") or name.endswith(".fits.gz")):
            continue

        m = rx.search(name)
        if not m:
            continue

        ev_key = m.group(1)
        if ev_key not in index:
            index[ev_key] = name
            if verbose:
                print(f"[gw_stat] S3 index: {ev_key} -> {name}")

    return index


def add_localization_area_from_s3(
    df: pd.DataFrame,
    *,
    catalog_key: str,
    cred: float = 0.9,
    column: str | None = None,
    bucket: str = "gwtc",
    base_prefix: str = "",   # set to "gwtc/" if you want: odahub/gwtc/...
    cache_index: bool = True,
    progress: bool = True,
    verbose: bool = False,
) -> pd.DataFrame:
    """
    Compute credible area by reading per-event skymap FITS from S3 (recursive scan).

    - Lists all .fits/.fits.gz under s3://bucket/<base_prefix>/<GWTC-*/...> recursively.
    - Builds an event->object index from filenames containing GWYYYYMM_DDHHMM.
    - Streams each FITS as bytes (no full download needed).
    """
    out = df.copy()
    if column is None:
        level = int(round(100 * cred))
        column = f"A{level}_deg2"
    out[column] = np.nan

    mc = _s3_client_from_env()
    cat_prefix = _catalog_to_s3_prefix(catalog_key)
    prefix = f"{base_prefix}{cat_prefix}"

    # build (or optionally cache) index
    index = _build_s3_skymap_index(mc, bucket=bucket, prefix=prefix, verbose=verbose)

    try:
        from tqdm.auto import tqdm  # type: ignore
    except Exception:
        tqdm = None  # type: ignore

    it = out.iterrows()
    if progress and tqdm is not None:
        it = tqdm(list(it), desc=f"A{int(cred*100)} from S3 {catalog_key}", total=len(out))

    found = ok = miss = errors = 0
    areas = [np.nan] * len(out)

    for idx, row in it:
        ev_id = row.get("event_id")
        ev = _normalize_event_name(ev_id)
        if ev is None:
            miss += 1
            continue

        m = re.search(r"(GW\d{6}_\d{6})", str(ev))
        if not m:
            miss += 1
            continue
        ev_key = m.group(1)

        obj_name = index.get(ev_key)
        if obj_name is None:
            miss += 1
            continue

        found += 1
        try:
            if verbose:
                print(f"[gw_stat] fetching {ev_id} from s3://{bucket}/{obj_name}")

            resp = mc.get_object(bucket, obj_name)
            try:
                data = resp.read()  # bytes
            finally:
                resp.close()
                resp.release_conn()

            prob = _read_sky_map_from_bytes(data)
            area = _credible_area_from_probmap(prob, cred=cred)
            areas[out.index.get_loc(idx)] = float(area)
            ok += 1

        except Exception as e:  # noqa: BLE001
            errors += 1
            print(f"[gw_stat] A{int(cred*100)} S3 error for {ev_id} using {obj_name}: {type(e).__name__}: {e}")

        if tqdm is not None and progress:
            try:
                it.set_postfix(found=found, ok=ok, miss=miss, errors=errors)  # type: ignore[attr-defined]
            except Exception:
                pass

    print(
        f"[gw_stat] A{int(cred*100)} S3 summary: "
        f"found={found} ok={ok} miss={miss} errors={errors}"
    )

    out[column] = areas
    return out


def add_localization_area_from_zenodo(
    df: pd.DataFrame,
    *,
    catalog_key: str,
    cred: float = 0.9,
    column: str | None = None,
    cache_dir: str = ".cache_gwosc",
    progress: bool = True,
    verbose: bool = False,
) -> pd.DataFrame:
    """
    Compute credible area for events by reading the official Zenodo skymap tarball for the given confident catalog.

    This bypasses the GWOSC v2 parameters endpoint entirely and should yield credible area for (almost) all events
    whose skymap FITS is included in the tarball.

    Returns a copy of df with `column` filled where possible.
    """
    tar_path = download_zenodo_skymaps_tarball(
        catalog_key, cache_dir=cache_dir, progress=progress, verbose=verbose
    )
    index = build_skymap_index_from_tar(tar_path, verbose=verbose)

    out = df.copy()
    if column is None:
        level = int(round(100 * cred))
        column = f"A{level}_deg2"

    areas = [np.nan] * len(out)

    try:
        from tqdm.auto import tqdm  # type: ignore
    except Exception:
        tqdm = None  # type: ignore

    it = out.iterrows()
    if progress and tqdm is not None:
        it = tqdm(list(it), desc=f"A{int(cred*100)} from Zenodo {catalog_key}", total=len(out))

    found = ok = miss = errors = 0

    with tarfile.open(tar_path, mode="r:gz") as tf:
        for idx, row in it:
            ev_id = row.get("event_id")
            ev = _normalize_event_name(ev_id)

            # fallback: try prefix match on GW id fragment
            if ev is None:
                miss += 1
                continue

            m = re.search(r"(GW\d{6}_\d{6})", str(ev))
            if not m:
                miss += 1
                continue
            ev_key = m.group(1)

            # Choose ONE member from the (event, approximant) index
            member = select_skymap_member(index, ev_key, prefer="Mixed")
            if member is None:
                miss += 1
                continue

            found += 1

            # Determine which approximant key was used (for logging)
            approx_used = None
            if verbose:
                for (e, a), mem in index.items():
                    if e == ev_key and mem == member:
                        approx_used = a
                        break

                if (ev_key, "Mixed") in index:
                    print(f"[gw_stat] A{int(cred*100)} selected Mixed for {ev_id}: {member}")
                else:
                    print(
                        f"[gw_stat] A{int(cred*100)} Mixed not found for {ev_id}; "
                        f"selected {approx_used}: {member}"
                    )

            try:
                # The exact file extracted (always printed if verbose)
                if verbose:
                    print(f"[gw_stat] extracting A{int(cred*100)} for {ev_id}: {member}")

                f = tf.extractfile(member)
                if f is None:
                    raise OSError(f"could not extract {member}")

                prob = _read_sky_map_from_bytes(f.read())
                area = _credible_area_from_probmap(prob, cred=cred)
                areas[out.index.get_loc(idx)] = float(area)
                ok += 1

            except Exception as e:  # noqa: BLE001
                errors += 1
                # Print at least one error line (otherwise NaNs are mysterious)
                print(f"[gw_stat] A{int(cred*100)} error for {ev_id} using {member}: {type(e).__name__}: {e}")

            if tqdm is not None and progress:
                try:
                    it.set_postfix(found=found, ok=ok, miss=miss, errors=errors)  # type: ignore[attr-defined]
                except Exception:
                    pass

    # Summary (correct scope: after the event loop)
    print(
        f"[gw_stat] A{int(cred*100)} summary: "
        f"found={found} ok={ok} miss={miss} errors={errors}"
    )

    out[column] = areas
    return out



if __name__ == "__main__":
    # tiny smoke test
    raw = fetch_gwtc_events(catalog="GWTC-3")
    df0 = events_to_dataframe(raw["events"])
    df = prepare_catalog_df(df0)
    print("rows:", len(df))
