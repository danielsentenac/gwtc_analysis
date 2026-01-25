from __future__ import annotations

"""
GW event parameter-estimation (PE) plotting pipeline.

This module provides a Galaxy/MMODA-friendly entry point:

    run_parameters_estimation(...) -> dict

Supported data repositories (data_repo):
  - "s3"     : MinIO bucket "gwtc"
  - "zenodo" : public Zenodo PEDataRelease records
  - "galaxy" : PE files staged by Galaxy collections under ./galaxy_inputs
              (expected collection directories: GWTC-2.1-PE, GWTC-3-PE, GWTC-4-PE)
  - "local"  : look for PEDataRelease file in plots_dir (kept for dev/debug)

Notes
-----
- The control-flow uses a notebook-like `go_next_cell` flag to preserve behaviour/logs.
- The last run outputs are stored in the module global `LAST_RESULT`.
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional

import html
import json
import re
import warnings

import requests


# ---------------------------------------------------------------------
# Result container
# ---------------------------------------------------------------------


@dataclass
class PEResult:
    """Collected outputs from a PE run."""
    fig_distribution: List[Any] = field(default_factory=list)
    fig_strain: List[Any] = field(default_factory=list)
    fig_psd: List[Any] = field(default_factory=list)
    fig_skymap: List[Any] = field(default_factory=list)
    tool_log: List[str] = field(default_factory=list)

    # raw filenames (always populated, even if ODA is unavailable)
    files_distribution: List[str] = field(default_factory=list)
    files_strain: List[str] = field(default_factory=list)
    files_psd: List[str] = field(default_factory=list)
    files_skymap: List[str] = field(default_factory=list)


LAST_RESULT: Optional[PEResult] = None


# ---------------------------------------------------------------------
# Zenodo PE helpers
# ---------------------------------------------------------------------

ZENODO_PE_RECORD_IDS = [6513631, 8177023, 17014085]
USER_AGENT = "gwtc_analysis (https://github.com/danielsentenac/gwtc_analysis)"

def _tqdm_or_none():
    try:
        from tqdm.auto import tqdm  # type: ignore
        return tqdm
    except Exception:
        return None


def _download_http_with_progress(
    url: str,
    out_path: Path,
    *,
    desc: str = "download",
    chunk_size: int = 1024 * 1024,
) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)

    headers = {"User-Agent": USER_AGENT, "Accept": "*/*"}
    tqdm = _tqdm_or_none()

    with requests.get(url, stream=True, headers=headers, allow_redirects=True, timeout=(10, 300)) as r:
        r.raise_for_status()
        total = r.headers.get("Content-Length")
        total_size = int(total) if total and total.isdigit() else None

        pbar = tqdm(total=total_size, unit="B", unit_scale=True, desc=desc) if tqdm else None
        try:
            with out_path.open("wb") as f:
                for chunk in r.iter_content(chunk_size=chunk_size):
                    if not chunk:
                        continue
                    f.write(chunk)
                    if pbar:
                        pbar.update(len(chunk))
        finally:
            if pbar:
                pbar.close()


def download_zenodo_pe_file(
    chosen: dict[str, Any],
    *,
    outdir: Path,
    progress_cb: Callable[[str, int, str], None] | None = None,
    log_cb: Callable[[str], None] | None = None,
) -> Path:
    """
    Download a PE file from Zenodo using the public web URL:

      https://zenodo.org/records/<rid>/files/<filename>?download=1

    `chosen` must contain: {'filename', 'url'}.
    """
    fname = str(chosen["filename"])
    url = str(chosen["url"])

    outdir.mkdir(parents=True, exist_ok=True)
    out_path = outdir / fname

    if out_path.exists() and out_path.stat().st_size > 0:
        if log_cb:
            log_cb(f"ℹ️ [CACHE] Using existing Zenodo PE file: {out_path}")
        return out_path

    if log_cb:
        log_cb(f"ℹ️ [DOWNLOAD] Zenodo PE: {fname}")
    if progress_cb:
        progress_cb("Download data", 15, f"Zenodo download: {fname}")

    _download_http_with_progress(url, out_path, desc=f"Zenodo: {fname}")
    return out_path


def _extract_event_id_from_filename(name: str) -> str | None:
    m = re.search(r"(GW\d{6}_\d{6})", name)
    return m.group(1) if m else None


def _zenodo_record_files(record_id: int) -> list[dict[str, Any]]:
    url = f"https://zenodo.org/api/records/{record_id}"
    r = requests.get(url, timeout=(10, 120), headers={"User-Agent": USER_AGENT})
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
    Build/return an index:
        index[event_id] = [{'record_id','filename','url'}, ...]

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
        try:
            files = _zenodo_record_files(rid)
        except Exception as e:
            # Zenodo can transiently fail (503, etc.). Do not fail the whole run.
            print(f"[pe][zenodo] WARN: could not fetch record {rid}: {type(e).__name__}: {e}")
            continue

        for f in files:
            fname = f.get("key") or f.get("filename") or ""
            ev = _extract_event_id_from_filename(fname)
            if not ev:
                continue

            url = f"https://zenodo.org/records/{rid}/files/{fname}?download=1"
            index.setdefault(ev, []).append({"record_id": rid, "filename": fname, "url": url})

    cache_path.write_text(json.dumps(index, indent=2), encoding="utf-8")
    return index


def choose_best_pe_file(candidates: list[dict[str, Any]]) -> dict[str, Any] | None:
    """
    Simple selection rule:
      - prefer 'mixed_cosmo' when exists
      - otherwise 'recombined'
      - otherwise prefer .hdf5/.h5
      - else first
    """
    if not candidates:
        return None

    def low(s: str) -> str:
        return (s or "").lower()

    for c in candidates:
        if "mixed_cosmo" in low(c.get("filename", "")):
            return c
    for c in candidates:
        if "recombined" in low(c.get("filename", "")):
            return c
    for ext in (".hdf5", ".h5"):
        for c in candidates:
            if low(c.get("filename", "")).endswith(ext):
                return c
    return candidates[0]


# ---------------------------------------------------------------------
# S3 helpers
# ---------------------------------------------------------------------


def _download_s3_with_progress(
    client,
    *,
    bucket: str,
    object_name: str,
    out_path: Path,
    desc: str = "s3 download",
    chunk_size: int = 1024 * 1024,
) -> None:
    """Download from MinIO/S3 with a tqdm byte progress bar."""
    out_path.parent.mkdir(parents=True, exist_ok=True)

    tqdm = _tqdm_or_none()

    total_size: int | None = None
    try:
        st = client.stat_object(bucket, object_name)
        total_size = int(getattr(st, "size", 0)) or None
    except Exception:
        total_size = None

    pbar = tqdm(total=total_size, unit="B", unit_scale=True, desc=desc) if tqdm else None
    resp = client.get_object(bucket, object_name)
    try:
        with out_path.open("wb") as f:
            while True:
                chunk = resp.read(chunk_size)
                if not chunk:
                    break
                f.write(chunk)
                if pbar:
                    pbar.update(len(chunk))
    finally:
        try:
            resp.close()
            resp.release_conn()
        except Exception:
            pass
        if pbar:
            pbar.close()


# ---------------------------------------------------------------------
# Plot ordering + report helpers
# ---------------------------------------------------------------------


def _is_skymap_plot(p: Path) -> bool:
    name = p.name.lower()
    return any(k in name for k in ("skymap", "localization", "allsky", "radec", "ra_dec", "sky"))


def _pe_plot_sort_key(p: Path) -> tuple[int, int, str]:
    """Sort plots in the HTML report into a stable, readable order."""
    name = p.name.lower()

    # 0) Posteriors / corners
    if any(k in name for k in ("corner", "posterior", "posteriors", "credible", "ci", "marginal")):
        group = 0
    # 1) PSD
    elif any(k in name for k in ("psd", "noise", "spectrum", "asd")):
        group = 1
    # 2) Strain / waveform / time-frequency
    elif any(k in name for k in ("strain", "waveform", "qtransform", "q-transform", "spectrogram", "timefreq", "tf")):
        group = 2
    # 3) Skymaps / localization (last)
    elif any(k in name for k in ("skymap", "localization", "allsky", "sky", "radec", "ra_dec")):
        group = 3
    else:
        group = 9

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
    posterior_keys_by_label: dict[str, list[str]] | None = None,
    missing_pe_vars_by_label: dict[str, list[str]] | None = None,
    requested_pe_vars: list[str] | None = None,
    posterior_pairs_missing_by_label: dict[str, list[str]] | None = None,
    requested_pe_pairs: list[str] | None = None,
    event_logs: list[str] | None = None,
) -> None:
    """Write a simple (non-embedded) HTML report referencing local PNGs."""
    out_path.parent.mkdir(parents=True, exist_ok=True)

    def rel(p: Path) -> str:
        try:
            return str(p.relative_to(out_path.parent))
        except Exception:
            return str(p)

    def _norm_key_for_pairing(stem: str) -> str:
        s = stem.lower()
        for tok in ("_ra", "_dec", "ra_", "dec_", "-ra", "-dec", "ra-", "dec-"):
            s = s.replace(tok, "_")
        parts = [p for p in s.replace("-", "_").split("_") if p and p not in ("ra", "dec")]
        return "_".join(parts) if parts else s

    def _is_ra_plot(p: Path) -> bool:
        n = p.name.lower()
        return ("_ra" in n) or n.startswith("ra_") or n.endswith("_ra.png") or (n == "ra.png")

    def _is_dec_plot(p: Path) -> bool:
        n = p.name.lower()
        return ("_dec" in n) or n.startswith("dec_") or n.endswith("_dec.png") or (n == "dec.png")

    # Optional: extract only the "Label report" section (fallback to full log)
    import re

    def _extract_label_report(log_text: str) -> str:
        m = re.search(
            r"(Label report:\n.*?Available labels in PE file.*?\n(?:\s+- .*?\n)+)",
            log_text,
            flags=re.S,
        )
        return m.group(1).strip() if m else ""

    lines: list[str] = []

    # --- HTML header ---
    lines.append("<!doctype html>")
    lines.append("<html>")
    lines.append("<head>")
    lines.append("""
        <meta charset="utf-8">
        <title>Parameters estimation</title>

        <style>
            .pe-label-report {
                font-size: 12px;
                line-height: 1.25;
                padding: 10px;
                background: #f6f8fa;
                border: 1px solid #d0d7de;
                border-radius: 6px;
                overflow-x: auto;
                white-space: pre;
                font-family: monospace;
            }
        </style>
    """)
    lines.append("</head>")
    lines.append("<body>")

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

    # Logs report
    if event_logs:
        full_log = "\n".join(event_logs)
        label_block = _extract_label_report(full_log)
        text_to_show = label_block if label_block else full_log

        lines.append("<h2>Log report</h2>")
        lines.append("<pre class='pe-label-report'>")
        lines.append(html.escape(text_to_show))
        lines.append("</pre>")

    # Output summary
    lines.append("<h2>Outputs</h2>")
    lines.append("<ul>")
    lines.append(f"<li>Plots: {len(plots)}</li>")
    lines.append(f"<li>Other files: {len(files)}</li>")
    lines.append("</ul>")

    if files:
        lines.append("<h3>Files</h3><ul>")
        for p in files:
            rp = html.escape(rel(p))
            lines.append(f"<li><a href='{rp}'><code>{rp}</code></a></li>")
        lines.append("</ul>")

    # Posterior keys section
    posterior_keys_by_label = posterior_keys_by_label or {}
    missing_pe_vars_by_label = missing_pe_vars_by_label or {}
    requested_pe_vars = requested_pe_vars or []
    posterior_pairs_missing_by_label = posterior_pairs_missing_by_label or {}
    requested_pe_pairs = requested_pe_pairs or []

    if posterior_keys_by_label:
        lines.append("<h2>Posterior samples</h2>")

        if requested_pe_vars:
            lines.append(
                "<p><b>Requested extra posterior variables:</b> "
                + html.escape(", ".join(requested_pe_vars))
                + "</p>"
            )
        if requested_pe_pairs:
            lines.append(
                "<p><b>Requested 2D posterior pairs:</b> "
                + html.escape(", ".join(requested_pe_pairs))
                + "</p>"
            )

        for label, keys in sorted(posterior_keys_by_label.items(), key=lambda kv: kv[0]):
            missing = missing_pe_vars_by_label.get(label, [])
            missing_pairs = posterior_pairs_missing_by_label.get(label, [])

            lines.append(f"<h3>{html.escape(label)}</h3>")
            lines.append(f"<p>Available posterior keys: <b>{len(keys)}</b></p>")

            if missing:
                lines.append(
                    "<p style='color:#b00020;'><b>Requested keys not found:</b> "
                    + html.escape(", ".join(missing))
                    + "</p>"
                )
            if missing_pairs:
                lines.append(
                    "<p style='color:#b00020;'><b>Requested pairs not plotted:</b> "
                    + html.escape(", ".join(missing_pairs))
                    + "</p>"
                )

            lines.append("<details>")
            lines.append("<summary>Show/hide full list</summary>")
            lines.append(
                "<pre style='white-space: pre-wrap; background:#f6f8fa; padding:10px; "
                "border-radius:8px; border:1px solid #ddd;'>"
                + html.escape("\n".join(keys))
                + "</pre>"
            )
            lines.append("</details>")

    if not plots:
        lines.append("<p><b>No plots found</b> in the plots directory.</p>")
        lines.append("</body></html>\n")
        out_path.write_text("\n".join(lines), encoding="utf-8")
        return

    # -----------------------
    # Plots
    # -----------------------
    lines.append("<h2>Plots</h2>")

    # 1) Pair RA/DEC plots on the same row
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

            if ra_p is not None:
                used.add(ra_p)
                rp = html.escape(rel(ra_p))
                lines.append(
                    "<td style='padding:10px; vertical-align:top;'>"
                    f"<a href='{rp}'><img src='{rp}' style='max-width:450px; height:auto;'/></a>"
                    f"<br/><code>{rp}</code></td>"
                )
            else:
                lines.append("<td style='padding:10px;'><i>missing RA plot</i></td>")

            if dec_p is not None:
                used.add(dec_p)
                dp = html.escape(rel(dec_p))
                lines.append(
                    "<td style='padding:10px; vertical-align:top;'>"
                    f"<a href='{dp}'><img src='{dp}' style='max-width:450px; height:auto;'/></a>"
                    f"<br/><code>{dp}</code></td>"
                )
            else:
                lines.append("<td style='padding:10px;'><i>missing DEC plot</i></td>")

            lines.append("</tr>")
        lines.append("</table>")

    remaining = [p for p in plots if p not in used]

    # Skymaps (other than RA/DEC) — render full-width
    skymaps = [p for p in remaining if _is_skymap_plot(p)]
    others = [p for p in remaining if p not in skymaps]

    for p in others:
        rp = html.escape(rel(p))
        lines.append(f"<p><a href='{rp}'><code>{rp}</code></a></p>")
        lines.append(f"<p><a href='{rp}'><img src='{rp}' style='max-width:900px; height:auto;'/></a></p>")

    if skymaps:
        lines.append("<h3>Sky localization</h3>")
        for p in skymaps:
            rp = html.escape(rel(p))
            lines.append(f"<p><a href='{rp}'><code>{rp}</code></a></p>")
            lines.append(f"<p><a href='{rp}'><img src='{rp}' style='max-width:900px; height:auto;'/></a></p>")

    lines.append("</body></html>\n")
    out_path.write_text("\n".join(lines), encoding="utf-8")


# ---------------------------------------------------------------------
# Galaxy PE discovery
# ---------------------------------------------------------------------


def _guess_galaxy_pe_dirs(catalog: str | None) -> list[Path]:
    """
    Galaxy stages input collections under ./galaxy_inputs/<collection_name>/...
    Your convention:
      - GWTC-2.1-PE
      - GWTC-3-PE
      - GWTC-4-PE
    """
    base = Path("galaxy_inputs")
    if not base.exists():
        return []
    if catalog:
        cand = base / f"{catalog}-PE"
        return [cand] if cand.exists() else []
    # any catalog: all *-PE dirs
    return sorted([p for p in base.iterdir() if p.is_dir() and p.name.endswith("-PE")])


def _find_local_pe_file(src_name: str, search_dirs: list[Path]) -> Path | None:
    """
    Find the first PEDataRelease file matching src_name under search_dirs.
    """
    pats = [
        f"*{src_name}*PEDataRelease*.h5",
        f"*{src_name}*PEDataRelease*.hdf5",
        f"*{src_name}*.h5",
        f"*{src_name}*.hdf5",
    ]
    for d in search_dirs:
        if not d.exists():
            continue
        for pat in pats:
            hits = sorted(d.rglob(pat))
            if hits:
                return hits[0]
    return None


# ---------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------


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
    data_repo: str = "s3",
    catalog: str | None = None,
    pe_vars: list[str] | None = None,
    pe_pairs: list[str] | None = None,
) -> dict[str, Any]:
    """
    Run the PE plotting pipeline and optionally write an HTML report.

    Returns a dict of produced filenames (always relative/absolute paths as written).
    """
    global LAST_RESULT
    
    import os

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
        plot_posterior_pairs,
        select_label,
        label_report,
        pe_log,
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

    posterior_keys_by_label: dict[str, list[str]] = {}
    missing_by_label: dict[str, list[str]] = {}
    posterior_pairs_missing_by_label: dict[str, list[str]] = {}

    go_next_cell = True

    def _progress(stage: str, progress: int, substage: str) -> None:
        if pr is not None:
            try:
                pr.report_progress(stage=stage, progress=progress, substage=substage)
            except Exception:
                pass

    # ---------------------------------------------------------------------
    # Resolve PE file (s3 / zenodo / galaxy / local)
    # ---------------------------------------------------------------------
    _progress("Download data", 10, "step 1")

    data = None
    label_waveform: Optional[str] = None
    local_pe_path: Optional[str] = None

    try:
        if data_repo == "zenodo":
            from difflib import get_close_matches

            index = build_zenodo_pe_index(cache_dir=".cache_gwosc", force_refresh=False)
            if src_name not in index:
                pe_log(f"⚠️ [WARN] Event {src_name} not found in cached Zenodo index. Refreshing index…", event_logs)
                index = build_zenodo_pe_index(cache_dir=".cache_gwosc", force_refresh=True)

            cands = index.get(src_name, [])
            chosen = choose_best_pe_file(cands)
            if chosen is None:
                keys = list(index.keys())
                sugg = get_close_matches(src_name, keys, n=5, cutoff=0.6)
                msg = (
                    f"No Zenodo PE file found for event {src_name}. "
                    f"Closest matches: {', '.join(sugg) if sugg else '(none)'}"
                )
                pe_log(f"❌ [ERROR] {msg}", event_logs)
                raise ValueError(msg)

            local_path = download_zenodo_pe_file(
                chosen, outdir=outdir, progress_cb=_progress if pr is not None else None, log_cb=lambda m: pe_log(m, event_logs),
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
                    pe_log(f"⚠️ [WARN] Could not parse S3_CREDENTIALS JSON: {e}", event_logs)
            else:
                credentials = {"endpoint": "minio-dev.odahub.fr", "secure": True}

            client = Minio(
                endpoint=credentials["endpoint"],
                secure=credentials.get("secure", True),
                access_key=credentials.get("access_key"),
                secret_key=credentials.get("secret_key"),
            )

            remote_file: Optional[str] = None
            for obj in client.list_objects("gwtc", recursive=True):
                file_name = obj.object_name
                if (
                    src_name in file_name
                    and "PEDataRelease" in file_name
                    and (file_name.endswith(".h5") or file_name.endswith(".hdf5"))
                ):
                    remote_file = file_name
                    pe_log(f"ℹ️ [INFO] Found remote PE file: {file_name}", event_logs)
                    break

            if remote_file is None:
                msg = f"❌ [ERROR] No PE file for event {src_name} found in S3 bucket 'gwtc'"
                pe_log(msg, event_logs)
                go_next_cell = False
            else:
                local_pe_path = remote_file  # keep same relative structure
                local_path = Path(local_pe_path)
                local_path.parent.mkdir(parents=True, exist_ok=True)

                if local_path.is_file() and local_path.stat().st_size > 0:
                    pe_log(f"ℹ️ [CACHE] Using existing local PE file: {local_path}", event_logs)
                else:
                    pe_log(f"ℹ️ [DOWNLOAD] Fetching from S3: {remote_file}", event_logs)
                    _download_s3_with_progress(
                        client,
                        bucket="gwtc",
                        object_name=remote_file,
                        out_path=local_path,
                        desc=f"S3: {Path(remote_file).name}",
                    )

        elif data_repo == "galaxy":
            pe_dirs = _guess_galaxy_pe_dirs(catalog)
            if not pe_dirs:
                pe_log("❌ [ERROR] Galaxy mode: ./galaxy_inputs/<CATALOG>-PE not found.", event_logs)
                go_next_cell = False
            else:
                found = _find_local_pe_file(src_name, pe_dirs)
                if found is None:
                    pe_log(
                        "❌ [ERROR] Galaxy mode: no PE file found for "
                        f"{src_name} under: {', '.join(str(p) for p in pe_dirs)}", event_logs
                    )
                    go_next_cell = False
                else:
                    local_pe_path = str(found)
                    pe_log(f"ℹ️ [INFO] Galaxy mode: using PE file {local_pe_path}", event_logs)

        else:
            # local/dev fallback
            candidates = sorted(outdir.glob(f"*{src_name}*PEDataRelease*.h5")) + sorted(
                outdir.glob(f"*{src_name}*PEDataRelease*.hdf5")
            )
            if not candidates:
                pe_log(
                    "❌ [ERROR] local mode: no PE file found in output directory. "
                    "Provide a PEDataRelease .h5/.hdf5 file or use --data-repo s3/zenodo/galaxy.", event_logs
                )
                go_next_cell = False
            else:
                local_pe_path = str(candidates[0])
                pe_log(f"ℹ️ [INFO] local mode: using PE file {local_pe_path}", event_logs)

        if go_next_cell and local_pe_path is not None:
            pe_log(f"ℹ️ [INFO] Reading PE data from: {local_pe_path}", event_logs)
            _progress("Read data", 25, "step 2")
            try:
                data = read(local_pe_path)
                label_report(data, sample_method=sample_method, strain_approximant=strain_approximant)
            except Exception as e:
                pe_log(f"❌ [ERROR] Failed reading PE data: {e}", event_logs)
                go_next_cell = False

    except Exception as e:
        pe_log(f"❌ [ERROR] Failed resolving PE input ({data_repo}): {e}", event_logs)
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
                require_psd=False,
                show_labels=True,
                context="samples",
                event_logs=event_logs,
            )

            if label is None:
                pe_log(f"❌ [ERROR] PE samples: no suitable label found; labels present: {labels}", event_logs)
                go_next_cell = False
            else:
                label_waveform = label
                pe_log(f"ℹ️ [INFO] PE samples: final label = {label}", event_logs)

                posterior_samples = samples_dict[label]
                approximant_for_title = label.split(":", 1)[1] if ":" in label else label

                posterior_info = plot_basic_posteriors(
                    posterior_samples,
                    src_name,
                    outdir,
                    approximant=approximant_for_title,
                    extra_params=pe_vars,
                )
                posterior_files: dict[str, str] = posterior_info.get("plots", {})
                posterior_keys_by_label[label] = posterior_info.get("available_keys", [])
                missing_by_label[label] = posterior_info.get("missing_requested", [])

                # Optional 2D pairs
                pair_info = plot_posterior_pairs(
                    posterior_samples,
                    src_name,
                    outdir,
                    pairs=pe_pairs,
                    approximant=approximant_for_title,
                )
                pair_files: dict[str, str] = pair_info.get("plots", {})
                posterior_pairs_missing_by_label[label] = pair_info.get("missing_pairs", [])

                for _tok, fname in {**pair_files, **posterior_files}.items():
                    result.files_distribution.append(fname)
                    if oda_available and PictureProduct is not None:
                        fig_distributionList.append(PictureProduct.from_file(fname))

        except Exception as e:
            pe_log(f"❌ [ERROR] Failed creating posterior plots: {e}", event_logs)
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
                require_psd=True,
                show_labels=False,
                context="strain",
                event_logs=event_logs,
            )

            if label_waveform is None:
                pe_log(
                    "❌ [ERROR] Waveform model: no suitable label found in the PE file; "
                    "unable to load strain / build waveform.", event_logs
                )
                go_next_cell = False
            else:
                pe_log(f"ℹ️ [INFO] Waveform model: using label_waveform = {label_waveform}", event_logs)

            if go_next_cell:
                posterior_samples_wave = samples_dict[label_waveform]

                try:
                    dets_available = list(data.psd[label_waveform].keys())
                except Exception:
                    dets_available = ["H1", "L1", "V1"]

                pe_log(f"ℹ️ [INFO] Detectors in PE file (label {label_waveform}): {dets_available}", event_logs)

                for det in dets_available:
                    det_time_key = f"{det}_time"
                    t0_det: Optional[float] = None
                    try:
                        t0_det = float(posterior_samples_wave.maxL[det_time_key][0])
                        pe_log(f"ℹ️ [INFO] Using maxL {det_time_key} = {t0_det} for {det}", event_logs)
                    except Exception:
                        try:
                            t0_det = float(posterior_samples_wave.maxL["geocent_time"][0])
                            pe_log(f"ℹ️ [INFO] Using geocent_time = {t0_det} for {det}", event_logs)
                        except Exception:
                            pe_log(f"⚠️ [WARN] No time information for {det} in maxL table", event_logs)
                            continue

                    try:
                        strain_data[det] = load_strain(src_name, t0_det, det)
                        merger_times[det] = t0_det
                    except Exception as e:
                        pe_log(f"⚠️ [WARN] Could not load strain for {det}: {e}", event_logs)

                if not strain_data:
                    pe_log("⚠️ [WARN] No strain data could be loaded for any detector.", event_logs)
                    go_next_cell = False

        except Exception as e:
            pe_log(f"❌ [ERROR] Failed loading strain data: {e}", event_logs)
            go_next_cell = False

    # Generate projected waveforms / q-transforms
    if go_next_cell and data is not None and label_waveform is not None:
        _progress("Generate projected waveforms", 90, "step 5")
        try:
            start_before = float(start)
            stop_after = float(stop)

            def _sanitize_for_engine(aprx: str | None) -> str | None:
                """Convert PE-style labels (e.g. IMRPhenomXPHM-SpinTaylor) into LAL-friendly names."""
                if not aprx:
                    return None
                a = str(aprx).strip().rstrip(",;")
                if ":" in a and a.lstrip().startswith("C"):
                    # If user passed a full PE label like 'C00:SEOBNRv5PHM', take RHS
                    a = a.split(":", 1)[1].strip()
                # Drop PE-specific suffixes that LAL/GWSignal doesn't understand
                if "-" in a:
                    a = a.split("-", 1)[0].strip()
                return a or None


            for det, strain in strain_data.items():
                t0 = merger_times[det]

                # Keep the PE label for reading PSD/samples/skymap
                pe_label = label_waveform

                # Derive the "RHS" of the PE label (e.g. "IMRPhenomXPHM-SpinTaylor")
                pe_rhs = pe_label.split(":", 1)[1] if ":" in pe_label else pe_label

                # Engine approximant to pass to maxL_td_waveform (LAL/GWSignal-friendly)
                requested_engine_aprx = (
                    _sanitize_for_engine(strain_approximant)
                    or _sanitize_for_engine(pe_rhs)
                )

                pe_log(
                    f"ℹ️ [INFO] Overlay config for {det}: "
                    f"PE label={pe_label}, PE rhs={pe_rhs}, engine requested={requested_engine_aprx}", event_logs
                )

                bp_cropped, crop_temp, used_aprx = generate_projected_waveform(
                    strain=strain,
                    event=src_name,
                    det=det,
                    t0=t0,
                    pedata=data,
                    label=pe_label,                          # PE label for PSD/maxL context
                    requested_approximant=requested_engine_aprx,  # engine approximant
                    allow_fallback=True,
                    freqrange=(fs_low, fs_high),
                    time_window=(start_before, stop_after),
                    event_logs=event_logs,
                )

                if bp_cropped is not None and crop_temp is not None:
                
                    if requested_engine_aprx == "SEOBNRv5PHM":
                        pe_log(
                            "ℹ️ [INFO] SEOBNRv5PHM exists in PE but may not be regenerable "
                            "for strain overlays with GWSignal; fallback may occur.", event_logs
                        )
                    fname_overlay = plot_whitened_overlay(
                        bp_cropped,
                        crop_temp,
                        src_name,
                        det,
                        outdir,
                        approximant=used_aprx,
                        t0=t0,
                        pe_label=pe_rhs,
                        engine_requested=requested_engine_aprx,   # what you tried first
                        engine_used=used_aprx,                    # what succeeded (may be fallback)
                    )
                    result.files_strain.append(fname_overlay)
                    if oda_available and PictureProduct is not None:
                        fig_strainList.append(PictureProduct.from_file(fname_overlay))
                else:
                    pe_log(f"⚠️ [WARN] Skipping overlay for {det} (no waveform).", event_logs)

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
                    if oda_available and PictureProduct is not None:
                        fig_strainList.append(PictureProduct.from_file(fname_q))
                except Exception as e_q:
                    pe_log(f"⚠️ [WARN] Could not create q-transform for {det}: {e_q}", event_logs)


        except Exception as e:
            pe_log(f"❌ [ERROR] Failed creating projected waveforms / q-transforms: {e}", event_logs)
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
            ax.set_title(f"PSD model (waveform: {label_waveform})", fontsize=10)
            plt.tight_layout()

            waveform = label_waveform.split(":", 1)[1] if ":" in label_waveform else label_waveform
            fname_psd = os.path.join(outdir, f"{src_name}_{waveform}_psd.png")
            fig.savefig(fname_psd, dpi=150)
            plt.close(fig)

            result.files_psd.append(fname_psd)
            if oda_available and PictureProduct is not None:
                fig_psdList.append(PictureProduct.from_file(fname_psd))

        except Exception as e:
            pe_log(f"❌ [ERROR] Failed creating PSD plot: {e}", event_logs)
            go_next_cell = False

    # ---------------------------------------------------------------------
    # Skymap plot
    # ---------------------------------------------------------------------
    if go_next_cell and data is not None and label_waveform is not None:
        try:
            fig, ax = data.skymap[label_waveform].plot(contour=[50, 90])
            
            # shrink the "50% area / 90% area" annotation
            for t in ax.texts + fig.texts:
                s = t.get_text()
                if "area" in s and "deg" in s:
                    t.set_fontsize(6)
            ax.set_title(f"Skymap (waveform: {label_waveform})", fontsize=8)
            ax.tick_params(axis="both", labelsize=7)
            ax.set_xlabel(ax.get_xlabel(), fontsize=7)
            ax.set_ylabel(ax.get_ylabel(), fontsize=7)

            try:
                from matplotlib.colorbar import Colorbar  # type: ignore
                for obj in fig.get_children():
                    if isinstance(obj, Colorbar):
                        obj.ax.tick_params(labelsize=6)
                        if obj.ax.yaxis.label is not None:
                            obj.ax.yaxis.label.set_fontsize(7)
            except Exception:
                pass

            waveform = label_waveform.split(":", 1)[1] if ":" in label_waveform else label_waveform
            fname_sky = os.path.join(outdir, f"{src_name}_{waveform}_skymap.png")
            fig.savefig(fname_sky, dpi=150)
            plt.close(fig)

            result.files_skymap.append(fname_sky)
            if oda_available and PictureProduct is not None:
                fig_skymapList.append(PictureProduct.from_file(fname_sky))

        except Exception as e:
            pe_log(f"❌ [ERROR] Failed creating skymap plot: {e}", event_logs)
            go_next_cell = False

    _progress("Finish", 100, "step 6")

    # ---------------------------------------------------------------------
    # HTML report
    # ---------------------------------------------------------------------
    if out_report_html:
        out_path = Path(out_report_html)
        out_path.parent.mkdir(parents=True, exist_ok=True)

        plot_files = sorted(outdir.glob("*.png"), key=_pe_plot_sort_key)
        other_files = sorted([p for p in outdir.iterdir() if p.is_file() and p.suffix.lower() != ".png"])

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
                catalog=catalog,
            ),
            posterior_keys_by_label=posterior_keys_by_label,
            missing_pe_vars_by_label=missing_by_label,
            requested_pe_vars=pe_vars or [],
            posterior_pairs_missing_by_label=posterior_pairs_missing_by_label,
            requested_pe_pairs=pe_pairs or [],
            event_logs=event_logs,
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
        last_msg = event_logs[-1].strip()
        raise ValueError(last_msg or "Parameter estimation failed before producing outputs.")

    return {
        "fig_distribution": result.files_distribution,
        "fig_strain": result.files_strain,
        "fig_psd": result.files_psd,
        "fig_skymap": result.files_skymap,
        "tool_log": event_logs,
    }
