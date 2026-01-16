from __future__ import annotations

import csv
import json
import math
import os
import re
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple


# ----------------------------
# Helpers: PE file resolution
# ----------------------------

def _is_hdf5_file(p: Path) -> bool:
    return p.is_file() and p.suffix.lower() in {".h5", ".hdf5", ".hdf"}


def _normalize_token(s: str) -> str:
    return re.sub(r"[^A-Za-z0-9_]+", "_", s).strip("_").lower()


def _find_pe_files_in_dir(root: Path) -> List[Path]:
    if not root.exists():
        return []
    return [p for p in root.rglob("*") if _is_hdf5_file(p)]


def _pick_best_match(event: str, candidates: List[Path]) -> Optional[Path]:
    if not candidates:
        return None

    ev = _normalize_token(event)

    def score(p: Path) -> Tuple[int, int, int]:
        name = _normalize_token(p.name)
        exact = 1 if ev in name.split("_") else 0
        substr = 1 if ev in name else 0
        return (exact, substr, -len(p.name))

    ranked = sorted(candidates, key=score, reverse=True)
    best = ranked[0]
    best_score = score(best)
    tied = [p for p in ranked if score(p) == best_score]
    uniq = {str(p.resolve()) for p in tied}
    if len(uniq) > 1:
        return None
    return best


def resolve_pe_input_from_collections(
    *,
    event: str,
    catalogs: List[str],
    collections: Dict[str, str],
    logs: List[str],
) -> Path:
    # Map catalogs -> preferred collections (tune if your catalog labels differ)
    search_order: List[str] = []
    for c in catalogs:
        cu = c.upper()
        if "2.1" in cu or "GWTC-2.1" in cu:
            search_order.append("GWTC-2.1-PE")
        elif "GWTC-3" in cu:
            search_order.append("GWTC-3-PE")
        elif "GWTC-4" in cu:
            search_order.append("GWTC-4-PE")

    if not search_order:
        search_order = ["GWTC-4-PE", "GWTC-3-PE", "GWTC-2.1-PE"]

    logs.append(f"PE search order: {search_order}")

    all_hits: List[Path] = []
    details: Dict[str, List[Path]] = {}

    evtok = _normalize_token(event)
    for col in search_order:
        col_dir = collections.get(col)
        if not col_dir:
            logs.append(f"Collection not provided: {col}")
            continue
        root = Path(col_dir).expanduser().resolve()
        files = _find_pe_files_in_dir(root)
        hits = [p for p in files if evtok in _normalize_token(p.name)]
        details[col] = hits
        logs.append(f"{col}: {len(files)} hdf5 files, {len(hits)} event-token matches")
        all_hits.extend(hits)

    best = _pick_best_match(event, all_hits)
    if best is None:
        any_hits = any(details[c] for c in details)
        if not any_hits:
            raise FileNotFoundError(
                f"No PE .h5/.hdf5 found for event='{event}' in provided collections."
            )
        # ambiguous
        msg = [f"Ambiguous PE match for event='{event}'. Candidates:"]
        for col, hits in details.items():
            if hits:
                msg.append(f"  {col}: " + ", ".join(p.name for p in hits[:10]) + (" ..." if len(hits) > 10 else ""))
        raise RuntimeError("\n".join(msg))

    logs.append(f"Selected PE file: {best}")
    return best


# ----------------------------
# Optional imports (best-effort)
# ----------------------------

def _safe_import_pesummary_read():
    try:
        from pesummary.io import read  # type: ignore
        return read
    except Exception as e:
        raise RuntimeError(
            "pesummary is required for parameters_estimation, but could not be imported."
        ) from e


def _extract_t0_by_detector_from_pesummary(data: Any, logs: List[str], label: Optional[str] = None) -> Dict[str, float]:
    """Best-effort extraction of detector-specific merger times from PE max-likelihood samples.

    The Renku notebook uses e.g. posterior_samples_wave.maxL['H1_time'][0].
    pesummary objects may expose ``maxL`` either as a dict of arrays, or nested by label.

    Returns a dict like: {"H1": t0, "L1": t0, "V1": t0}. Missing detectors are omitted.
    """
    out: Dict[str, float] = {}
    try:
        maxl = getattr(data, "maxL", None)
        if maxl is None:
            logs.append("INFO: PE file has no maxL attribute -> cannot derive detector merger times.")
            return {}

        # If nested by label, select sub-dict
        if isinstance(maxl, dict) and label and label in maxl and isinstance(maxl[label], dict):
            maxl_dict = maxl[label]
        else:
            maxl_dict = maxl

        for det in ("H1", "L1", "V1"):
            key = f"{det}_time"
            try:
                if isinstance(maxl_dict, dict) and key in maxl_dict:
                    t0 = float(maxl_dict[key][0])
                    out[det] = t0
                    logs.append(f"ℹ️ [INFO] Using maxL {key} = {t0} for {det}")
            except Exception:
                # ignore per-detector failures
                pass

        if not out:
            logs.append("INFO: No detector *_time keys found in PE maxL samples.")
        return out
    except Exception as e:
        logs.append(f"WARNING: Failed to extract detector merger times from PE: {e}")
        return {}


def _maybe_import_numpy():
    try:
        import numpy as np  # type: ignore
        return np
    except Exception:
        return None


def _maybe_import_matplotlib():
    try:
        import matplotlib.pyplot as plt  # type: ignore
        return plt
    except Exception:
        return None


def _maybe_import_gwpy():
    """
    Used for strain/PSD quicklooks via GWOSC open data.
    """
    try:
        from gwpy.timeseries import TimeSeries  # type: ignore
        return TimeSeries
    except Exception:
        return None


# ----------------------------
# Reporting
# ----------------------------

def _html_escape(s: str) -> str:
    return (
        s.replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
        .replace("'", "&#039;")
    )


def _write_html_report(
    out_html: Path,
    *,
    title: str,
    summary_kv: List[Tuple[str, str]],
    plots: List[Path],
    logs: List[str],
) -> None:
    out_html.parent.mkdir(parents=True, exist_ok=True)

    rows = "\n".join(
        f"<tr><td style='padding:6px 10px; border:1px solid #ddd;'><b>{_html_escape(k)}</b></td>"
        f"<td style='padding:6px 10px; border:1px solid #ddd;'>{_html_escape(v)}</td></tr>"
        for k, v in summary_kv
    )

    plot_items = "\n".join(
        f"<div style='margin:16px 0;'>"
        f"<div><b>{_html_escape(p.name)}</b></div>"
        f"<img src='{_html_escape(p.name)}' style='max-width:100%; border:1px solid #ddd; border-radius:6px;'/>"
        f"</div>"
        for p in plots
    )

    logs_block = "\n".join(
        f"<li><pre style='white-space:pre-wrap; margin:0;'>{_html_escape(x)}</pre></li>"
        for x in logs
    )

    html = f"""<!doctype html>
<html>
<head>
  <meta charset="utf-8"/>
  <title>{_html_escape(title)}</title>
</head>
<body style="font-family: system-ui, -apple-system, Segoe UI, Roboto, Helvetica, Arial, sans-serif; margin: 24px;">
  <h1 style="margin-top:0;">{_html_escape(title)}</h1>

  <h2>Summary</h2>
  <table style="border-collapse:collapse;">{rows}</table>

  <h2>Plots</h2>
  {plot_items if plot_items else "<p>No plots produced.</p>"}

  <h2>Log</h2>
  <ul>{logs_block}</ul>
</body>
</html>
"""
    out_html.write_text(html, encoding="utf-8")


# ----------------------------
# Plotting
# ----------------------------

def _ensure_plots_dir(plots_dir: Optional[str], out_report_html: str) -> Path:
    if plots_dir:
        p = Path(plots_dir)
    else:
        # like other modes: default next to report
        p = Path(out_report_html).resolve().parent / "plots"
    p.mkdir(parents=True, exist_ok=True)
    return p


def _plot_1d_histograms(
    samples: Dict[str, Any],
    out_dir: Path,
    *,
    label: str,
    params: Iterable[str],
    logs: List[str],
) -> List[Path]:
    plt = _maybe_import_matplotlib()
    np = _maybe_import_numpy()
    if plt is None:
        logs.append("WARNING: matplotlib not available -> skipping posterior histograms.")
        return []

    out: List[Path] = []
    for par in params:
        if par not in samples:
            continue
        try:
            arr = samples[par]
            if np is not None:
                arr = np.asarray(arr).reshape(-1)
            else:
                arr = list(arr)

            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.hist(arr, bins=60)
            ax.set_title(f"{label}: {par}")
            ax.set_xlabel(par)
            ax.set_ylabel("count")

            fp = out_dir / f"{label}_{par}_hist.png"
            fig.savefig(fp, dpi=150, bbox_inches="tight")
            plt.close(fig)

            out.append(fp)
            logs.append(f"Saved histogram: {fp.name}")
        except Exception as e:
            logs.append(f"WARNING: failed histogram {par}: {e}")
    return out


def _plot_skymap_from_radec(
    samples: Dict[str, Any],
    out_dir: Path,
    *,
    label: str,
    logs: List[str],
) -> List[Path]:
    """
    Simple RA/Dec density plot (mollweide) without healpy/ligo.skymap dependency.
    Requires samples['ra'] and samples['dec'] in radians (common in PE outputs).
    If in degrees, it still produces something reasonable but may look wrong.
    """
    plt = _maybe_import_matplotlib()
    np = _maybe_import_numpy()
    if plt is None or np is None:
        logs.append("WARNING: matplotlib/numpy not available -> skipping skymap plot.")
        return []

    if "ra" not in samples or "dec" not in samples:
        logs.append("INFO: ra/dec not found in samples -> skipping skymap plot.")
        return []

    try:
        ra = np.asarray(samples["ra"]).reshape(-1)
        dec = np.asarray(samples["dec"]).reshape(-1)

        # Mollweide expects lon in [-pi, +pi]. If ra is [0, 2pi], shift.
        ra_wrapped = ((ra + np.pi) % (2 * np.pi)) - np.pi
        # Flip RA for astronomical convention (optional but common)
        ra_plot = -ra_wrapped

        fig = plt.figure()
        ax = fig.add_subplot(111, projection="mollweide")
        ax.grid(True)
        ax.scatter(ra_plot, dec, s=1, alpha=0.15)
        ax.set_title(f"{label}: RA/Dec posterior (mollweide)")

        fp = out_dir / f"{label}_skymap_mollweide.png"
        fig.savefig(fp, dpi=150, bbox_inches="tight")
        plt.close(fig)

        logs.append(f"Saved skymap: {fp.name}")
        return [fp]
    except Exception as e:
        logs.append(f"WARNING: failed skymap plot: {e}")
        return []


def _load_event_time_from_events_json(event: str, events_json: str, logs: List[str]) -> Optional[float]:
    """
    Try to pull a GPS time from events_json. We keep it permissive:
    - if entry has keys like 'gps', 'gpstime', 'geocent_time', etc.
    """
    try:
        d = json.loads(Path(events_json).read_text(encoding="utf-8"))
    except Exception as e:
        logs.append(f"WARNING: cannot read events_json: {e}")
        return None

    # common structures: dict[event]=..., or list of dicts with name field
    entry = None
    if isinstance(d, dict):
        entry = d.get(event) or d.get(event.upper()) or d.get(event.lower())
    if entry is None and isinstance(d, list):
        for x in d:
            if isinstance(x, dict) and (x.get("event") == event or x.get("name") == event):
                entry = x
                break

    if not isinstance(entry, dict):
        logs.append("INFO: events_json did not contain a usable entry -> skipping strain/PSD.")
        return None

    for k in ["gps", "gpstime", "geocent_time", "geocent_time_gps", "time"]:
        if k in entry:
            try:
                return float(entry[k])
            except Exception:
                pass

    logs.append("INFO: could not find GPS time field in events_json -> skipping strain/PSD.")
    return None
    
def _plot_strain_and_psd(
    *,
    gps: Optional[float] = None,
    t0_by_det: Optional[Dict[str, float]] = None,
    out_dir: Path,
    label: str,
    logs: List[str],
    start: float,
    stop: float,
    fs_low: float,
    fs_high: float,

    # PSD controls
    psd_mode: str,  # "estimate" | "file"
    psd_duration: float,
    psd_offset: float,
    psd_files: Dict[str, Optional[str]],

    detectors: Optional[List[str]] = None,
    make_strain: bool = True,
    make_psd: bool = True,
    # optional future:
    make_overlay: bool = False,
    strain_approximant: str = "IMRPhenomXPHM",
) -> List[Path]:
    """
    Fetch open data from GWOSC via gwpy (best-effort). If network access isn't available,
    this will fail gracefully (log + return []).

    PSD modes:
      - estimate: estimate PSD from an off-source segment using psd_offset and psd_duration
      - file: require a detector PSD file path (psd_files[det]) and load it
    """
    TimeSeries = _maybe_import_gwpy()
    plt = _maybe_import_matplotlib()
    np = _maybe_import_numpy()
    if TimeSeries is None or plt is None:
        logs.append("WARNING: gwpy/matplotlib not available -> skipping strain/PSD plots.")
        return []

    psd_mode_norm = (psd_mode or "").strip().lower()
    if psd_mode_norm not in {"estimate", "file"}:
        raise ValueError(f"psd_mode must be 'estimate' or 'file' (got: {psd_mode})")

    dets = detectors or ["H1", "L1"]  # default guess
    out: List[Path] = []
    # We can use either:
    #   - detector-specific merger times (t0_by_det) derived from PE maxL samples, or
    #   - a single gps time (gps) as a fallback.
    # Windows are computed per-detector in the loop below.

    # ---- local helper: load PSD from file (simple 2-column text) ----
    def _load_psd_file_text(psd_path: Path) -> Tuple[Any, Any]:
        """
        Load PSD from a text-like file with at least 2 columns: frequency[Hz]  psd[strain^2/Hz]
        Returns (freqs, psd_values) as numpy arrays if numpy is available, else Python lists.
        """
        freqs: List[float] = []
        vals: List[float] = []
        with psd_path.open("r", encoding="utf-8") as f:
            for line in f:
                s = line.strip()
                if not s or s.startswith("#"):
                    continue
                parts = s.split()
                if len(parts) < 2:
                    continue
                try:
                    freqs.append(float(parts[0]))
                    vals.append(float(parts[1]))
                except Exception:
                    continue
        if not freqs:
            raise ValueError(f"No PSD data parsed from file: {psd_path}")
        if np is not None:
            return np.asarray(freqs), np.asarray(vals)
        return freqs, vals

    for det in dets:
        # Choose reference time for this detector
        if t0_by_det and det in t0_by_det:
            t0 = float(t0_by_det[det])
        elif gps is not None:
            t0 = float(gps)
        else:
            logs.append(f"WARNING: No t0 available for {det} (need PE maxL {det}_time or gps fallback) -> skipping")
            continue

        # Event window [t0 - start, t0 + stop]
        t1 = t0 - float(start)
        t2 = t0 + float(stop)

        # PSD off-source window, by default BEFORE the event:
        # [t0 - psd_offset - psd_duration, t0 - psd_offset]
        psd_t2 = t0 - float(psd_offset)
        psd_t1 = psd_t2 - float(psd_duration)
        try:
            # 1) Fetch event-window strain
            logs.append(f"Fetching open data (event window): {det} [{t1}, {t2}]")
            ts = TimeSeries.fetch_open_data(det, t1, t2)

            # Bandpass for plotting/PSD convenience
            try:
                ts_bp = ts.bandpass(fs_low, fs_high)
            except Exception:
                ts_bp = ts

            # 2) Build/obtain PSD depending on mode (only if needed)
            psd_obj = None          # gwpy FrequencySeries if estimate
            psd_freqs = None        # arrays/lists if file
            psd_vals = None         # arrays/lists if file

            if make_psd or make_overlay:
                if psd_mode_norm == "estimate":
                    try:
                        logs.append(f"Fetching open data (PSD window): {det} [{psd_t1}, {psd_t2}]")
                        ts_for_psd = TimeSeries.fetch_open_data(det, psd_t1, psd_t2)
                        try:
                            ts_for_psd = ts_for_psd.bandpass(fs_low, fs_high)
                        except Exception:
                            pass

                        # 4-second FFT segments is a decent default; you can expose later if needed
                        psd_obj = ts_for_psd.psd(4)
                        logs.append(f"Estimated PSD for {det} using {psd_duration}s off-source window (offset={psd_offset}s).")
                    except Exception as e:
                        logs.append(f"WARNING: PSD estimate failed for {det}: {e}")
                        psd_obj = None

                elif psd_mode_norm == "file":
                    psd_path_str = psd_files.get(det)
                    if not psd_path_str:
                        logs.append(f"WARNING: psd_mode=file but no PSD file provided for {det} (psd_files['{det}'] missing).")
                    else:
                        psd_path = Path(psd_path_str).expanduser().resolve()
                        if not psd_path.exists():
                            logs.append(f"WARNING: PSD file does not exist for {det}: {psd_path}")
                        else:
                            try:
                                psd_freqs, psd_vals = _load_psd_file_text(psd_path)
                                logs.append(f"Loaded PSD from file for {det}: {psd_path.name}")
                            except Exception as e:
                                logs.append(f"WARNING: Failed to load PSD file for {det}: {e}")

            # 3) Strain plot
            if make_strain:
                fig = plt.figure()
                ax = fig.add_subplot(111)
                ax.plot(ts_bp.times.value - t0, ts_bp.value)
                ax.set_title(f"{label}: {det} strain (bandpass {fs_low}-{fs_high} Hz)")
                ax.set_xlabel("t - t0 (s)")
                ax.set_ylabel("strain")
                fp = out_dir / f"{label}_{det}_strain.png"
                fig.savefig(fp, dpi=150, bbox_inches="tight")
                plt.close(fig)
                out.append(fp)
                logs.append(f"Saved strain plot: {fp.name}")

            # 4) PSD plot
            if make_psd:
                try:
                    fig = plt.figure()
                    ax = fig.add_subplot(111)

                    if psd_mode_norm == "estimate":
                        if psd_obj is None:
                            logs.append(f"WARNING: No PSD available for {det} (estimate mode failed) -> skipping PSD plot.")
                        else:
                            ax.loglog(psd_obj.frequencies.value, psd_obj.value)
                            ax.set_title(f"{label}: {det} PSD (estimated off-source)")
                            ax.set_xlabel("frequency (Hz)")
                            ax.set_ylabel("strain^2 / Hz")
                            fp = out_dir / f"{label}_{det}_psd.png"
                            fig.savefig(fp, dpi=150, bbox_inches="tight")
                            out.append(fp)
                            logs.append(f"Saved PSD plot: {fp.name}")

                    elif psd_mode_norm == "file":
                        if psd_freqs is None or psd_vals is None:
                            logs.append(f"WARNING: No PSD available for {det} (file mode missing/failed) -> skipping PSD plot.")
                        else:
                            ax.loglog(psd_freqs, psd_vals)
                            ax.set_title(f"{label}: {det} PSD (from file)")
                            ax.set_xlabel("frequency (Hz)")
                            ax.set_ylabel("strain^2 / Hz")
                            fp = out_dir / f"{label}_{det}_psd.png"
                            fig.savefig(fp, dpi=150, bbox_inches="tight")
                            out.append(fp)
                            logs.append(f"Saved PSD plot: {fp.name}")

                    plt.close(fig)
                except Exception as e:
                    logs.append(f"WARNING: PSD plot failed for {det}: {e}")

            # 5) (Future) overlay hook — keep inputs propagated, but don’t break current flow
            if make_overlay:
                logs.append(
                    f"INFO: overlay requested for {det} with strain_approximant={strain_approximant}, "
                    "but overlay generation is not implemented yet."
                )

        except Exception as e:
            logs.append(f"WARNING: strain fetch/plot failed for {det}: {e}")

    return out


# ----------------------------
# TSV outputs
# ----------------------------

def _quantiles(arr: Any, qs: Tuple[float, float, float]) -> Optional[Tuple[float, float, float]]:
    np = _maybe_import_numpy()
    try:
        if np is None:
            return None
        a = np.asarray(arr).reshape(-1)
        return tuple(float(x) for x in np.quantile(a, qs))
    except Exception:
        return None

def _write_out_events_tsv(
        out_path: str,
        *,
        event: str,
        catalogs: List[str],
        pe_file: str,
        label: str,
        sample_method: str,
        strain_approximant: str,
        psd_mode: str,
        psd_duration: float,
        psd_offset: float,
        psd_h1: Optional[str],
        psd_l1: Optional[str],
        psd_v1: Optional[str],
        samples: Dict[str, Any],
        logs: List[str],
) -> None:

    # Keep consistent with other modes: one TSV, one row
    # Include a few common parameters if present
   
    fields = [
        "event",
        "catalogs",
        "pe_file",
        "label",
        "sample_method",
        "strain_approximant",
        "psd_mode",
        "psd_duration",
        "psd_offset",
        "psd_h1",
        "psd_l1",
        "psd_v1",
        "n_samples",
        "mchirp_p05", "mchirp_p50", "mchirp_p95",
        "dl_p05", "dl_p50", "dl_p95",
        "chi_eff_p05", "chi_eff_p50", "chi_eff_p95",
    ]


    n_samples = None
    try:
        # find any param to estimate length
        for k, v in samples.items():
            q = _quantiles(v, (0.05, 0.5, 0.95))
            if q is not None:
                # if quantiles worked, we can estimate length via numpy
                np = _maybe_import_numpy()
                if np is not None:
                    n_samples = int(np.asarray(v).reshape(-1).shape[0])
                break
    except Exception:
        pass

    def q_or_blank(param: str) -> Tuple[str, str, str]:
        if param not in samples:
            return ("", "", "")
        q = _quantiles(samples[param], (0.05, 0.5, 0.95))
        if q is None:
            return ("", "", "")
        return (f"{q[0]:.6g}", f"{q[1]:.6g}", f"{q[2]:.6g}")

    mchirp = q_or_blank("chirp_mass_source")
    dl = q_or_blank("luminosity_distance")
    chi = q_or_blank("chi_eff")

    row = {
        "event": event,
        "catalogs": ",".join(catalogs),
        "pe_file": pe_file,
        "label": label,

        "sample_method": sample_method,
        "strain_approximant": strain_approximant,
        "psd_mode": psd_mode,
        "psd_duration": str(psd_duration),
        "psd_offset": str(psd_offset),
        "psd_h1": psd_h1 or "",
        "psd_l1": psd_l1 or "",
        "psd_v1": psd_v1 or "",
        "n_samples": "" if n_samples is None else str(n_samples),
        "mchirp_p05": mchirp[0], "mchirp_p50": mchirp[1], "mchirp_p95": mchirp[2],
        "dl_p05": dl[0], "dl_p50": dl[1], "dl_p95": dl[2],
        "chi_eff_p05": chi[0], "chi_eff_p50": chi[1], "chi_eff_p95": chi[2],
    }

    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerow(row)

    logs.append(f"Wrote out-events TSV: {out_path}")


# ----------------------------
# Main entry point (Galaxy pattern)
# ----------------------------

@dataclass
class PEConfig:
    catalogs: List[str]
    event: str
    out_events_tsv: str
    out_report_html: str
    events_json: Optional[str]
    # PE collections
    pe_collection_gwtc21: Optional[str]
    pe_collection_gwtc3: Optional[str]
    pe_collection_gwtc4: Optional[str]
    # plot controls
    make_plots: str
    plots_dir: Optional[str]
    # strain window + bandpass
    start: float
    stop: float
    fs_low: float
    fs_high: float

def run_parameters_estimation(
    *,
    catalogs: List[str],
    out_events_tsv: str,
    out_report_html: str,
    event: str,
    events_json: Optional[str] = None,
    make_plots: str = "posteriors",
    plots_dir: Optional[str] = None,

    # Posterior selection + waveform label
    sample_method: str = "Mixed",
    strain_approximant: str = "IMRPhenomXPHM",

    # Strain window + bandpass
    start: float = 0.5,
    stop: float = 0.1,
    fs_low: float = 20.0,
    fs_high: float = 300.0,

    # PSD controls
    psd_mode: str = "estimate",      # estimate|file
    psd_duration: float = 32.0,
    psd_offset: float = 64.0,
    psd_h1: Optional[str] = None,
    psd_l1: Optional[str] = None,
    psd_v1: Optional[str] = None,

    # NEW: preferred Galaxy input (expanded from collections)
    pe_files: Optional[List[str]] = None,

    # OPTIONAL fallback (only if you still want directory scanning)
    pe_collection_gwtc21: Optional[str] = None,
    pe_collection_gwtc3: Optional[str] = None,
    pe_collection_gwtc4: Optional[str] = None,
) -> int:
    """
    Galaxy-friendly PE mode: writes out_events_tsv + out_report_html.

    make_plots choices:
      - none
      - posteriors
      - skymap
      - strain
      - psd
      - all
    """
    logs: List[str] = []
    logs.append(f"event={event}")
    logs.append(f"catalogs={catalogs}")
    logs.append(f"make_plots={make_plots}")
    logs.append(f"sample_method={sample_method}")
    logs.append(f"strain_approximant={strain_approximant}")
    logs.append(f"psd_mode={psd_mode} psd_duration={psd_duration} psd_offset={psd_offset}")
    logs.append(f"psd_files: H1={psd_h1} L1={psd_l1} V1={psd_v1}")
    logs.append(f"pe_files count={0 if not pe_files else len(pe_files)}")

    # Output locations
    out_report_path = Path(out_report_html).resolve()
    out_dir = out_report_path.parent
    out_dir.mkdir(parents=True, exist_ok=True)

        # ------------------------------------------------------------
    # 1) Resolve PE file (prefer explicit pe_files from Galaxy)
    # ------------------------------------------------------------
    pe_file: Optional[Path] = None

    candidates: List[Path] = []
    if pe_files:
        for x in pe_files:
            if not x:
                continue
            p = Path(x).expanduser().resolve()
            if p.exists() and p.is_file():
                candidates.append(p)

        logs.append(f"pe_files provided={len(pe_files)} existing_files={len(candidates)}")

        # In Galaxy, staged datasets typically end in .dat even if the original was .h5.
        # So do NOT filter by suffix here.

        if len(candidates) == 1:
            pe_file = candidates[0]
            logs.append(f"Selected the only PE file provided: {pe_file}")
        elif len(candidates) > 1:
            # Try best-effort filename token match (may fail in Galaxy because filenames are dataset_*.dat)
            pe_file = _pick_best_match(event, candidates)
            if pe_file is None:
                sample_names = ", ".join(p.name for p in candidates[:10])
                raise FileNotFoundError(
                    f"Multiple PE files were provided but none matched event='{event}' by filename token. "
                    f"(Galaxy often stages files as dataset_*.dat.) "
                    f"Candidates (first 10): {sample_names}. "
                    "Tip: provide only one PE file for the event, or extend the tool to accept "
                    "identifier:path pairs."
                )
            logs.append(f"Selected PE file by best-match: {pe_file}")

    # Optional fallback: directory scanning (only if pe_file not resolved)
    if pe_file is None:
        collections = {
            "GWTC-2.1-PE": pe_collection_gwtc21 or "",
            "GWTC-3-PE": pe_collection_gwtc3 or "",
            "GWTC-4-PE": pe_collection_gwtc4 or "",
        }
        provided = {k: v for k, v in collections.items() if v}
        if not provided:
            raise SystemExit(
                "No PE files provided. Pass --pe-files (recommended, from Galaxy collections), "
                "or provide at least one of --pe-collection-gwtc21 / --pe-collection-gwtc3 / --pe-collection-gwtc4."
            )

        pe_file = resolve_pe_input_from_collections(
            event=event,
            catalogs=catalogs,
            collections=provided,
            logs=logs,
        )
        logs.append(f"Selected PE file from collection directory scan: {pe_file}")

    assert pe_file is not None


    # Copy PE file next to outputs for provenance
    pe_copy = out_dir / pe_file.name
    if pe_copy.resolve() != pe_file.resolve():
        shutil.copy2(pe_file, pe_copy)
    pe_file_used = pe_copy
    logs.append(f"PE file used (copied): {pe_file_used}")

    # ------------------------------------------------------------
    # 2) Read with pesummary
    # ------------------------------------------------------------
    read = _safe_import_pesummary_read()

    def _read_pe_file(path: str):
        # HDF5 magic number starts with 0x89 'H' 'D' 'F'
        with open(path, "rb") as f:
            magic = f.read(4)
        if magic == b"\x89HDF":
            return read(path, file_format="hdf5")
        return read(path)

    data = _read_pe_file(str(pe_file_used))
    labels = list(getattr(data, "labels", []) or [])
    label = labels[0] if labels else "analysis"

    # Prefer a label matching sample_method
    sm = (sample_method or "").strip().lower()
    if labels and sm:
        preferred = [lab for lab in labels if sm in str(lab).lower()]
        if preferred:
            label = preferred[0]

    logs\.append\(f"pesummary labels=\{labels\} selected_label=\{label\}"\)
    t0_by_det = _extract_t0_by_detector_from_pesummary(data, logs, label=label)

    # Extract samples dict robustly
    samples_dict: Dict[str, Any] = {}
    try:
        sd = getattr(data, "samples_dict", None) or {}
        if isinstance(sd, dict) and label in sd and isinstance(sd[label], dict):
            samples_dict = sd[label]
        elif isinstance(sd, dict):
            samples_dict = sd
    except Exception:
        samples_dict = {}

    logs.append(f"samples keys={list(samples_dict.keys())[:20]}")

    # ------------------------------------------------------------
    # 3) Write out-events TSV
    # ------------------------------------------------------------
    _write_out_events_tsv(
        out_events_tsv,
        event=event,
        catalogs=catalogs,
        pe_file=str(pe_file_used.name),
        label=label,
        sample_method=sample_method,
        strain_approximant=strain_approximant,
        psd_mode=psd_mode,
        psd_duration=psd_duration,
        psd_offset=psd_offset,
        psd_h1=psd_h1,
        psd_l1=psd_l1,
        psd_v1=psd_v1,
        samples=samples_dict,
        logs=logs,
    )

    # ------------------------------------------------------------
    # 4) Plots (explicit make_plots options)
    # ------------------------------------------------------------
    plots: List[Path] = []
    plots_root = _ensure_plots_dir(plots_dir, out_report_html)
    label_tag = f"{event}_{label}"

    def wants(kind: str) -> bool:
        return make_plots == "all" or make_plots == kind

    if make_plots != "none":
        if wants("posteriors"):
            default_params = [
                "mass_1_source",
                "mass_2_source",
                "chirp_mass_source",
                "luminosity_distance",
                "chi_eff",
                "ra",
                "dec",
                "theta_jn",
            ]
            plots.extend(
                _plot_1d_histograms(
                    samples_dict,
                    plots_root,
                    label=label_tag,
                    params=default_params,
                    logs=logs,
                )
            )

        if wants("skymap"):
            plots.extend(_plot_skymap_from_radec(samples_dict, plots_root, label=label_tag, logs=logs))

        if wants("strain") or wants("psd"):
            gps = None
            if events_json:
                gps = _load_event_time_from_events_json(event, events_json, logs)
            psd_files_map = {"H1": psd_h1, "L1": psd_l1, "V1": psd_v1}
            plots.extend(
                _plot_strain_and_psd(
                            gps=gps,
                            t0_by_det=t0_by_det,
                            out_dir=plots_root,
                            label=label_tag,
                            logs=logs,
                            start=start,
                            stop=stop,
                            fs_low=fs_low,
                            fs_high=fs_high,
                            psd_mode=psd_mode,
                            psd_duration=psd_duration,
                            psd_offset=psd_offset,
                            psd_files=psd_files_map,
                            detectors=None,
                            make_strain=wants("strain"),
                            make_psd=wants("psd"),
                            make_overlay=False,
                            strain_approximant=strain_approximant,
                        )
                    )

    # ------------------------------------------------------------
    # 5) Write HTML report
    # ------------------------------------------------------------
    report_dir = out_report_path.parent
    report_dir.mkdir(parents=True, exist_ok=True)
    report_plots: List[Path] = []
    for p in plots:
        dst = report_dir / p.name
        if dst.resolve() != p.resolve():
            shutil.copy2(p, dst)
        report_plots.append(dst)

    summary_kv = [
        ("Event", event),
        ("Catalogs", ",".join(catalogs)),
        ("PE file", pe_file_used.name),
        ("Label", label),
        ("sample_method", sample_method),
        ("strain_approximant", strain_approximant),
        ("psd_mode", psd_mode),
        ("psd_duration (s)", str(psd_duration)),
        ("psd_offset (s)", str(psd_offset)),
        ("psd_h1", psd_h1 or ""),
        ("psd_l1", psd_l1 or ""),
        ("psd_v1", psd_v1 or ""),
        ("Plots produced", str(len(report_plots))),
        ("Plots dir", str(plots_root)),
    ]

    _write_html_report(
        out_report_path,
        title=f"Parameter estimation: {event}",
        summary_kv=summary_kv,
        plots=report_plots,
        logs=logs,
    )

    return 0


