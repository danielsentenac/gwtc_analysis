"""
gwpe_utils.py

Utility functions for:
- caching downloaded files / strain
- loading GW open data strain via GWPy
- generating projected waveforms from PE posteriors (using pesummary)
- producing whitened overlay plots and q-transform time–frequency maps
"""

import os
import re
from typing import Tuple, Optional, List, Dict

import numpy as np
import matplotlib.pyplot as plt

from gwpy.timeseries import TimeSeries
from gwpy.frequencyseries import FrequencySeries
from scipy.interpolate import interp1d
from matplotlib.ticker import MaxNLocator
from pathlib import Path

__all__ = [
    "get_cache_dir",
    "cached_path",
    "get_cached_or_download",
    "file_already_exists",
    "local_pe_path",
    "load_strain",
    "ensure_outdir",
    "generate_projected_waveform",
    "plot_whitened_overlay",
    "plot_time_frequency",
    "compare_spectrogram_vs_qtransform",
    "plot_basic_posteriors",
    "local_pe_path",
    "label_has_psd",
    "select_label",
    "label_report",
]

# ---------------------------------------------------------------------
# Select label helpers
# ---------------------------------------------------------------------

def find_label_for_approximant(labels, aprx_str):
    """
    Return the label matching a given approximant name, e.g. 'SEOBNRv4PHM'
    → 'C01:SEOBNRv4PHM'. Returns None if not found.
    """
    pat = re.compile(rf"C\d{{2}}:{re.escape(aprx_str)}")
    for lab in labels:
        if pat.fullmatch(lab):
            return lab
    return None


def label_has_psd(pedata, label):
    """
    Return True if the PE dataset contains PSD info for the given label.
    """
    try:
        return (
            hasattr(pedata, "psd")
            and label in pedata.psd
            and pedata.psd[label] is not None
            and len(pedata.psd[label]) > 0
        )
    except Exception:
        return False

def select_label(
    pedata,
    requested_approximant,
    event_logs=None,
    require_psd: bool = True,
    show_labels: bool = True,
):
    """
    Select the best label for a given context (PE samples or strain), with:

      1. Requested approximant/method (if present),
      2. Other available approximants found in pedata.labels,
      3. Optionally a 'Mixed' label (useful for PE samples),
      4. First label as last resort,

    and *optionally* enforcing that the final choice has PSD available.

    Parameters
    ----------
    pedata : pesummary.gw.file.File
        The PE data returned by pesummary.read().
    requested_approximant : str or None
        Desired waveform approximant / method, e.g. "IMRPhenomXPHM", "Mixed".
    event_logs : list or None
        If provided, log strings are appended to event_logs[-1].
    require_psd : bool, optional
        If True, try very hard to return a label with PSD
        (use this for STRAIN / WAVEFORM generation).
        If False, return the best label by approximant/method logic only
        (use this for PE SAMPLES / POSTERIORS).
    show_labels : bool, optional
        If True, print a pretty, sorted list of available labels/waveforms.

    Returns
    -------
    str or None
        The chosen label, or None if no suitable label was found.
    """

    # Sorted labels for nicer logging
    try:
        labels = sorted(list(pedata.labels))
    except Exception:
        labels = []

    # Helper logger with emoji-ish tags
    def log(msg: str):
        print(msg)
        if event_logs is not None:
            event_logs[-1] += f"\n{msg}\n"

    # ---- Pretty, multi-line, sorted label listing ----
    if show_labels:
        if labels:
            bullet_lines = "\n".join(f"  - {lab}" for lab in labels)
            log(
                "ℹ️ [INFO] Available labels/waveforms in PE file (sorted):\n"
                f"{bullet_lines}"
            )
        else:
            log("⚠️ [WARN] No labels/waveforms found in PE file.")

    # If there are no labels at all, nothing more to do
    if not labels:
        log("❌ [ERROR] No labels found – cannot select any label.")
        return None

    # ------------------------------------------------------------------
    # Build a dynamic list of candidate approximants from the labels
    # Example: "C01:IMRPhenomXPHM" → "IMRPhenomXPHM"
    # ------------------------------------------------------------------
    available_aprx: List[str] = []
    for lab in labels:
        try:
            aprx = lab.split(":", 1)[1]
        except Exception:
            continue
        if aprx not in available_aprx:
            available_aprx.append(aprx)

    # Dynamic candidate list:
    #   - requested_approximant (if any) first
    #   - then all approximants actually present in pedata.labels
    candidate_aprx: List[str] = []
    if requested_approximant is not None:
        candidate_aprx.append(requested_approximant)

    for aprx in available_aprx:
        if aprx not in candidate_aprx:
            candidate_aprx.append(aprx)

    # 1a) Warn if requested approximant is not present
    if requested_approximant is not None:
        if find_label_for_approximant(labels, requested_approximant) is None:
            log(
                "⚠️ [WARN] Requested approximant/method "
                f"'{requested_approximant}' not in labels {labels}; "
                "will try other available approximants derived from labels."
            )

    # ------------------------------------------------------------------
    # 1b) Choose label based on candidate approximant list (no PSD check)
    # ------------------------------------------------------------------
    label = None
    for aprx_try in candidate_aprx:
        lab = find_label_for_approximant(labels, aprx_try)
        if lab is not None:
            label = lab
            log(f"✅ [INFO] Requested using label {label} (from method/approximant {aprx_try})")
            break

    # ------------------------------------------------------------------
    # 1c) Try 'Mixed' label as a special case (useful for PE samples)
    #     We do this only if no label chosen yet. We don't treat 'Mixed'
    #     as a normal approximant in the candidate_aprx logic above.
    # ------------------------------------------------------------------
    if label is None:
        mixed_label = next((lab for lab in labels if lab.endswith(":Mixed")), None)
        if mixed_label:
            label = mixed_label
            log(f"ℹ️ [INFO] Falling back to Mixed label: {label}")

    # ------------------------------------------------------------------
    # 1d) Final fallback: first label
    # ------------------------------------------------------------------
    if label is None and labels:
        label = labels[0]
        log(
            "⚠️ [WARN] No label matched any approximant/method; "
            f"using first label {label}"
        )

    if label is None:
        log("❌ [ERROR] No suitable label found in PE dataset.")
        return None

    # -----------------------------------------------------------------
    # If we *don’t* require PSD (PE samples / method choice), we’re done
    # -----------------------------------------------------------------
    if not require_psd:
        return label

    # -----------------------------------------------------------------
    # 2) PSD-enforcing branch (for strain / waveform generation)
    # -----------------------------------------------------------------
    if not label_has_psd(pedata, label):
        log(
            "⚠️ [WARN] Label "
            f"{label} does not provide PSD; searching for PSD-capable labels."
        )

        # Re-try with the same dynamic approximant list but requiring PSD.
        # Here we *intentionally* skip 'Mixed' as a primary PSD candidate
        # and treat it separately later.
        for aprx_try in candidate_aprx:
            if aprx_try == "Mixed":
                continue
            lab = find_label_for_approximant(labels, aprx_try)
            if lab is not None and label_has_psd(pedata, lab):
                log(f"✅ [INFO] Switching to {lab} which has PSD.")
                return lab

        # Try Mixed-with-PSD explicitly, if present
        mixed_label = next((lab2 for lab2 in labels if lab2.endswith(":Mixed")), None)
        if mixed_label and label_has_psd(pedata, mixed_label):
            log(f"✅ [INFO] Using Mixed label {mixed_label} which has PSD.")
            return mixed_label

        # Try ANY label that has PSD
        for lab in labels:
            if label_has_psd(pedata, lab):
                log(f"✅ [INFO] Using {lab} because it provides PSD.")
                return lab

        # No PSD anywhere
        log("❌ [ERROR] No label in PE file provides PSD. Cannot load strain.")
        return None

    # Label has PSD → good for strain
    return label

def label_report(
    pedata,
    sample_method: Optional[str] = None,
    strain_approximant: Optional[str] = None,
):
    """
    Print a compact summary of the PE labels / approximants, including:

      - Available labels (Cxx:Approximant)
      - Parsed approximant / method names
      - PSD availability per label
      - Whether a label is 'Mixed'
      - Which label would be chosen by select_label() for:
          * PE SAMPLES  (require_psd = False, using sample_method)
          * STRAIN      (require_psd = True,  using strain_approximant)

    Parameters
    ----------
    pedata : pesummary.gw.file.File
        The PE data returned by pesummary.read().
    sample_method : str or None
        Requested method / approximant for PE samples (posteriors),
        e.g. "Mixed", "IMRPhenomXPHM". Passed to select_label with
        require_psd=False.
    strain_approximant : str or None
        Requested approximant for strain / waveform generation,
        e.g. "IMRPhenomXPHM". Passed to select_label with require_psd=True.

    Returns
    -------
    dict
        {
          "labels": [...],
          "approximants": [...],
          "has_psd": {...},
          "is_mixed": {...},
          "best_for_sample": <label or None>,
          "best_for_strain": <label or None>,
        }
    """
    try:
        labels = sorted(list(pedata.labels))
    except Exception:
        labels = []

    # Parse approximant names from labels
    approximants: List[str] = []
    for lab in labels:
        try:
            aprx = lab.split(":", 1)[1]
        except Exception:
            aprx = "?"
        approximants.append(aprx)

    # PSD availability and Mixed flags
    has_psd_map = {lab: label_has_psd(pedata, lab) for lab in labels}
    is_mixed_map = {lab: lab.endswith(":Mixed") for lab in labels}

    # What would select_label() pick in the two key modes?
    best_for_sample = select_label(
        pedata,
        requested_approximant=sample_method,
        event_logs=None,
        require_psd=False,
        show_labels=False,
    )
    best_for_strain = select_label(
        pedata,
        requested_approximant=strain_approximant,
        event_logs=None,
        require_psd=True,
        show_labels=False,
    )

    # Header
    print("===")
    print("Label report:")
    print(f"- requested sample_method      (PE samples) : {sample_method}")
    print(f"- requested strain_approximant  (strain)     : {strain_approximant}")
    print(f"- best_for_sample (no PSD requirement) : {best_for_sample}")
    print(f"- best_for_strain  (require PSD=True)   : {best_for_strain}")
    print("===")

    if not labels:
        print("No labels found in PE file.")
        return {
            "labels": [],
            "approximants": [],
            "has_psd": {},
            "is_mixed": {},
            "best_for_sample": None,
            "best_for_strain": None,
        }

    # Compact table header
    header = (
        f"{'idx':>3}  {'label':<24}  {'aprx':<18}  "
        f"{'PSD':<3}  {'Mixed':<5}  {'best_samples':<13}  {'best_strain':<11}"
    )
    print(header)
    print("-" * len(header))

    for i, (lab, aprx) in enumerate(zip(labels, approximants)):
        psd_flag = "✔" if has_psd_map[lab] else "·"
        mixed_flag = "✔" if is_mixed_map[lab] else "·"
        best_s = "✅" if lab == best_for_sample else ""
        best_t = "✅" if lab == best_for_strain else ""
        print(
            f"{i:3d}  {lab:<24}  {aprx:<18}  "
            f"{psd_flag:<3}  {mixed_flag:<5}  {best_s:<13}  {best_t:<11}"
        )

    print("===")

    return {
        "labels": labels,
        "approximants": approximants,
        "has_psd": has_psd_map,
        "is_mixed": is_mixed_map,
        "best_for_sample": best_for_sample,
        "best_for_strain": best_for_strain,
    }

# ---------------------------------------------------------------------
# PE FILE LOCAL CACHE HELPERS (mirror S3 directory structure)
# ---------------------------------------------------------------------


def local_pe_path(object_name: str) -> str:
    """
    Return the local path where the PE file should be stored,
    mirroring the S3 directory structure, but inside the working directory.

    Example:
        object_name = 'GWTC-4/PE/...hdf5'
        → './GWTC-4/PE/...hdf5'
    """
    full_path = os.path.join(".", object_name)
    local_dir = os.path.dirname(full_path)

    if local_dir and not os.path.isdir(local_dir):
        os.makedirs(local_dir, exist_ok=True)

    return full_path


def get_cache_dir() -> str:
    """Return the directory where files will be cached."""
    cache_dir = os.path.expanduser("~/.gwcache")
    os.makedirs(cache_dir, exist_ok=True)
    return cache_dir


def cached_path(filename: str) -> str:
    """Return a full path into the cache directory."""
    return os.path.join(get_cache_dir(), filename)


def file_already_exists(filename: str) -> bool:
    """Check whether a cached file exists."""
    return os.path.isfile(cached_path(filename))


def get_cached_or_download(filename: str, download_fn):
    """
    Generic cache helper:
    - If filename exists in ~/.gwcache → return that path.
    - Otherwise call download_fn(path) to save it, then return path.

    Parameters
    ----------
    filename : str
        File name to store in cache, e.g. 'GW150914_PSD.hdf5'
    download_fn : function(path)
        A function that accepts a local path and writes the file there.

    Returns
    -------
    str
        Full local path to the cached file.
    """
    cache_path = cached_path(filename)

    if os.path.isfile(cache_path):
        print(f"ℹ️ [CACHE] Using cached file: {cache_path}")
        return cache_path

    print(f"ℹ️ [DOWNLOAD] No cached file found → downloading to {cache_path}")
    download_fn(cache_path)
    print(f"ℹ️ [OK] Cached file saved: {cache_path}")
    return cache_path


# ---------------------------------------------------------------------
# Strain helpers
# ---------------------------------------------------------------------
def label_has_psd(pedata, label: str) -> bool:
    """
    Return True if the PE dataset provides PSD for the given label.
    A valid PSD entry must:
      - exist in pedata.psd
      - contain at least one detector (H1, L1, etc.)
      - have non-empty PSD data

    Parameters
    ----------
    pedata : pesummary.gw.file.File
        The PE data returned by pesummary's read().
    label : str
        A Cxx:Approximant label stored in pedata.labels

    Returns
    -------
    bool
        True if PSD exists and is non-empty, otherwise False.
    """
    try:
        return (
            hasattr(pedata, "psd")
            and label in pedata.psd
            and pedata.psd[label] is not None
            and len(pedata.psd[label]) > 0
        )
    except Exception:
        return False


def load_strain(
    event: str,
    t0: float,
    detector: str,
    window: float = 14.0,
    cache: bool = True,
) -> TimeSeries:
    """
    Fetch ~2*window seconds of open data around t0 (GPS) for one detector.

    Parameters
    ----------
    event : str
        Event name, used only for cache file naming (e.g. "GW150914").
    t0 : float
        Reference GPS time (merger or detector-specific).
    detector : str
        Detector name (e.g. "H1", "L1", "V1").
    window : float
        Half-length of the segment around t0 to download, in seconds.
    cache : bool
        If True, cached copies in ~/.gwcache are used/created.

    Returns
    -------
    TimeSeries
        Strain time series.
    """
    start = float(t0) - float(window)
    end = float(t0) + float(window)

    cache_name = f"{event}_{detector}_{int(start)}_{int(end)}_strain.hdf5"
    cache_file = cached_path(cache_name)

    if cache and file_already_exists(cache_name):
        print(f"ℹ️ [CACHE] Using cached strain for {detector}: {cache_name}")
        return TimeSeries.read(cache_file)

    print(f"ℹ️ [INFO] Fetching open strain for {detector} from {start} to {end} ...")
    strain = TimeSeries.fetch_open_data(detector, start, end, cache=False)

    if cache:
        strain.write(cache_file, format="hdf5")
        print(f"ℹ️ [INFO] Cached strain saved to {cache_file}")

    return strain


# ---------------------------------------------------------------------
# Plot / filesystem helpers
# ---------------------------------------------------------------------


def ensure_outdir(outdir: str) -> None:
    """Create output directory if it does not exist."""
    if outdir and not os.path.isdir(outdir):
        os.makedirs(outdir, exist_ok=True)


# ---------------------------------------------------------------------
# Projected waveform utilities
# ---------------------------------------------------------------------


def _get_approximant_for_label(pedata, label: Optional[str]) -> str:
    """Internal helper to pick an approximant string."""
    try:
        labels = list(pedata.labels)
    except Exception:
        labels = []

    if label is not None and label in labels:
        try:
            idx = labels.index(label)
            return pedata.approximant[idx]
        except Exception:
            pass

    if hasattr(pedata, "approximant"):
        approx = pedata.approximant
        if isinstance(approx, str):
            return approx

    # Fallback
    return "IMRPhenomXPHM"


def generate_projected_waveform(
    strain: TimeSeries,
    event: str,
    det: str,
    t0: float,
    pedata,
    label: Optional[str] = None,
    freqrange: Tuple[float, float] = (30.0, 400.0),
    time_window: Tuple[float, float] = (0.2, 0.2),
):
    """
    Build projected waveform on detector `det`, whiten + band-pass + crop.

    Parameters
    ----------
    time_window : (float, float)
        (start_before, stop_after) in seconds, where:
        - start_before > 0 → seconds BEFORE merger (t0 - start_before)
        - stop_after  > 0 → seconds AFTER  merger (t0 + stop_after)
    """

    """

    Returns
    -------
    bp_cropped : TimeSeries or None
        Whitened, band-passed, cropped detector data.
    crop_temp : TimeSeries or None
        Whitened, band-passed, cropped projected waveform.
    used_aprx : str or None
        Name of the approximant actually used to generate the waveform.
    """
    print(f"ℹ️ [INFO] Building projected waveform for {det} ...")

    samples_dict = pedata.samples_dict
    all_labels = list(samples_dict.keys())
    if label is None or label not in all_labels:
        if not all_labels:
            print("⚠️ [WARN] No samples_dict labels found.")
            return None, None, None
        label = all_labels[0]

    posterior_samples = samples_dict[label]

    # Approximant from metadata
    aprx = _get_approximant_for_label(pedata, label)
    print(f"ℹ️ [INFO] Initial approximant guess {aprx} and label {label}")

    # Reference frequency
    fref = 20.0
    try:
        fref = float(pedata.config[label]["engine"]["fref"])
    except Exception:
        try:
            fref = float(pedata.config[label]["config"]["reference-frequency"])
        except Exception:
            pass

    # Choose f_low based on chirp mass
    try:
        loglike = posterior_samples["log_likelihood"]
        maxl_index = loglike.argmax()
        chirp_mass = posterior_samples["chirp_mass"][maxl_index]
        f_low = 60.0 if chirp_mass < 10 else 20.0
    except Exception:
        f_low = 20.0

    # PSD for this label
    try:
        psd_dict = pedata.psd[label]  # dict keyed by 'H1', 'L1', 'V1', ...
    except Exception:
        print("⚠️ [WARN] No PSD in PE file – cannot build projected waveform.")
        return None, None, None

    if det not in psd_dict:
        print(f"⚠️ [WARN] No PSD for detector {det} – skipping waveform.")
        return None, None, None

    zippedpsd = psd_dict[det]
    psdfreq, psdvalue = zip(*zippedpsd)

    # Build ASD to match data sampling
    fs = float(strain.sample_rate.value)
    duration = len(strain) * strain.dt.value
    target_frequencies = np.linspace(0.0, fs / 2.0, int(duration * fs / 2.0), endpoint=False)

    asdsquare = FrequencySeries(
        interp1d(psdfreq, psdvalue, bounds_error=False, fill_value=np.inf)(
            target_frequencies
        ),
        frequencies=target_frequencies,
    )
    asd = np.sqrt(asdsquare)

    # Whiten, band-pass, crop the detector data
    start_before, stop_after = time_window
    cropstart = float(t0) - float(start_before)
    cropend = float(t0) + float(stop_after)

    white_data = strain.whiten(asd=asd)
    bp_data = white_data.bandpass(freqrange[0], freqrange[1])
    bp_cropped = bp_data.crop(cropstart, cropend)

    # Build max-L template projected on det
    # Try the PE approximant first; if it fails (e.g. missing pyseob_wf for SEOBNRv5PHM),
    # fall back to standard LAL waveforms that are usually available.
    hp = None
    used_aprx: Optional[str] = None
    tried: List[str] = []

    # candidate approximants, in order of preference
    candidate_aprx = [aprx, "IMRPhenomXPHM", "IMRPhenomPv2"]

    for aprx_try in candidate_aprx:
        if aprx_try in tried or aprx_try is None:
            continue
        tried.append(aprx_try)

        try:
            print(f"ℹ️ [INFO] Trying maxL_td_waveform with approximant {aprx_try}")
            hp_dict = posterior_samples.maxL_td_waveform(
                aprx_try,
                delta_t=1.0 / fs,
                f_low=f_low,
                f_ref=fref,
                project=det,
            )

            if isinstance(hp_dict, dict):
                hp = hp_dict.get("h_plus", None)
                if hp is None and len(hp_dict):
                    # fallback to the first entry in the dict
                    hp = list(hp_dict.values())[0]
            else:
                hp = hp_dict

            if hp is not None:
                used_aprx = aprx_try
                print(f"[OK] Built waveform with {aprx_try}")
                break

        except NameError as e:
            # This is where "pyseob_wf is not defined" will land
            print(f"⚠️ [WARN] Failed to build waveform with {aprx_try}: {e}")
        except Exception as e:
            print(f"⚠️ [WARN] Failed to build waveform with {aprx_try}: {e}")

    if hp is None:
        print("⚠️ [WARN] Could not build waveform with any approximant; skipping.")
        return None, None, None

    hp = hp.taper()
    hp = hp.pad(int(60 * fs))

    white_temp = hp.whiten(asd=asd)
    bp_temp = white_temp.bandpass(freqrange[0], freqrange[1])
    crop_temp = bp_temp.crop(cropstart, cropend)

    print(f"ℹ️ [OK] Projected waveform for {det} built (approximant {used_aprx}).")
    return bp_cropped, crop_temp, used_aprx


# ---------------------------------------------------------------------
# Plotting utilities
# ---------------------------------------------------------------------

def plot_whitened_overlay(
    bp_cropped: TimeSeries,
    crop_temp: TimeSeries,
    event: str,
    det: str,
    outdir: str| Path,
    approximant: Optional[str] = None,
    t0: Optional[float] = None,
) -> str:

    """
    Save overlay plot: whitened data + projected waveform.
    Shows time relative to referenced merger t0, and puts GPS t0 in the title.
    """
    outdir = Path(outdir)  # normalize
    ensure_outdir(outdir)
    print(f"ℹ️ [INFO] Plotting whitened overlay for {det} at t0={t0}...")

    # Extract absolute GPS times
    t_abs = bp_cropped.times.value

    # If t0 is not provided, approximate with segment midpoint
    if t0 is None:
        t0 = t_abs[len(t_abs)//2]

    # Convert to time relative to t0 (s)
    t_rel_data = t_abs - float(t0)
    t_rel_temp = crop_temp.times.value - float(t0)


    # Compact fonts
    title_font = 7
    label_font = 7
    tick_font = 6
    legend_font = 6

    fig, ax = plt.subplots(figsize=(6, 3.5))

    ax.plot(t_rel_data, bp_cropped.value,
            label="Whitened data", alpha=0.6, linewidth=0.8)
    ax.plot(t_rel_temp, crop_temp.value,
            label="Projected waveform", alpha=0.85, linewidth=0.9)

    ax.set_xlabel("Time relative to merger (s)", fontsize=label_font)
    ax.set_ylabel("Whitened strain", fontsize=label_font)

    # Title
    title = f"{event} – {det} whitened data + projected waveform"
    if approximant:
        title += f"\nApproximant: {approximant},  t0 = {float(t0):.3f} s (GPS)"
    else:
        title += f"\nt0 = {float(t0):.3f} s (GPS)"

    ax.set_title(title, fontsize=title_font, pad=2, linespacing=0.9)

    # Tick formatting: plain numbers, no scientific notation
    ax.ticklabel_format(style="plain", axis="x", useOffset=False)
    ax.tick_params(axis="both", which="major", labelsize=tick_font)

    ax.legend(fontsize=legend_font, loc="upper right", frameon=False)

    fig.tight_layout(pad=1.0)

    fname = os.path.join(outdir, f"{event}_{det}_whitened_waveform.png")
    fig.savefig(fname, dpi=150)
    plt.close(fig)

    print(f"ℹ️ [OK] Saved {fname}")
    return fname


def plot_time_frequency(
    strain: TimeSeries,
    t0: float,
    event: str,
    det: str,
    *,
    outdir: str | Path,
    outseg: Tuple[float, float] = (-2.0, 2.0),
    frange: Tuple[float, float] = (20.0, 512.0),
    approximant: Optional[str] = None,
) -> str:
    """
    Make a q-transform around t0 and save as PNG, using pure matplotlib.

    - X-axis: time relative to t0 (seconds), centered around 0.
    - GPS reference t0 is shown in the title.
    """
    outdir = Path(outdir)  # normalize
    ensure_outdir(outdir)
    print(f"ℹ️ [INFO] Computing q-transform for {det} ...")

    # Time window around t0 in absolute GPS
    seg = (float(t0) + outseg[0], float(t0) + outseg[1])

    # Compute q-transform
    q = strain.q_transform(outseg=seg, frange=frange)

    # Absolute times (GPS) and frequencies
    t_abs = q.times.value
    f = q.frequencies.value
    z = np.abs(q.value)

    # Δt axis (seconds relative to merger)
    t_rel = t_abs - float(t0)

    # Log10 energy
    eps = 1e-24
    z_log = np.log10(np.maximum(z, eps)).T
    vmin, vmax = np.percentile(z_log, [5, 99])

    fig, ax = plt.subplots(figsize=(6, 4))
    im = ax.imshow(
        z_log,
        extent=[t_rel[0], t_rel[-1], f[0], f[-1]],
        origin="lower",
        aspect="auto",
        vmin=vmin,
        vmax=vmax,
    )

    small = 8
    smaller = 7

    # Axes labels and ticks
    ax.set_xlabel("Time relative to merger (s)", fontsize=small)
    ax.set_ylabel("Frequency (Hz)", fontsize=small)
    ax.set_ylim(frange)
    ax.tick_params(axis="both", which="major", labelsize=smaller)

    # Nice number of ticks on Δt axis, plain formatting
    ax.xaxis.set_major_locator(MaxNLocator(6))
    ax.ticklabel_format(style="plain", axis="x", useOffset=False)

    # Title includes GPS reference if available
    
    title = (f"{event} – {det} q-transform\n  t0 = {float(t0):.0f} s (GPS)")
    ax.set_title(title, fontsize=small)

    # Colorbar
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("log10(energy)", fontsize=small)
    cbar.ax.tick_params(labelsize=smaller)

    plt.tight_layout()

    fname = os.path.join(outdir, f"{event}_{det}_qtransform.png")
    fig.savefig(fname, dpi=150)
    plt.close(fig)

    print(f"ℹ️ [OK] Saved {fname}")
    return fname

def plot_basic_posteriors(
    posterior_samples,
    src_name: str,
    outdir: str | Path,
    plot_specs: Optional[list] = None,
    approximant: Optional[str] = None,
) -> Dict[str, str]:
    """
    Produce simple 1D posterior histograms for a set of parameters using matplotlib,
    with median ± 68% CI shown in a small inset box inside the plot.
    """
    outdir = Path(outdir)  # normalize
    ensure_outdir(outdir)

    # Default parameters to plot if not provided
    if plot_specs is None:
        plot_specs = [
            ("mass_1_source",       "mass1",              r"$m_1^{\mathrm{source}}\ [M_\odot]$"),
            ("mass_2_source",       "mass2",              r"$m_2^{\mathrm{source}}\ [M_\odot]$"),
            ("final_mass_source",   "finalmass",          r"$M_f^{\mathrm{source}}\ [M_\odot]$"),
            ("luminosity_distance", "luminositydistance", r"$D_L\ [\mathrm{Mpc}]$"),
            ("ra",                  "RA",                 r"Right ascension [deg]"),
            ("dec",                 "Dec",                r"Declination [deg]"),
        ]

    generated: Dict[str, str] = {}
    available_keys = set(posterior_samples.keys())

    for par_name, basename, xlabel in plot_specs:
        if par_name not in available_keys:
            print(f"⚠️ [WARN] plot_basic_posteriors: parameter '{par_name}' not found; skipping.")
            continue

        samples = np.array(posterior_samples[par_name])

        # Convert RA/Dec from radians → degrees
        if par_name in ["ra", "dec"]:
            samples = np.degrees(samples)

        # ---- Summary statistics: median and 68% CI ----
        median = np.median(samples)
        low, high = np.percentile(samples, [16, 84])
        err_minus = median - low
        err_plus = high - median

        # Text for the inset box
        # (keep it short; 3 sig figs is usually enough)
        stats_text = (
            f"median = {median:.3g}\n"
            f"68% CI: +{err_plus:.3g} / -{err_minus:.3g}"
        )

        # ---- Make the plot ----
        fig, ax = plt.subplots()
        ax.hist(samples, bins=50, density=True, histtype="step")
        ax.set_xlabel(xlabel)
        ax.set_ylabel("Posterior density")

        # Clean title: parameter + event (+ approximant if given)
        title = f"{par_name} – {src_name}"
        if approximant:
            title += f"\nApproximant: {approximant}"
        ax.set_title(title)

        # Inset box with stats in the top-right corner
        ax.text(
            0.97, 0.97,
            stats_text,
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=8,
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.8, edgecolor="none"),
        )

        plt.tight_layout()

        fname = os.path.join(outdir, f"{basename}_{src_name}.png")
        fig.savefig(fname)
        plt.close(fig)

        print(f"ℹ️ [OK] Saved posterior histogram for {par_name} → {fname}")
        generated[par_name] = fname

    return generated




def compare_spectrogram_vs_qtransform(
    strain: TimeSeries,
    t0: float,
    event: str,
    det: str,
    *,
    outdir: str | Path,
    outseg: Tuple[float, float] = (-2.0, 2.0),
    frange: Tuple[float, float] = (20.0, 512.0),
    
) -> str:
    """
    Compare:
      - left: spectrogram of whitened data
      - right: q-transform,
    on the same time–frequency window.

    (You said you don't need this in the notebook anymore, but the function
    is kept here for possible debugging or future use.)
    """
    outdir = Path(outdir)  # normalize
    ensure_outdir(outdir)

    t_start = float(t0) + outseg[0]
    t_end = float(t0) + outseg[1]

    seg = strain.crop(t_start, t_end)

    print(f"ℹ️ [INFO] Whitening segment for spectrogram ({det}) ...")
    white_seg = seg.whiten()

    print(f"ℹ️ [INFO] Computing spectrogram ({det}) ...")
    spec = white_seg.spectrogram(0.1, 0.05)

    print(f"ℹ️ [INFO] Computing q-transform ({det}) ...")
    q = strain.q_transform(outseg=(t_start, t_end), frange=frange)

    spec_t = spec.times.value
    spec_f = spec.frequencies.value
    spec_v = spec.value
    q_t = q.times.value
    q_f = q.frequencies.value
    q_v = q.value

    eps = 1e-24
    spec_log = np.log10(np.maximum(np.abs(spec_v), eps)).T
    q_log = np.log10(np.maximum(np.abs(q_v), eps)).T

    spec_vmin, spec_vmax = np.percentile(spec_log, [5, 99])
    q_vmin, q_vmax = np.percentile(q_log, [5, 99])

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Spectrogram
    im0 = axes[0].imshow(
        spec_log,
        extent=[spec_t[0], spec_t[-1], spec_f[0], spec_f[-1]],
        origin="lower",
        aspect="auto",
        cmap="viridis",
        vmin=spec_vmin,
        vmax=spec_vmax,
    )
    axes[0].set_title("Spectrogram (whitened)")
    axes[0].set_xlabel("Time (s)")
    axes[0].set_ylabel("Frequency (Hz)")
    axes[0].set_ylim(frange)
    fig.colorbar(im0, ax=axes[0], label="log10(power)")

    # Q-transform
    im1 = axes[1].imshow(
        q_log,
        extent=[q_t[0], q_t[-1], q_f[0], q_f[-1]],
        origin="lower",
        aspect="auto",
        cmap="viridis",
        vmin=q_vmin,
        vmax=q_vmax,
    )
    axes[1].set_title("Q-transform")
    axes[1].set_xlabel("Time (s)")
    axes[1].set_ylabel("Frequency (Hz)")
    axes[1].set_ylim(frange)
    fig.colorbar(im1, ax=axes[1], label="log10(energy)")

    fig.suptitle(f"{event} – {det}: Spectrogram (whitened) vs Q-transform")
    plt.tight_layout(rect=[0, 0, 1, 0.95])

    fname = os.path.join(outdir, f"{event}_{det}_spec_vs_q.png")
    plt.savefig(fname, dpi=150)
    plt.close()
    print(f"ℹ️ [OK] Saved {fname}")
    return fname

