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


def pe_log(msg: str, event_logs: list[str] | None = None) -> None:
    """
    Global logger for PE pipeline.
    - ALWAYS prints to stdout
    - ALWAYS appends to event_logs (for HTML report)
    """
    print(msg)
    if event_logs is not None:
        event_logs[-1] += f"\n{msg}\n"
        
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


def find_closest_label_by_token(
    pedata,
    labels: list[str],
    token: str,
    *,
    require_psd: bool,
):
    """Return the closest PE label whose RHS contains the given token.

    This is a pure string-based match (no hardcoded waveform names).

    Scoring heuristics (generic, waveform-agnostic):
      - exact RHS match > prefix match > substring match
      - shorter RHS preferred (smaller suffix)
      - when require_psd=True, only PSD-capable labels are considered
    """

    def _norm(s: str) -> str:
        return re.sub(r"\s+", "", str(s)).lower()

    if not token:
        return None

    tok = _norm(token)
    best = None
    best_score = None

    for lab in labels:
        if require_psd and not label_has_psd(pedata, lab):
            continue

        if ":" not in lab:
            continue

        rhs_raw = lab.split(":", 1)[1]
        rhs = _norm(rhs_raw)

        if tok not in rhs:
            continue

        score = 0
        if rhs == tok:
            score += 1000
        if rhs.startswith(tok):
            score += 500
        score += 200  # substring baseline
        score -= len(rhs_raw)  # prefer shorter strings

        if best_score is None or score > best_score:
            best = lab
            best_score = score

    return best


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
    require_psd: bool = True,
    show_labels: bool = True,
    context: str = "samples",  # "samples" or "strain"
    event_logs: list[str] | None = None,
):
    """
    Select the best PE label for a given context (PE samples or strain/PSD), with:

      1. Direct PE label match if user provided e.g. "C00:SEOBNRv5PHM",
      2. Requested approximant/method (matched against label RHS),
      3. Other available approximants found in pedata.labels,
      4. Optionally a 'Mixed' label (useful for PE samples),
      5. First label as last resort,

    and *optionally* enforcing that the final choice has PSD available.

    Notes
    -----
    - In context="samples", requested_approximant is treated as a PE method/label RHS
      (e.g. "Mixed", "SEOBNRv5PHM") and warnings will mention PE labels.
    - In context="strain", requested_approximant may be a waveform *engine* name
      (e.g. "IMRPhenomPv2") which may not exist as a PE label RHS; in this case
      we do NOT emit the confusing "not in labels" warning and instead explain
      that we will pick the closest PSD-capable PE label.

    Returns
    -------
    str or None
        The chosen PE label, or None if no suitable label was found.
    """
    import re
    from typing import List

    # Sorted labels for nicer logging
    try:
        labels = sorted(list(pedata.labels))
    except Exception:
        labels = []

    def _normalize_requested(s: str | None) -> str | None:
        """Normalize a user-requested approximant/method/label string.

        - strips whitespace and trailing punctuation from CLI copy/paste
        - if a full PE label like 'C00:SEOBNRv5PHM' is provided, keeps only the RHS
        """
        if s is None:
            return None
        s2 = str(s).strip().rstrip(",;")
        if re.match(r"^C\d{2}:", s2):
            s2 = s2.split(":", 1)[1].strip()
        return s2 or None

    requested_raw = requested_approximant
    requested_norm = _normalize_requested(requested_approximant)

    # ------------------------------------------------------------------
    # 0) If the user directly provided a full PE label, honor it.
    # ------------------------------------------------------------------
    requested_label = None
    if isinstance(requested_raw, str):
        requested_label = str(requested_raw).strip().rstrip(",;")

    if requested_label and requested_label in labels:
        pe_log(f"✅ [INFO] Requested using label {requested_label} (direct label match)",event_logs)
        if (not require_psd) or label_has_psd(pedata, requested_label):
            return requested_label
        pe_log(f"⚠️ [WARN] Requested label {requested_label} has no PSD; will select another label.",event_logs)

    # ------------------------------------------------------------------
    # 1) Pretty, multi-line, sorted label listing
    # ------------------------------------------------------------------
    if show_labels:
        if labels:
            bullet_lines = "\n".join(f"  - {lab}" for lab in labels)
            pe_log(
                "ℹ️ [INFO] Available labels in PE file (sorted):\n"
                f"{bullet_lines}",event_logs
            )
        else:
            pe_log("⚠️ [WARN] No labels found in PE file.",event_logs)

    if not labels:
        pe_log("❌ [ERROR] No labels found – cannot select any label.",event_logs)
        return None

    # ------------------------------------------------------------------
    # 2) Build list of available approximants from label RHS.
    #     Example: "C01:IMRPhenomXPHM" → "IMRPhenomXPHM"
    # ------------------------------------------------------------------
    available_aprx: List[str] = []
    for lab in labels:
        try:
            rhs = lab.split(":", 1)[1]
        except Exception:
            continue
        if rhs not in available_aprx:
            available_aprx.append(rhs)

    # ------------------------------------------------------------------
    # 3) Dynamic candidate list:
    #    - requested_norm (if any) first
    #    - then all RHS values present in pedata.labels
    # ------------------------------------------------------------------
    candidate_aprx: List[str] = []
    if requested_norm is not None:
        candidate_aprx.append(requested_norm)
    for rhs in available_aprx:
        if rhs not in candidate_aprx:
            candidate_aprx.append(rhs)

    # ------------------------------------------------------------------
    # 3a) If a request was provided, first try a *substring* match on RHS.
    #     This yields (e.g.) IMRPhenomXPHM -> IMRPhenomXPHM-SpinTaylor
    #     without hardcoding waveform names.
    # ------------------------------------------------------------------
    label = None
    if requested_norm is not None:
        closest = find_closest_label_by_token(
            pedata,
            labels,
            requested_norm,
            require_psd=require_psd,
        )
        if closest is not None:
            label = closest
            rhs = label.split(":", 1)[1] if ":" in label else "?"
            pe_log(
                f"✅ [INFO] Selected closest label {label} (matched token '{requested_raw}' in '{rhs}')",
                event_logs,
            )
        else:
            # Context-aware message when no substring match exists
            if context == "samples":
                pe_log(
                    "⚠️ [WARN] Requested method/label "
                    f"'{requested_raw}' not found in PE labels {labels}; "
                    "will try other labels derived from the PE file.",
                    event_logs,
                )
            else:
                pe_log(
                    f"ℹ️ [INFO] Requested strain engine '{requested_raw}' does not match any PE label RHS; "
                    "selecting the closest PSD-capable label from the PE file.",
                    event_logs,
                )

    # ------------------------------------------------------------------
    # 4) If no substring match was found, try exact RHS matches from the
    #    candidate list (no PSD check yet; PSD enforcement happens later).
    # ------------------------------------------------------------------
    if label is None:
        for rhs_try in candidate_aprx:
            lab = find_label_for_approximant(labels, rhs_try)
            if lab is not None:
                label = lab
                pe_log(f"✅ [INFO] Selected label {label} (from method {rhs_try})", event_logs)
                break

    # ------------------------------------------------------------------
    # 5) Try 'Mixed' label as a special case (useful for PE samples)
    #    IMPORTANT: only do this when the user did not request a specific
    #    approximant/engine. If the user requested an approximant, "Mixed"
    #    should NOT silently override that intent.
    # ------------------------------------------------------------------
    if label is None and requested_norm is None:
        mixed_label = next((lab for lab in labels if lab.endswith(":Mixed")), None)
        if mixed_label:
            label = mixed_label
            pe_log(f"ℹ️ [INFO] Falling back to Mixed label: {label}", event_logs)

    # ------------------------------------------------------------------
    # 6) Final fallback: first label
    # ------------------------------------------------------------------
    if label is None:
        label = labels[0]
        pe_log(
            "⚠️ [WARN] No label matched any method/approximant; "
            f"using first label {label}", event_logs
        )

    if label is None:
        pe_log("❌ [ERROR] No suitable label found in PE dataset.", event_logs)
        return None

    # -----------------------------------------------------------------
    # If we *don’t* require PSD (PE samples / method choice), we’re done
    # -----------------------------------------------------------------
    if not require_psd:
        return label

    # -----------------------------------------------------------------
    # 7) PSD-enforcing branch (for strain / waveform generation)
    # -----------------------------------------------------------------
    if not label_has_psd(pedata, label):
        pe_log(
            "⚠️ [WARN] Selected label "
            f"{label} does not provide PSD; searching for PSD-capable labels.", event_logs
        )

        # If the user requested a specific approximant/engine, prefer a PSD-capable
        # label whose RHS *contains* that token (substring match), before any other
        # fallback. This keeps PSD/maxL inputs consistent with the user's intent
        # without hardcoding waveform names.
        if requested_norm is not None:
            closest_psd = find_closest_label_by_token(
                pedata,
                labels,
                requested_norm,
                require_psd=True,
            )
            if closest_psd is not None:
                pe_log(f"✅ [INFO] Switching to {closest_psd} which matches '{requested_raw}' and has PSD.", event_logs)
                return closest_psd

        # Re-try with the same RHS candidate list but requiring PSD.
        # Intentionally skip 'Mixed' as a primary PSD candidate.
        for rhs_try in candidate_aprx:
            if rhs_try == "Mixed":
                continue
            lab = find_label_for_approximant(labels, rhs_try)
            if lab is not None and label_has_psd(pedata, lab):
                pe_log(f"✅ [INFO] Switching to {lab} which has PSD.", event_logs)
                return lab

        # Try Mixed-with-PSD explicitly, if present
        mixed_label = next((lab2 for lab2 in labels if lab2.endswith(":Mixed")), None)
        if mixed_label and label_has_psd(pedata, mixed_label):
            pe_log(f"✅ [INFO] Using Mixed label {mixed_label} which has PSD.", event_logs)
            return mixed_label

        # Try ANY label that has PSD
        for lab in labels:
            if label_has_psd(pedata, lab):
                pe_log(f"✅ [INFO] Using {lab} because it provides PSD.", event_logs)
                return lab

        # No PSD anywhere
        pe_log("❌ [ERROR] No label in PE file provides PSD. Cannot load strain.", event_logs)
        return None

    # Label has PSD → good for strain
    return label


def label_report(
    pedata,
    pe_label: Optional[str] = None,
    waveform_engine: Optional[str] = None,
):
    """
    print a compact summary of the PE labels / approximants, including:

      - Available labels (Cxx:Approximant)
      - Parsed approximant / method names
      - PSD availability per label
      - Whether a label is 'Mixed'
      - Which label would be chosen by select_label() for:
          * PE SAMPLES  (require_psd = False, using pe_label)
          * STRAIN      (require_psd = True,  using waveform_engine)

    Parameters
    ----------
    pedata : pesummary.gw.file.File
        The PE data returned by pesummary.read().
    pe_label : str or None
        Requested method / approximant for PE samples (posteriors),
        e.g. "Mixed", "IMRPhenomXPHM". Passed to select_label with
        require_psd=False.
    waveform_engine : str or None
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
        requested_approximant=pe_label,
        event_logs=None,
        require_psd=False,
        show_labels=False,
    )
    best_for_strain = select_label(
        pedata,
        requested_approximant=waveform_engine,
        event_logs=None,
        require_psd=True,
        show_labels=False,
    )

    # Header
    print("===")
    print("Label report:")
    print(f"- requested pe_label      (PE samples) : {pe_label}")
    print(f"- requested waveform_engine  (strain)     : {waveform_engine}")
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
    """Return the approximant associated with a PE label, without inventing one.

    Preference order:
    1) pedata.approximant[idx] if label is in pedata.labels
    2) parse from label string like 'C00:SEOBNRv5PHM' -> 'SEOBNRv5PHM'
    3) if pedata.approximant is a single string, use it
    4) otherwise return empty string (unknown)
    """
    # 1) indexed lookup via pedata.labels / pedata.approximant
    labels: list[str] = []
    try:
        labels = list(pedata.labels)
    except Exception:
        labels = []

    if label and labels and label in labels:
        try:
            idx = labels.index(label)
            aprx = pedata.approximant[idx]
            if isinstance(aprx, str) and aprx.strip():
                return aprx.strip()
        except Exception:
            pass

    # 2) derive from label string
    if label and isinstance(label, str) and ":" in label:
        rhs = label.split(":", 1)[1].strip()
        if rhs:
            return rhs

    # 3) sometimes pedata.approximant is a single string
    if hasattr(pedata, "approximant"):
        approx = pedata.approximant
        if isinstance(approx, str) and approx.strip():
            return approx.strip()

    # 4) unknown
    return ""
        
from typing import Callable, Optional
            
def generate_projected_waveform(
    strain: TimeSeries,
    event: str,
    det: str,
    t0: float,
    pedata,
    label: Optional[str] = None,
    freqrange: Tuple[float, float] = (30.0, 400.0),
    time_window: Tuple[float, float] = (0.2, 0.2),
    requested_approximant: Optional[str] = None,
    allow_fallback: bool = True,
    event_logs: list[str] | None = None,
):
    """
    Build projected waveform on detector `det`, whiten + band-pass + crop.

    Parameters
    ----------
    time_window : (float, float)
        (start_before, stop_after) in seconds, where:
        - start_before > 0 → seconds BEFORE merger (t0 - start_before)
        - stop_after  > 0 → seconds AFTER  merger (t0 + stop_after)

    requested_approximant : str or None
        If provided, try this approximant first (and optionally only this one).

    allow_fallback : bool
        If False, do not fall back to generic LAL waveforms if the requested/label
        approximant fails. This is useful to avoid silently producing overlays with
        a different waveform than the user asked for.

    Returns
    -------
    bp_cropped : TimeSeries or None
        Whitened, band-passed, cropped detector data.
    crop_temp : TimeSeries or None
        Whitened, band-passed, cropped projected waveform.
    used_aprx : str or None
        Name of the approximant actually used to generate the waveform.
    """
    import builtins

    pe_log(f"ℹ️ [INFO] Building projected waveform for {det} ...",event_logs)

    samples_dict = pedata.samples_dict
    all_labels = list(samples_dict.keys())
    if label is None or label not in all_labels:
        if not all_labels:
            pe_log("⚠️ [WARN] No samples_dict labels found.",event_logs)
            return None, None, None
        label = all_labels[0]

    posterior_samples = samples_dict[label]

    # Approximant from metadata
    aprx = _get_approximant_for_label(pedata, label)
    pe_log(f"ℹ️ [INFO] Initial approximant guess {aprx} and label {label}",event_logs)

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
        pe_log("⚠️ [WARN] No PSD in PE file – cannot build projected waveform.",event_logs)
        return None, None, None

    if det not in psd_dict:
        pe_log(f"⚠️ [WARN] No PSD for detector {det} – skipping waveform.",event_logs)
        return None, None, None

    zippedpsd = psd_dict[det]
    psdfreq, psdvalue = zip(*zippedpsd)

    # Build ASD to match data sampling
    fs = float(strain.sample_rate.value)
    duration = len(strain) * strain.dt.value
    target_frequencies = np.linspace(0.0, fs / 2.0, int(duration * fs / 2.0), endpoint=False)

    asdsquare = FrequencySeries(
        interp1d(psdfreq, psdvalue, bounds_error=False, fill_value=np.inf)(target_frequencies),
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
    hp = None
    last_err: Optional[str] = None
    used_aprx: Optional[str] = None
    tried: List[str] = []

    # Candidate approximants, in order of preference:
    # 1) user-requested approximant (if any)
    # 2) PE-derived approximant from metadata
    # 3) generic fallbacks (only if allow_fallback=True)
    candidate_aprx: List[str] = []

    if requested_approximant:
        candidate_aprx.append(requested_approximant)

    if aprx and aprx not in candidate_aprx:
        candidate_aprx.append(aprx)

    if allow_fallback:
        for fb in ("IMRPhenomXPHM", "IMRPhenomPv2"):
            if fb not in candidate_aprx:
                candidate_aprx.append(fb)

    pe_log(f"ℹ️ [INFO] Approximants to try (in order): {candidate_aprx}",event_logs)

    for aprx_try in candidate_aprx:
        if aprx_try in tried or aprx_try is None:
            continue
        tried.append(aprx_try)

        try:
            pe_log(f"ℹ️ [INFO] Trying maxL_td_waveform with approximant {aprx_try}",event_logs)
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
                    hp = list(hp_dict.values())[0]
            else:
                hp = hp_dict

            if hp is not None:
                used_aprx = aprx_try
                pe_log(f"[OK] Built waveform with {aprx_try}",event_logs)
                break

        except NameError as e:
            last_err = str(e)
            pe_log(f"⚠️ [WARN] Failed to build waveform with {aprx_try}: {e}",event_logs)
        except Exception as e:
            last_err = str(e)
            pe_log(f"⚠️ [WARN] Failed to build waveform with {aprx_try}: {e}",event_logs)

    if hp is None:
        if requested_approximant and not allow_fallback:
            pe_log(
                f"⚠️ [WARN] Could not build waveform with requested approximant "
                f"{requested_approximant}; fallback disabled, skipping.",event_logs
            )
        else:
            pe_log("⚠️ [WARN] Could not build waveform with any approximant; skipping.",event_logs)
        return None, None, None

    if requested_approximant and used_aprx != requested_approximant:
        pe_log(
            f"⚠️ [WARN] Requested approximant {requested_approximant} but used {used_aprx} fallback.",event_logs
        )

    hp = hp.taper()
    hp = hp.pad(int(60 * fs))

    white_temp = hp.whiten(asd=asd)
    bp_temp = white_temp.bandpass(freqrange[0], freqrange[1])
    crop_temp = bp_temp.crop(cropstart, cropend)

    pe_log(f"ℹ️ [OK] Projected waveform for {det} built (approximant {used_aprx}).",event_logs)
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
    pe_label: str | None = None,
    engine_requested: str | None = None,
    engine_used: str | None = None,
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
    title1 = f"{event} – {det} whitened data + projected waveform"

    # Build a truthful second line
    parts = []
    if engine_requested and engine_used and engine_requested != engine_used:
        parts.append(f"Engine: {engine_used} (fallback from {engine_requested})")
    elif engine_used:
        parts.append(f"Engine: {engine_used}")
    elif engine_requested:
        parts.append(f"Engine: {engine_requested}")

    parts.append(f"t0 = {t0:.3f} s (GPS)")
    title2 = " | ".join(parts)

    ax.set_title(title1 + "\n" + title2, fontsize=8)

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

from typing import Dict, Optional, Tuple, List, Any

def plot_posterior_pairs(
    posterior_samples: dict[str, Any],
    src_name: str,
    outdir: str | Path,
    *,
    pairs: list[str] | None = None,
    approximant: str | None = None,
    bins: int = 60,
    max_points_scatter: int = 30000,
) -> Dict[str, Any]:
    """
    Plot 2D posterior pairs given as tokens 'x:y'.

    Returns:
      - 'plots': dict 'x:y' -> filepath
      - 'available_keys': sorted list of keys
      - 'missing_pairs': list of tokens that couldn't be plotted (missing key or bad format)
      - 'plotted_pairs': list of tokens successfully plotted
    """
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    pairs = pairs or []
    available = set(posterior_samples.keys())
    available_list = sorted(map(str, available))

    plots: dict[str, str] = {}
    missing_pairs: list[str] = []
    plotted_pairs: list[str] = []

    def _get_arr(key: str) -> np.ndarray:
        arr = np.array(posterior_samples[key])
        if key in ("ra", "dec"):
            arr = np.degrees(arr)
        return arr

    for tok in pairs:
        if ":" not in tok:
            missing_pairs.append(tok)
            continue

        xk, yk = tok.split(":", 1)
        xk = xk.strip()
        yk = yk.strip()

        if (xk not in available) or (yk not in available):
            missing_pairs.append(tok)
            continue

        x = _get_arr(xk)
        y = _get_arr(yk)

        # Drop non-finite
        m = np.isfinite(x) & np.isfinite(y)
        x = x[m]
        y = y[m]
        if len(x) < 10:
            missing_pairs.append(tok)
            continue

        # If huge, sub-sample scatter overlay
        do_scatter = len(x) <= max_points_scatter
        if not do_scatter:
            # sub-sample for scatter if you still want it
            idx = np.random.choice(len(x), size=max_points_scatter, replace=False)
            xs = x[idx]
            ys = y[idx]
        else:
            xs, ys = x, y

        fig, ax = plt.subplots(figsize=(6, 5))

        # 2D histogram density
        h = ax.hist2d(x, y, bins=bins)

        # optional scatter overlay (helps when bins are coarse)
        ax.scatter(xs, ys, s=3, alpha=0.15)

        ax.set_xlabel(xk if xk not in ("ra", "dec") else f"{xk} [deg]")
        ax.set_ylabel(yk if yk not in ("ra", "dec") else f"{yk} [deg]")

        title = f"{xk} vs {yk} – {src_name}"
        if approximant:
            title += f"\nApproximant: {approximant}"
        ax.set_title(title, fontsize=10)

        plt.tight_layout()

        safe = f"{xk}_vs_{yk}".replace("/", "_").replace(" ", "_")
        fname = os.path.join(outdir, f"pair_{safe}_{src_name}.png")
        fig.savefig(fname, dpi=150)
        plt.close(fig)

        print(f"ℹ️ [OK] Saved posterior 2D plot {tok} → {fname}")
        plots[tok] = fname
        plotted_pairs.append(tok)

    return {
        "plots": plots,
        "available_keys": available_list,
        "missing_pairs": sorted(set(missing_pairs)),
        "plotted_pairs": plotted_pairs,
    }

def plot_basic_posteriors(
    posterior_samples,
    src_name: str,
    outdir: str | Path,
    plot_specs: Optional[list[tuple[str, str, str]]] = None,
    approximant: Optional[str] = None,
    *,
    extra_params: Optional[list[str]] = None,
) -> Dict[str, Any]:
    """
    Produce 1D posterior histograms for a set of parameters using matplotlib,
    with median ± 68% CI shown in a small inset box.

    Returns a dict with:
      - 'plots': dict par_name -> filepath
      - 'available_keys': sorted list of keys
      - 'missing_requested': list of requested keys not found
      - 'plotted': list of plotted keys
    """
    outdir = Path(outdir)
    ensure_outdir(outdir)

    # What is available in the posterior object?
    available_keys = set(posterior_samples.keys())
    available_list = sorted(map(str, available_keys))

    # Default parameters to plot if not provided
    default_specs: list[tuple[str, str, str]] = [
        ("mass_1_source",       "mass1",              r"$m_1^{\mathrm{source}}\ [M_\odot]$"),
        ("mass_2_source",       "mass2",              r"$m_2^{\mathrm{source}}\ [M_\odot]$"),
        ("final_mass_source",   "finalmass",          r"$M_f^{\mathrm{source}}\ [M_\odot]$"),
        ("luminosity_distance", "luminositydistance", r"$D_L\ [\mathrm{Mpc}]$"),
        ("ra",                  "RA",                 r"Right ascension [deg]"),
        ("dec",                 "Dec",                r"Declination [deg]"),
    ]

    # If caller provided explicit plot_specs, respect it exactly.
    # Otherwise use defaults + any extra_params requested by user.
    if plot_specs is None:
        plot_specs = list(default_specs)

        # Add user-requested variables (if any), using generic labels
        # basename must be filesystem-safe and stable
        extra_params = extra_params or []
        seen = {p for p, _, _ in plot_specs}

        for par_name in extra_params:
            if par_name in seen:
                continue
            seen.add(par_name)
            basename = par_name.replace("/", "_").replace(" ", "_")
            xlabel = par_name  # generic label; can be improved later with a mapping
            plot_specs.append((par_name, basename, xlabel))

    generated: Dict[str, str] = {}
    missing_requested: list[str] = []
    plotted: list[str] = []

    # Track missing only for user-requested keys (not defaults)
    requested_set = set(extra_params or [])

    for par_name, basename, xlabel in plot_specs:
        if par_name not in available_keys:
            if par_name in requested_set:
                missing_requested.append(par_name)
            print(f"⚠️ [WARN] plot_basic_posteriors: parameter '{par_name}' not found; skipping.")
            continue

        samples = np.array(posterior_samples[par_name])

        # Convert RA/Dec from radians → degrees
        if par_name in ["ra", "dec"]:
            samples = np.degrees(samples)

        # ---- Summary statistics ----
        median = np.median(samples)
        low, high = np.percentile(samples, [16, 84])
        err_minus = median - low
        err_plus = high - median

        stats_text = (
            f"median = {median:.3g}\n"
            f"68% CI: +{err_plus:.3g} / -{err_minus:.3g}"
        )

        # ---- Plot ----
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.hist(samples, bins=50, density=True, histtype="step")
        ax.set_xlabel(xlabel)
        ax.set_ylabel("Posterior density")

        title = f"{par_name} – {src_name}"
        if approximant:
            title += f"\nApproximant: {approximant}"
        ax.set_title(title)

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
        fig.savefig(fname, dpi=150)
        plt.close(fig)

        print(f"ℹ️ [OK] Saved posterior histogram for {par_name} → {fname}")
        generated[par_name] = fname
        plotted.append(par_name)

    return {
        "plots": generated,
        "available_keys": available_list,
        "missing_requested": sorted(set(missing_requested)),
        "plotted": plotted,
    }




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

