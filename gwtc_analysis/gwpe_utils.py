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
import traceback

# Remember last requested waveform engine to support callers that
# do not pass it through to sample-label selection.
_LAST_REQUESTED_WAVEFORM_ENGINE: str | None = None

import matplotlib.pyplot as plt

from gwpy.frequencyseries import FrequencySeries
from scipy.interpolate import interp1d
from matplotlib.ticker import MaxNLocator
from pathlib import Path
from gwpy.timeseries import TimeSeries as GWpyTimeSeries
from pycbc.types import TimeSeries as PyCBCTimeSeries

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


# ---------------------------------------------------------------------
# Debug helpers
# ---------------------------------------------------------------------
# Enable with:  GWPE_DEBUG=1 python -m gwtc_analysis.cli ...
GWPE_DEBUG = os.environ.get("GWPE_DEBUG", "0") not in ("0", "false", "False", "")
GWPE_ALWAYS_TRACE = os.environ.get("GWPE_ALWAYS_TRACE", "0") not in ("0", "false", "False", "")

GWPE_DEBUG_BANNER = (GWPE_DEBUG or GWPE_ALWAYS_TRACE)
if GWPE_DEBUG_BANNER:
    print(f"🐛 [GWPE_DEBUG] gwpe_utils imported from: {__file__}")

def debug_wrap(name: str):
    """Decorator to print a full traceback when the CLI catches exceptions.

    Enable with:
      - GWPE_DEBUG=1 (prints tracebacks for any exception)
      - GWPE_ALWAYS_TRACE=1 (prints tracebacks even if GWPE_DEBUG is off)

    Also: if a TypeError contains 'NoneType' and 'iterable', we print a traceback
    automatically because the CLI otherwise hides the real line number.
    """
    def _decorator(fn):
        def _wrapped(*args, **kwargs):
            try:
                return fn(*args, **kwargs)
            except Exception as e:
                msg = str(e)
                auto = (
                    isinstance(e, TypeError)
                    and ("NoneType" in msg)
                    and ("iterable" in msg or "not iterable" in msg)
                )
                if GWPE_DEBUG or GWPE_ALWAYS_TRACE or auto:
                    print(f"\n🐛 [GWPE_DEBUG] Exception in {name} ({fn.__module__}.{fn.__name__})")
                    print(f"🐛 [GWPE_DEBUG] gwpe_utils file = {__file__}")
                    print(f"🐛 [GWPE_DEBUG] GWPE_DEBUG={GWPE_DEBUG} GWPE_ALWAYS_TRACE={GWPE_ALWAYS_TRACE} auto={auto}")
                    traceback.print_exc()
                    print("🐛 [GWPE_DEBUG] End traceback\n")
                raise
        return _wrapped
    return _decorator
        
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


@debug_wrap('label_has_psd')
def label_has_psd(pedata, label: str) -> bool:
    """
    Robust PSD availability check for a PE label.

    `pesummary` releases sometimes expose `pedata.psd` as None, or
    `pedata.psd[label]` as None, or store PSDs in different containers.
    This function avoids membership tests on None and handles dict/array-like entries.
    """
    psd_container = getattr(pedata, "psd", None)
    if psd_container is None:
        return False

    # Fetch entry safely
    try:
        if isinstance(psd_container, dict):
            psd_entry = psd_container.get(label, None)
        elif hasattr(psd_container, "get"):
            psd_entry = psd_container.get(label, None)
        else:
            psd_entry = psd_container[label]
    except Exception:
        return False

    if psd_entry is None:
        return False

    if isinstance(psd_entry, dict):
        for v in psd_entry.values():
            if v is None:
                continue
            try:
                if len(v) > 0:
                    return True
            except Exception:
                return True
        return False

    try:
        return len(psd_entry) > 0
    except Exception:
        return True


@debug_wrap('select_label')
def select_label(
    pedata,
    pe_label: str | None = None,
    waveform_engine: str | None = None,
    require_psd: bool = False,
    event_logs: list[str] | None = None,
    requested_approximant: str | None = None,
    return_reason: bool = False,
    context: str = "samples",
    show_labels: bool = True,
    **_ignored_kwargs,
):
    """Select which PE label/run to use from a PEDataRelease.

    Most callers expect a *string label*. Set return_reason=True to get (label, reason).

    Compatibility:
      - accepts 'requested_approximant' as an alias for 'waveform_engine'
      - accepts 'context' and 'show_labels' used by older call sites

    Behavior:
      - If waveform_engine is provided (or requested_approximant), select the matching PE label.
      - If not provided, but a previous call provided it, reuse that remembered engine.
        This fixes cases where the pipeline selects posterior-sample labels without
        threading through the CLI waveform-engine option.
    """
    global _LAST_REQUESTED_WAVEFORM_ENGINE

    # Back-compat alias
    if (not waveform_engine) and requested_approximant:
        waveform_engine = requested_approximant

    # Remember the engine when we learn it (typically from strain selection)
    if waveform_engine:
        _LAST_REQUESTED_WAVEFORM_ENGINE = str(waveform_engine).strip()
    elif _LAST_REQUESTED_WAVEFORM_ENGINE:
        # Reuse for callers that don't pass it (e.g. posterior plot selection)
        waveform_engine = _LAST_REQUESTED_WAVEFORM_ENGINE

    samples_dict = getattr(pedata, "samples_dict", {}) or {}
    labels = sorted(list(samples_dict.keys()))
    if not labels:
        raise RuntimeError("No labels found in PE file (samples_dict is empty).")

    if show_labels:
        pe_log("ℹ️ [INFO] Available labels in PE file (sorted):", event_logs)
        for lab in labels:
            pe_log(f"  - {lab}", event_logs)

    def _label_aprx(lab: str) -> str:
        try:
            return str(_get_approximant_for_label(pedata, lab) or "")
        except Exception:
            return ""

    def _has_psd(lab: str) -> bool:
        return label_has_psd(pedata, lab)

    def _ret(lab: str, reason: str):
        pe_log(f"✅ [INFO] Selected label {lab} ({reason})", event_logs)
        return (lab, reason) if return_reason else lab

    # 1) Explicit pe_label (only if provided)
    if pe_label:
        if pe_label in labels:
            return _ret(pe_label, "explicit pe_label")
        needle = str(pe_label).lower()
        for lab in labels:
            if needle in lab.lower():
                return _ret(lab, f"closest match for '{pe_label}'")
        pe_log(f"⚠️ [WARN] Requested pe_label '{pe_label}' not found; falling back.", event_logs)

    # 2) Waveform engine match (preferred when available)
    if waveform_engine:
        eng = str(waveform_engine).strip()
        eng_low = eng.lower()

        candidates = []
        for lab in labels:
            aprx = _label_aprx(lab)
            if eng_low in lab.lower() or (aprx and eng_low in aprx.lower()):
                if (not require_psd) or _has_psd(lab):
                    candidates.append(lab)

        if not candidates:
            token = re.sub(r"[^a-zA-Z0-9]+", "", eng_low)
            for lab in labels:
                aprx = _label_aprx(lab)
                lab_tok = re.sub(r"[^a-zA-Z0-9]+", "", lab.lower())
                aprx_tok = re.sub(r"[^a-zA-Z0-9]+", "", (aprx or "").lower())
                if token and (token in lab_tok or token in aprx_tok):
                    if (not require_psd) or _has_psd(lab):
                        candidates.append(lab)

        if candidates:
            exact = [
                lab for lab in candidates
                if _label_aprx(lab).lower() == eng_low
                or lab.split(":", 1)[-1].lower() == eng_low
            ]
            sel = exact[0] if exact else candidates[0]
            return _ret(sel, f"matched waveform_engine {eng}")

    # 3) Fallback
    if require_psd:
        psd_labels = [lab for lab in labels if _has_psd(lab)]
        if psd_labels:
            return _ret(psd_labels[0], "fallback: first PSD-capable label")

    return _ret(labels[0], "fallback: first label")


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
def load_strain(
    event: str,
    t0: float,
    detector: str,
    window: float = 14.0,
    cache: bool = True,
) -> GWpyTimeSeries:
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
    GWpyTimeSeries
        Strain time series.
    """
    start = float(t0) - float(window)
    end = float(t0) + float(window)

    cache_name = f"{event}_{detector}_{int(start)}_{int(end)}_strain.hdf5"
    cache_file = cached_path(cache_name)

    if cache and file_already_exists(cache_name):
        print(f"ℹ️ [CACHE] Using cached strain for {detector}: {cache_name}")
        return GWpyTimeSeries.read(cache_file)

    print(f"ℹ️ [INFO] Fetching open strain for {detector} from {start} to {end} ...")
    strain = GWpyTimeSeries.fetch_open_data(detector, start, end, cache=False)

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
            
from typing import Optional, Dict, Any, Tuple

def _first_present(d: Dict[str, Any], keys: list[str], default=None):
    for k in keys:
        if k in d:
            return d[k]
    return default

def _maxl_row(posterior_samples, maxl_index: int) -> Dict[str, float]:
    """
    Extract one maxL sample as a flat dict: key -> float.
    posterior_samples is typically a pesummary SamplesDict / pandas-like object.
    """
    row = {}
    # posterior_samples behaves like a dict of arrays
    for k in posterior_samples.keys():
        try:
            v = posterior_samples[k][maxl_index]
            # Convert numpy scalars cleanly
            row[k] = float(v) if np.ndim(v) == 0 else v
        except Exception:
            pass
    return row
    
import astropy.units as u
from lalsimulation import gwsignal

def _replace_key(d: dict, old: str, new: str):
    d[new] = d.pop(old)


def _try_generate_with_alias_rotation(gen, base_params: dict, groups: dict, event_logs=None):
    """
    groups: dict[group_name] = (aliases, value)  where aliases is a list[str] (in preference order)
    """
    # Start by using first alias for each group
    active = dict(base_params)
    active_groups = {k: list(v[0]) for k, v in groups.items()}
    for gname, (aliases, value) in groups.items():
        active[aliases[0]] = value

    for _ in range(20):
        try:
            return gen.generate_td_waveform(**active)
        except Exception as e:
            msg = str(e)
            if "Parameter " in msg and " not in accepted list" in msg:
                bad = msg.split("Parameter ", 1)[1].split(" ", 1)[0].strip()

                # Find which group owns this alias and rotate to next alias
                rotated = False
                for gname, aliases in active_groups.items():
                    if bad in aliases:
                        if len(aliases) <= 1:
                            break
                        # rotate: drop bad alias, switch to next
                        aliases.remove(bad)
                        new_key = aliases[0]
                        val = groups[gname][1]
                        if bad in active:
                            _replace_key(active, bad, new_key)
                        else:
                            active[new_key] = val
                        pe_log(f"⚠️ [WARN] GWSignal rejected '{bad}', retrying with '{new_key}'", event_logs)
                        rotated = True
                        break

                if rotated:
                    continue

            raise  # not an alias issue; propagate

    raise RuntimeError("GWSignal alias rotation exhausted (too many rejections).")

def _debug_gwsignal_param_sources(gen, event_logs=None):
    import inspect

    pe_log(f"ℹ️ [INFO] GWSignal gen type: {type(gen)}", event_logs)
    pe_log(f"ℹ️ [INFO] GWSignal gen dir has metadata: {hasattr(gen, 'metadata')}", event_logs)

    if hasattr(gen, "metadata"):
        try:
            pe_log(f"ℹ️ [INFO] GWSignal metadata keys: {list(gen.metadata.keys())}", event_logs)
            pe_log(f"ℹ️ [INFO] GWSignal metadata:parameters: {gen.metadata.get('parameters')}", event_logs)
        except Exception as e:
            pe_log(f"⚠️ [WARN] Could not read gen.metadata: {e}", event_logs)

    # Try common internal names where parameter lists often live
    candidates = [
        "parameters", "parameter_list", "accepted_parameters", "allowed_parameters",
        "_parameters", "_parameter_list", "_accepted_parameters", "_allowed_parameters",
        "PARAMETERS", "VALID_PARAMS",
    ]
    for name in candidates:
        if hasattr(gen, name):
            try:
                v = getattr(gen, name)
                pe_log(f"ℹ️ [INFO] gen.{name} = {v}", event_logs)
            except Exception as e:
                pe_log(f"⚠️ [WARN] gen.{name} exists but unreadable: {e}", event_logs)

    # Signature (sometimes useless, but cheap)
    try:
        pe_log(f"ℹ️ [INFO] signature(generate_td_waveform) = {inspect.signature(gen.generate_td_waveform)}", event_logs)
    except Exception as e:
        pe_log(f"⚠️ [WARN] Could not inspect signature: {e}", event_logs)

def _to_pycbc_timeseries(x, delta_t: float, epoch: float = 0.0) -> PyCBCTimeSeries:
    if isinstance(x, PyCBCTimeSeries):
        return x

    if hasattr(x, "value") and hasattr(x, "dt") and hasattr(x, "t0"):
        data = np.asarray(x.value)
        dt = float(getattr(x.dt, "value", x.dt))
        t0 = float(getattr(x.t0, "value", x.t0))
        return PyCBCTimeSeries(data, delta_t=dt, epoch=t0)

    data = np.asarray(getattr(x, "value", x))
    dt = float(getattr(x, "delta_t", delta_t))
    t0 = float(getattr(x, "epoch", epoch))
    
    return PyCBCTimeSeries(data, delta_t=dt, epoch=t0)

def _gwsignal_td_waveform_scalar_spins(
    approximant: str,
    params: Dict[str, Any],
    delta_t: float,
    f_low: float,
    f_ref: float,
    event_logs: list[str] | None = None,
):
    """
    Build (hp, hc) at the source (not projected) using LALSimulation GWSignal,
    forcing scalar aligned spins (chi1z/chi2z).

    Returns
    -------
    hp, hc : TimeSeries-like
        GWSignal / PyCBC-compatible plus and cross polarizations.
    """
    # ------------------------------------------------------------------
    # Required physical parameters
    # ------------------------------------------------------------------
    m1 = _first_present(params, ["mass_1", "mass1", "m1", "mass_1_detector", "mass_1_source"])
    m2 = _first_present(params, ["mass_2", "mass2", "m2", "mass_2_detector", "mass_2_source"])
    dist = _first_present(params, ["luminosity_distance", "distance", "d_l", "dl"])
    inc = _first_present(params, ["inclination", "iota"], default=0.0)
    phi0 = _first_present(params, ["phase", "coa_phase", "phi_ref", "phi0"], default=0.0)

    chi1z = _first_present(params, ["chi_1z", "chi1z", "spin_1z", "s1z", "a_1z"])
    chi2z = _first_present(params, ["chi_2z", "chi2z", "spin_2z", "s2z", "a_2z"])

    if m1 is None or m2 is None or dist is None or chi1z is None or chi2z is None:
        missing = [k for k in ("m1", "m2", "dist", "chi1z", "chi2z") if locals().get(k) is None]
        raise RuntimeError(f"Cannot build scalar-spin GWSignal params; missing: {missing}")

    # ------------------------------------------------------------------
    # Frequencies
    # ------------------------------------------------------------------
    f_low = float(f_low)
    f_ref = float(f_ref)
    if f_ref < f_low:
        f_ref = f_low

    # ------------------------------------------------------------------
    # Get GWSignal generator (class or instance depending on version)
    # ------------------------------------------------------------------
    try:
        gen_or_cls = gwsignal.gwsignal_get_waveform_generator(approximant)
    except Exception as e:
        raise RuntimeError(f"Failed to get GWSignal generator for {approximant}: {e}")

    gen = gen_or_cls() if isinstance(gen_or_cls, type) else gen_or_cls

    # ------------------------------------------------------------------
    # Build GWSignal kwargs: ONE alias per physical quantity
    # (avoid passing multiple aliases that trigger strict "accepted list" errors)
    # ------------------------------------------------------------------
    
    base_params = {
        "mass1": float(m1) * u.solMass,
        "mass2": float(m2) * u.solMass,
        "spin1z": float(chi1z) * u.dimensionless_unscaled,
        "spin2z": float(chi2z) * u.dimensionless_unscaled,
        "distance": float(dist) * u.Mpc,
        "inclination": float(inc) * u.rad,
    }
    
    alias_groups = {
        "phase": (
            ["phi_ref", "phase", "reference_phase"],
            float(phi0) * u.rad,
        ),
        "flow": (
            ["f22_start", "f_start","f_low", "f_min", "f_lower"],
            float(f_low) * u.Hz,
        ),
        "dt": (
            ["deltaT", "delta_t", "dt"],
            float(delta_t) * u.s,
        ),
    }

    # ------------------------------------------------------------------
    # Generate time-domain waveform
    # ------------------------------------------------------------------
    
    # Generate with alias rotation
    try:
        out = _try_generate_with_alias_rotation(gen, base_params, alias_groups, event_logs=event_logs)
    except Exception as e:
        raise RuntimeError(f"GWSignal generate_td_waveform failed for {approximant}: {e}")
        try:
            out = gen.generate_td_waveform(**gws_params)
        except Exception as e:
            raise RuntimeError(f"GWSignal generate_td_waveform failed for {approximant}: {e}")

    # ------------------------------------------------------------------
    # Extract h+ / h×
    # ------------------------------------------------------------------
    hp = None
    hc = None

    if isinstance(out, dict):
        hp = out.get("h_plus") or out.get("hp") or out.get("plus") or out.get("hplus")
        hc = out.get("h_cross") or out.get("hc") or out.get("cross") or out.get("hcross")
    else:
        try:
            hp, hc = out
        except Exception:
            pass

    if hp is None or hc is None:
        keys_or_type = list(out.keys()) if isinstance(out, dict) else type(out)
        raise RuntimeError(
            "GWSignal td waveform output not understood "
            f"(keys/type={keys_or_type}); cannot find hp/hc."
        )

    # Ensure PyCBC TimeSeries for downstream projection/time shifting
    hp = _to_pycbc_timeseries(hp, delta_t=float(delta_t), epoch=0.0)
    hc = _to_pycbc_timeseries(hc, delta_t=float(delta_t), epoch=0.0)
    
    return hp, hc






def _project_hp_hc_to_detector(
    hp, hc,
    det: str,
    ra: float,
    dec: float,
    psi: float,
    geocent_time: float,
):
    """
    Project (hp, hc) onto a detector using antenna patterns.
    Returns a PyCBC TimeSeries-like strain for that detector.
    """
    try:
        from pycbc.detector import Detector
    except Exception as e:
        raise RuntimeError(f"pycbc is required for detector projection but not available: {e}")

    ifo = Detector(det)

    # Antenna pattern
    fp, fc = ifo.antenna_pattern(ra, dec, psi, geocent_time)

    # Time delay from geocenter to detector
    dt = ifo.time_delay_from_earth_center(ra, dec, geocent_time)

    # Combine
    ht = fp * hp + fc * hc

    # Shift by dt (PyCBC TimeSeries has .cyclic_time_shift or .time_shift depending on type)
    try:
        ht = ht.time_shift(dt)
    except Exception:
        try:
            ht = ht.cyclic_time_shift(dt)
        except Exception:
            # If ht is not a PyCBC type, you’ll need to convert it first.
            raise RuntimeError("Projected waveform object does not support time shifting; convert to PyCBC TimeSeries.")

    return ht
    
# ------------------------------------------------------------------
# Scalar / time helpers (PEViewer-compatible)
# ------------------------------------------------------------------
def _sec(x) -> float:
    """Convert Quantity / LIGOTimeGPS / numeric to float seconds."""
    if x is None:
        raise TypeError("Cannot convert None to seconds")
    try:
        return float(x.to_value("s"))
    except Exception:
        pass
    try:
        return float(x)
    except Exception:
        pass
    return float(getattr(x, "value"))

def _scalar(x) -> float:
    """Convert Quantity / numpy scalar / python numeric to plain float."""
    if x is None:
        raise TypeError("Cannot convert None to float")
    if hasattr(x, "to_value"):
        try:
            return float(x.to_value(x.unit))
        except Exception:
            pass
    if hasattr(x, "value") and not isinstance(x.value, (str, bytes)):
        try:
            return float(x.value)
        except Exception:
            pass
    return float(x)




def _embed_template_on_strain_grid(template: GWpyTimeSeries, strain: GWpyTimeSeries, fs: float) -> GWpyTimeSeries:
    """Embed a (possibly short) template onto the strain time grid with zero padding.

    Used for SEOBNRv5 scalar-spin adapter outputs so whitening uses a consistent grid.
    """
    n = int(len(strain))
    out = np.zeros(n, dtype=float)

    try:
        t0_strain = _sec(strain.t0.value)
    except Exception:
        t0_strain = _sec(getattr(strain, "t0", 0.0))

    try:
        t0_tpl = _sec(template.t0.value)
    except Exception:
        t0_tpl = _sec(getattr(template, "t0", 0.0))

    start_idx = int(round((t0_tpl - t0_strain) * fs))
    tpl_vals = np.asarray(template.value, dtype=float)

    i0 = max(0, start_idx)
    j0 = max(0, -start_idx)
    ncopy = min(n - i0, len(tpl_vals) - j0)
    if ncopy > 0:
        out[i0:i0+ncopy] = tpl_vals[j0:j0+ncopy]

    return GWpyTimeSeries(out, sample_rate=fs, t0=t0_strain)

def _as_gwpy_timeseries(ts, fs: float) -> GWpyTimeSeries:
    """Convert a PyCBC (or array-like) TimeSeries to GWpy TimeSeries, preserving epoch when possible."""
    if isinstance(ts, GWpyTimeSeries):
        return ts
    # PyCBC TimeSeries
    if isinstance(ts, PyCBCTimeSeries):
        data = np.asarray(ts.numpy())
        try:
            t0 = _sec(getattr(ts, "start_time"))
        except Exception:
            try:
                t0 = _sec(getattr(ts, "epoch"))
            except Exception:
                t0 = 0.0
        return GWpyTimeSeries(data, sample_rate=fs, t0=t0)
    # Something with .value and .t0
    if hasattr(ts, "value") and hasattr(ts, "t0"):
        data = np.asarray(ts.value)
        try:
            t0 = _sec(getattr(ts, "t0"))
        except Exception:
            t0 = 0.0
        return GWpyTimeSeries(data, sample_rate=fs, t0=t0)
    # Fallback: raw array, unknown epoch
    data = np.asarray(ts)
    return GWpyTimeSeries(data, sample_rate=fs, t0=0.0)
    
@debug_wrap('generate_projected_waveform')
def generate_projected_waveform(
    strain: GWpyTimeSeries,
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
    """PEViewer-style projected waveform + whitened overlay inputs.

    This reproduces PEViewer's logic:
      - crop time window is centered on the event geocenter GPS (GWOSC datasets.event_gps),
        not on a detector-specific maxL time
      - ASD is built by interpolating the PE PSD onto a target frequency grid with
        fill_value=np.inf (so whitening ignores frequencies outside PSD support)
      - template comes directly from posterior_samples.maxL_td_waveform(..., project=det)
        and is then tapered + padded, whitened, bandpassed, cropped

    SEOBNRv5 special case:
      - if requested_approximant starts with 'SEOBNRv5' and maxL_td_waveform fails due
        to invalid vector spins, we fall back to a scalar aligned-spin adapter using
        GWSignal + PyCBC detector projection.
    """
    # Normalize numeric inputs (Quantities -> floats)
    f_lo, f_hi = _scalar(freqrange[0]), _scalar(freqrange[1])
    dt_before, dt_after = _sec(time_window[0]), _sec(time_window[1])

    pe_log(f"ℹ️ [INFO] Building projected waveform for {det} ...", event_logs)

    # Posterior samples label
    samples_dict = pedata.samples_dict
    all_labels = list(samples_dict.keys())
    if label is None or label not in all_labels:
        if not all_labels:
            pe_log("⚠️ [WARN] No samples_dict labels found.", event_logs)
            return None, None, None
        label = all_labels[0]
    posterior_samples = samples_dict[label]

    # Approximant guess from label
    aprx_guess = _get_approximant_for_label(pedata, label)
    pe_log(f"ℹ️ [INFO] Initial approximant guess {aprx_guess} and label {label}", event_logs)

    # Reference frequency
    fref = 20.0
    try:
        fref = float(pedata.config[label]["engine"]["fref"])
    except Exception:
        try:
            fref = float(pedata.config[label]["config"]["reference-frequency"])
        except Exception:
            pass

    # f_low heuristic (as in PEViewer)
    try:
        loglike = posterior_samples["log_likelihood"]
        maxl_index = int(np.argmax(loglike))
        chirp_mass = posterior_samples["chirp_mass"][maxl_index]
        f_low = 60.0 if chirp_mass < 10 else 20.0
    except Exception:
        maxl_index = 0
        f_low = 20.0

    # --- PSD / ASD from PE for this label+det ---
    # PSDs in PE files are not always available (or may be stored as None / missing detector).
    # In that case, fall back to GWPy's internal PSD estimate via strain.whiten() (no asd=...).
    asd = None
    fs = int(_scalar(strain.sample_rate.value))
    duration = len(strain) * _sec(strain.dt.value)

    try:
        psd_container = getattr(pedata, "psd", None)
        psd_entry = psd_container.get(label, None) if hasattr(psd_container, "get") else (psd_container[label] if psd_container is not None else None)
    except Exception:
        psd_entry = None

    if isinstance(psd_entry, dict):
        zippedpsd = psd_entry.get(det, None)
        if zippedpsd is None:
            pe_log(f"⚠️ [WARN] No PE PSD for detector {det} under label {label}; whitening will use GWPy PSD estimate.", event_logs)
        else:
            # Accept: list-of-(f, psd), Nx2 array, or (freqs, psdvals)
            try:
                arr = np.asarray(zippedpsd)
                if arr.ndim == 2 and arr.shape[1] == 2:
                    psdfreq = np.asarray(arr[:, 0], dtype=float)
                    psdvalue = np.asarray(arr[:, 1], dtype=float)
                elif isinstance(zippedpsd, (list, tuple)) and len(zippedpsd) == 2:
                    psdfreq = np.asarray(zippedpsd[0], dtype=float)
                    psdvalue = np.asarray(zippedpsd[1], dtype=float)
                else:
                    # Fallback: try unpacking iterable of pairs
                    psdfreq, psdvalue = zip(*zippedpsd)
                    psdfreq = np.asarray(psdfreq, dtype=float)
                    psdvalue = np.asarray(psdvalue, dtype=float)

                target_frequencies = np.linspace(0, fs / 2, int(duration * fs / 2), endpoint=False)
                asdsquare = FrequencySeries(
                    interp1d(psdfreq, psdvalue, bounds_error=False, fill_value=np.inf)(target_frequencies),
                    frequencies=target_frequencies,
                )
                asd = np.sqrt(asdsquare)
            except Exception as e:
                pe_log(f"⚠️ [WARN] Failed to parse/interpolate PE PSD for {det} ({label}): {e}; using GWPy PSD estimate.", event_logs)
                asd = None
    else:
        pe_log(f"⚠️ [WARN] No PE PSD dict for label {label}; whitening will use GWPy PSD estimate.", event_logs)
# Crop window is centered on GWOSC event geocenter time (PEViewer)
    try:
        from gwosc import datasets
        t0_center = _sec(datasets.event_gps(event))
    except Exception:
        # fallback to provided t0 if GWOSC lookup unavailable
        t0_center = _sec(t0)

    cropstart = t0_center - dt_before
    cropend = t0_center + dt_after

    # Whiten + bandpass + crop strain
    white_data = strain.whiten(asd=asd) if asd is not None else strain.whiten()
    bp_data = white_data.bandpass(f_lo, f_hi)
    bp_cropped = bp_data.crop(cropstart, cropend)

    # --- Build projected template (PEViewer) ---
    tries: List[str] = []
    if requested_approximant:
        tries.append(str(requested_approximant))
    if aprx_guess and aprx_guess not in tries:
        tries.append(aprx_guess)
    if allow_fallback:
        for fb in ("IMRPhenomXPHM", "IMRPhenomPv2"):
            if fb not in tries:
                tries.append(fb)

    pe_log(f"ℹ️ [INFO] Approximants to try (in order): {tries}", event_logs)

    hp = None
    used_aprx = None

    for aprx_try in tries:
        try:
            pe_log(f"ℹ️ [INFO] Trying maxL_td_waveform with approximant {aprx_try}", event_logs)
            hp = posterior_samples.maxL_td_waveform(
                aprx_try,
                delta_t=1 / fs,
                f_low=f_low,
                f_ref=fref,
                project=det,
            )
            used_aprx = aprx_try
            pe_log(f"[OK] Built waveform with {aprx_try}", event_logs)
            break
        except Exception as e:
            pe_log(f"⚠️ [WARN] Failed to build waveform with {aprx_try}: {e}", event_logs)

            # Special case: SEOBNRv5 vector-spin failure -> scalar-spin adapter fallback
            if requested_approximant and str(requested_approximant).upper().startswith("SEOBNRV5"):
                emsg = str(e).lower()
                spin_hint = ("chi1 has to be in" in emsg) or ("chi_1" in emsg) or ("seobnrv5" in emsg)
                if spin_hint:
                    pe_log("ℹ️ [INFO] Retrying SEOBNRv5 using scalar aligned spins (chi1z/chi2z) adapter...", event_logs)
                    try:
                        row = _maxl_row(posterior_samples, maxl_index)
                        ra = _scalar(_first_present(row, ["ra", "right_ascension"], default=0.0))
                        dec = _scalar(_first_present(row, ["dec", "declination"], default=0.0))
                        psi = _scalar(_first_present(row, ["psi", "polarization"], default=0.0))
                        tgps = _sec(_first_present(row, ["geocent_time", "gps_time", "tc"], default=t0_center))

                        hp_src, hc_src = _gwsignal_td_waveform_scalar_spins(
                            approximant=str(requested_approximant),
                            params=row,
                            delta_t=1 / fs,
                            f_low=f_low,
                            f_ref=fref,
                            event_logs=event_logs,
                        )
                        hp_proj = _project_hp_hc_to_detector(hp_src, hc_src, det, ra, dec, psi, tgps)

                        # Place adapter waveform in absolute time by aligning its peak to crop center
                        hp_proj = _to_pycbc_timeseries(hp_proj, delta_t=1 / fs, epoch=0.0)
                        arr = np.asarray(hp_proj.numpy())
                        if arr.size > 0:
                            peak = int(np.argmax(np.abs(arr)))
                            hp_start = t0_center - peak * (1 / fs)
                            hp_proj._epoch = hp_start

                        hp = hp_proj
                        used_aprx = str(requested_approximant) + " (scalar-spin adapter)"
                        pe_log(f"[OK] Built waveform with {used_aprx}", event_logs)
                        break
                    except Exception as ee:
                        pe_log(f"⚠️ [WARN] Scalar-spin adapter failed: {ee}", event_logs)
                        hp = None
                        used_aprx = None
                        continue

    if hp is None:
        pe_log("⚠️ [WARN] Could not build waveform with any approximant; skipping.", event_logs)
        return None, None, None

    if requested_approximant and used_aprx and used_aprx != requested_approximant:
        pe_log(f"⚠️ [WARN] Requested approximant {requested_approximant} but used {used_aprx} fallback.", event_logs)

    # Template processing exactly like PEViewer
    try:
        hp = hp.taper()
    except Exception:
        pass
    try:
        hp = hp.pad(60 * fs)
    except Exception:
        pass

    white_temp = hp_gw = _as_gwpy_timeseries(hp, fs)
    if used_aprx and "scalar-spin adapter" in str(used_aprx):
        hp_gw = _embed_template_on_strain_grid(hp_gw, strain, fs)
    white_temp = hp_gw.whiten(asd=asd) if asd is not None else hp_gw.whiten()
    bp_temp = white_temp.bandpass(f_lo, f_hi)
    crop_temp = bp_temp.crop(cropstart, cropend)

    pe_log(f"ℹ️ [OK] Projected waveform for {det} built (approximant {used_aprx}).", event_logs)
    return bp_cropped, crop_temp, used_aprx

def _sec(x) -> float:
    """Convert Quantity/LIGOTimeGPS/float-like to float seconds."""
    if x is None:
        raise TypeError("Cannot convert None to seconds")
    # astropy Quantity
    try:
        return float(x.to_value("s"))  # Quantity -> seconds
    except Exception:
        pass
    # LIGOTimeGPS / numpy scalar / python float
    try:
        return float(x)
    except Exception:
        pass
    # last resort: objects with .value
    return float(getattr(x, "value"))

# ---------------------------------------------------------------------
# Plotting utilities
# ---------------------------------------------------------------------

@debug_wrap('plot_whitened_overlay')
def plot_whitened_overlay(
    bp_cropped: GWpyTimeSeries,
    crop_temp: GWpyTimeSeries,
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
    t0s = _sec(t0)
    t_rel_data = t_abs - t0s
    t_rel_temp = crop_temp.times.value - t0s


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


@debug_wrap('plot_time_frequency')
def plot_time_frequency(
    strain: GWpyTimeSeries,
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
    strain: GWpyTimeSeries,
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

