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


def run_parameters_estimation(args: Any) -> int:
    """Entry point: run the PE plotting pipeline.

    Parameters
    ----------
    args:
        Object/dict/Namespace with fields:
        src_name, sample_method, strain_approximant, start, stop, fs_low, fs_high

    Returns
    -------
    int
        0 on success, non-zero on failure.
    """

    global LAST_RESULT

    # --- read inputs (exactly those defined at the top of the notebook) ---
    src_name = _get_arg(args, "src_name")
    sample_method = _get_arg(args, "sample_method")
    strain_approximant = _get_arg(args, "strain_approximant")
    start = float(_get_arg(args, "start"))
    stop = float(_get_arg(args, "stop"))
    fs_low = float(_get_arg(args, "fs_low"))
    fs_high = float(_get_arg(args, "fs_high"))

    # --- imports kept inside to make module import cheap and Galaxy-friendly ---
    import os
    import json
    import warnings

    warnings.filterwarnings("ignore", category=UserWarning, append=True)
    try:
        from astropy.wcs import FITSFixedWarning

        warnings.simplefilter("ignore", category=FITSFixedWarning)
    except Exception:
        pass

    warnings.filterwarnings("ignore", "Wswiglal-redir-stdio")

    import lal

    lal.swig_redirect_standard_output_error(False)

    import matplotlib

    # Headless backends for Galaxy batch environments
    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt
    from matplotlib.colorbar import Colorbar

    from minio import Minio
    from pesummary.io import read

    # gwpe_utils are assumed to be provided by your package environment
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
        from oda_api.data_products import PictureProduct
        from oda_api.api import ProgressReporter

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
        # Mirror notebook behaviour: append to current log chunk
        event_logs[-1] += f"\n{msg}\n"

    def _progress(stage: str, progress: int, substage: str) -> None:
        if pr is not None:
            try:
                pr.report_progress(stage=stage, progress=progress, substage=substage)
            except Exception:
                pass

    # ---------------------------------------------------------------------
    # Download/read PE file from S3
    # ---------------------------------------------------------------------
    _progress("Download data", 10, "step 1")

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

    data = None
    label_waveform: Optional[str] = None

    try:
        found_file = False
        remote_file: Optional[str] = None

        for obj in client.list_objects("gwtc", recursive=True):
            file_name = obj.object_name
            if (
                src_name in file_name
                and "PEDataRelease" in file_name
                and (file_name.endswith(".h5") or file_name.endswith(".hdf5"))
            ):
                remote_file = file_name
                found_file = True
                _log(f"ℹ️ [INFO] Found remote PE file: {file_name}")
                break

        if not found_file:
            msg = f"❌ [ERROR] Failed reading data. No such event name {src_name} found in catalogs"
            _log(msg)
            go_next_cell = False

        if found_file and go_next_cell and remote_file is not None:
            local_pe_path = remote_file  # keep same relative structure
            local_dir = os.path.dirname(local_pe_path)
            if local_dir and not os.path.isdir(local_dir):
                os.makedirs(local_dir, exist_ok=True)

            if os.path.isfile(local_pe_path):
                _log(f"ℹ️ [CACHE] Using existing local PE file: {local_pe_path}")
            else:
                _log(f"ℹ️ [DOWNLOAD] Fetching from S3: {remote_file}")
                client.fget_object("gwtc", remote_file, local_pe_path)

            _log(f"ℹ️ [INFO] Reading PE data from: {local_pe_path}")
            _progress("Read data", 25, "step 2")

            try:
                data = read(local_pe_path)
                label_report(data, sample_method=sample_method, strain_approximant=strain_approximant)
            except Exception as e:
                msg = f"❌ [ERROR] Failed reading PE data: {e}"
                _log(msg)
                go_next_cell = False

    except Exception as e:
        msg = f"❌ [ERROR] Failed reading data. Unable to access S3 repository: {e}"
        _log(msg)
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
                msg = f"❌ [ERROR] PE samples: no suitable label found; labels present: {labels}"
                _log(msg)
                go_next_cell = False
            else:
                label_waveform = label
                _log(f"ℹ️ [INFO] PE samples: final label = {label}")

                posterior_samples = samples_dict[label]

                # Derive approximant name from the label, if possible
                approximant_for_title = label.split(":", 1)[1] if ":" in label else label

                posterior_files = plot_basic_posteriors(
                    posterior_samples,
                    src_name,
                    outdir=".",
                    approximant=approximant_for_title,
                )

                for _par_name, fname in posterior_files.items():
                    result.files_distribution.append(fname)
                    if oda_available:
                        fig_distributionList.append(PictureProduct.from_file(fname))

        except Exception as e:
            msg = f"❌ [ERROR] Failed creating posterior samples: {e}"
            _log(msg)
            # keep going? notebook stops only in some cases; keep conservative:
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
                msg = (
                    "❌ [ERROR] Waveform model: no suitable label found in the PE file; "
                    "unable to load strain / build waveform."
                )
                _log(msg)
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
                        outdir=".",
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
                        outseg=(-start_before, stop_after),
                        frange=(fs_low, fs_high),
                        outdir=".",
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
            fname_psd = os.path.join(".", f"{src_name}_{waveform}_psd.png")
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
            fname_sky = os.path.join(".", f"{src_name}_{waveform}_skymap.png")
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
    # Finalize outputs
    # ---------------------------------------------------------------------
    result.tool_log = event_logs
    result.fig_distribution = fig_distributionList
    result.fig_strain = fig_strainList
    result.fig_psd = fig_psdList
    result.fig_skymap = fig_skymapList

    LAST_RESULT = result

    # Attach outputs back to args (Galaxy/MMODA pattern)
    _set_output(args, "fig_distribution", fig_distributionList if oda_available else result.files_distribution)
    _set_output(args, "fig_strain", fig_strainList if oda_available else result.files_strain)
    _set_output(args, "fig_psd", fig_psdList if oda_available else result.files_psd)
    _set_output(args, "fig_skymap", fig_skymapList if oda_available else result.files_skymap)
    _set_output(args, "tool_log", event_logs)

    return 0 if go_next_cell else 1
