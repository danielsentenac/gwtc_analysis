from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Callable


@dataclass(frozen=True)
class UnofficialPEAnalysisSpec:
    dataset_name: str
    label: str
    approximant: str
    lalinference_samples_path: Path | None = None


@dataclass(frozen=True)
class UnofficialPEBundleSpec:
    event_id: str
    raw_samples_path: Path | None
    psd_path: Path | None
    skymap_path: Path
    output_filename: str
    analyses: tuple[UnofficialPEAnalysisSpec, ...]
    gps_time: float
    asd_paths: tuple[tuple[str, Path], ...] = ()


_PE_SAMPLES_REPO = Path.home() / "GW170817"

GW170817_SPEC = UnofficialPEBundleSpec(
    event_id="GW170817",
    raw_samples_path=Path.home() / ".gwcache" / "GW170817_GWTC-1.hdf5",
    psd_path=Path.home() / ".gwcache" / "GWTC1_GW170817_PSDs.dat",
    skymap_path=Path.home() / ".gwcache" / "GW170817_skymap.fits.gz",
    output_filename="Unofficial-GWTC1-GW170817_PEDataRelease.h5",
    analyses=(
        UnofficialPEAnalysisSpec(
            dataset_name="IMRPhenomPv2NRT_lowSpin_posterior",
            label="C02:IMRPhenomPv2_NRTidal-LowSpin",
            approximant="IMRPhenomPv2_NRTidal",
            lalinference_samples_path=_PE_SAMPLES_REPO / "MCMC_IMRPPNRT_LowSpin_post.dat",
        ),
        UnofficialPEAnalysisSpec(
            dataset_name="IMRPhenomPv2NRT_highSpin_posterior",
            label="C02:IMRPhenomPv2_NRTidal-HighSpin",
            approximant="IMRPhenomPv2_NRTidal",
            lalinference_samples_path=_PE_SAMPLES_REPO / "MCMC_IMRPPNRT_HighSpin_post.dat",
        ),
    ),
    gps_time=1187008882.429464,
    asd_paths=(
        ("H1", _PE_SAMPLES_REPO / "psd" / "BayesWave_PSD_H1_IFO0_asd_median.dat"),
        ("L1", _PE_SAMPLES_REPO / "psd" / "BayesWave_PSD_L1_IFO0_asd_median.dat"),
        ("V1", _PE_SAMPLES_REPO / "psd" / "BayesWave_PSD_V1_IFO0_asd_median.dat"),
    ),
)


UNOFFICIAL_PE_BUNDLES: dict[str, UnofficialPEBundleSpec] = {
    GW170817_SPEC.event_id: GW170817_SPEC,
}


def list_unofficial_pe_specs() -> list[str]:
    return sorted(UNOFFICIAL_PE_BUNDLES.keys())


def get_unofficial_pe_spec(src_name: str) -> UnofficialPEBundleSpec | None:
    return UNOFFICIAL_PE_BUNDLES.get(str(src_name).strip())


def build_unofficial_pe_bundle(
    src_name: str,
    *,
    cache_dir: str | Path = ".cache_gwosc",
    log_cb: Callable[[str], None] | None = None,
    force_rebuild: bool = False,
) -> Path | None:
    spec = get_unofficial_pe_spec(src_name)
    if spec is None:
        return None

    sources: list[Path] = [spec.skymap_path.expanduser()]
    if spec.asd_paths:
        sources.extend(Path(p).expanduser() for _, p in spec.asd_paths)
    elif spec.psd_path is not None:
        sources.append(spec.psd_path.expanduser())

    needs_raw_hdf5 = any(a.lalinference_samples_path is None for a in spec.analyses)
    if needs_raw_hdf5:
        if spec.raw_samples_path is None:
            if log_cb:
                log_cb(
                    "⚠️ [WARN] Unofficial PE bundle for "
                    f"{spec.event_id} has analyses without lalinference_samples_path "
                    "but no raw_samples_path is configured."
                )
            return None
        sources.append(spec.raw_samples_path.expanduser())
    for analysis in spec.analyses:
        if analysis.lalinference_samples_path is not None:
            sources.append(analysis.lalinference_samples_path.expanduser())

    missing = [str(p) for p in sources if not p.exists()]
    if missing:
        if log_cb:
            log_cb(
                "⚠️ [WARN] Unofficial PE bundle for "
                f"{spec.event_id} cannot be built; missing source files: {', '.join(missing)}"
            )
        return None

    target_dir = Path(cache_dir) / "unofficial_pe"
    target_dir.mkdir(parents=True, exist_ok=True)
    target = target_dir / spec.output_filename

    if not force_rebuild and target.exists() and target.stat().st_size > 0:
        target_mtime = target.stat().st_mtime
        if all(p.stat().st_mtime <= target_mtime for p in sources):
            if log_cb:
                log_cb(f"ℹ️ [CACHE] Using unofficial PE bundle: {target}")
            return target

    if log_cb:
        if force_rebuild and target.exists():
            log_cb(f"ℹ️ [BUILD] Rebuilding unofficial PE bundle for {spec.event_id}: {target}")
        log_cb(f"ℹ️ [BUILD] Building unofficial PE bundle for {spec.event_id} from local cache files")

    _write_unofficial_pesummary_bundle(spec, target, log_cb=log_cb)

    if log_cb:
        log_cb(f"ℹ️ [OK] Built unofficial PE bundle: {target}")
    return target


def _write_unofficial_pesummary_bundle(
    spec: UnofficialPEBundleSpec,
    out_path: Path,
    *,
    log_cb: Callable[[str], None] | None = None,
) -> None:
    import h5py
    import numpy as np
    from astropy.io import fits
    from pesummary.gw.file.formats.pesummary import write_pesummary
    from pesummary.gw.file.skymap import SkyMap
    from pesummary.utils.samples_dict import MultiAnalysisSamplesDict, SamplesDict

    samples_by_label: dict[str, SamplesDict] = {}

    h5f = None
    needs_raw_hdf5 = any(a.lalinference_samples_path is None for a in spec.analyses)
    try:
        if needs_raw_hdf5:
            assert spec.raw_samples_path is not None  # validated earlier
            h5f = h5py.File(spec.raw_samples_path.expanduser(), "r")
        for analysis in spec.analyses:
            if analysis.lalinference_samples_path is not None:
                if log_cb:
                    log_cb(
                        f"ℹ️ [INFO] {analysis.label}: using full LALInference posterior "
                        f"({analysis.lalinference_samples_path.name}) — psi/phase/time/logl are real."
                    )
                samples_by_label[analysis.label] = _build_samples_from_lalinference(
                    analysis.lalinference_samples_path.expanduser(),
                    gps_time=spec.gps_time,
                )
            else:
                if log_cb:
                    log_cb(
                        f"⚠️ [INFO] {analysis.label}: using public {spec.raw_samples_path.name} dataset "
                        f"{analysis.dataset_name!r} — psi/phase from priors, log_likelihood synthetic."
                    )
                if analysis.dataset_name not in h5f:
                    raise KeyError(
                        f"Dataset {analysis.dataset_name!r} not found in {spec.raw_samples_path}"
                    )
                dataset = h5f[analysis.dataset_name][:]
                samples_by_label[analysis.label] = _build_samples_dict(dataset, gps_time=spec.gps_time)
    finally:
        if h5f is not None:
            h5f.close()

    if spec.asd_paths:
        if log_cb:
            log_cb(
                "ℹ️ [INFO] Using per-detector BayesWave-median ASDs (squared to PSD) "
                f"for {len(spec.asd_paths)} detector(s)."
            )
        psds = _read_per_detector_asds(spec.asd_paths)
    elif spec.psd_path is not None:
        psds = _read_multidetector_psd(spec.psd_path.expanduser())
    else:
        psds = {}
    skymap = _read_skymap(spec.skymap_path.expanduser(), event_id=spec.event_id, gps_time=spec.gps_time)

    labels = [analysis.label for analysis in spec.analyses]
    approximant = {analysis.label: analysis.approximant for analysis in spec.analyses}
    file_versions = {label: "GWTC-1 unofficial bundle" for label in labels}
    file_kwargs = {label: {"sampler": {}, "meta_data": {}} for label in labels}
    psd_by_label = {label: psds for label in labels}
    skymap_by_label = {label: skymap for label in labels}

    if out_path.exists():
        out_path.unlink()

    write_pesummary(
        MultiAnalysisSamplesDict(samples_by_label),
        outdir=str(out_path.parent),
        filename=out_path.name,
        file_versions=file_versions,
        file_kwargs=file_kwargs,
        approximant=approximant,
        psd=psd_by_label,
        skymap=skymap_by_label,
    )


def _build_samples_dict(dataset, *, gps_time: float):
    import numpy as np
    from pesummary.utils.samples_dict import SamplesDict

    base = {name: np.asarray(dataset[name], dtype=float) for name in dataset.dtype.names}
    return _finalize_samples_from_base(base, gps_time=gps_time)


def _build_samples_from_lalinference(path: Path, *, gps_time: float):
    """Read a LALInference posterior_samples.dat and build a PESummary SamplesDict.

    LALInference column names differ from PESummary's standard set; in particular
    ``time`` is the geocenter coalescence time, ``{h1,l1,v1}_end_time`` are the
    detector arrival times, and ``logl``/``phi12``/``tilt1``/``tilt2``/``costilt1``/
    ``costilt2``/``cosiota``/``costheta_jn``/``a1z``/``a2z``/``lam_tilde``/
    ``dlam_tilde``/``mc``/``distance``/``m1``/``m2``/``mtotal``/``mf``/``af``
    have analogous renames. ``standardize_parameter_names`` covers most; the
    remaining few are renamed explicitly below.
    """
    import numpy as np

    arr = np.genfromtxt(path, names=True)
    base = {name: np.asarray(arr[name], dtype=float) for name in arr.dtype.names}

    rename_pairs = (
        ("time", "geocent_time"),
        ("h1_end_time", "H1_time"),
        ("l1_end_time", "L1_time"),
        ("v1_end_time", "V1_time"),
        ("phi12", "phi_12"),
        ("logl", "log_likelihood"),
        ("logprior", "log_prior"),
        ("logpost", "log_posterior"),
    )
    for src, dst in rename_pairs:
        if src in base and dst not in base:
            base[dst] = base.pop(src)

    return _finalize_samples_from_base(base, gps_time=gps_time)


def _finalize_samples_from_base(base: dict, *, gps_time: float):
    import numpy as np
    from pesummary.utils.samples_dict import SamplesDict

    samples = SamplesDict(base).standardize_parameter_names()
    nsamples = len(next(iter(samples.values())))

    # Parameters not exposed by the public GWTC-1 release. Sample from the
    # LALInference priors so the bundle is self-consistent: maxL projection
    # via posterior_samples.maxL_td_waveform() relies on (psi, phase) for the
    # antenna-pattern factor F+·cos(2ψ) + F×·sin(2ψ); a constant ψ=0 silently
    # drops the F× contribution and biases the projected amplitude. Seed from
    # gps_time so rebuilds are reproducible.
    rng = np.random.default_rng(seed=int(gps_time) & 0xFFFFFFFF)

    if "geocent_time" not in samples:
        samples["geocent_time"] = np.full(nsamples, float(gps_time), dtype=float)
    if "phase" not in samples:
        samples["phase"] = rng.uniform(0.0, 2.0 * np.pi, size=nsamples)
    if "psi" not in samples:
        samples["psi"] = rng.uniform(0.0, np.pi, size=nsamples)
    if "phi_jl" not in samples:
        samples["phi_jl"] = rng.uniform(0.0, 2.0 * np.pi, size=nsamples)
    if "phi_12" not in samples:
        samples["phi_12"] = rng.uniform(0.0, 2.0 * np.pi, size=nsamples)

    _ensure_compatibility_columns(samples)
    _ensure_detector_arrival_times(samples, gps_time=gps_time)
    if "log_likelihood" not in samples:
        samples["log_likelihood"] = _synthetic_log_likelihood(samples)

    try:
        samples.generate_all_posterior_samples()
    except Exception:
        # The bundle only needs a minimal compatible set; keep explicitly-added
        # columns if PESummary cannot derive the full conversion set.
        pass

    _ensure_compatibility_columns(samples)
    _ensure_detector_arrival_times(samples, gps_time=gps_time)
    return samples


def _ensure_detector_arrival_times(samples, *, gps_time: float) -> None:
    """Add H1/L1/V1 detector arrival times derived from geocenter time and sky position."""
    import numpy as np

    if "ra" not in samples or "dec" not in samples:
        return

    try:
        import lal
    except Exception:
        return

    try:
        nsamples = len(next(iter(samples.values())))
    except Exception:
        return

    geocent_times = np.asarray(
        samples.get("geocent_time", np.full(nsamples, float(gps_time), dtype=float)),
        dtype=float,
    )
    if geocent_times.ndim == 0:
        geocent_times = np.full(nsamples, float(geocent_times), dtype=float)

    ra = np.asarray(samples["ra"], dtype=float)
    dec = np.asarray(samples["dec"], dtype=float)

    detectors = {
        "H1": lal.LALDetectorIndexLHODIFF,
        "L1": lal.LALDetectorIndexLLODIFF,
        "V1": lal.LALDetectorIndexVIRGODIFF,
    }

    for det, detector_index in detectors.items():
        key = f"{det}_time"
        if key in samples:
            continue

        detector = lal.CachedDetectors[detector_index]
        arrivals = np.empty(nsamples, dtype=float)
        for idx in range(nsamples):
            tc = float(geocent_times[idx])
            if not (np.isfinite(tc) and np.isfinite(ra[idx]) and np.isfinite(dec[idx])):
                arrivals[idx] = tc if np.isfinite(tc) else float(gps_time)
                continue
            delay = lal.TimeDelayFromEarthCenter(
                detector.location,
                float(ra[idx]),
                float(dec[idx]),
                lal.LIGOTimeGPS(tc),
            )
            arrivals[idx] = tc + float(delay)

        samples[key] = arrivals


def _ensure_compatibility_columns(samples) -> None:
    import numpy as np

    if "theta_jn" not in samples and "cos_theta_jn" in samples:
        samples["theta_jn"] = np.arccos(np.clip(np.asarray(samples["cos_theta_jn"], dtype=float), -1.0, 1.0))
    if "tilt_1" not in samples and "cos_tilt_1" in samples:
        samples["tilt_1"] = np.arccos(np.clip(np.asarray(samples["cos_tilt_1"], dtype=float), -1.0, 1.0))
    if "tilt_2" not in samples and "cos_tilt_2" in samples:
        samples["tilt_2"] = np.arccos(np.clip(np.asarray(samples["cos_tilt_2"], dtype=float), -1.0, 1.0))
    if "iota" not in samples and "theta_jn" in samples:
        samples["iota"] = np.asarray(samples["theta_jn"], dtype=float)
    if "spin_1z" not in samples and "a_1" in samples and "cos_tilt_1" in samples:
        samples["spin_1z"] = np.asarray(samples["a_1"], dtype=float) * np.asarray(samples["cos_tilt_1"], dtype=float)
    if "spin_2z" not in samples and "a_2" in samples and "cos_tilt_2" in samples:
        samples["spin_2z"] = np.asarray(samples["a_2"], dtype=float) * np.asarray(samples["cos_tilt_2"], dtype=float)


def _synthetic_log_likelihood(samples) -> "object":
    import numpy as np

    anchor_params = [
        "mass_1",
        "mass_2",
        "luminosity_distance",
        "ra",
        "dec",
        "a_1",
        "a_2",
        "cos_theta_jn",
        "cos_tilt_1",
        "cos_tilt_2",
        "lambda_1",
        "lambda_2",
    ]
    present = [p for p in anchor_params if p in samples]
    if not present:
        n = len(next(iter(samples.values())))
        return np.zeros(n, dtype=float)

    score = np.zeros(len(np.asarray(samples[present[0]], dtype=float)), dtype=float)
    for param in present:
        vals = np.asarray(samples[param], dtype=float)
        if vals.size == 0:
            continue
        center = float(np.nanmedian(vals))
        spread = float(np.nanstd(vals))
        if not np.isfinite(spread) or spread <= 0:
            spread = 1.0
        score -= ((vals - center) / spread) ** 2
    return score


def _read_multidetector_psd(path: Path) -> dict[str, "object"]:
    import numpy as np

    data = np.genfromtxt(path, comments="#")
    if data.ndim != 2 or data.shape[1] < 2:
        raise ValueError(f"PSD file {path} does not contain at least two columns")

    detector_columns = ("H1", "L1", "V1")
    out: dict[str, object] = {}
    for idx, detector in enumerate(detector_columns, start=1):
        if idx >= data.shape[1]:
            break
        out[detector] = np.column_stack([data[:, 0], data[:, idx]])

    if not out:
        raise ValueError(f"Could not extract detector PSD columns from {path}")
    return out


def _read_per_detector_asds(asd_paths: tuple[tuple[str, Path], ...]) -> dict[str, "object"]:
    """Read per-detector ASD files (freq, asd) and return PSD dict (asd squared)."""
    import numpy as np

    out: dict[str, object] = {}
    for detector, path in asd_paths:
        data = np.genfromtxt(Path(path).expanduser(), comments="#")
        if data.ndim != 2 or data.shape[1] < 2:
            raise ValueError(f"ASD file {path} does not contain (freq, asd) columns")
        freqs = data[:, 0]
        psd = data[:, 1] ** 2
        out[detector] = np.column_stack([freqs, psd])
    return out


def _read_skymap(path: Path, *, event_id: str, gps_time: float):
    import numpy as np
    from astropy.io import fits
    from pesummary.gw.file.skymap import SkyMap

    with fits.open(path) as hdus:
        if len(hdus) < 2:
            raise ValueError(f"Skymap file {path} does not contain a FITS table extension")
        table = hdus[1].data
        header = hdus[1].header
        meta = {
            "nest": str(header.get("ORDERING", "")).upper() == "NESTED",
            "objid": header.get("OBJECT", event_id),
            "gps_time": float(gps_time),
            "creator": header.get("CREATOR", "unofficial bundle"),
            "origin": header.get("ORIGIN", "unknown"),
            "distmean": float(header.get("DISTMEAN", 0.0)),
            "diststd": float(header.get("DISTSTD", 0.0)),
        }
        return SkyMap(np.asarray(table["PROB"], dtype=float), meta_data=meta)
