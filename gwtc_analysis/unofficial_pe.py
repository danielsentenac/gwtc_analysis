from __future__ import annotations

import hashlib
import json
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
class PublicSource:
    """A public file the bundle is built from, downloaded when missing.

    If ``member`` is set, ``url`` is a .tar.gz archive and ``path`` is the
    local copy of that archive member.
    """
    url: str
    path: Path
    member: str | None = None


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
    calibration_paths: tuple[tuple[str, Path], ...] = ()
    downloads: tuple[PublicSource, ...] = ()
    # Fit (geocent_time, phase, psi) of the maxL sample to the GWOSC strain, for
    # samples released without them (the overlay would otherwise be incoherent).
    fit_extrinsic: bool = False


# Bumped when the bundle layout or the build procedure changes, so cached
# bundles built by an older recipe are rebuilt.
_BUNDLE_RECIPE_VERSION = 2

_GWCACHE = Path.home() / ".gwcache"
_DCC = "https://dcc.ligo.org/public"

# GW170817 from public GWTC-1 releases only:
#   samples      LIGO-P1800370 (GWTC-1 PE sample release)
#   PSDs         LIGO-P1900011 (GWTC-1 PSD release)
#   calibration  LIGO-P1900040 (GWTC-1 calibration uncertainty envelopes)
#   skymap       LIGO-P1800381 (GWTC-1 skymap release)
# The public samples have no psi/phase/coalescence time/likelihood; they are
# drawn from their priors, and fitted to the strain for the maxL sample.
_GW170817_CALENV = _GWCACHE / "GWTC1_GW170817_CalEnv"
GW170817_SPEC = UnofficialPEBundleSpec(
    event_id="GW170817",
    raw_samples_path=_GWCACHE / "GW170817_GWTC-1.hdf5",
    psd_path=_GWCACHE / "GWTC1_GW170817_PSDs.dat",
    skymap_path=_GWCACHE / "GW170817_skymap.fits.gz",
    output_filename="Unofficial-GWTC1-GW170817_PEDataRelease.h5",
    analyses=(
        UnofficialPEAnalysisSpec(
            dataset_name="IMRPhenomPv2NRT_lowSpin_posterior",
            label="C02:IMRPhenomPv2_NRTidal-LowSpin",
            approximant="IMRPhenomPv2_NRTidal",
        ),
        UnofficialPEAnalysisSpec(
            dataset_name="IMRPhenomPv2NRT_highSpin_posterior",
            label="C02:IMRPhenomPv2_NRTidal-HighSpin",
            approximant="IMRPhenomPv2_NRTidal",
        ),
    ),
    gps_time=1187008882.429464,
    calibration_paths=tuple(
        (det, _GW170817_CALENV / f"GWTC1_GW170817_{det[0]}_CalEnv.txt") for det in ("H1", "L1", "V1")
    ),
    downloads=(
        PublicSource(f"{_DCC}/0157/P1800370/005/GW170817_GWTC-1.hdf5", _GWCACHE / "GW170817_GWTC-1.hdf5"),
        PublicSource(f"{_DCC}/0158/P1900011/001/GWTC1_GW170817_PSDs.dat", _GWCACHE / "GWTC1_GW170817_PSDs.dat"),
        PublicSource(f"{_DCC}/0157/P1800381/007/GW170817_skymap.fits.gz", _GWCACHE / "GW170817_skymap.fits.gz"),
        *(
            PublicSource(
                f"{_DCC}/0158/P1900040/001/GWTC1_GW170817_CalEnv.tar.gz",
                _GW170817_CALENV / f"GWTC1_GW170817_{d}_CalEnv.txt",
                member=f"GWTC1_GW170817_CalEnv/GWTC1_GW170817_{d}_CalEnv.txt",
            )
            for d in ("H", "L", "V")
        ),
    ),
    fit_extrinsic=True,
)


UNOFFICIAL_PE_BUNDLES: dict[str, UnofficialPEBundleSpec] = {
    GW170817_SPEC.event_id: GW170817_SPEC,
}


def list_unofficial_pe_specs() -> list[str]:
    return sorted(UNOFFICIAL_PE_BUNDLES.keys())


def get_unofficial_pe_spec(src_name: str) -> UnofficialPEBundleSpec | None:
    return UNOFFICIAL_PE_BUNDLES.get(str(src_name).strip())


def _download_public_sources(spec: UnofficialPEBundleSpec, log_cb: Callable[[str], None] | None) -> None:
    """Download the spec's public source files that are not on disk yet."""
    import io
    import tarfile

    import requests

    archives: dict[str, bytes] = {}
    for src in spec.downloads:
        path = src.path.expanduser()
        if path.exists() and path.stat().st_size > 0:
            continue
        try:
            if src.url not in archives:
                if log_cb:
                    log_cb(f"ℹ️ [DOWNLOAD] {src.url}")
                r = requests.get(src.url, timeout=(10, 600))
                r.raise_for_status()
                archives[src.url] = r.content
            data = archives[src.url]
            if src.member is not None:
                with tarfile.open(fileobj=io.BytesIO(data), mode="r:gz") as tf:
                    member = next(m for m in tf.getmembers() if m.name.lstrip("./") == src.member)
                    data = tf.extractfile(member).read()
            path.parent.mkdir(parents=True, exist_ok=True)
            tmp = path.with_name(path.name + ".part")
            tmp.write_bytes(data)
            tmp.replace(path)
        except Exception as e:  # reported as a missing source by the caller
            if log_cb:
                log_cb(f"⚠️ [WARN] Could not download {src.url}: {type(e).__name__}: {e}")


def _recipe_fingerprint(spec: UnofficialPEBundleSpec) -> str:
    return hashlib.sha256(f"{_BUNDLE_RECIPE_VERSION}|{spec!r}".encode()).hexdigest()[:16]


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

    _download_public_sources(spec, log_cb)

    sources: list[Path] = [spec.skymap_path.expanduser()]
    if spec.asd_paths:
        sources.extend(Path(p).expanduser() for _, p in spec.asd_paths)
    elif spec.psd_path is not None:
        sources.append(spec.psd_path.expanduser())
    sources.extend(Path(p).expanduser() for _, p in spec.calibration_paths)

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
    recipe_file = target.with_name(target.name + ".recipe.json")
    fingerprint = _recipe_fingerprint(spec)

    if not force_rebuild and target.exists() and target.stat().st_size > 0:
        try:
            cached_fp = json.loads(recipe_file.read_text(encoding="utf-8")).get("fingerprint")
        except (OSError, ValueError):
            cached_fp = None  # bundle from an older recipe: rebuild it
        target_mtime = target.stat().st_mtime
        if cached_fp == fingerprint and all(p.stat().st_mtime <= target_mtime for p in sources):
            if log_cb:
                log_cb(f"ℹ️ [CACHE] Using unofficial PE bundle: {target}")
            return target

    if log_cb:
        if target.exists():
            log_cb(f"ℹ️ [BUILD] Rebuilding unofficial PE bundle for {spec.event_id}: {target}")
        log_cb(f"ℹ️ [BUILD] Building unofficial PE bundle for {spec.event_id} from local cache files")

    _write_unofficial_pesummary_bundle(spec, target, log_cb=log_cb)
    recipe_file.write_text(
        json.dumps({"fingerprint": fingerprint, "sources": [str(p) for p in sources]}, indent=2),
        encoding="utf-8",
    )

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
                        f"{analysis.dataset_name!r} — psi/phase/time from priors, log_likelihood synthetic"
                        + (" (maxL sample fitted to the strain below)." if spec.fit_extrinsic else ".")
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
    calibration = {det: np.genfromtxt(Path(p).expanduser()) for det, p in spec.calibration_paths}

    if spec.fit_extrinsic:
        for analysis in spec.analyses:
            if analysis.lalinference_samples_path is not None:
                continue  # full LALInference samples already carry the real values
            try:
                _fit_maxl_extrinsics(
                    samples_by_label[analysis.label],
                    approximant=analysis.approximant,
                    psds=psds,
                    gps_time=spec.gps_time,
                    label=analysis.label,
                    log_cb=log_cb,
                )
            except Exception as e:
                if log_cb:
                    log_cb(
                        f"⚠️ [WARN] {analysis.label}: could not fit time/phase/polarization to the strain "
                        f"({type(e).__name__}: {e}); the maxL waveform overlay will not be coherent."
                    )

    labels = [analysis.label for analysis in spec.analyses]
    approximant = {analysis.label: analysis.approximant for analysis in spec.analyses}
    file_versions = {label: "GWTC-1 unofficial bundle" for label in labels}
    file_kwargs = {label: {"sampler": {}, "meta_data": {}} for label in labels}
    psd_by_label = {label: psds for label in labels}
    skymap_by_label = {label: skymap for label in labels}
    calibration_kwargs = {"calibration": {label: calibration for label in labels}} if calibration else {}

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
        **calibration_kwargs,
    )


def _fit_maxl_extrinsics(
    samples,
    *,
    approximant: str,
    psds: dict,
    gps_time: float,
    label: str,
    log_cb: Callable[[str], None] | None = None,
    f_low: float = 23.0,
    f_ref: float = 20.0,
    duration: int = 256,
    sample_rate: int = 4096,
) -> dict[str, float]:
    """Fit geocent_time, phase and psi of the maxL sample to the GWOSC strain.

    Public sample releases (e.g. GWTC-1) omit the coalescence time, phase and
    polarization, which the builder then draws from their priors: the maxL
    waveform used for the strain overlay is then incoherent with the data (or
    even sign-flipped). Keeping the sample's intrinsic parameters and sky
    position, this maximizes the coherent network log-likelihood ratio

        ln L = |sum_d F+_d Z+_d(t_d) + Fx_d Zx_d(t_d)| - 1/2 sum_d <h_d|h_d>

    on a (geocent time, psi) grid, where Z are the complex correlations of the
    data with h+ and hx; the phase maximization is analytic. The fitted values
    replace psi, phase, geocent_time and the detector times of that sample only.
    ``f_ref`` must match the one used for the overlay (20 Hz without a config).
    """
    import numpy as np
    from gwpy.timeseries import TimeSeries
    from pesummary.gw.waveform import fd_waveform
    from pycbc.detector import Detector

    dets = [d for d in ("H1", "L1", "V1") if d in psds]
    i = int(np.argmax(np.asarray(samples["log_likelihood"], dtype=float)))
    row = {k: np.array([float(samples[k][i])]) for k in samples.keys()}
    ra, dec = float(row["ra"][0]), float(row["dec"][0])

    start = int(gps_time) - duration + 16
    n = duration * sample_rate
    df = 1.0 / duration
    freqs = np.arange(n // 2 + 1) * df
    band = (freqs >= f_low) & (freqs <= sample_rate / 2)
    taper = np.ones(n)
    nt = sample_rate // 2
    taper[:nt] = 0.5 * (1 - np.cos(np.pi * np.arange(nt) / nt))
    taper[-nt:] = taper[:nt][::-1]

    dtilde, invpsd = {}, {}
    for det in dets:
        ts = TimeSeries.fetch_open_data(det, start - 8, start + duration + 8, cache=True)
        if int(round(ts.sample_rate.value)) != sample_rate:
            ts = ts.resample(sample_rate)
        x = ts.crop(start, start + duration).value[:n]
        dtilde[det] = np.fft.rfft(x * taper) / sample_rate
        psd = np.asarray(psds[det], dtype=float)
        s = np.interp(freqs, psd[:, 0], psd[:, 1], left=np.inf, right=np.inf)
        invpsd[det] = np.where(band & np.isfinite(s) & (s > 0), 1.0 / s, 0.0)

    def polarizations(phase):
        r = dict(row)
        r["phase"] = np.array([phase])
        h = fd_waveform(r, approximant, df, f_low, sample_rate / 2, f_ref=f_ref)
        hp = np.zeros(len(freqs), complex)
        hc = np.zeros(len(freqs), complex)
        m = min(len(freqs), len(h["h_plus"]))
        hp[:m] = h["h_plus"].value[:m]
        hc[:m] = h["h_cross"].value[:m]
        return hp, hc

    def inner(a, b, det):  # <a|b> = 4 Re sum(a b* / S) df
        return 4 * df * float(np.real(np.sum(a * np.conj(b) * invpsd[det])))

    up = 4  # correlation time resolution: 1 / (up * sample_rate)

    def correlation(h, det):
        full = np.zeros(n * up, complex)
        full[: len(freqs)] = 4 * df * dtilde[det] * np.conj(h) * invpsd[det]
        return np.fft.ifft(full) * (n * up)

    phase0 = float(row["phase"][0])
    hp, hc = polarizations(phase0)
    # How the overall frequency-domain phase moves with 'phase' (2 for the 22 mode).
    hp_shift, _ = polarizations(phase0 + 0.1)
    weight = np.abs(hp) ** 2 * invpsd[dets[0]]
    kfac = float(np.angle(np.sum(hp_shift * np.conj(hp) * weight)) / 0.1)
    if abs(kfac) < 0.5:
        raise RuntimeError(f"waveform phase does not respond to 'phase' (d arg h / d phase = {kfac:.2f})")

    corr = {d: (correlation(hp, d), correlation(hc, d)) for d in dets}
    norms = {d: (inner(hp, hp, d), inner(hc, hc, d), inner(hp, hc, d)) for d in dets}
    detectors = {d: Detector(d) for d in dets}
    delays = {d: detectors[d].time_delay_from_earth_center(ra, dec, gps_time) for d in dets}
    dt_up = 1.0 / (sample_rate * up)
    tc_grid = gps_time + np.arange(-0.1, 0.1, dt_up)

    best = (-np.inf, gps_time, 0.0, 0.0)
    for psi in np.linspace(0.0, np.pi, 181)[:-1]:
        coh = np.zeros(len(tc_grid), complex)
        hh = 0.0
        for d in dets:
            fp, fx = detectors[d].antenna_pattern(ra, dec, psi, gps_time)
            idx = np.round((tc_grid + delays[d] - start) / dt_up).astype(int)
            coh += fp * corr[d][0][idx] + fx * corr[d][1][idx]
            a, b, c = norms[d]
            hh += fp * fp * a + fx * fx * b + 2 * fp * fx * c
        lnl = np.abs(coh) - 0.5 * hh
        j = int(np.argmax(lnl))
        if lnl[j] > best[0]:
            best = (float(lnl[j]), float(tc_grid[j]), float(psi), float(np.angle(coh[j])))
    _, tc, psi, alpha = best

    def network(phase):
        hp_, hc_ = polarizations(phase)
        lnl, snrs = 0.0, {}
        for d in dets:
            fp, fx = detectors[d].antenna_pattern(ra, dec, psi, tc)
            shift = tc + detectors[d].time_delay_from_earth_center(ra, dec, tc) - start
            h = (fp * hp_ + fx * hc_) * np.exp(-2j * np.pi * freqs * shift)
            dh, hh_ = inner(dtilde[d], h, d), inner(h, h, d)
            lnl += dh - 0.5 * hh_
            snrs[d] = dh / np.sqrt(hh_) if hh_ > 0 else 0.0
        return lnl, snrs

    # 'phase' is defined modulo 2 pi / |kfac|: evaluate each branch with the full waveform.
    branches = [
        (phase0 + alpha / kfac + m * 2 * np.pi / abs(kfac)) % (2 * np.pi)
        for m in range(max(1, int(round(abs(kfac)))))
    ]
    lnl, snrs, phase = max((network(p) + (p,) for p in branches), key=lambda x: x[0])

    samples["geocent_time"][i] = tc
    samples["phase"][i] = phase
    samples["psi"][i] = psi
    for d in ("H1", "L1", "V1"):
        key = f"{d}_time"
        if key in samples:
            samples[key][i] = tc + Detector(d).time_delay_from_earth_center(ra, dec, tc)

    if log_cb:
        per_det = ", ".join(f"{d} {v:+.1f}" for d, v in snrs.items())
        log_cb(
            f"ℹ️ [FIT] {label}: maxL sample fitted to GWOSC strain: geocent_time={tc:.4f}, "
            f"phase={phase:.3f}, psi={psi:.3f}; matched-filter SNR {per_det} "
            f"(network {np.sqrt(max(2 * lnl, 0.0)):.1f})"
        )
    return {"geocent_time": tc, "phase": phase, "psi": psi, "log_likelihood_ratio": lnl, **snrs}


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
