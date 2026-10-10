"""Binary black holes in a scalar-field cloud (``scalar_cloud`` mode): the IMRPhenomXPHM_Scalar waveform of
Roy, Vicente, Aurrekoetxea, Clough & Ferreira (arXiv:2510.17967) and a best-fit comparison with the vacuum model.

An ultralight scalar of mass m_phi around a binary of total mass M forms a monopole structure of radius
r_phi = M / alpha^2 (alpha = M m_phi, G = c = hbar = 1). In the late inspiral (beta = R / r_phi << 1) the binary
ionizes it, and loses angular momentum at the rate (Eq. 2 and Eq. 6 of the paper)

    dL_phi / dt = 4 pi rho M^3 I(v, alpha, q) / v^4,
    I = sum_{l, m > 0} pi^2 beta^(1/4 + l/2) [2 (1 + q2) / m]^(1/2 - l) [1 + (-1)^l q2^(l - 1)]^2
                       [Y_lm(pi/2, 0) / Gamma(l + 3/2)]^2,

with rho the mean scalar density within r_phi, v = (M Omega)^(1/3), beta = alpha^2 / v^2, q2 = m2 / M, and the
sum up to l = 4. Against the Newtonian GW loss dL_GW/dt = (32/5) eta^2 M v^7, the orbital frequency then evolves
faster by the factor 1 + delta(v), delta = (5 pi / 8) rho M^2 I / (eta^2 v^11). In the stationary-phase
approximation this changes the phase of the (2, 2) mode by

    Delta Psi_22(f) = 2 pi int_f^{f_c} (f' - f) Delta(dt/df)(f') df',   Delta(dt/df) = (dt/df)_GW [1/(1 + delta) - 1],

zero with its derivative at the end of the inspiral f_c (M f_c = 0.018), so that only the inspiral changes, and
of the (l, m) modes by (|m| / 2) Delta Psi_22(2 f / |m|). The masses in delta are source-frame (rho is a physical
density); the kinematics use the detector-frame masses.

The waveform is the IMRPhenomXPHM polarizations of bilby (the catalog model) plus that correction, computed on the
modes of SimInspiralChooseFDModes (J frame, seen at theta_jn) and rotated into the polarization frame of
SimInspiralChooseFDWaveform by an angle fitted on the vacuum modes. With rho = 0 it is exactly IMRPhenomXPHM. The
paper applies the correction in the co-precessing frame; here it is applied per J-frame mode, which is exact for
aligned spins and an approximation for precessing ones.
"""
from __future__ import annotations

from math import gamma, pi

import numpy as np

G_SI = 6.67430e-11
C_SI = 299792458.0
MSUN_SI = 1.988409870698051e30
HBAR_EV_S = 6.582119569e-16
MF_INSPIRAL_END = 0.018            # end of the inspiral of the (2, 2) mode, in units of 1 / M (PhenomX)
L_MAX = 4
_trapz = getattr(np, "trapezoid", None) or np.trapz


def _log(msg: str) -> None:
    print(f"[scalar_cloud] {msg}", flush=True)


# ---------------------------------------------------------------------------
# physics
# ---------------------------------------------------------------------------
def _ylm_equator(l: int, m: int) -> float:
    """Y_lm(theta = pi/2, phi = 0) = N_lm P_l^m(0) (only its square enters, so the phase convention does not)."""
    from math import factorial

    from scipy.special import lpmv

    return float(np.sqrt((2 * l + 1) / (4 * pi) * factorial(l - m) / factorial(l + m)) * lpmv(m, l, 0.0))


_YLM = {(l, m): _ylm_equator(l, m) for l in range(1, L_MAX + 1) for m in range(1, l + 1)}


def torque_coefficients(q2: float) -> list[tuple[float, float]]:
    """(power of beta, coefficient) of the terms of I(v, alpha, q) for q2 = m2 / M (Eq. 6)."""
    out = []
    for (l, m), y in _YLM.items():
        c = (pi ** 2 * (2 * (1 + q2) / m) ** (0.5 - l) * (1 + (-1) ** l * q2 ** (l - 1)) ** 2
             * (y / gamma(l + 1.5)) ** 2)
        if c > 0:
            out.append((0.25 + 0.5 * l, c))
    return out


def torque_integral(v: np.ndarray, alpha: float, q2: float) -> np.ndarray:
    """I(v, alpha, q) of Eq. (6), with beta = alpha^2 / v^2."""
    beta = alpha ** 2 / np.asarray(v, float) ** 2
    return sum(c * beta ** p for p, c in torque_coefficients(q2))


def density_geometric(rho_cgs: float) -> float:
    """A mass density in g/cm^3 in geometric units (1/m^2)."""
    return rho_cgs * 1e3 * G_SI / C_SI ** 2


def scalar_mass_ev(alpha: float, m_source_msun: float) -> float:
    """m_phi (eV) for alpha = G M m_phi / (hbar c^3), M the source-frame total mass."""
    return alpha * HBAR_EV_S / (G_SI * m_source_msun * MSUN_SI / C_SI ** 3)


def cloud_mass_fraction(rho_cgs: float, alpha: float, m_source_msun: float) -> float:
    """rho r_phi^3 / M, the mass of the scalar structure against the binary (the paper's prior keeps it < 0.1)."""
    m_geo = G_SI * m_source_msun * MSUN_SI / C_SI ** 2
    return density_geometric(rho_cgs) * (m_geo / alpha ** 2) ** 3 / m_geo


def dephasing_22(f: np.ndarray, m_det_msun: float, m_src_msun: float, q: float, rho_cgs: float, alpha: float,
                 n_grid: int = 4000) -> np.ndarray:
    """Delta Psi_22 (rad) at the frequencies f (Hz, > 0) of the (2, 2) mode; 0 at and above the end of the
    inspiral. q = m2 / m1 <= 1."""
    f = np.asarray(f, float)
    out = np.zeros_like(f)
    if rho_cgs <= 0:
        return out
    m_s = G_SI * m_det_msun * MSUN_SI / C_SI ** 3                     # detector-frame total mass, seconds
    eta = q / (1 + q) ** 2
    q2 = q / (1 + q)
    rho_m2 = density_geometric(rho_cgs) * (G_SI * m_src_msun * MSUN_SI / C_SI ** 2) ** 2
    f_c = MF_INSPIRAL_END / m_s
    lo = max(min(f[f > 0].min(), f_c) * 0.999, 1e-3) if np.any(f > 0) else f_c
    fg = np.geomspace(lo, f_c, n_grid)
    v = (pi * m_s * fg) ** (1 / 3)
    delta = (5 * pi / 8) * rho_m2 * torque_integral(v, alpha, q2) / (eta ** 2 * v ** 11)
    dtdf_gw = (5 * m_s / (32 * eta * v ** 9)) * (pi * m_s / (3 * v ** 2))
    g = dtdf_gw * (1 / (1 + delta) - 1)                                # Delta(dt/df) <= 0: a faster inspiral
    # Delta Psi(f) = 2 pi int_f^fc (f' - f) g(f') df' = 2 pi [A(f) - f B(f)], A = int f' g, B = int g (from f to fc)
    cum = lambda y: np.concatenate([[0.0], np.cumsum(0.5 * (y[1:] + y[:-1]) * np.diff(fg))])
    a_tot, b_tot = cum(fg * g), cum(g)
    a = a_tot[-1] - a_tot
    b = b_tot[-1] - b_tot
    psi = 2 * pi * (a - fg * b)
    inside = (f > 0) & (f < f_c)
    out[inside] = np.interp(f[inside], fg, psi)
    return out


# ---------------------------------------------------------------------------
# waveform
# ---------------------------------------------------------------------------
_Z_GRID = np.linspace(0, 3, 3001)
_DL_GRID = None


def _source_frame_total_mass(m_det_total: float, luminosity_distance: float) -> float:
    """Detector- to source-frame mass with the redshift of the distance (Planck15, as the GWTC PE files)."""
    global _DL_GRID
    if _DL_GRID is None:
        from astropy.cosmology import Planck15

        _DL_GRID = Planck15.luminosity_distance(_Z_GRID).value
    return m_det_total / (1 + float(np.interp(luminosity_distance, _DL_GRID, _Z_GRID)))


def _mode_sum(modes: dict, theta_jn: float, n_pos: int, factors: dict | None = None):
    """(plus, cross) on the positive frequencies from two-sided modes h_lm, seen at (theta_jn, 0), each mode
    multiplied by factors[(l, m)] (two-sided arrays) if given."""
    import lal

    total = None
    for (l, m), h in modes.items():
        term = h * lal.SpinWeightedSphericalHarmonic(theta_jn, 0.0, -2, l, m)
        if factors is not None:
            term = term * factors[(l, m)]
        total = term if total is None else total + term
    c = len(total) // 2
    pos = total[c:c + n_pos]
    neg = np.conj(total[c::-1][:n_pos])
    return 0.5 * (pos + neg), 0.5j * (pos - neg)


def scalar_cloud_binary_black_hole(frequency_array, mass_1, mass_2, luminosity_distance, a_1, tilt_1, phi_12, a_2,
                                   tilt_2, phi_jl, theta_jn, phase, rho_phi=0.0, alpha_phi=0.1, **kwargs):
    """bilby frequency-domain source model: IMRPhenomXPHM in a scalar cloud of mean density rho_phi (g/cm^3) and
    alpha_phi = M m_phi (source-frame M). rho_phi = 0 is bilby's lal_binary_black_hole, unchanged."""
    import lal
    import lalsimulation as ls
    from bilby.gw.conversion import bilby_to_lalsimulation_spins
    from bilby.gw.source import lal_binary_black_hole

    wf_kw = dict(waveform_approximant="IMRPhenomXPHM", reference_frequency=20.0, minimum_frequency=20.0,
                 maximum_frequency=frequency_array[-1])
    wf_kw.update({k: v for k, v in kwargs.items() if k in wf_kw})
    vac = lal_binary_black_hole(frequency_array, mass_1, mass_2, luminosity_distance, a_1, tilt_1, phi_12, a_2,
                                tilt_2, phi_jl, theta_jn, phase, **wf_kw)
    if rho_phi <= 0 or vac is None:
        return vac
    f_ref, f_min = wf_kw["reference_frequency"], wf_kw["minimum_frequency"]
    df = frequency_array[1] - frequency_array[0]
    q = min(mass_1, mass_2) / max(mass_1, mass_2)
    m_det = mass_1 + mass_2
    m_src = _source_frame_total_mass(m_det, luminosity_distance)
    # the correction vanishes above m f_c / 2 <= 2 f_c (l = |m| <= 4): the modes are needed only up to there
    f_c = MF_INSPIRAL_END / (G_SI * m_det * MSUN_SI / C_SI ** 3)
    f_max = min(frequency_array[-1], 2.05 * f_c)
    iota, s1x, s1y, s1z, s2x, s2y, s2z = bilby_to_lalsimulation_spins(
        theta_jn=theta_jn, phi_jl=phi_jl, tilt_1=tilt_1, tilt_2=tilt_2, phi_12=phi_12, a_1=a_1, a_2=a_2,
        mass_1=mass_1 * MSUN_SI, mass_2=mass_2 * MSUN_SI, reference_frequency=f_ref, phase=phase)
    hlm = ls.SimInspiralChooseFDModes(mass_1 * MSUN_SI, mass_2 * MSUN_SI, s1x, s1y, s1z, s2x, s2y, s2z, df, f_min,
                                      f_max, f_ref, phase, luminosity_distance * 1e6 * lal.PC_SI, iota,
                                      lal.CreateDict(), ls.IMRPhenomXPHM)
    modes = {}
    while hlm is not None:
        modes[(hlm.l, hlm.m)] = hlm.mode.data.data.copy()
        hlm = hlm.next
    n_two = len(next(iter(modes.values())))
    f_two = df * (np.arange(n_two) - (n_two - 1) // 2)            # LAL rounds the band up to a power of two
    n_pos = min((n_two - 1) // 2 + 1, len(frequency_array))
    # the polarization frame of ChooseFDWaveform: a rotation of the J-frame mode sum, fitted on the vacuum modes
    p0, x0 = _mode_sum(modes, theta_jn, n_pos)
    band = (frequency_array[:n_pos] >= f_min) & (frequency_array[:n_pos] <= f_max) & (np.abs(vac["plus"][:n_pos]) > 0)
    a_mat = np.concatenate([np.stack([p0[band], x0[band]], 1), np.stack([x0[band], -p0[band]], 1)])
    rhs = np.concatenate([vac["plus"][:n_pos][band], vac["cross"][:n_pos][band]])
    a_mat = np.concatenate([a_mat.real, a_mat.imag])
    rhs = np.concatenate([rhs.real, rhs.imag])
    (cz, sz), *_ = np.linalg.lstsq(a_mat, rhs, rcond=None)
    af = np.abs(f_two)
    factors = {}
    for (l, m) in modes:
        if m == 0:
            factors[(l, m)] = 0.0
            continue
        dpsi = (abs(m) / 2) * dephasing_22(2 * af / abs(m), m_det, m_src, q, rho_phi, alpha_phi)
        # modes with m > 0 live at f < 0 in LAL's two-sided convention: the phase shift changes sign with f
        factors[(l, m)] = np.exp(-1j * np.sign(f_two) * np.sign(m) * dpsi * SIGN_CONVENTION) - 1.0
    dp, dx = _mode_sum(modes, theta_jn, n_pos, factors)
    plus, cross = vac["plus"].copy(), vac["cross"].copy()
    plus[:n_pos] += cz * dp + sz * dx
    cross[:n_pos] += cz * dx - sz * dp
    return {"plus": plus, "cross": cross}


SIGN_CONVENTION = 1.0            # fixed by the test that the scalar inspiral is shorter (see tests)


# ---------------------------------------------------------------------------
# data: everything the fits need, extracted once into the work directory
# ---------------------------------------------------------------------------
# The work directory holds event.json (settings), psd_<det>.txt (the noise spectra of the PE file),
# catalog_samples.npz (the highest-likelihood catalog samples: starting points and the metric of the fits) and the
# strain: strain_<det>.npy for the event, offsource/<k>_<det>.npy for the background segments.
FIT_KEYS = ["chirp_mass", "mass_ratio", "a_1", "a_2", "tilt_1", "tilt_2", "phi_12", "phi_jl", "theta_jn",
            "luminosity_distance", "phase", "psi", "geocent_time"]
SKY_KEYS = ["ra", "dec"]
N_CATALOG = 2000                   # top catalog samples kept: the covariance that whitens the fit coordinates


def _read_scalar(ds):
    v = ds[()]
    v = v[0] if hasattr(v, "__len__") and not isinstance(v, (bytes, str)) and len(v) else v
    return v.decode() if isinstance(v, bytes) else v


def extract_pe(pe_file: str, label: str, workdir) -> dict:
    """Settings, noise spectra and the top catalog samples of `label`, written into `workdir`."""
    import json
    from pathlib import Path

    import h5py

    wd = Path(workdir)
    wd.mkdir(parents=True, exist_ok=True)
    with h5py.File(pe_file, "r") as h:
        if label not in h:
            raise ValueError(f"{label} not in {pe_file}: {[k for k in h if 'posterior_samples' in h[k]]}")
        g = h[label]
        md = g["meta_data/meta_data"]
        ps = g["posterior_samples"][()]
        settings = dict(
            label=label, pe_file=str(pe_file),
            duration=float(_read_scalar(md["segment_length"])),
            segment_start=float(_read_scalar(md["segment_start"])),
            f_low=float(_read_scalar(md["f_low"])), f_ref=float(_read_scalar(md["f_ref"])),
            pe_sample_rate=float(_read_scalar(md["sample_rate"])),
            detectors=sorted(g["psds"].keys()),
        )
        try:
            other = g["meta_data/other"]
            settings["pe_max_log_likelihood_ratio"] = float(
                _read_scalar(other["log_max_likelihood"]) - _read_scalar(other["log_noise_evidence"]))
        except KeyError:
            pass
        for det in settings["detectors"]:
            np.savetxt(wd / f"psd_{det}.txt", g["psds"][det][()])
    order = np.argsort(ps["log_likelihood"])[::-1][:N_CATALOG]
    cols = FIT_KEYS + SKY_KEYS + [k for k in ("log_likelihood", "chi_eff", "chi_p", "redshift") if k in ps.dtype.names]
    np.savez(wd / "catalog_samples.npz", **{k: np.asarray(ps[k][order], float) for k in cols})
    (wd / "event.json").write_text(json.dumps(settings, indent=1))
    return settings


def _gwosc_coverage(det: str, start: float, end: float) -> list[tuple[float, float]]:
    from gwosc.timeline import get_segments

    return [(float(a), float(b)) for a, b in get_segments(f"{det}_DATA", int(start), int(end) + 1)]


def fetch_strain(workdir, n_offsource: int = 0, sampling_frequency: float = 2048.0, window: float = 3600.0,
                 avoid: float = 64.0) -> int:
    """The event's strain, and n_offsource segments of the same length within +- window s of it, at least `avoid`
    s away and inside the GWOSC data of every detector, resampled and written as text. Returns the number of
    off-source segments written."""
    import json
    from pathlib import Path

    from gwpy.timeseries import TimeSeries

    wd = Path(workdir)
    s = json.loads((wd / "event.json").read_text())
    dur, t0 = s["duration"], s["segment_start"]
    s["sampling_frequency"] = sampling_frequency
    (wd / "event.json").write_text(json.dumps(s, indent=1))
    lo, hi = (t0 - window, t0 + dur + window) if n_offsource else (t0, t0 + dur)
    ok = None
    for det in s["detectors"]:
        cov = _gwosc_coverage(det, lo, hi) if n_offsource else [(lo, hi)]
        ok = cov if ok is None else [(max(a, c), min(b, d)) for a, b in ok for c, d in cov if min(b, d) > max(a, c)]
    starts = []
    if n_offsource:
        k = 1
        while len(starts) < n_offsource and k * dur < window:
            for t in (t0 - avoid - k * dur, t0 + dur + avoid + (k - 1) * dur):     # before, after the event
                if len(starts) < n_offsource and any(a + 8 <= t and t + dur + 8 <= b for a, b in ok):
                    starts.append(t)
            k += 1
    (wd / "offsource").mkdir(exist_ok=True)
    for det in s["detectors"]:
        span = [t0] + starts
        ts = TimeSeries.fetch_open_data(det, min(span) - 8, max(span) + dur + 8, sample_rate=4096, cache=True)
        ts = ts.resample(sampling_frequency)
        for k, t in enumerate([t0] + starts):
            seg = ts.crop(t, t + dur)
            if not np.isfinite(seg.value).all():
                raise ValueError(f"{det}: missing data in [{t}, {t + dur}]")
            out = wd / (f"strain_{det}.npy" if k == 0 else f"offsource/{k - 1:03d}_{det}.npy")
            np.save(out, np.column_stack([seg.times.value, seg.value]))
    (wd / "offsource" / "starts.json").write_text(json.dumps(starts))
    _log(f"strain written: the event and {len(starts)} off-source segments ({sampling_frequency:g} Hz)")
    return len(starts)


class Analysis:
    """The likelihood of one data segment (the event or an off-source one), and the fits on it."""

    def __init__(self, workdir, segment: str | None = None):
        import json
        from pathlib import Path

        import bilby

        bilby.core.utils.logger.setLevel("ERROR")
        self.wd = Path(workdir)
        self.s = json.loads((self.wd / "event.json").read_text())
        self.catalog = dict(np.load(self.wd / "catalog_samples.npz"))
        fs, dur = self.s["sampling_frequency"], self.s["duration"]
        ifos = bilby.gw.detector.InterferometerList([])
        for det in self.s["detectors"]:
            path = self.wd / (f"strain_{det}.npy" if segment is None else f"offsource/{segment}_{det}.npy")
            t, x = np.load(path).T
            ifo = bilby.gw.detector.get_empty_interferometer(det)
            ifo.strain_data.set_from_time_domain_strain(x, sampling_frequency=fs, duration=dur, start_time=t[0])
            ifo.minimum_frequency, ifo.maximum_frequency = self.s["f_low"], fs / 2
            psd = np.loadtxt(self.wd / f"psd_{det}.txt")
            ifo.power_spectral_density = bilby.gw.detector.PowerSpectralDensity(frequency_array=psd[:, 0],
                                                                                psd_array=psd[:, 1])
            ifos.append(ifo)
        self.ifos = ifos
        self.start_time = float(ifos[0].strain_data.start_time)
        self.waveform_generator = bilby.gw.WaveformGenerator(
            duration=dur, sampling_frequency=fs, start_time=self.start_time,
            frequency_domain_source_model=scalar_cloud_binary_black_hole,
            parameter_conversion=bilby.gw.conversion.convert_to_lal_binary_black_hole_parameters,
            waveform_arguments=dict(waveform_approximant="IMRPhenomXPHM", reference_frequency=self.s["f_ref"],
                                    minimum_frequency=self.s["f_low"]))
        self.likelihood = bilby.gw.likelihood.GravitationalWaveTransient(ifos, self.waveform_generator)
        x = np.column_stack([self.catalog[k] for k in FIT_KEYS])
        self.chol = np.linalg.cholesky(np.cov(x.T) + 1e-12 * np.eye(len(FIT_KEYS)))
        self.n_eval = 0

    def inject(self, params: dict, time_shift: float) -> None:
        """Add the signal of `params`, moved by `time_shift` s, to the strain (off-source segments)."""
        p = dict(params, geocent_time=params["geocent_time"] + time_shift)
        for ifo in self.ifos:
            ifo.inject_signal(parameters=p, waveform_generator=self.waveform_generator, raise_error=False)

    def lnl(self, p: dict) -> float:
        self.n_eval += 1
        self.likelihood.parameters.update(p)
        return float(self.likelihood.log_likelihood_ratio())

    def catalog_starts(self) -> list[dict]:
        """The best catalog sample, and the best one in the other mode of the mass ratio when there is one."""
        c, keys = self.catalog, FIT_KEYS + SKY_KEYS
        out = [{k: float(c[k][0]) for k in keys}]
        q0 = c["mass_ratio"][0]
        other = np.where(np.abs(c["mass_ratio"] - q0) > 0.25)[0]
        if len(other):
            out.append({k: float(c[k][other[0]]) for k in keys})
        return out

    # -- fits ---------------------------------------------------------------------------------------------------
    @staticmethod
    def _valid(p: dict) -> bool:
        return (0.05 <= p["mass_ratio"] <= 1 and 0 <= p["a_1"] <= 0.99 and 0 <= p["a_2"] <= 0.99
                and 0 <= p["tilt_1"] <= pi and 0 <= p["tilt_2"] <= pi and 0 <= p["theta_jn"] <= pi
                and p["luminosity_distance"] > 10 and p["chirp_mass"] > 0)

    def fit(self, start: dict, rho: float = 0.0, alpha: float = 0.1, free_cloud: bool = False,
            maxfev: int = 1500, step: float = 0.3) -> tuple[dict, float]:
        """Nelder-Mead in coordinates whitened by the catalog covariance (the binary parameters), from `start`;
        the cloud (log10 rho, alpha) fixed, or free with `free_cloud`. Returns (parameters, ln L ratio)."""
        from scipy.optimize import minimize

        th0 = np.array([start[k] for k in FIT_KEYS])
        n = len(FIT_KEYS) + (2 if free_cloud else 0)
        l0, a0 = (np.log10(rho) if rho > 0 else 3.0), alpha

        def unpack(y):
            p = dict(start)
            p.update(zip(FIT_KEYS, th0 + self.chol @ y[:len(FIT_KEYS)]))
            if free_cloud:
                p["rho_phi"], p["alpha_phi"] = 10 ** (l0 + 0.5 * y[-2]), a0 + 0.05 * y[-1]
            else:
                p["rho_phi"], p["alpha_phi"] = rho, alpha
            return p

        def obj(y):
            p = unpack(y)
            if not self._valid(p) or not (0.005 <= p["alpha_phi"] <= 0.5 * (1 + p["mass_ratio"])):
                return 1e9
            if free_cloud and not (2 <= np.log10(p["rho_phi"]) <= 11):
                return 1e9
            return -self.lnl(p)

        sim = np.vstack([np.zeros(n)] + [step * np.eye(n)[i] for i in range(n)])
        r = minimize(obj, np.zeros(n), method="Nelder-Mead",
                     options=dict(initial_simplex=sim, maxfev=maxfev, xatol=1e-3, fatol=1e-3, adaptive=True))
        return unpack(r.x), -float(r.fun)

    def density_for_dephasing(self, p: dict, alpha: float, target: float) -> float:
        """The density (g/cm^3) for which |Delta Psi_22(f_low)| = target rad, at the binary of `p`."""
        from bilby.gw.conversion import chirp_mass_and_mass_ratio_to_total_mass

        m_det = chirp_mass_and_mass_ratio_to_total_mass(p["chirp_mass"], p["mass_ratio"])
        m_src = _source_frame_total_mass(m_det, p["luminosity_distance"])
        f = np.array([self.s["f_low"]])
        g = lambda lr: abs(dephasing_22(f, m_det, m_src, p["mass_ratio"], 10 ** lr, alpha)[0]) - target
        lo, hi = 0.0, 14.0
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            lo, hi = (mid, hi) if g(mid) < 0 else (lo, mid)
        return 10 ** (0.5 * (lo + hi))


# the fixed procedure that defines the statistic, identical on the event and on the background
ALPHAS = (0.02, 0.04, 0.07, 0.1, 0.15, 0.2, 0.3)
DEPHASING_TARGETS = (0.5, 2.0, 8.0, 30.0)          # |Delta Psi_22(f_low)| in rad: the densities of the profile
CHIRP_MASS_SHIFTS = (0.0, -0.02, -0.04, -0.06)    # starts of the profile fits, relative to the vacuum best fit
VACUUM_MAXFEV, PROFILE_MAXFEV, POLISH_MAXFEV, N_POLISH = 2500, 1200, 3000, 2


def fit_vacuum(an: Analysis, starts: list[dict]) -> dict:
    best = None
    for st in starts:
        p, l = an.fit(dict(st, rho_phi=0.0, alpha_phi=0.1), maxfev=VACUUM_MAXFEV)
        _log(f"vacuum fit from q = {st['mass_ratio']:.2f}: ln L = {l:.2f}")
        if best is None or l > best[1]:
            best = (p, l)
    return dict(params=best[0], log_likelihood_ratio=best[1])


def profile_point(an: Analysis, vacuum: dict, i: int, alpha: float | None = None,
                  target: float | None = None) -> dict:
    """Grid point i of the profile (or the point (alpha, target) given): the cloud fixed, the binary re-fitted
    from the vacuum best fit."""
    alpha = ALPHAS[i // len(DEPHASING_TARGETS)] if alpha is None else alpha
    target = DEPHASING_TARGETS[i % len(DEPHASING_TARGETS)] if target is None else target
    rho = an.density_for_dephasing(vacuum["params"], alpha, target)
    best = None
    # a strong dephasing at low frequency is partly compensated by a lower chirp mass (Roy et al., Fig. 8): the
    # binary is re-fitted from the vacuum best fit and from starts with the chirp mass lowered and the spins
    # in the orbital plane
    for shift in CHIRP_MASS_SHIFTS:
        start = dict(vacuum["params"])
        if shift:
            start["chirp_mass"] *= 1 + shift
            start["tilt_1"] = start["tilt_2"] = pi / 2
        p, l = an.fit(start, rho=rho, alpha=alpha, maxfev=PROFILE_MAXFEV, step=0.5 if shift else 0.3)
        if best is None or l > best[1]:
            best = (p, l, shift)
    return dict(index=i, alpha=alpha, target_dephasing=target, rho=rho, params=best[0],
                log_likelihood_ratio=best[1], chirp_mass_shift_of_start=best[2])


def polish(an: Analysis, points: list[dict], vacuum: dict) -> dict:
    """Free fit (binary and cloud) from the best profile points; the statistic Lambda = 2 (ln L_cloud - ln L_vac)."""
    best = None
    for pt in sorted(points, key=lambda d: -d["log_likelihood_ratio"])[:N_POLISH]:
        p, l = an.fit(pt["params"], rho=pt["rho"], alpha=pt["alpha"], free_cloud=True, maxfev=POLISH_MAXFEV)
        if best is None or l > best[1]:
            best = (p, l)
    if best[1] < vacuum["log_likelihood_ratio"]:
        best = (dict(vacuum["params"]), vacuum["log_likelihood_ratio"])
    lam = 2 * (best[1] - vacuum["log_likelihood_ratio"])
    return dict(params=best[0], log_likelihood_ratio=best[1], Lambda=lam,
                vacuum_log_likelihood_ratio=vacuum["log_likelihood_ratio"])


def statistic(an: Analysis, starts: list[dict]) -> dict:
    """The whole procedure on one data segment: vacuum fit, profile, polish."""
    vac = fit_vacuum(an, starts)
    points = [profile_point(an, vac, i) for i in range(len(ALPHAS) * len(DEPHASING_TARGETS))]
    res = polish(an, points, vac)
    res.update(vacuum=vac, profile=points)
    return res


N_SEEDS, SEEDED_MAXFEV, SEEDED_STEPS = 4, 2000, (0.1, 0.05)


def cloud_seeds(an: Analysis, posterior, n_top: int = 25, n_seeds: int = N_SEEDS) -> list[dict]:
    """Starting points in the cloud mode: the highest-likelihood samples of a nested-sampling posterior of the
    event (time marginalized: the merger time of each is set by a scan of the fixed-time likelihood), the n_seeds
    best at fixed time."""
    keys = FIT_KEYS + SKY_KEYS + ["rho_phi", "alpha_phi"]
    t0 = float(an.catalog["geocent_time"][0])
    out = []
    for _, row in posterior.sort_values("log_likelihood", ascending=False).head(n_top).iterrows():
        p = {k: float(row[k]) for k in keys if k in row}
        best = (-np.inf, t0)
        for grid in (lambda c: c + np.arange(-0.02, 0.02, 2.5e-4), lambda c: c + np.arange(-3e-4, 3e-4, 1e-5)):
            for tc in grid(best[1]):
                l = an.lnl(dict(p, geocent_time=tc))
                if l > best[0]:
                    best = (l, tc)
        p["geocent_time"] = best[1]
        out.append(dict(params=p, log_likelihood_ratio=best[0]))
    return sorted(out, key=lambda d: -d["log_likelihood_ratio"])[:n_seeds]


def seeded_statistic(an: Analysis, vacuum_start: dict, seeds: list[dict], time_shift: float = 0.0) -> dict:
    """The statistic with the cloud fits started in the cloud mode: vacuum fit from vacuum_start, then free fits
    (binary and cloud, Nelder-Mead twice) from each seed moved by time_shift s. Identical on the event
    (time_shift 0) and on the injections."""
    vac = fit_vacuum(an, [vacuum_start])
    best = (dict(vac["params"]), vac["log_likelihood_ratio"])
    fits = []
    for s in seeds:
        p = dict(s["params"], geocent_time=s["params"]["geocent_time"] + time_shift)
        l = an.lnl(p)
        for step in SEEDED_STEPS:
            p, l = an.fit(p, rho=p["rho_phi"], alpha=p["alpha_phi"], free_cloud=True, maxfev=SEEDED_MAXFEV, step=step)
        _log(f"cloud fit from seed: ln L = {l:.2f} (rho {p['rho_phi']:.3g}, alpha {p['alpha_phi']:.3f})")
        fits.append(dict(params=p, log_likelihood_ratio=l))
        if l > best[1]:
            best = (p, l)
    return dict(params=best[0], log_likelihood_ratio=best[1], Lambda=2 * (best[1] - vac["log_likelihood_ratio"]),
                vacuum_log_likelihood_ratio=vac["log_likelihood_ratio"], vacuum=vac, fits=fits)


# ---------------------------------------------------------------------------
# nested sampling: the evidences of vacuum and of a cloud (the paper's Bayes factor)
# ---------------------------------------------------------------------------
RHO_MAX_AGNOSTIC = 1e9             # g/cm^3: the paper's "broad uniform prior" (its upper end is not given)
SAMPLER_SETTINGS = dict(nlive=1000, naccept=60, sample="acceptance-walk", dlogz=0.1)   # as the paper


def superradiance_window(m1_source: float, tau_d_yr: float) -> tuple[float, float]:
    """Bounds on alpha / (1 + q) from tau_SR <= tau_d <= tau_ann (Eq. 9 of the paper's SM)."""
    x = (m1_source / 10.0) * (1e6 / tau_d_yr)
    return 0.02 * x ** (1 / 9), 0.11 * x ** (1 / 15)


def superradiance_density_max(alpha: float, q: float, m1_source: float) -> float:
    """Upper bound on rho (g/cm^3) from the saturation mass of the cloud (Eq. 10 of the paper's SM)."""
    return 1.3e8 * (1 + q) ** 3 / 8 * ((alpha / (1 + q)) / 0.07) ** 7 * (10.0 / m1_source) ** 2


def _cloud_conversion(tau_d_yr):
    def convert(p):
        p = dict(p)
        from bilby.gw.conversion import convert_to_lal_binary_black_hole_parameters

        q_in = p.copy()
        out, _ = convert_to_lal_binary_black_hole_parameters(q_in)
        q = out["mass_ratio"]
        p["mass_1"], p["mass_2"] = out["mass_1"], out["mass_2"]
        if "alpha_phi" in p:
            m_det = np.asarray(out["mass_1"]) + np.asarray(out["mass_2"])
            dl = np.asarray(p["luminosity_distance"], float)
            ratio = np.array([_source_frame_total_mass(1.0, d) for d in np.atleast_1d(dl)]).reshape(np.shape(dl))
            m_src = m_det * ratio
            p["alpha_bound"] = p["alpha_phi"] / (0.5 * (1 + q))
            m_geo = G_SI * m_src * MSUN_SI / C_SI ** 2
            p["cloud_mass_fraction"] = density_geometric(p["rho_phi"]) * (m_geo / p["alpha_phi"] ** 2) ** 3 / m_geo
            if tau_d_yr:
                m1s = np.asarray(out["mass_1"]) * ratio
                lo, hi = superradiance_window(m1s, tau_d_yr)
                p["superradiance_alpha"] = (p["alpha_phi"] / (1 + q) - lo) / (hi - lo)
                p["superradiance_density"] = p["rho_phi"] / superradiance_density_max(p["alpha_phi"], q, m1s)
        return p
    return convert


def sampling_priors(an: Analysis, model: str, tau_d_yr: float | None = None):
    """The binary priors of the GWTC-2.1 LALInference analysis, and those of the cloud: none (vacuum), agnostic
    (rho uniform, alpha <= (1 + q) / 2, cloud mass < 0.1 M) or superradiance (Eqs. 9-10 for tau_d)."""
    import bilby
    from bilby.core.prior import Constraint, Cosine, PowerLaw, Sine, Uniform
    from bilby.gw.prior import UniformInComponentsChirpMass, UniformInComponentsMassRatio

    tc = float(an.catalog["geocent_time"][0])
    pri = bilby.core.prior.PriorDict(conversion_function=_cloud_conversion(tau_d_yr if model == "superradiance"
                                                                           else None))
    pri.update(dict(
        chirp_mass=UniformInComponentsChirpMass(7.83981847070932, 12.811796079336329, name="chirp_mass"),
        mass_ratio=UniformInComponentsMassRatio(0.05, 1.0, name="mass_ratio"),
        mass_1=Constraint(1.0, 1000.0, name="mass_1"), mass_2=Constraint(1.0, 1000.0, name="mass_2"),
        a_1=Uniform(0, 0.99, name="a_1"), a_2=Uniform(0, 0.99, name="a_2"),
        tilt_1=Sine(name="tilt_1"), tilt_2=Sine(name="tilt_2"),
        phi_12=Uniform(0, 2 * pi, name="phi_12", boundary="periodic"),
        phi_jl=Uniform(0, 2 * pi, name="phi_jl", boundary="periodic"),
        luminosity_distance=PowerLaw(2, 10.0, 2000.0, name="luminosity_distance"),
        dec=Cosine(name="dec"), ra=Uniform(0, 2 * pi, name="ra", boundary="periodic"),
        theta_jn=Sine(name="theta_jn"), psi=Uniform(0, pi, name="psi", boundary="periodic"),
        phase=Uniform(0, 2 * pi, name="phase", boundary="periodic"),
        geocent_time=Uniform(tc - 0.1, tc + 0.1, name="geocent_time"),
    ))
    if model == "vacuum":
        pri["rho_phi"], pri["alpha_phi"] = 0.0, 0.1
    else:
        pri["rho_phi"] = Uniform(0.0, RHO_MAX_AGNOSTIC, name="rho_phi")
        pri["alpha_phi"] = Uniform(0.008, 1.0, name="alpha_phi")
        pri["alpha_bound"] = Constraint(0.0, 1.0, name="alpha_bound")
        pri["cloud_mass_fraction"] = Constraint(0.0, 0.1, name="cloud_mass_fraction")
        if model == "superradiance":
            pri["superradiance_alpha"] = Constraint(0.0, 1.0, name="superradiance_alpha")
            pri["superradiance_density"] = Constraint(0.0, 1.0, name="superradiance_density")
    return pri


def sample(workdir, model: str, tau_d_yr: float | None = None, npool: int = 1, nlive: int | None = None) -> None:
    """Nested sampling (bilby, dynesty, the paper's settings) of the event for one model; resumable."""
    import bilby

    an = Analysis(workdir)
    pri = sampling_priors(an, model, tau_d_yr)
    like = bilby.gw.likelihood.GravitationalWaveTransient(an.ifos, an.waveform_generator, priors=pri,
                                                         time_marginalization=True, jitter_time=True)
    label = model if model != "superradiance" else f"superradiance_tau{tau_d_yr:.0e}".replace("+", "")
    kw = dict(SAMPLER_SETTINGS, nlive=nlive or SAMPLER_SETTINGS["nlive"])
    bilby.run_sampler(likelihood=like, priors=pri, sampler="dynesty", outdir=str(an.wd / "sampling"),
                      label=label, npool=npool, resume=True, check_point_delta_t=1800, save="hdf5",
                      conversion_function=None, plot=False, **kw)


# ---------------------------------------------------------------------------
# command line of the jobs
# ---------------------------------------------------------------------------
def _dump(path, obj) -> None:
    import json

    def conv(o):
        if isinstance(o, dict):
            return {k: conv(v) for k, v in o.items()}
        if isinstance(o, (list, tuple)):
            return [conv(v) for v in o]
        if isinstance(o, (np.floating, np.integer)):
            return o.item()
        return o

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(conv(obj), indent=1))


def main(argv=None) -> int:
    import argparse
    import json
    import time
    from pathlib import Path

    ap = argparse.ArgumentParser(prog="scalar_cloud.py")
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name in ("vacuum", "profile", "polish", "background", "seeded"):
        sp = sub.add_parser(name)
        sp.add_argument("--workdir", required=True)
        sp.add_argument("--task", type=int, default=0)
    sub.choices["seeded"].add_argument("--inject", choices=("vacuum", "cloud"), default="vacuum",
                                       help="signal of the injections: the vacuum best fit (background) or the "
                                            "event's cloud best fit (sensitivity), output in seeded_cloud/")
    sp = sub.add_parser("seeds", help="cloud-mode starting points from a sampling result of the event")
    sp.add_argument("--workdir", required=True)
    sp.add_argument("--result", required=True)
    sp = sub.add_parser("extra", help="profile points beyond the grid of the statistic (for the figures)")
    sp.add_argument("--workdir", required=True)
    sp.add_argument("--alpha", type=float, required=True)
    sp.add_argument("--target", type=float, required=True, help="|Delta Psi_22(f_low)| in rad")
    sp = sub.add_parser("sample")
    sp.add_argument("--workdir", required=True)
    sp.add_argument("--model", choices=("vacuum", "agnostic", "superradiance"), required=True)
    sp.add_argument("--tau-d", type=float, default=1e6, help="superradiance: years between BH formation and merger")
    sp.add_argument("--npool", type=int, default=1)
    sp.add_argument("--nlive", type=int, default=None)
    args = ap.parse_args(argv)
    wd = Path(args.workdir)
    t = time.time()
    if args.cmd == "vacuum":
        an = Analysis(wd)
        _dump(wd / "vacuum.json", fit_vacuum(an, an.catalog_starts()))
    elif args.cmd == "profile":
        an = Analysis(wd)
        vac = json.loads((wd / "vacuum.json").read_text())
        _dump(wd / "profile" / f"point_{args.task:03d}.json", profile_point(an, vac, args.task))
    elif args.cmd == "polish":
        an = Analysis(wd)
        vac = json.loads((wd / "vacuum.json").read_text())
        points = [json.loads(p.read_text()) for p in sorted((wd / "profile").glob("point_*.json"))]
        res = polish(an, points, vac)
        res.update(vacuum=vac, profile=points)
        _dump(wd / "scalar.json", res)
        _log(f"Lambda = {res['Lambda']:.2f} (ln L {res['vacuum_log_likelihood_ratio']:.2f} -> "
             f"{res['log_likelihood_ratio']:.2f})")
    elif args.cmd == "extra":
        an = Analysis(wd)
        vac = json.loads((wd / "vacuum.json").read_text())
        res = profile_point(an, vac, -1, alpha=args.alpha, target=args.target)
        _dump(wd / "profile_extra" / f"point_a{args.alpha:g}_t{args.target:g}.json", res)
    elif args.cmd == "sample":
        sample(wd, args.model, args.tau_d if args.model == "superradiance" else None, npool=args.npool,
               nlive=args.nlive)
    elif args.cmd == "seeds":
        import bilby

        an = Analysis(wd)
        post = bilby.core.result.read_in_result(args.result).posterior
        _dump(wd / "cloud_seeds.json", cloud_seeds(an, post))
    elif args.cmd == "seeded":                     # task -1: the event; k >= 0: off-source injection k
        seeds = json.loads((wd / "cloud_seeds.json").read_text())
        vac = json.loads((wd / "vacuum.json").read_text())
        if args.task < 0:
            an, shift, st, name = Analysis(wd), 0.0, vac["params"], "event"
        else:
            name = f"inj_{args.task:03d}"
            starts = json.loads((wd / "offsource" / "starts.json").read_text())
            an = Analysis(wd, segment=f"{args.task:03d}")
            inj = dict(vac["params"], rho_phi=0.0)
            if args.inject == "cloud":
                inj = json.loads((wd / "seeded" / "event.json").read_text())["params"]
            shift = starts[args.task] - an.s["segment_start"]
            an.inject(inj, shift)
            st = dict(inj, geocent_time=inj["geocent_time"] + shift)
        res = seeded_statistic(an, st, seeds, shift)
        res.update(segment=name, time_shift=shift, injection_model=args.inject)
        _dump(wd / ("seeded" if args.inject == "vacuum" else "seeded_cloud") / f"{name}.json", res)
        _log(f"seeded {name}: Lambda = {res['Lambda']:.2f}")
    elif args.cmd == "background":
        k = f"{args.task:03d}"
        starts = json.loads((wd / "offsource" / "starts.json").read_text())
        an = Analysis(wd, segment=k)
        vac = json.loads((wd / "vacuum.json").read_text())
        inj = dict(vac["params"], rho_phi=0.0)
        shift = starts[args.task] - an.s["segment_start"]
        an.inject(inj, shift)
        st = dict(inj, geocent_time=inj["geocent_time"] + shift)
        res = statistic(an, [st])
        res.update(segment=k, segment_start=starts[args.task], injection=st)
        _dump(wd / "background" / f"inj_{k}.json", res)
        _log(f"background {k}: Lambda = {res['Lambda']:.2f}")
    _log(f"{args.cmd} done in {time.time() - t:.0f} s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
