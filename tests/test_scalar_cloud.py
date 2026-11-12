"""scalar_cloud: the torque of Eq. (6), the dephasing, and the IMRPhenomXPHM_Scalar waveform."""
from __future__ import annotations

import numpy as np
import pytest

from gwtc_analysis import scalar_cloud as sc

F = np.arange(0, 1024 + 1 / 16, 1 / 16)
WF = dict(waveform_approximant="IMRPhenomXPHM", reference_frequency=20.0, minimum_frequency=20.0,
          maximum_frequency=1024.0)
ALIGNED = dict(mass_1=14.0, mass_2=9.5, luminosity_distance=850.0, a_1=0.3, tilt_1=0.0, phi_12=0.0, a_2=0.2,
               tilt_2=0.0, phi_jl=0.0, theta_jn=0.0, phase=0.4)
PRECESSING = dict(ALIGNED, tilt_1=0.9, tilt_2=1.2, phi_12=1.0, phi_jl=0.5, theta_jn=1.0)


def test_torque_terms():
    """Y_22(pi/2, 0)^2 = 15 / (32 pi); the dipole and the odd l + m terms vanish; l = 2 dominates."""
    assert sc._YLM[(2, 2)] ** 2 == pytest.approx(15 / (32 * np.pi))
    terms = sc.torque_coefficients(0.375)
    powers = [p for p, _ in terms]
    assert powers == [1.25, 1.75, 1.75, 2.25, 2.25]               # (2,2), (3,1), (3,3), (4,2), (4,4)
    assert terms[0][1] > 10 * sum(c for _, c in terms[1:])
    v = np.array([0.2, 0.3])
    assert sc.torque_integral(v, 0.1, 0.375) == pytest.approx(
        sum(c * (0.01 / v ** 2) ** p for p, c in terms))


def test_units():
    """alpha = 0.1 for M = 20.6 Msun is m_phi = 6.5e-13 eV; the cloud mass fraction rhobar r_phi^3 / M."""
    assert sc.scalar_mass_ev(0.1, 20.6) == pytest.approx(6.49e-13, rel=0.01)
    assert sc.cloud_mass_fraction(1e8, 0.1, 20.6) == pytest.approx(0.0687, rel=0.01)
    assert sc.cloud_mass_fraction(2e8, 0.1, 20.6) == pytest.approx(2 * sc.cloud_mass_fraction(1e8, 0.1, 20.6))


def test_dephasing():
    """Zero without a cloud and above the end of the inspiral; negative (a faster inspiral), growing towards low
    frequencies, linear in the density while the torque is small against the GW loss."""
    f = np.array([20.0, 30.0, 50.0, 100.0, 149.0, 200.0])
    args = (23.5, 23.5 / 1.18, 9.5 / 14)
    assert not sc.dephasing_22(f, *args, 0.0, 0.1).any()
    d = sc.dephasing_22(f, *args, 1e6, 0.1)
    assert (d[:-1] < 0).all() and d[-1] == 0.0
    assert (np.diff(d[:-1]) > 0).all()
    assert sc.dephasing_22(f, *args, 1e4, 0.1) == pytest.approx(sc.dephasing_22(f, *args, 1e3, 0.1) * 10, rel=1e-3)


def _duration(h, f1=25.0, f2=100.0):
    """Time from f1 to f2 of the inspiral, from the frequency-domain phase: t(f) = -dphi/df / 2 pi."""
    k = (F >= 21) & (F <= 300)
    ph = np.unwrap(np.angle(h["plus"][k] + 1j * h["cross"][k]))
    t = -np.gradient(ph, F[k]) / (2 * np.pi)
    return np.interp(f2, F[k], t) - np.interp(f1, F[k], t)


def test_waveform_vacuum_limit_and_faster_inspiral():
    """rhobar = 0 is IMRPhenomXPHM exactly (aligned and precessing); the cloud shortens the inspiral by the time the
    dephasing predicts."""
    from bilby.gw.source import lal_binary_black_hole

    for base in (ALIGNED, PRECESSING):
        vac = lal_binary_black_hole(F, **base, **WF)
        for rhobar in (0.0, 1e-30):
            h = sc.scalar_cloud_binary_black_hole(F, **base, rhobar_phi=rhobar, alpha_phi=0.1)
            assert np.allclose(h["plus"], vac["plus"], rtol=0, atol=1e-35)
    vac = lal_binary_black_hole(F, **ALIGNED, **WF)
    h = sc.scalar_cloud_binary_black_hole(F, **ALIGNED, rhobar_phi=1e6, alpha_phi=0.1)
    ff = np.linspace(25, 100, 4000)
    d = sc.dephasing_22(ff, 23.5, sc._source_frame_total_mass(23.5, 850.0), 9.5 / 14, 1e6, 0.1)
    expected = (np.gradient(d, ff)[-1] - np.gradient(d, ff)[0]) / (2 * np.pi)        # change of t(100) - t(25)
    assert _duration(h) - _duration(vac) == pytest.approx(expected, rel=0.05)
    assert _duration(h) < _duration(vac)


def test_superradiance_bounds():
    """Eq. (9): alpha / (1 + q) in [0.02, 0.11] for m1 = 10 Msun and tau_d = 1 Myr, narrower for longer delays;
    Eq. (10): 1.3e8 g/cm^3 at alpha / (1 + q) = 0.07, q = 1, m1 = 10 Msun."""
    assert sc.superradiance_window(10.0, 1e6) == pytest.approx((0.02, 0.11))
    lo8, hi8 = sc.superradiance_window(10.0, 1e8)
    assert lo8 < 0.02 and hi8 < 0.11 and lo8 < hi8
    assert sc.superradiance_density_max(0.14, 1.0, 10.0) == pytest.approx(1.3e8)


def test_background_p_value():
    from gwtc_analysis.scalar_cloud_mode import background_p_value

    lams = np.arange(10.0)
    assert background_p_value(lams, 8.5) == (pytest.approx(2 / 11), 1)
    assert background_p_value(lams, 100.0) == (pytest.approx(1 / 11), 0)


def test_slurm_chains(tmp_path):
    """fit then background in one chain, the sampling models as one array job in another."""
    from gwtc_analysis.scalar_cloud_mode import write_slurm_chains

    scripts = write_slurm_chains(tmp_path, "/opt/env/bin/python", ["fit", "background", "sample", "seeded"], 20,
                                 ["vacuum", "agnostic"], 1e6, 16, ["--partition=htc"])
    assert [s.parent.name for s in scripts] == ["slurm", "slurm_sample", "slurm_seeded"]
    names = sorted(p.name for p in (tmp_path / "slurm").glob("0*.sh"))
    assert names == ["00_vacuum.sh", "01_profile.sh", "02_polish.sh", "03_background.sh"]
    prof = (tmp_path / "slurm" / "01_profile.sh").read_text()
    assert f"--array=0-{len(sc.ALPHAS) * len(sc.DEPHASING_TARGETS) - 1}" in prof
    assert "--array=0-19" in (tmp_path / "slurm" / "03_background.sh").read_text()
    samp = (tmp_path / "slurm_sample" / "00_sample.sh").read_text()
    assert "--array=0-1" in samp and "MODELS=(vacuum agnostic)" in samp and "--cpus-per-task=16" in samp
    seeded = (tmp_path / "slurm_seeded" / "01_seeded.sh").read_text()
    assert "--array=0-20" in seeded and "--task $((SLURM_ARRAY_TASK_ID-1))" in seeded
    assert "superradiance_tau1e06_result.hdf5" in (tmp_path / "slurm_seeded" / "00_seeds.sh").read_text()


# ---------------------------------------------------------------------------
# helpers of the matched-filter plot
# ---------------------------------------------------------------------------
class _FakeIfo:
    """The few attributes of a bilby Interferometer the plot helpers read."""

    def __init__(self, fs=1024.0, dur=8.0, psd=1.0):
        self.time_array = np.arange(int(fs * dur)) / fs
        self.frequency_array = np.fft.rfftfreq(len(self.time_array), 1 / fs)
        self.power_spectral_density_array = np.full(len(self.frequency_array), psd)
        self.frequency_mask = np.ones(len(self.frequency_array), bool)
        self.frequency_domain_strain = np.zeros(len(self.frequency_array), complex)


def test_noise_sigma_ignores_the_signal():
    from gwtc_analysis import scalar_cloud_mode as m

    rng = np.random.default_rng(1)
    t = np.arange(16 * 1024) / 1024.0
    x = rng.normal(0, 1, len(t))
    x[(t > 10) & (t < 12)] += 50 * np.sin(2 * np.pi * 40 * t[(t > 10) & (t < 12)])     # a loud signal
    assert m._noise_sigma(x, t, t_merger=12.0, duration=3.0) == pytest.approx(1.0, rel=0.05)


def test_frequency_track_of_a_chirp():
    """A linear chirp from 30 to 90 Hz: the track follows f(t) = 30 + 15 t away from the edges."""
    from gwtc_analysis import scalar_cloud_mode as m

    ifo = _FakeIfo()
    t = ifo.time_array
    x = np.where((t > 2) & (t < 6), np.sin(2 * np.pi * (30 * (t - 2) + 7.5 * (t - 2) ** 2)), 0.0)
    tr = m._frequency_track(ifo, np.fft.rfft(x), f_lo=10.0)
    for tx in (3.0, 4.0, 5.0):
        assert tr[np.argmin(abs(t - tx))] == pytest.approx(30 + 15 * (tx - 2), abs=1.5)


def test_snr_band_and_band_snr():
    """White noise, flat |h|: the 99% frequency and the band SNR follow from the bin counts."""
    from gwtc_analysis import scalar_cloud_mode as m

    ifo = _FakeIfo(psd=4.0)
    f = ifo.frequency_array
    h = np.where((f >= 20) & (f <= 220), 1.0, 0.0).astype(complex)
    assert m._snr_band(ifo, h, 20.0, 0.99) == pytest.approx(218.0, abs=1.0)
    df = f[1] - f[0]
    n = ((f >= 50) & (f <= 100)).sum()
    assert m._band_snr(ifo, h, 50.0, 100.0) == pytest.approx(np.sqrt(4 * df * n / 4.0))


def test_matched_filter_peaks_at_the_shift_with_the_template_norm():
    """Data = template moved 25 ms later, no noise: |z| peaks at +25 ms with the value sqrt(<h, h>)."""
    from gwtc_analysis import scalar_cloud_mode as m

    ifo = _FakeIfo(psd=2.0)
    f = ifo.frequency_array
    h = np.where((f >= 20) & (f <= 300), (f / 100.0) ** (-7 / 6), 0.0) * np.exp(-2j * np.pi * f * 4.0)
    ifo.frequency_domain_strain = h * np.exp(-2j * np.pi * f * 0.025)
    lags, z = m._matched_filter(ifo, h, 20.0, 300.0)
    norm = np.sqrt(4 * (f[1] - f[0]) * np.sum(np.abs(h) ** 2 / 2.0))
    lag, rho = m._peak(lags, z, 0.1)
    assert lag == pytest.approx(0.025, abs=1.5 / 1024) and rho == pytest.approx(norm, rel=0.02)
    assert abs(z[np.argmin(abs(lags))]) < rho


def test_gwpy_import_keeps_the_matplotlib_axes():
    """gwpy replaces the default Axes class (thick legend lines in every later plot): the helper restores it."""
    pytest.importorskip("gwpy")
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.axes import Axes

    from gwtc_analysis import scalar_cloud as sc

    handlelength = matplotlib.rcParams["legend.handlelength"]
    sc.gwpy_timeseries()
    fig, ax = plt.subplots()
    assert type(ax) is Axes and matplotlib.rcParams["legend.handlelength"] == handlelength
    plt.close(fig)


@pytest.mark.parametrize("names", [("segment_length", "segment_start", "sample_rate"),
                                   ("duration", "start_time", "sampling_frequency")])      # GWTC-3 files: bilby's
def test_extract_pe_reads_both_namings_of_the_settings(tmp_path, names):
    import h5py

    from gwtc_analysis import scalar_cloud as sc

    cols = sc.FIT_KEYS + sc.SKY_KEYS + ["log_likelihood"]
    ps = np.zeros(50, dtype=[(c, float) for c in cols])
    ps["log_likelihood"] = np.arange(50)
    f = tmp_path / "pe.h5"
    with h5py.File(f, "w") as h:
        g = h.create_group("C01:IMRPhenomXPHM")
        g["posterior_samples"] = ps
        md = g.create_group("meta_data/meta_data")
        for k, v in zip(names + ("f_low", "f_ref"), (4.0, 1257296851.2, 1024.0, 20.0, 20.0)):
            md[k] = [v]
        g.create_group("psds")["H1"] = np.ones((10, 2))
        h.create_group("C01:Mixed")["posterior_samples"] = ps
        h["C01:Mixed"].create_group("meta_data/meta_data")["f_low"] = [20.0]
        h["C01:Mixed"].create_group("psds")
    s = sc.extract_pe(str(f), "C01:IMRPhenomXPHM", tmp_path / "wd")
    assert (s["duration"], s["segment_start"], s["pe_sample_rate"]) == (4.0, 1257296851.2, 1024.0)
    with pytest.raises(ValueError, match="does not record its segment_length"):
        sc.extract_pe(str(f), "C01:Mixed", tmp_path / "wd2")
