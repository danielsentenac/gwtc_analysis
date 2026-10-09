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
    """alpha = 0.1 for M = 20.6 Msun is m_phi = 6.5e-13 eV; the cloud mass fraction rho r_phi^3 / M."""
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
    """rho = 0 is IMRPhenomXPHM exactly (aligned and precessing); the cloud shortens the inspiral by the time the
    dephasing predicts."""
    from bilby.gw.source import lal_binary_black_hole

    for base in (ALIGNED, PRECESSING):
        vac = lal_binary_black_hole(F, **base, **WF)
        for rho in (0.0, 1e-30):
            h = sc.scalar_cloud_binary_black_hole(F, **base, rho_phi=rho, alpha_phi=0.1)
            assert np.allclose(h["plus"], vac["plus"], rtol=0, atol=1e-35)
    vac = lal_binary_black_hole(F, **ALIGNED, **WF)
    h = sc.scalar_cloud_binary_black_hole(F, **ALIGNED, rho_phi=1e6, alpha_phi=0.1)
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

    scripts = write_slurm_chains(tmp_path, "/opt/env/bin/python", ["fit", "background", "sample"], 20,
                                 ["vacuum", "agnostic"], 1e6, 16, ["--partition=htc"])
    assert [s.parent.name for s in scripts] == ["slurm", "slurm_sample"]
    names = sorted(p.name for p in (tmp_path / "slurm").glob("0*.sh"))
    assert names == ["00_vacuum.sh", "01_profile.sh", "02_polish.sh", "03_background.sh"]
    prof = (tmp_path / "slurm" / "01_profile.sh").read_text()
    assert f"--array=0-{len(sc.ALPHAS) * len(sc.DEPHASING_TARGETS) - 1}" in prof
    assert "--array=0-19" in (tmp_path / "slurm" / "03_background.sh").read_text()
    samp = (tmp_path / "slurm_sample" / "00_sample.sh").read_text()
    assert "--array=0-1" in samp and "MODELS=(vacuum agnostic)" in samp and "--cpus-per-task=16" in samp
