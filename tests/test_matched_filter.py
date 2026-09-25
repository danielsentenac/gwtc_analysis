"""Offline tests for ``compute_matched_filter_snr`` on injected signals.

A chirp template is injected at a known time into white Gaussian noise; the
matched filter must peak at that time with the injected SNR, both for a short
(BBH-like) template and for a long (BNS-like) one that does not fit in the
overlay strain segment and needs the longer fetch.
"""
from __future__ import annotations

import numpy as np
import pytest
from gwpy.timeseries import TimeSeries

from gwtc_analysis import gwpe_utils as gu

FS = 4096
T0 = 1_000_000_000.0          # merger arrival time in the detector
SIGMA = 1e-21                 # white noise standard deviation per sample
TARGET_SNR = 20.0


def _chirp(duration: float, delta_t: float) -> np.ndarray:
    """Chirp ending at t = 0, frequency rising to 200 Hz, tapered at the start."""
    t = np.arange(-duration, 0.0, delta_t)
    tau = np.maximum(-t, 1e-3)
    f = np.minimum(200.0, 30.0 * (tau / duration) ** (-3.0 / 8.0))
    phase = 2 * np.pi * np.cumsum(f) * delta_t
    h = np.cos(phase)
    ramp = min(len(h), int(1.0 / delta_t))
    h[:ramp] *= 0.5 * (1 - np.cos(np.pi * np.arange(ramp) / ramp))
    return h


class _Posterior:
    def __init__(self, duration):
        self.duration = duration

    def maxL_td_waveform(self, approximant, delta_t, f_low, f_ref, project):
        h = _chirp(self.duration, delta_t)
        return TimeSeries(h, t0=T0 - self.duration, dt=delta_t)


class _PEData:
    def __init__(self, duration):
        self.labels = ["C00:Fake"]
        self.approximant = ["Fake"]
        self.samples_dict = {"C00:Fake": _Posterior(duration)}
        self.config = {}


def _strain(duration_tmpl: float, before: float, after: float, seed: int = 1) -> TimeSeries:
    """White noise with the template injected at T0, scaled to TARGET_SNR."""
    rng = np.random.default_rng(seed)
    start = T0 - before
    n = int((before + after) * FS)
    x = rng.normal(0.0, SIGMA, n)
    h = _chirp(duration_tmpl, 1.0 / FS)
    band = 1.0 / FS
    # optimal SNR of h in white noise restricted to 20-1024 Hz (the filter band)
    hf = np.fft.rfft(h) / FS
    fr = np.fft.rfftfreq(len(h), 1.0 / FS)
    psd = 2 * SIGMA ** 2 * band
    m = (fr >= 20) & (fr <= 1024)
    snr1 = np.sqrt(4 * np.sum(np.abs(hf[m]) ** 2) / psd * (fr[1] - fr[0]))
    k0 = int(round((T0 - duration_tmpl - start) * FS))
    lo, hi = max(0, k0), min(n, k0 + len(h))
    x[lo:hi] += (TARGET_SNR / snr1) * h[lo - k0:hi - k0]
    return TimeSeries(x, t0=start, dt=1.0 / FS)


def _run(pedata, strain, **kw):
    logs = [""]
    t, rho, _ = gu.compute_matched_filter_snr(
        strain=strain, pedata=pedata, det="H1", label="C00:Fake", t0=T0,
        fmin=20.0, fmax=1000.0, event_logs=logs, **kw,
    )
    return t, rho, logs[-1]


def _peak(t, rho):
    a = np.abs(rho)
    k = int(np.argmax(a))
    return float(a[k]), float(t[k] - T0)


def test_short_template_peaks_at_the_merger():
    """BBH-like 2 s template in 28 s of strain: peak at T0 with the injected SNR."""
    t, rho, logs = _run(_PEData(2.0), _strain(2.0, 14, 14))
    val, dt = _peak(t, rho)
    assert abs(dt) < 1e-3, logs
    assert 17 < val < 23, logs


def test_long_template_fetches_longer_strain(monkeypatch):
    """BNS-like 60 s template: longer strain is fetched and the peak is recovered."""
    calls = []

    def fake_load_strain(event, t0, detector, window=14.0, cache=True):
        calls.append(window)
        return _strain(60.0, window, window)

    monkeypatch.setattr(gu, "load_strain", fake_load_strain)
    t, rho, logs = _run(_PEData(60.0), _strain(60.0, 14, 14), event="GWFAKE")
    assert calls and calls[0] >= 60, logs
    val, dt = _peak(t, rho)
    assert abs(dt) < 1e-3, logs
    assert 17 < val < 23, logs


def test_long_template_without_event_is_skipped():
    """Without an event name to fetch more strain, a long template is skipped cleanly."""
    t, rho, logs = _run(_PEData(60.0), _strain(60.0, 14, 14))
    assert t is None and rho is None
    assert "Skipping matched-filter SNR" in logs


def _rotated_strain(alpha: float, shift: float, seed: int = 3) -> TimeSeries:
    """Noise plus the 2 s chirp shifted by `shift` s and phase-rotated by `alpha`."""
    from scipy.signal import hilbert

    rng = np.random.default_rng(seed)
    start = T0 - 14
    n = 28 * FS
    x = rng.normal(0.0, SIGMA, n)
    h = np.real(hilbert(_chirp(2.0, 1.0 / FS)) * np.exp(1j * alpha))
    k0 = int(round((T0 + shift - 2.0 - start) * FS))  # the injection lands on the sample grid
    x[k0:k0 + len(h)] += 6e-22 * h / np.std(h)
    return TimeSeries(x, t0=start, dt=1.0 / FS)


def _injected_shift(shift: float) -> float:
    start = T0 - 14
    return start + round((T0 + shift - 2.0 - start) * FS) / FS + 2.0 - T0


class _ShiftedPosterior(_Posterior):
    """Posterior whose maxL template is the aligned one."""

    def __init__(self, duration, align):
        super().__init__(duration)
        self.align = align

    def maxL_td_waveform(self, approximant, delta_t, f_low, f_ref, project):
        return gu._align_template(super().maxL_td_waveform(approximant, delta_t, f_low, f_ref, project), *self.align)


@pytest.mark.parametrize("alpha, shift", [(0.9, 0.004), (-2.0, -0.0063), (2.8, 0.0)])
def test_alignment_moves_template_onto_the_signal(alpha, shift):
    """The peak time/phase recovered from the matched filter align the template with the data."""
    strain = _rotated_strain(alpha, shift)
    t, rho, _ = _run(_PEData(2.0), strain)
    dt, dphi, val = gu.matched_filter_alignment(t, rho, T0)
    # recovered within the noise scatter of the injection (a few tenths of a ms at |rho|~50)
    assert abs(dt - _injected_shift(shift)) < 5e-4
    assert abs(np.angle(np.exp(1j * (dphi - alpha)))) < 0.15

    aligned = _PEData(2.0)
    aligned.samples_dict["C00:Fake"] = _ShiftedPosterior(2.0, (dt, dphi))
    t2, rho2, _ = _run(aligned, strain)
    dt2, dphi2, val2 = gu.matched_filter_alignment(t2, rho2, T0)
    # the aligned template sits on the data's best fit: no residual shift/phase, full SNR in phase
    assert abs(dt2) < 1e-4 and abs(dphi2) < 0.05
    assert val2 == pytest.approx(val, rel=0.01)
    assert rho2[int(np.argmin(np.abs(t2 - T0)))].real == pytest.approx(val, rel=0.01)


def test_alignment_ignores_noise():
    """Below the SNR threshold no alignment is returned."""
    noise = TimeSeries(np.random.default_rng(9).normal(0.0, SIGMA, 28 * FS), t0=T0 - 14, dt=1.0 / FS)
    t, rho, _ = _run(_PEData(2.0), noise)
    assert gu.matched_filter_alignment(t, rho, T0) is None
