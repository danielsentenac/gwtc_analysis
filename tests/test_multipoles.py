"""Offline tests of the higher-multipole and precession summary of parameters_estimation."""
from __future__ import annotations

import numpy as np
import pytest

from gwtc_analysis import multipoles as mp


def _samples(rng, n=4000, shift=0.0, fields=tuple(mp.FIELDS)):
    """Rayleigh-distributed SNRs (noise only), offset by `shift` (a real multipole)."""
    return {f: np.hypot(rng.normal(shift, 1, n), rng.normal(0, 1, n)) for f in fields} | {"mass_1": np.ones(n)}


def test_noise_scale():
    """P(ρ > x) = exp(-x²/2) for the noise-only Rayleigh distribution: 11% at 2.1, 1% at 3."""
    rng = np.random.default_rng(1)
    rho = np.hypot(rng.normal(size=200_000), rng.normal(size=200_000))
    assert np.mean(rho > mp.HINT) == pytest.approx(np.exp(-mp.HINT ** 2 / 2), abs=0.003)
    assert np.exp(-mp.HINT ** 2 / 2) == pytest.approx(0.11, abs=0.005)
    assert np.exp(-mp.CLEAR ** 2 / 2) == pytest.approx(0.011, abs=0.001)
    assert [mp.evidence(x) for x in (1.0, 2.5, 3.2, np.nan)] == ["none", "hint", "clear", ""]


def test_summary_and_report(tmp_path):
    rng = np.random.default_rng(2)
    sd = {"C00:Noise": _samples(rng), "C00:Loud": _samples(rng, shift=5.0)}
    table, files, html = mp.analyse(sd, "C00:Loud", "GWTEST", tmp_path)
    assert len(table) == 8 and set(table["quantity"]) == set(mp.FIELDS.values())
    loud = table[table["label"] == "C00:Loud"]
    assert (loud["evidence"] == "clear").all() and (table[table["label"] == "C00:Noise"]["median"] < 1.5).all()
    assert loud["frac_above_3"].min() > 0.9
    assert [p.name for p in files] == ["GWTEST_multipoles_precession.tsv", "multipoles_precession_GWTEST.png"]
    assert all(p.stat().st_size > 0 for p in files)
    assert "clear evidence for ρ₃₃" in html and "waveform models disagree" in html


def test_files_without_the_snrs(tmp_path):
    """GWTC-1 to GWTC-3 files: no table, no files, and the report says so."""
    rng = np.random.default_rng(3)
    table, files, html = mp.analyse({"C01:Mixed": {"mass_1": np.ones(10)}}, "C01:Mixed", "GWOLD", tmp_path)
    assert table.empty and files == [] and "stores no multipole" in html
    partial = {"C00:A": _samples(rng, fields=("network_precessing_snr",))}
    t, f, _ = mp.analyse(partial, "C00:A", "GWP", tmp_path)
    assert t["quantity"].tolist() == ["ρ_p"] and len(f) == 2
