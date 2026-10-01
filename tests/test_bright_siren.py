"""Offline tests for the `bright_siren` mode (likelihood, sky conditioning, combination, report and CLI)."""
from __future__ import annotations

import h5py
import numpy as np
import pandas as pd
import pytest

from gwtc_analysis import bright_siren as bs
from gwtc_analysis.cli import build_parser

CP = bs.COUNTERPARTS["GW170817"]


def test_gw170817_hubble_flow_velocity():
    """v_H = 3327 - 310 = 3017 ± 166 km/s, as in LVK 2017."""
    assert CP.v_recession - CP.v_peculiar == 3017
    assert np.hypot(CP.sigma_recession, CP.sigma_peculiar) == pytest.approx(166, abs=0.5)


def test_likelihood_peaks_at_v_over_d():
    """Narrow distances around 43 Mpc and v_H = 3010 km/s give H0 = 70."""
    rng = np.random.default_rng(1)
    h0 = np.linspace(10, 200, 3801)
    p = bs.h0_likelihood(h0, rng.normal(43.0, 0.2, 5000), 3010.0, 10.0)
    assert h0[np.argmax(p)] == pytest.approx(70.0, abs=0.3)


def test_summarize_gaussian():
    """MAP, 68.3% HPD and 90% interval of a Gaussian."""
    h0 = np.linspace(10, 200, 19001)
    s = bs.summarize(h0, np.exp(-0.5 * ((h0 - 70) / 5) ** 2))
    assert s["map"] == pytest.approx(70, abs=0.02)
    assert s["hpd68_low"] == pytest.approx(65, abs=0.1) and s["hpd68_high"] == pytest.approx(75, abs=0.1)
    assert s["median"] == pytest.approx(70, abs=0.02)
    assert s["low_90"] == pytest.approx(70 - 1.645 * 5, abs=0.05)


def test_sky_conditioning():
    """Samples fixed to the counterpart are all kept; otherwise only those near it, or an error."""
    n = 1000
    ra0, dec0 = np.radians(CP.ra_deg), np.radians(CP.dec_deg)
    d = np.linspace(20, 50, n)
    fixed = dict(luminosity_distance=d, ra=np.full(n, ra0), dec=np.full(n, dec0))
    out, how = bs.sky_conditioned_distances(fixed, CP, 3.0)
    assert len(out) == n and "fixed" in how

    spread = dict(luminosity_distance=d, ra=ra0 + np.radians(np.linspace(-10, 10, n)), dec=np.full(n, dec0))
    out, how = bs.sky_conditioned_distances(spread, CP, 3.0)
    assert 250 < len(out) < 330 and "within 3.0°" in how
    with pytest.raises(ValueError, match="increase --sky-radius"):
        bs.sky_conditioned_distances(spread, CP, 0.5)


def test_spectral_density_and_reader(tmp_path):
    """A spectral posterior is read from a TSV or a work directory, and its density integrates to 1 in the prior."""
    rng = np.random.default_rng(2)
    samples = np.clip(rng.normal(90, 30, 4000), 10.5, 199.5)
    pd.DataFrame({"H0": samples, "alpha": 1.0}).to_csv(tmp_path / "posterior.tsv", sep="\t", index=False)
    assert bs.read_spectral_posterior(tmp_path).shape == (4000,)
    h0 = np.linspace(10, 200, 3801)
    p = bs.spectral_density(h0, bs.read_spectral_posterior(tmp_path / "posterior.tsv"))
    assert bs._trapz(p, h0) == pytest.approx(1, abs=0.02)
    (tmp_path / "empty").mkdir()
    with pytest.raises(ValueError, match="no posterior"):
        bs.read_spectral_posterior(tmp_path / "empty")


def _pe_file(path, labels=("C02:Test-HighSpin", "C02:Test-LowSpin")):
    rng = np.random.default_rng(3)
    with h5py.File(path, "w") as h:
        for lab in labels:
            n = 3000
            arr = np.zeros(n, dtype=[("luminosity_distance", "f8"), ("ra", "f8"), ("dec", "f8")])
            arr["luminosity_distance"] = rng.normal(43.0, 3.0, n)
            arr["ra"], arr["dec"] = np.radians(CP.ra_deg), np.radians(CP.dec_deg)
            h.create_group(lab).create_dataset("posterior_samples", data=arr)
    return path


def test_run_bright_siren_end_to_end(tmp_path):
    """Both labels (LowSpin first), the combination, the TSVs and the report."""
    pe =_pe_file(tmp_path / "pe.h5")
    rng = np.random.default_rng(4)
    pd.DataFrame({"H0": np.clip(rng.normal(90, 30, 4000), 10.5, 199.5)}).to_csv(tmp_path / "spec.tsv", sep="\t",
                                                                                  index=False)
    table = bs.run_bright_siren(pe_file=pe, spectral_posterior=tmp_path / "spec.tsv",
                                out_report_html=tmp_path / "r.html", out_summary_tsv=tmp_path / "s.tsv",
                                plots_dir=tmp_path / "plots")
    assert list(table["analysis"])[:2] == ["bright siren, C02:Test-LowSpin", "bright siren, C02:Test-HighSpin"]
    assert table["map"].iloc[0] == pytest.approx(3017 / 43.0, rel=0.03)
    assert "spectral" in table["analysis"].iloc[-1]
    grid = pd.read_csv(tmp_path / "s.posterior.tsv", sep="\t")
    assert {"H0", "p_spectral", "p_combined"} <= set(grid.columns)
    assert (tmp_path / "plots" / "h0_bright_siren_GW170817.png").exists()
    assert (tmp_path / "plots" / "h0_combined.png").exists()
    assert "NGC 4993" in (tmp_path / "r.html").read_text()

    with pytest.raises(ValueError, match="not in pe.h5"):
        bs.run_bright_siren(pe_file=pe, pe_labels=["C02:Nope"], out_report_html=None, out_summary_tsv=None)


def test_cli_parses_bright_siren():
    args = build_parser().parse_args(["bright_siren", "--v-peculiar", "300", "100", "--pe-label", "A", "B"])
    assert args.src_name == "GW170817" and args.v_peculiar == [300.0, 100.0] and args.pe_label == ["A", "B"]
    assert args.h0_range == [10.0, 200.0]
