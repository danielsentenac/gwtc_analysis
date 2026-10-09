"""Offline tests for the `counterpart` mode: registry, sky searched probability, viewing angle and its constraint,
the report, the CLI."""
from __future__ import annotations

import h5py
import numpy as np
import pandas as pd
import pytest

from gwtc_analysis import counterpart as cpm
from gwtc_analysis.cli import build_parser

CP = cpm.COUNTERPARTS["GW170817"]
Z170817 = 3017 / cpm.C_KMS


def test_registry_and_position_override():
    assert cpm.get_counterpart("GW170817") is CP
    other = cpm.get_counterpart("GW170817", 10.0, -5.0)
    assert (other.ra_deg, other.dec_deg, other.host) == (10.0, -5.0, "NGC 4993") and other.published is None
    with pytest.raises(ValueError, match="no counterpart registered"):
        cpm.get_counterpart("GW150914")
    with pytest.raises(ValueError, match="both"):
        cpm.get_counterpart("GW170817", 10.0)


def test_sky_conditioning():
    rng = np.random.default_rng(1)
    n = 1000
    d = rng.normal(40, 5, n)
    ra0, dec0 = np.radians(CP.ra_deg), np.radians(CP.dec_deg)
    keep, how = cpm.sky_conditioned(dict(luminosity_distance=d, ra=np.full(n, ra0), dec=np.full(n, dec0)), CP, 3.0)
    assert keep.all() and "fixed" in how
    spread = dict(luminosity_distance=d, ra=ra0 + np.radians(rng.normal(0, 3, n)), dec=np.full(n, dec0))
    keep, how = cpm.sky_conditioned(spread, CP, 3.0)
    assert 0.4 < keep.mean() < 0.95 and "within 3.0°" in how
    with pytest.raises(ValueError, match="increase --sky-radius"):
        cpm.sky_conditioned(spread, CP, 0.5)


def _fisher_samples(rng, n, ra0, dec0, sigma_deg):
    """Samples of a narrow 2D Gaussian on the sky (tangent plane) around (ra0, dec0), in radians."""
    x, y = np.radians(rng.normal(0, sigma_deg, (2, n)))
    dec = np.radians(dec0) + y
    return np.radians(ra0) + x / np.cos(np.radians(dec0)), dec


@pytest.mark.parametrize("offset_sigma", [0.5, 1.0, 2.0])
def test_sky_credible_level_of_a_gaussian(offset_sigma):
    """For a 2D Gaussian, the searched probability at r sigma is 1 - exp(-r^2 / 2)."""
    rng = np.random.default_rng(7)
    ra, dec = _fisher_samples(rng, 4000, 120.0, 30.0, 2.0)
    level = cpm.sky_credible_level(ra, dec, 120.0, 30.0 + 2.0 * offset_sigma)
    assert level == pytest.approx(1 - np.exp(-offset_sigma ** 2 / 2), abs=0.06)


def test_weighted_quantiles_and_viewing_angle_weights():
    x = np.arange(1001.0)
    assert cpm.weighted_quantiles(x, (0.05, 0.5, 0.95)) == pytest.approx([50, 500, 950], abs=1)
    w = cpm.viewing_angle_weights(np.array([20.0, 25.0, 50.0]), (20.0, 5.0))
    assert w == pytest.approx([1.0, np.exp(-0.5), np.exp(-18)])


def test_viewing_angle_and_degeneracy_table():
    """theta_jn folded to 0-90 degrees; H0 implied by the host redshift; the table by viewing angle."""
    v = cpm.viewing_angle(dict(theta_jn=np.radians([10.0, 170.0, 90.0])))
    assert v == pytest.approx([10.0, 10.0, 90.0])
    assert cpm.viewing_angle(dict(iota=np.radians([30.0]))) == pytest.approx([30.0])
    assert cpm.viewing_angle(dict(luminosity_distance=[1.0])) is None
    assert cpm.implied_h0(np.array([43.0]), Z170817)[0] == pytest.approx(3017 / 43.0 * 1.0077, rel=0.002)
    d = np.array([45.0, 44.0, 36.0, 23.0])
    tab = cpm.degeneracy_table(d, np.array([10.0, 20.0, 45.0, 70.0]), Z170817)
    assert tab["fraction"].tolist() == [0.5, 0.25, 0.25] and tab["h0_median"].is_monotonic_increasing


def _pe_file(path, fixed_sky=True):
    rng = np.random.default_rng(3)
    n = 4000
    view = rng.uniform(0, 80, n)
    arr = np.zeros(n, dtype=[("luminosity_distance", "f8"), ("ra", "f8"), ("dec", "f8"), ("theta_jn", "f8")])
    arr["luminosity_distance"] = 47.0 - 0.3 * view + rng.normal(0, 1.0, n)
    arr["theta_jn"] = np.radians(view)
    if fixed_sky:
        arr["ra"], arr["dec"] = np.radians(CP.ra_deg), np.radians(CP.dec_deg)
    else:
        arr["ra"], arr["dec"] = _fisher_samples(rng, n, CP.ra_deg, CP.dec_deg + 2.0, 2.0)
    with h5py.File(path, "w") as h:
        h.create_group("C02:Test-LowSpin").create_dataset("posterior_samples", data=arr)
    return path


def test_run_counterpart_end_to_end(tmp_path):
    table = cpm.run_counterpart(pe_file=_pe_file(tmp_path / "pe.h5"), viewing_angle_constraint=(20.0, 3.0),
                                out_report_html=tmp_path / "c.html", out_summary_tsv=tmp_path / "c.tsv",
                                plots_dir=tmp_path / "plots")
    r = table.iloc[0]
    assert np.isnan(r["sky_credible_level"])                       # sky fixed to the counterpart
    assert r["distance_Planck"] == pytest.approx(3017 / 67.4 * 1.0077, rel=0.003)
    assert 0 < r["percentile_SH0ES"] < r["percentile_Planck"] < 1
    assert r["distance_constrained_median"] == pytest.approx(41.0, abs=1.0)
    assert r["viewing_angle_constrained_median"] == pytest.approx(20.0, abs=1.5)
    assert r["distance_constrained_high_90"] - r["distance_constrained_low_90"] < r["distance_high_90"] - r["distance_low_90"]
    assert pd.read_csv(tmp_path / "c.tsv", sep="\t").shape[0] == 1
    assert (tmp_path / "plots" / "distance_GW170817.png").exists()
    assert (tmp_path / "plots" / "distance_inclination_GW170817.png").exists()
    html = (tmp_path / "c.html").read_text()
    assert "degenerate with the inclination" in html and "--method bright" in html


def test_run_counterpart_sky_level(tmp_path):
    """Samples spread on the sky, centred 1 sigma away from the counterpart: searched probability ≈ 39%."""
    table = cpm.run_counterpart(pe_file=_pe_file(tmp_path / "pe.h5", fixed_sky=False), out_report_html=None,
                                out_summary_tsv=None)
    assert table.iloc[0]["sky_credible_level"] == pytest.approx(1 - np.exp(-0.5), abs=0.06)
    # a position far outside the sky posterior: no samples along its line of sight, reported without a distance
    far = cpm.run_counterpart(pe_file=tmp_path / "pe.h5", ra_deg=CP.ra_deg + 60, dec_deg=CP.dec_deg,
                              out_report_html=tmp_path / "far.html", out_summary_tsv=None)
    assert far.iloc[0]["sky_credible_level"] > 0.99 and "distance_median" not in far.columns
    assert "not estimated" in (tmp_path / "far.html").read_text()


def test_cli_parses_counterpart():
    args = build_parser().parse_args(["counterpart", "--event", "GW190521", "--viewing-angle", "30", "10",
                                      "--ra", "192.4", "--dec", "34.8"])
    assert args.event == "GW190521" and args.viewing_angle == [30.0, 10.0] and (args.ra, args.dec) == (192.4, 34.8)
    assert build_parser().parse_args(["counterpart"]).out_report == "counterpart.html"
    with pytest.raises(SystemExit):
        build_parser().parse_args(["bright_siren"])
