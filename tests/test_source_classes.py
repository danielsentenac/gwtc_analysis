"""Offline tests of the remnant estimates and the source-class presets (source_classes, event_selection)."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from gwtc_analysis import gw_stat as gw
from gwtc_analysis import source_classes as sc


def test_final_spin_fit_limits():
    """0.686 for equal masses without spin, the spin itself for a test particle, ~0.83 for equal aligned 0.5."""
    assert sc.final_spin_aligned(0.25, 0.0) == pytest.approx(0.6865, abs=1e-3)
    assert sc.final_spin_aligned(1e-6, 0.4) == pytest.approx(0.4, abs=1e-4)
    assert sc.final_spin_aligned(0.25, 0.5) == pytest.approx(0.832, abs=2e-3)
    assert sc.final_spin_aligned(0.25, -0.5) < sc.final_spin_aligned(0.25, 0.0)


def _raw(name, **kw):
    v = dict(commonName=name, GPS=1.0, version=1)
    for k, (med, lo, hi) in kw.items():
        v.update({k: med, f"{k}_lower": lo, f"{k}_upper": hi})
    return v


def test_bounds_are_absolute_and_swapped_with_the_masses():
    """GWOSC offsets become absolute 90% bounds; they follow the masses when m2 > m1 is swapped."""
    raw = {"A": _raw("A", mass_1_source=(10.0, -2.0, 3.0), mass_2_source=(30.0, -5.0, 4.0), chi_eff=(-0.2, -0.1, 0.15))}
    df = gw.events_to_dataframe(raw)
    assert df.loc[0, "mass_1_source_lo90"] == 8.0 and df.loc[0, "mass_1_source_hi90"] == 13.0
    assert df.loc[0, "chi_eff_hi90"] == pytest.approx(-0.05)
    out = gw.prepare_catalog_df(df.assign(chirp_mass_source=np.nan))
    assert out.loc[0, "mass_1_source"] == 30.0
    assert (out.loc[0, "mass_1_source_lo90"], out.loc[0, "mass_1_source_hi90"]) == (25.0, 34.0)
    assert (out.loc[0, "mass_2_source_lo90"], out.loc[0, "mass_2_source_hi90"]) == (8.0, 13.0)


def test_remnant_columns():
    """E_rad = M_total - M_final (GWOSC total, else m1 + m2); no final-spin estimate with a neutron star."""
    df = pd.DataFrame(dict(mass_1_source=[36.0, 30.0, 1.5], mass_2_source=[29.0, 30.0, 1.3],
                           total_mass_source=[65.0, 60.0, 2.8], total_mass_source_gwosc=[64.6, np.nan, 2.8],
                           final_mass_source=[61.5, 57.1, np.nan], chi_eff=[-0.04, 0.0, 0.0],
                           binary_type=["BBH", "BBH", "NS-NS"]))
    out = sc.add_remnant_columns(df)
    assert out["radiated_energy_msun"].tolist()[:2] == pytest.approx([3.1, 2.9])
    assert out.loc[0, "radiated_energy_erg"] == pytest.approx(3.1 * 1.7877e54)
    assert out.loc[1, "radiated_fraction"] == pytest.approx(2.9 / 60.0)
    assert out.loc[1, "final_spin_estimate"] == pytest.approx(0.6865, abs=1e-3)
    assert np.isnan(out.loc[2, "radiated_energy_msun"]) and np.isnan(out.loc[2, "final_spin_estimate"])


def _table():
    return pd.DataFrame(dict(
        event_id=["BNS", "NSBH", "GAP", "GW190814", "HEAVY", "NEG", "PLAIN"],
        mass_1_source=[1.5, 8.0, 3.7, 23.0, 98.0, 16.0, 30.0], mass_1_source_lo90=[1.4, 7.0, 2.5, 21.0, 77.0, 14.0, 25.0],
        mass_1_source_hi90=[1.6, 9.0, 4.5, 25.0, 110.0, 18.0, 35.0],
        mass_2_source=[1.3, 1.4, 1.4, 2.6, 57.0, 8.0, 25.0], mass_2_source_lo90=[1.2, 1.3, 1.2, 2.5, 40.0, 7.0, 20.0],
        mass_2_source_hi90=[1.4, 1.6, 2.0, 2.7, 70.0, 9.0, 30.0],
        chi_eff=[0.0, 0.0, 0.0, 0.0, -0.1, -0.3, 0.1], chi_eff_lo90=[-0.1] * 4 + [-0.4, -0.5, 0.0],
        chi_eff_hi90=[0.1] * 4 + [0.2, -0.08, 0.2]))


def test_presets():
    df = _table()
    mask, extra, _ = sc.apply_preset(df, "neutron-stars")
    assert df.event_id[mask].tolist() == ["BNS", "NSBH", "GAP", "GW190814"]
    assert extra["class"][mask].tolist() == ["BNS", "NSBH", "NSBH", "NSBH"]
    mask, extra, _ = sc.apply_preset(df, "mass-gap")
    assert df.event_id[mask].tolist() == ["GAP"] and extra.loc[2, "gap_component"] == "primary"
    assert not extra.loc[2, "gap_confident"]                  # its 90% interval leaves the gap
    assert sc.apply_preset(df, "mass-gap", mass_gap=(2.0, 5.0))[1].loc[2, "gap_confident"]
    mask, extra, _ = sc.apply_preset(df, "hierarchical")
    assert df.event_id[mask].tolist() == ["HEAVY", "NEG"]
    assert extra.loc[4, "pisn_gap_confident"] and extra.loc[5, "negative_chi_eff"]
    with pytest.raises(ValueError, match="unknown preset"):
        sc.apply_preset(df, "exotic")


def test_event_selection_preset(tmp_path, monkeypatch):
    """--preset adds its columns to the TSV and combines with the cuts; the plot is written."""
    from gwtc_analysis import event_selection as es

    monkeypatch.setattr(es.gw, "fetch_gwtc_events", lambda catalog: {"events": {}})
    monkeypatch.setattr(es.gw, "events_to_dataframe", lambda ev: _table().assign(luminosity_distance=500.0,
                                                                                 redshift=0.1))
    es.run_event_selection(catalogs=["GWTC-4"], out_tsv=tmp_path / "h.tsv", preset="hierarchical",
                           out_plot=tmp_path / "h.png")
    out = pd.read_csv(tmp_path / "h.tsv", sep="\t")
    assert out["event_id"].tolist() == ["HEAVY", "NEG"] and {"pisn_gap", "negative_chi_eff", "chi_eff"} <= set(out)
    assert (tmp_path / "h.png").stat().st_size > 1000
    es.run_event_selection(catalogs=["GWTC-4"], out_tsv=tmp_path / "n.tsv", preset="neutron-stars", m1_min=5)
    assert pd.read_csv(tmp_path / "n.tsv", sep="\t")["event_id"].tolist() == ["GW190814", "NSBH"]
    with pytest.raises(ValueError, match="unknown preset"):
        es.run_event_selection(catalogs=["GWTC-4"], out_tsv=tmp_path / "x.tsv", preset="exotic")


def test_cli_presets():
    from gwtc_analysis.cli import build_parser

    a = build_parser().parse_args(["event_selection", "--catalogs", "ALL", "--preset", "mass-gap", "--mass-gap", "2.5", "5"])
    assert a.preset == "mass-gap" and a.mass_gap == [2.5, 5.0] and a.pisn_gap_min == sc.PISN_GAP_MIN
