"""search_skymaps: GWTC-1 read from the GWTC-2.1 skymaps, restricted to its events."""
from __future__ import annotations

from gwtc_analysis import gw_stat
from gwtc_analysis import search_skymaps as ss


def test_gwtc1_reads_the_gwtc21_skymaps_with_its_own_events(monkeypatch):
    lists = {"GWTC-1-confident": {"events": {"GW150914-v3": {}, "GW170818-v1": {}}}}
    monkeypatch.setattr(gw_stat, "fetch_gwtc_events", lambda name: lists[name])
    assert ss._skymap_sources(["GWTC-1"]) == (["GWTC-2.1"], {"GWTC-2.1": {"GW150914-v3", "GW170818-v1"}})
    # named together with GWTC-2.1, the whole release is read once
    assert ss._skymap_sources(["GWTC-1", "GWTC-2.1", "GWTC-3"]) == (["GWTC-2.1", "GWTC-3"], {})


def test_event_selected_matches_short_and_full_names():
    sel = {"GW170818-v1", "GW200311_115853-v1"}
    assert ss._event_selected("GW170818_022509_PEDataRelease_cosmo_reweight_C01", sel)
    assert ss._event_selected("GW200311_115853_PEDataRelease", sel)
    assert not ss._event_selected("GW200311_000000_PEDataRelease", sel)
    assert not ss._event_selected("GW190521_030229_PEDataRelease", sel)
