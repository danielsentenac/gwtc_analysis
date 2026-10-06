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


def test_skymap_waveform_and_choice():
    """GWTC-4.0/5.0 names have no ':': each waveform must get its own index key (they all fell under
    'Unknown' and the first map of the archive was kept), and GWTC-5.0, without Mixed maps, falls back to
    IMRPhenomXPHM_SpinTaylor."""
    from gwtc_analysis.gw_stat import choose_skymap, select_skymap_member, skymap_waveform

    g5 = "IGWN-GWTC5p0-29ebe06b7_25-GW240615_113620-{}_Skymap_PEDataRelease.fits.gz"
    assert skymap_waveform("parameter_estimation/skymaps/" + g5.format("IMRPhenomXPHM_SpinTaylor")) == "IMRPhenomXPHM_SpinTaylor"
    assert skymap_waveform("IGWN-GWTC3p0-v3-GW191219_163120_PEDataRelease_mixed_cosmo_reweight_C01:"
                           "SEOBNRv4_ROM_NRTidalv2_NSBH:HighSpin.fits") == "SEOBNRv4_ROM_NRTidalv2_NSBH:HighSpin"
    assert skymap_waveform("IGWN-GWTC2p1-v2-GW150914_095045_PEDataRelease_cosmo_reweight_C01_Mixed.fits") == "Mixed"
    index = {("GW240615_113620", w): g5.format(w) for w in ("SEOBNRv5PHM", "NRSur7dq4", "IMRPhenomXPHM_SpinTaylor")}
    assert select_skymap_member(index, "GW240615_113620", prefer="Mixed") == g5.format("IMRPhenomXPHM_SpinTaylor")
    assert select_skymap_member(index, "GW240615_113620", prefer="C00:NRSur7dq4") == g5.format("NRSur7dq4")
    assert select_skymap_member(index, "GW240615_160735") is None
    maps = {"Mixed:NSBH:HighSpin": 1, "Mixed": 2, "IMRPhenomXPHM": 3}
    assert choose_skymap(maps, "Mixed") == ("Mixed", 2)
    assert choose_skymap({"IMRPhenomXAS:HighSpin": 1, "IMRPhenomNSBH:LowSpin": 2}, "Mixed") == ("IMRPhenomNSBH:LowSpin", 2)

