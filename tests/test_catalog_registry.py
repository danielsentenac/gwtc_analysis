"""The catalog registry: consistency, and a new catalog propagating to every derived table."""
from __future__ import annotations

import dataclasses

import pytest

from gwtc_analysis import catalog_registry as reg


def test_registry_is_consistent():
    """Every run, catalog and release refers to known entries, and the runs are ordered without overlap."""
    runs = list(reg.OBSERVING_RUNS.values())
    for a, b in zip(runs, runs[1:]):
        assert a.start_gps < a.end_gps <= b.start_gps
    for c in reg.CATALOGS.values():
        assert c.runs and all(r in reg.OBSERVING_RUNS for r in c.runs)
        assert c.zenodo or c.products_from in reg.CATALOGS, c.key       # PE files and skymaps come from somewhere
    for s in reg.SENSITIVITY_RELEASES.values():
        assert all(r in reg.OBSERVING_RUNS for r in s.runs)
    assert reg.DEFAULT_RATES_RELEASE in reg.SENSITIVITY_RELEASES and reg.DEFAULT_H0_RELEASE in reg.SENSITIVITY_RELEASES
    # the historical tables are reproduced
    assert reg.allowed_catalogs() == ("GWTC-1", "GWTC-2.1", "GWTC-3", "GWTC-4", "GWTC-5", "ALL")
    assert list(reg.zenodo_releases()) == ["GWTC-5", "GWTC-4", "GWTC-3", "GWTC-2.1"]     # newest first
    assert reg.products_catalog("GWTC-1") == "GWTC-2.1" and reg.products_catalog("GWTC-4") == "GWTC-4"
    assert reg.release_catalogs("gwtc4") == ("GWTC-1", "GWTC-2.1", "GWTC-3", "GWTC-4")


def test_a_new_catalog_is_one_registry_entry(monkeypatch):
    """Adding a catalog (here a mock GWTC-6 for an O4c run) updates every derived view."""
    runs = dict(reg.OBSERVING_RUNS, O4c=reg.ObservingRun("O4c", 1422118818, 1447430418))
    cats = dict(reg.CATALOGS, **{"GWTC-6": reg.Catalog("GWTC-6", "GWTC-6.0", ("O4c",),
                                                       zenodo=(reg.ZenodoRelease("99999999", "GWTC6_Skymaps.tar.gz"),),
                                                       s3_prefix="GWTC-6/")})
    rel = dict(reg.SENSITIVITY_RELEASES, gwtc6=dataclasses.replace(
        reg.SENSITIVITY_RELEASES["gwtc5"], key="gwtc6", record="88888888", runs=tuple(runs)))
    monkeypatch.setattr(reg, "OBSERVING_RUNS", runs)
    monkeypatch.setattr(reg, "CATALOGS", cats)
    monkeypatch.setattr(reg, "SENSITIVITY_RELEASES", rel)
    assert reg.allowed_catalogs()[-2:] == ("GWTC-6", "ALL")
    assert reg.gwosc_aliases()["GWTC-6"] == "GWTC-6.0" and reg.gwosc_lists()[-1] == "GWTC-6.0"
    assert list(reg.zenodo_releases())[0] == "GWTC-6"
    assert reg.skymap_catalogs()[-1] == "GWTC-6" and reg.s3_prefix("GWTC-6") == "GWTC-6/"
    assert reg.catalog_runs_map()["GWTC-6"] == ("O4c",) and "O4c" in reg.observing_runs()
    assert "GWTC-6" in reg.release_catalogs("gwtc6") and "GWTC-6" not in reg.release_catalogs("gwtc5")
    assert "GWTC-6: O4c" in reg.catalog_runs_help()
    with pytest.raises(KeyError):
        reg.release_catalogs("gwtc7")
