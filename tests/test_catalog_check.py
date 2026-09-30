"""check_catalogs against a fake GWOSC and Zenodo (no network)."""
from __future__ import annotations

from gwtc_analysis import catalog_check as cc
from gwtc_analysis import catalog_registry as reg

G, Z = cc.GWOSC, cc.ZENODO


def _fake_get(extra_lists=(), extra_runs=(), gwtc5_latest="20348005"):
    lists = [{"name": n} for n in ("GWTC", "GWTC-1-confident", "GWTC-1-marginal", "GWTC-2", "GWTC-2.1-auxiliary",
                                   "GWTC-2.1-confident", "GWTC-2.1-marginal", "GWTC-3-confident", "GWTC-3-marginal",
                                   "GWTC-4.0", "GWTC-4.1", "GWTC-5.0", "O3_Discovery_Papers", "IAS-O3a")]
    lists += [{"name": n} for n in extra_lists]
    runs = [{"name": r.name, "gps_start": r.start_gps, "gps_end": r.end_gps} for r in reg.OBSERVING_RUNS.values()]
    runs += [{"name": n, "gps_start": a, "gps_end": b} for n, a, b in extra_runs]
    table = {
        f"{G}/api/v2/catalogs": {"results": lists, "next": None},
        f"{G}/api/v2/runs": {"results": runs, "next": None},
        # a new catalog of the O4c run, with its PE files on Zenodo record 777 (concept 776)
        f"{G}/eventapi/jsonfull/GWTC-6.0/": {"events": {
            "GW250301_000000-v1": {"GPS": 1424822418.0}, "GW250601_000000-v1": {"GPS": 1432771218.0}}},
        f"{G}/api/v2/event-versions/GW250301_000000-v1/parameters": {"results": [
            {"pipeline_type": "pe", "data_url": f"{Z}/api/records/777/files/GW250301_PEDataRelease.hdf5/content"},
            {"pipeline_type": "search", "data_url": ""}], "next": None},
        f"{G}/api/v2/event-versions/GW250601_000000-v1/parameters": {"results": [
            {"pipeline_type": "pe", "data_url": f"{Z}/api/records/777/files/GW250601_PEDataRelease.hdf5/content"}],
            "next": None},
        f"{Z}/api/records/777": {"id": 777, "conceptrecid": "776", "metadata": {"title": "GWTC-6.0: PE data release"}},
        f"{Z}/api/records/777/versions?size=25&sort=version": {"hits": {"hits": [
            {"id": 777, "metadata": {"title": "GWTC-6.0: PE data release", "publication_date": "2026-12-15"},
             "files": [{"key": "GW250301_PEDataRelease.hdf5"}, {"key": "IGWN-GWTC6p0-Archived_Skymaps.tar.gz"}]}]},
            "links": {}},
    }

    def get(url):
        if url in table:
            return table[url]
        for part in (p for ps in reg.zenodo_releases().values() for p in ps):      # the registry's own records
            rid = part.record_id
            if url == f"{Z}/api/records/{rid}":
                return {"id": int(rid), "conceptrecid": "1", "metadata": {"title": "t"}}
            if url == f"{Z}/api/records/{rid}/versions?size=25&sort=version":
                latest = gwtc5_latest if rid == "20348005" else rid
                return {"hits": {"hits": [{"id": int(latest), "metadata": {"publication_date": "2026-01-01"},
                                           "files": []}]}, "links": {}}
        raise KeyError(url)

    return get


def test_nothing_to_report_when_registry_is_complete():
    """The registry already describes every GWTC list and run (GWTC-4.1 included)."""
    r = cc.check_catalogs(get=_fake_get())
    assert r["new_lists"] == [] and r["new_runs"] == {} and r["missing_lists"] == [] and r["newer_versions"] == []
    assert "describes every GWTC event list" in cc.format_report(r)


def test_new_catalog_and_run_are_reported_with_a_draft_entry():
    """A new GWOSC list and run: the draft entry has the Zenodo record, concept, tarball and the new run note."""
    r = cc.check_catalogs(get=_fake_get(extra_lists=["GWTC-6.0"], extra_runs=[("O4c", 1422118818, 1447430418),
                                                                              ("O4c1DiscC00", 1422962688, 1422966784)]))
    assert [d["name"] for d in r["new_lists"]] == ["GWTC-6.0"]
    assert list(r["new_runs"]) == ["O4c"]                         # data-release runs are not observing runs
    d = r["new_lists"][0]
    assert d["key"] == "GWTC-6" and d["n_after_known_runs"] == 2
    z = d["zenodo"][0]
    assert (z["record"], z["concept"], z["skymap_tarball"]) == ("777", "776", "IGWN-GWTC6p0-Archived_Skymaps.tar.gz")
    draft = cc.draft_entry(d)
    assert 'Catalog("GWTC-6", "GWTC-6.0"' in draft and 'ZenodoRelease("777", "IGWN-GWTC6p0-Archived_Skymaps.tar.gz")' in draft
    assert "add the new run to OBSERVING_RUNS" in draft and "update_of" not in draft
    assert "engineering run" in cc.format_report(r)


def test_reanalysis_is_drafted_as_an_update_and_newer_versions_reported():
    """A list over the runs of an existing catalog is drafted with update_of; a newer Zenodo version is reported."""
    d = dict(name="GWTC-5.1", key="GWTC-5.1", n_events=3, gps_range=(1.4e9, 1.41e9), runs={"O4b": 3},
             n_after_known_runs=0, zenodo=[])
    assert 'update_of="GWTC-5"' in cc.draft_entry(d)
    r = cc.check_catalogs(get=_fake_get(gwtc5_latest="30000000"))
    assert r["newer_versions"] == [dict(catalog="GWTC-5", registry="20348005", latest="30000000", n_versions=1)]


def test_suggested_keys():
    assert cc.suggested_key("GWTC-6.0") == "GWTC-6" and cc.suggested_key("GWTC-4.1") == "GWTC-4.1"
    assert cc.suggested_key("GWTC-3-confident") == "GWTC-3"
