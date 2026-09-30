"""check_catalogs: compare what GWOSC and Zenodo publish with the catalog registry.

Reports the GWOSC event lists and observing runs that `catalog_registry` does not describe, and for each new
list proposes a draft registry entry: its observing runs (from the event times), and its Zenodo records,
found from the PE data links of its events in the GWOSC v2 API, with their concept IDs, latest versions and
skymap tarballs. It also reports registry records whose Zenodo release has a newer version (informational:
the latest version is used at run time anyway).
"""
from __future__ import annotations

import json
import re
import urllib.request
from collections import Counter
from typing import Callable, Optional

from . import catalog_registry as reg

GWOSC = "https://gwosc.org"
ZENODO = "https://zenodo.org"
_SKYMAP_TARBALL = re.compile(r"(skymap|skylocalization).*\.tar(\.gz)?$", re.I)
_RECORD_IN_URL = re.compile(r"zenodo\.org/(?:api/)?records/(\d+)/")
_OBSERVING_RUN = re.compile(r"^O\d+[a-z]?$")      # O4c, not data-release runs such as O4c1DiscC00


def _get_json(url: str):
    req = urllib.request.Request(url, headers={"User-Agent": "gwtc_analysis check_catalogs"})
    with urllib.request.urlopen(req, timeout=60) as r:
        return json.load(r)


def _paged(get: Callable, url: str) -> list:
    """All the results of a paged GWOSC ("results", "next") or Zenodo ("hits", "links.next") listing."""
    out = []
    while url:
        d = get(url)
        if isinstance(d, list):
            return out + d
        if "hits" in d:
            out += d["hits"].get("hits", [])
            url = (d.get("links") or {}).get("next")
        else:
            out += d.get("results", [])
            url = d.get("next")
    return out


def suggested_key(gwosc_list: str) -> str:
    """Catalog key for a GWOSC list: 'GWTC-6.0' -> 'GWTC-6', 'GWTC-4.1' -> 'GWTC-4.1'."""
    key = re.sub(r"-confident$", "", gwosc_list)
    return re.sub(r"\.0$", "", key)


def _zenodo_record(get: Callable, rid: str) -> dict:
    d = get(f"{ZENODO}/api/records/{rid}")
    versions = _paged(get, f"{ZENODO}/api/records/{rid}/versions?size=25&sort=version")
    latest = max(versions, key=lambda h: (h.get("metadata", {}).get("publication_date", ""), int(h["id"]))) if versions else d
    files = [f.get("key") for f in latest.get("files", [])]
    return dict(record=str(latest["id"]), concept=str(d.get("conceptrecid")), n_versions=len(versions) or 1,
                title=latest.get("metadata", {}).get("title", ""), n_files=len(files),
                skymap_tarball=next((f for f in files if f and _SKYMAP_TARBALL.search(f)), None))


def _discover_list(get: Callable, name: str, sample_events: int) -> dict:
    events = get(f"{GWOSC}/eventapi/jsonfull/{name}/").get("events", {})
    gps = [float(v["GPS"]) for v in events.values() if v.get("GPS") is not None]
    runs = Counter(reg.run_of(g) for g in gps)
    records = Counter()
    for key in list(events)[:sample_events]:
        try:
            params = _paged(get, f"{GWOSC}/api/v2/event-versions/{key}/parameters")
        except Exception:
            continue
        for p in params:
            if p.get("pipeline_type") == "pe" and p.get("data_url"):
                m = _RECORD_IN_URL.search(p["data_url"])
                if m:
                    records[m.group(1)] += 1
    concepts: dict[str, dict] = {}
    for rid in records:
        try:
            info = _zenodo_record(get, rid)
        except Exception as e:
            info = dict(record=rid, concept=f"? ({type(e).__name__})", n_versions=0, title="", n_files=0, skymap_tarball=None)
        concepts.setdefault(info["concept"], info)          # one entry per dataset: its latest version
    last_end = max(r.end_gps for r in reg.OBSERVING_RUNS.values())
    return dict(name=name, key=suggested_key(name), n_events=len(events),
                gps_range=(min(gps), max(gps)) if gps else None,
                n_after_known_runs=sum(g > last_end for g in gps),
                runs={str(k): v for k, v in runs.items()}, zenodo=list(concepts.values()))


def draft_entry(d: dict) -> str:
    """Draft `Catalog(...)` line for the registry."""
    known = [r for r in d["runs"] if r != "None"]
    runs = ", ".join(repr(r) for r in known) + ("," if len(known) == 1 else "")
    parts = []
    for z in d["zenodo"]:
        tar = f', "{z["skymap_tarball"]}"' if z.get("skymap_tarball") else ""
        parts.append(f'ZenodoRelease("{z["record"]}"{tar}),   # concept {z["concept"]}, {z["n_versions"]} version(s): '
                     f'{z["title"][:60]}')
    zen = ("\n            zenodo=(" + "\n                    ".join(parts) + "),") if parts else ""
    # a list whose runs are all those of an existing catalog re-analyses it
    same = [c.key for c in reg.CATALOGS.values() if known and set(known) <= set(c.runs) and not c.update_of]
    upd = f'\n            update_of="{same[0]}",   # same runs as {same[0]}: a re-analysis' if same else ""
    note = ("   # events after the last known run: add the new run to OBSERVING_RUNS and to this entry"
            if d.get("n_after_known_runs") else "")
    return f'Catalog("{d["key"]}", "{d["name"]}", ({runs}),{upd}{zen}\n            s3_prefix=None),{note}'


def check_catalogs(get: Callable = _get_json, sample_events: int = 3) -> dict:
    """What GWOSC and Zenodo publish that the registry does not describe."""
    lists = [c["name"] for c in _paged(get, f"{GWOSC}/api/v2/catalogs")]
    gwtc = [n for n in lists if re.match(r"^GWTC-\d", n)]
    new_lists = [n for n in gwtc if n not in reg.known_gwosc_lists()]
    missing_known = sorted(n for n in reg.known_gwosc_lists() if n not in lists and n not in reg.IGNORED_GWOSC_LISTS)
    gwosc_runs = {r["name"]: (r.get("gps_start"), r.get("gps_end")) for r in _paged(get, f"{GWOSC}/api/v2/runs")}
    new_runs = {k: v for k, v in gwosc_runs.items() if _OBSERVING_RUN.match(k) and k not in reg.OBSERVING_RUNS}
    newer = []
    for key, parts in reg.zenodo_releases().items():
        for part in parts:
            try:
                latest = _zenodo_record(get, part.record_id)
            except Exception:
                continue
            if latest["record"] != part.record_id:
                newer.append(dict(catalog=key, registry=part.record_id, latest=latest["record"],
                                  n_versions=latest["n_versions"]))
    return dict(new_lists=[_discover_list(get, n, sample_events) for n in new_lists], new_runs=new_runs,
                missing_lists=missing_known, newer_versions=newer)


def format_report(r: dict) -> str:
    lines = []
    if not (r["new_lists"] or r["new_runs"] or r["missing_lists"]):
        lines.append("The registry describes every GWTC event list and observing run published by GWOSC.")
    for d in r["new_lists"]:
        lo, hi = d["gps_range"] or (None, None)
        runs = ", ".join(f"{k} ({v})" for k, v in d["runs"].items() if k != "None")
        outside = d["runs"].get("None", 0)
        lines += [f"New GWOSC list {d['name']}: {d['n_events']} events, GPS {lo:.0f}-{hi:.0f}, runs {runs}"
                  + (f"; {outside} outside the known runs ({d['n_after_known_runs']} after the last one, the others in "
                     "engineering runs or gaps)" if outside else ""), "  draft registry entry:",
                  "    " + draft_entry(d).replace("\n", "\n    ")]
    for k, (a, b) in r["new_runs"].items():
        lines.append(f"New GWOSC run {k}: GPS {a}-{b}. The GWOSC run limits may include the engineering run before it: "
                     "use the official start of observing in OBSERVING_RUNS.")
    for n in r["missing_lists"]:
        lines.append(f"Registry list {n} is not published by GWOSC (renamed or removed?)")
    for v in r["newer_versions"]:
        lines.append(f"{v['catalog']}: the registry names record {v['registry']}, the latest of its "
                     f"{v['n_versions']} versions is {v['latest']} (used at run time; update the fallback if wished)")
    return "\n".join(lines)


def run_check_catalogs(out_json: Optional[str] = None, sample_events: int = 3) -> dict:
    r = check_catalogs(sample_events=sample_events)
    print(format_report(r))
    if out_json:
        with open(out_json, "w") as f:
            json.dump(r, f, indent=1)
    return r
