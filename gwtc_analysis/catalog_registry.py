"""The gravitational-wave transient catalogs, observing runs and sensitivity releases known to gwtc_analysis.

This is the single place to edit when the LVK publishes a new catalog: add its observing run(s) to
``OBSERVING_RUNS``, an entry to ``CATALOGS`` and, if it comes with injections, one to
``SENSITIVITY_RELEASES``. Every other table of the package (catalog keys, GWOSC list names, Zenodo
records, S3 prefixes, run mappings, the `ALL` expansions, the rates and hubble_constant releases) is
derived from these three.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Iterable, Optional


@dataclass(frozen=True)
class ZenodoRelease:
    """One Zenodo record series (all the versions of a data release).

    ``record_id`` may be any version of the series: it is used to list the
    versions through the Zenodo API, and as the fallback record when zenodo.org
    is unreachable and no version listing is cached. ``skymap_filename`` is the
    skymap tarball of that fallback record (None if the series has no skymaps).
    """
    record_id: str
    skymap_filename: Optional[str] = None


@dataclass(frozen=True)
class ObservingRun:
    """An observing run: GWOSC boundaries (GPS) and how its injections detect signals."""
    name: str
    start_gps: float
    end_gps: float
    semi_analytic: bool = False     # injections found from the semi-analytic SNR (O1, O2), not from searches


@dataclass(frozen=True)
class Catalog:
    """A catalog key of the command line.

    ``gwosc_list``: the GWOSC event list of its confident events; ``marginal_list``: its marginal candidates,
    if GWOSC publishes them; ``runs``: the observing runs of its new events; ``zenodo``: the Zenodo record
    series of its PE files and skymaps (several parts for a split release); ``products_from``: the catalog
    whose releases hold its PE files and skymaps when it has none of its own (GWTC-1 → GWTC-2.1);
    ``s3_prefix``: its folder in the S3 bucket; ``update_of``: the catalog it re-analyses (GWTC-4.1 →
    GWTC-4). An update covers the same events again, so it is used only when named explicitly: it is left
    out of ALL, of the default event lists and releases, and of the default PE index, which keeps the
    results of the original catalog unchanged.
    """
    key: str
    gwosc_list: str
    runs: tuple[str, ...]
    marginal_list: Optional[str] = None
    zenodo: tuple[ZenodoRelease, ...] = ()
    products_from: Optional[str] = None
    s3_prefix: Optional[str] = None
    update_of: Optional[str] = None

    @property
    def has_skymaps(self) -> bool:
        return any(p.skymap_filename for p in self.zenodo)


@dataclass(frozen=True)
class SensitivityRelease:
    """A cumulative LVK search-sensitivity release (injections) on Zenodo.

    ``real_file_re``: its real-injection mixture (O3 onward), used by `rates`; ``semi_file_re``: its mixture
    with the semi-analytic O1+O2 injections, used by `rates` with GWTC-1 and by `hubble_constant`;
    ``published``: spectral-siren H0 values of the matching cosmology paper, per mass model.
    """
    key: str
    record: str
    label: str
    semi_label: str
    runs: tuple[str, ...]
    real_file_re: str
    semi_file_re: str
    published: dict = field(default_factory=dict)


# ---------------------------------------------------------------------------
# The registry
# ---------------------------------------------------------------------------

OBSERVING_RUNS: dict[str, ObservingRun] = {r.name: r for r in (
    ObservingRun("O1", 1126051217, 1137254417, semi_analytic=True),
    ObservingRun("O2", 1164556817, 1187733618, semi_analytic=True),
    ObservingRun("O3a", 1238166018, 1253977218),
    ObservingRun("O3b", 1256655618, 1269363618),
    ObservingRun("O4a", 1368975618, 1389456018),
    ObservingRun("O4b", 1396796418, 1422118818),
)}

# in time order; the Zenodo record IDs are starting points, the version read is resolved at run time
CATALOGS: dict[str, Catalog] = {c.key: c for c in (
    Catalog("GWTC-1", "GWTC-1-confident", ("O1", "O2"), products_from="GWTC-2.1"),
    Catalog("GWTC-2.1", "GWTC-2.1-confident", ("O3a",), marginal_list="GWTC-2.1-marginal",
            zenodo=(ZenodoRelease("6513631", "IGWN-GWTC2p1-v2-PESkyMaps.tar.gz"),),   # concept 5117702
            s3_prefix="GWTC-2.1/"),
    Catalog("GWTC-3", "GWTC-3-confident", ("O3b",), marginal_list="GWTC-3-marginal",
            zenodo=(ZenodoRelease("22685054", "IGWN-GWTC3p0-v3-PESkyLocalizations.tar.gz"),),   # concept 5546662
            s3_prefix="GWTC-3/"),
    Catalog("GWTC-4", "GWTC-4.0", ("O4a",),
            zenodo=(ZenodoRelease("17602505", "IGWN-GWTC4p0-38214bd95_724-Archived_Skymaps.tar.gz"),),  # concept 16053483
            s3_prefix="GWTC-4/"),
    Catalog("GWTC-4.1", "GWTC-4.1", ("O4a",), update_of="GWTC-4",      # O4a re-analysed with the GWTC-5.0 methods
            zenodo=(ZenodoRelease("20275769", "IGWN-GWTC4p1-18965dda8_5-Archived_Skymaps.tar.gz"),)),  # concept 20275768
    Catalog("GWTC-5", "GWTC-5.0", ("O4b",),
            zenodo=(ZenodoRelease("20348005", "IGWN-GWTC5p0-29ebe06b7_25-Archived_Skymaps.tar.gz"),  # part 1 of 2 (concept 20276105)
                    ZenodoRelease("20348006")),                                                     # part 2 of 2 (concept 20291739)
            s3_prefix="GWTC-5/"),
)}

SENSITIVITY_RELEASES: dict[str, SensitivityRelease] = {s.key: s for s in (
    SensitivityRelease(
        "gwtc5", "19500052", "GWTC-5.0 cumulative, real O3 + O4a + O4b injections (~900 MB)",
        "GWTC-5.0 cumulative, semi-analytic O1+O2 + real O3+O4a+O4b injections",
        ("O1", "O2", "O3a", "O3b", "O4a", "O4b"),
        real_file_re=r"^mixture-real_.*cartesian_spins.*\.hdf5?$",
        semi_file_re=r"^mixture-semi_o1_o2-real_o3_o4a_o4b-cartesian_spins.*\.hdf5?$",
        # spectral sirens, GWTC-5.0 cosmology paper (arXiv:2605.27227), Table 11; it has no PLP analysis
        published={
            "mltp": dict(ref="GWTC-5.0 cosmology, arXiv:2605.27227 (MLTP)", median=71.0, plus=21.0, minus=17.5,
                         lo90=44.5, hi90=107.0),
        }),
    SensitivityRelease(
        "gwtc4", "16740128", "GWTC-4.0 cumulative, real O3 + O4a injections (~400 MB)",
        "GWTC-4.0 cumulative, semi-analytic O1+O2 + real O3+O4a injections",
        ("O1", "O2", "O3a", "O3b", "O4a"),
        real_file_re=r"^mixture-real_.*cartesian_spins.*\.hdf5?$",
        semi_file_re=r"^mixture-semi_o1_o2-real_o3_o4a-cartesian_spins.*\.hdf5?$",
        # spectral sirens, GWTC-4.0 cosmology paper v3 (the published version); 90% bounds from the quoted intervals
        published={
            "plp": dict(ref="GWTC-4.0 cosmology, arXiv:2509.04348 (PLP)", median=105.5, plus=46.4, minus=35.8,
                        lo90=50.5, hi90=176.1),
            "mltp": dict(ref="GWTC-4.0 cosmology, arXiv:2509.04348 (MLTP)", median=72.3, plus=42.5, minus=25.6,
                         lo90=34.2, hi90=154.1),
        }),
)}

# GWOSC event lists deliberately not used, with the reason (check_catalogs does not report them)
IGNORED_GWOSC_LISTS: dict[str, str] = {
    "GWTC": "cumulative list of the confident events of all catalogs",
    "GWTC-2": "superseded by GWTC-2.1",
    "GWTC-2.1-auxiliary": "GWTC-2 candidates below the GWTC-2.1 thresholds",
    "GWTC-1-marginal": "marginal GWTC-1 candidates, none with FAR below 1 per year",
}

DEFAULT_RATES_RELEASE = "gwtc5"       # the latest release
DEFAULT_H0_RELEASE = "gwtc4"          # the release of the reproduced, published analysis

# Date of the last `gwtc_analysis check_catalogs` finding nothing new: update it with each registry check
REGISTRY_CHECKED = "2026-09-30"


# ---------------------------------------------------------------------------
# Derived views
# ---------------------------------------------------------------------------

def catalog_keys() -> tuple[str, ...]:
    return tuple(CATALOGS)


def default_catalog_keys() -> tuple[str, ...]:
    """The catalogs ALL stands for: every catalog except the updates of another one."""
    return tuple(c.key for c in CATALOGS.values() if not c.update_of)


def expand_all(catalogs: Iterable[str]) -> list[str]:
    """Catalog keys with ALL replaced by the default catalogs, duplicates removed, order kept."""
    out: list[str] = []
    for c in catalogs or []:
        for k in (default_catalog_keys() if c == "ALL" else (c,)):
            if k not in out:
                out.append(k)
    return out


def update_catalogs(catalogs: Iterable[str]) -> tuple[str, ...]:
    """The update catalogs among some catalog keys."""
    return tuple(k for k in catalogs or () if k in CATALOGS and CATALOGS[k].update_of)


def allowed_catalogs() -> tuple[str, ...]:
    """The catalog keys of the command line, and ALL."""
    return catalog_keys() + ("ALL",)


def catalog_help(examples: Optional[Iterable[str]] = None) -> str:
    return " ".join(examples or catalog_keys())


def update_help() -> str:
    """'GWTC-4.1, update of GWTC-4' for the help texts."""
    return "; ".join(f"{c.key}, update of {c.update_of}" for c in CATALOGS.values() if c.update_of)


def catalog_name(key: str) -> str:
    """Published name of a catalog: 'GWTC-4' -> 'GWTC-4.0', 'GWTC-1' -> 'GWTC-1'."""
    return CATALOGS[key].gwosc_list.removesuffix("-confident")


def latest_catalog() -> Catalog:
    """The most recent catalog (not an update) of the registry."""
    return [c for c in CATALOGS.values() if not c.update_of][-1]


def coverage_text(version: Optional[str] = None, markdown: bool = False) -> str:
    """The catalogs covered by this version of the package, for --version, the README and the docs."""
    b = (lambda t: f"**{t}**") if markdown else (lambda t: t)
    keys = default_catalog_keys()
    updates = [c.key for c in CATALOGS.values() if c.update_of]
    runs = list(OBSERVING_RUNS)
    last = latest_catalog()
    names = ", ".join(catalog_name(k) for k in keys) + "".join(
        f", and the update {catalog_name(u)} (of {catalog_name(CATALOGS[u].update_of)})" for u in updates)
    head = f"Catalogs in version {version}" if version else "Catalogs covered"
    return (f"{b(head)}: {names}, observing runs {runs[0]} to {runs[-1]}. "
            f"{b(f'Latest catalog: {catalog_name(last.key)} ({chr(43).join(last.runs)})')}. "
            f"Registry checked against GWOSC and Zenodo on {REGISTRY_CHECKED}; catalogs published later need a "
            f"newer version of the package (`gwtc_analysis check_catalogs` tells whether GWOSC has published one).")


def gwosc_list(key: str) -> str:
    """GWOSC event list of a catalog key's confident events (keys already naming a list are returned as is)."""
    return CATALOGS[key].gwosc_list if key in CATALOGS else key


def gwosc_aliases() -> dict[str, str]:
    """Catalog key → GWOSC confident list, newest first."""
    return {c.key: c.gwosc_list for c in reversed(list(CATALOGS.values()))}


def gwosc_lists(catalogs: Optional[Iterable[str]] = None, marginal: bool = True) -> tuple[str, ...]:
    """GWOSC lists (confident, and marginal ones if published) of some catalogs (default: the default catalogs,
    without the updates), in time order."""
    keys = set(default_catalog_keys() if catalogs is None else catalogs)
    out = []
    for c in CATALOGS.values():
        if c.key not in keys:
            continue
        out.append(c.gwosc_list)
        if marginal and c.marginal_list:
            out.append(c.marginal_list)
    return tuple(out)


def zenodo_releases() -> dict[str, tuple[ZenodoRelease, ...]]:
    """Catalog key → Zenodo record series of its PE files and skymaps, newest first (updates included)."""
    return {c.key: c.zenodo for c in reversed(list(CATALOGS.values())) if c.zenodo}


def is_update(key: str) -> bool:
    return key in CATALOGS and bool(CATALOGS[key].update_of)


def products_catalog(key: str) -> str:
    """The catalog whose releases hold a catalog's PE files and skymaps (GWTC-1 → GWTC-2.1)."""
    c = CATALOGS.get((key or "").strip())
    return c.products_from if c and c.products_from else (key or "").strip()


def skymap_catalogs() -> list[str]:
    """Default catalogs with skymap releases, in time order (what ALL means for skymap searches)."""
    return [c.key for c in CATALOGS.values() if c.has_skymaps and not c.update_of]


def s3_prefix(key: str) -> str:
    for c in CATALOGS.values():
        if key.startswith(c.key) and c.s3_prefix:
            return c.s3_prefix
    return f"{key}/"


def observing_runs() -> dict[str, tuple[float, float]]:
    return {r.name: (r.start_gps, r.end_gps) for r in OBSERVING_RUNS.values()}


def catalog_runs_map() -> dict[str, tuple[str, ...]]:
    return {c.key: c.runs for c in CATALOGS.values()}


def semi_analytic_runs() -> tuple[str, ...]:
    return tuple(r.name for r in OBSERVING_RUNS.values() if r.semi_analytic)


def release_runs() -> dict[str, tuple[str, ...]]:
    return {s.key: s.runs for s in SENSITIVITY_RELEASES.values()}


def release_catalogs(release: str) -> tuple[str, ...]:
    """Catalog keys whose runs are all covered by a sensitivity release."""
    runs = set(SENSITIVITY_RELEASES[release].runs)
    return tuple(c.key for c in CATALOGS.values() if set(c.runs) <= runs and not c.update_of)


def _span(runs: Iterable[str]) -> str:
    runs = list(runs)
    return runs[0] if len(runs) == 1 else f"{runs[0]}-{runs[-1]}"


def catalog_runs_help() -> str:
    """'GWTC-1: O1-O2, GWTC-2.1: O3a, ...' for the help texts."""
    return ", ".join(f"{c.key}: {_span(c.runs)}" for c in CATALOGS.values())


def release_runs_help() -> str:
    """'gwtc4: O1-O4a, gwtc5: O1-O4b' for the help texts."""
    return "; ".join(f"{k}: {_span(r.runs)}" for k, r in sorted(SENSITIVITY_RELEASES.items()))


def known_gwosc_lists() -> set[str]:
    """GWOSC lists described by the registry or deliberately ignored."""
    out = set(IGNORED_GWOSC_LISTS)
    for c in CATALOGS.values():
        out.add(c.gwosc_list)
        if c.marginal_list:
            out.add(c.marginal_list)
    return out


def run_of(gps: float) -> Optional[str]:
    """Observing run containing a GPS time, or None (engineering runs, gaps)."""
    return next((r.name for r in OBSERVING_RUNS.values() if r.start_gps <= float(gps) <= r.end_gps), None)
