from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Optional, Tuple

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
class RepoConfig:
    galaxy_base: Path
    zenodo_cache: Path
    # catalog key -> release parts (GWTC-5.0 is split across two record series)
    zenodo_releases: Dict[str, Tuple[ZenodoRelease, ...]]
    s3_bucket: str
    s3_prefix: str
    # Public usegalaxy.org published history used as a network data source
    # (download links resolved via the Galaxy REST API).
    galaxy_server_url: str
    galaxy_history_id: str

DEFAULT_REPO_CONFIG = RepoConfig(
    galaxy_base=Path("/data/gwtc_analysis"),

    zenodo_cache=Path.home() / ".cache_gwtc_analysis" / "zenodo",

    # The version actually read is resolved at run time: the latest one by
    # default, or the one requested with --zenodo-version (see data_repo.py).
    zenodo_releases={
        "GWTC-5": (
            ZenodoRelease(  # Part 1 of 2 (concept 20276105), carries the skymaps
                record_id="20348005",
                skymap_filename="IGWN-GWTC5p0-29ebe06b7_25-Archived_Skymaps.tar.gz",
            ),
            ZenodoRelease(record_id="20348006"),  # Part 2 of 2 (concept 20291739)
        ),
        "GWTC-4": (
            ZenodoRelease(  # GWTC-4.0 (concept 16053483)
                record_id="17602505",
                skymap_filename="IGWN-GWTC4p0-38214bd95_724-Archived_Skymaps.tar.gz",
            ),
        ),
        "GWTC-3": (
            ZenodoRelease(  # concept 5546662
                record_id="22685054",
                skymap_filename="IGWN-GWTC3p0-v3-PESkyLocalizations.tar.gz",
            ),
        ),
        "GWTC-2.1": (
            ZenodoRelease(  # concept 5117702
                record_id="6513631",
                skymap_filename="IGWN-GWTC2p1-v2-PESkyMaps.tar.gz",
            ),
        ),
    },

    s3_bucket="gwtc",
    s3_prefix="PEDataRelease",

    galaxy_server_url="https://usegalaxy.org",
    galaxy_history_id="bbd44e69cb8906b53c9eaf6f96c3950e",
)
