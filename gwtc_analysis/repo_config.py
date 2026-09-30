from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Optional, Tuple

from .catalog_registry import ZenodoRelease, zenodo_releases  # noqa: E402  (ZenodoRelease re-exported)

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
    # catalog key -> Zenodo record series, from catalog_registry (newest first)
    zenodo_releases=zenodo_releases(),

    s3_bucket="gwtc",
    s3_prefix="PEDataRelease",

    galaxy_server_url="https://usegalaxy.org",
    galaxy_history_id="bbd44e69cb8906b53c9eaf6f96c3950e",
)
