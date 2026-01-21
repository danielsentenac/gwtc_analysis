from dataclasses import dataclass
from pathlib import Path
from typing import Dict

@dataclass(frozen=True)
class ZenodoSkymapSpec:
    record_id: str
    filename: str

@dataclass(frozen=True)
class RepoConfig:
    galaxy_base: Path
    zenodo_cache: Path
    zenodo_skymaps: Dict[str, ZenodoSkymapSpec]
    s3_bucket: str
    s3_prefix: str

DEFAULT_REPO_CONFIG = RepoConfig(
    galaxy_base=Path("/data/gwtc_analysis"),

    zenodo_cache=Path.home() / ".cache_gwtc_analysis" / "zenodo",

    zenodo_skymaps={
        # Former ZENODO_SKYMAP_TARBALLS entries (lossless refactor)
        "GWTC-4": ZenodoSkymapSpec(
            record_id="17014085",
            filename="IGWN-GWTC4p0-1a206db3d_721-Archived_Skymaps.tar.gz",
        ),
        "GWTC-3": ZenodoSkymapSpec(
            record_id="8177023",
            filename="IGWN-GWTC3p0-v2-PESkyLocalizations.tar.gz",
        ),
        "GWTC-2.1": ZenodoSkymapSpec(
            record_id="6513631",
            filename="IGWN-GWTC2p1-v2-PESkyMaps.tar.gz",
        ),
    },

    s3_bucket="gwtc",
    s3_prefix="PEDataRelease",
)

