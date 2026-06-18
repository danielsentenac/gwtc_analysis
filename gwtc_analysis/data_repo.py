from __future__ import annotations

from pathlib import Path
from .repo_config import DEFAULT_REPO_CONFIG

def zenodo_skymap_spec(catalog_key: str):
    try:
        return DEFAULT_REPO_CONFIG.zenodo_skymaps[catalog_key]
    except KeyError as e:
        raise ValueError(
            f"No Zenodo skymap configured for catalog '{catalog_key}'"
        ) from e


def zenodo_skymap_url(catalog_key: str) -> str:
    spec = zenodo_skymap_spec(catalog_key)
    return (
        f"https://zenodo.org/records/{spec.record_id}/files/"
        f"{spec.filename}?download=1"
    )


def zenodo_cache_dir() -> Path:
    return DEFAULT_REPO_CONFIG.zenodo_cache


# ---------------------------------------------------------------------
# Public usegalaxy.org history (network data source)
# ---------------------------------------------------------------------

def galaxy_server_url() -> str:
    return DEFAULT_REPO_CONFIG.galaxy_server_url


def galaxy_history_id() -> str:
    return DEFAULT_REPO_CONFIG.galaxy_history_id


def galaxy_contents_url(history_id: str | None = None) -> str:
    """REST endpoint listing the datasets of a (public) Galaxy history."""
    hid = history_id or galaxy_history_id()
    return f"{galaxy_server_url()}/api/histories/{hid}/contents"


def galaxy_dataset_download_url(hda_id: str, ext: str = "data") -> str:
    """Stable, anonymous download URL for a single Galaxy dataset (HDA id)."""
    return f"{galaxy_server_url()}/api/datasets/{hda_id}/display?to_ext={ext}"

