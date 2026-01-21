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

