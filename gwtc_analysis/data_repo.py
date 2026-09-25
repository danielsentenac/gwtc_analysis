from __future__ import annotations

import json
import re
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Optional

import requests

from .repo_config import DEFAULT_REPO_CONFIG, ZenodoRelease

ZENODO_API = "https://zenodo.org/api/records"
ZENODO_USER_AGENT = "gwtc_analysis (https://github.com/danielsentenac/gwtc_analysis)"
# Version listings are re-fetched after this delay, so a new release is picked up.
ZENODO_VERSIONS_TTL_S = 24 * 3600
# Skymap tarballs across releases: *-PESkyMaps, *-PESkyLocalizations,
# *-Archived_Skymaps, and a plain skymaps.tar.gz in GWTC-3 v1.
_SKYMAP_TARBALL_RE = re.compile(r"sky(maps|localizations)[^/]*\.tar\.gz$", re.IGNORECASE)


def zenodo_cache_dir() -> Path:
    return DEFAULT_REPO_CONFIG.zenodo_cache


# ---------------------------------------------------------------------
# Zenodo release versions
# ---------------------------------------------------------------------

@dataclass(frozen=True)
class ZenodoRecord:
    """One resolved version of a catalog's Zenodo release part."""
    catalog_key: str
    record_id: str
    version: Optional[int]      # 1 = oldest; None when unknown (offline fallback)
    n_versions: Optional[int]
    publication_date: str
    files: tuple[dict[str, Any], ...]   # {"key", "size", "checksum"}

    @property
    def label(self) -> str:
        if self.version is None:
            return f"record {self.record_id}"
        latest = " (latest)" if self.version == self.n_versions else ""
        return f"v{self.version}{latest}, record {self.record_id}"

    def file_url(self, filename: str) -> str:
        return f"https://zenodo.org/records/{self.record_id}/files/{filename}?download=1"


def zenodo_catalogs() -> list[str]:
    return list(DEFAULT_REPO_CONFIG.zenodo_releases)


def zenodo_release_parts(catalog_key: str) -> tuple[ZenodoRelease, ...]:
    try:
        return DEFAULT_REPO_CONFIG.zenodo_releases[catalog_key]
    except KeyError as e:
        raise ValueError(
            f"No Zenodo release configured for catalog '{catalog_key}'. "
            f"Catalogs with Zenodo releases: {', '.join(zenodo_catalogs())}"
        ) from e


def parse_zenodo_version(version: str | int | None) -> Optional[int]:
    """'latest'/None -> None; 'v2', '2' or 2 -> 2."""
    if version is None:
        return None
    s = str(version).strip().lower()
    if s in ("", "latest"):
        return None
    m = re.fullmatch(r"v?(\d+)", s)
    if not m or int(m.group(1)) < 1:
        raise ValueError(f"Invalid Zenodo version {version!r}: use 'latest', or v1, v2, ...")
    return int(m.group(1))


def _fetch_versions(record_id: str) -> list[dict[str, Any]]:
    """All versions of the series containing `record_id`, oldest first."""
    hits: list[dict[str, Any]] = []
    page = 1
    while True:
        r = requests.get(
            f"{ZENODO_API}/{record_id}/versions",
            params={"size": 25, "sort": "version", "page": page},  # 25 = anonymous page limit
            headers={"User-Agent": ZENODO_USER_AGENT},
            timeout=(10, 120),
        )
        r.raise_for_status()
        j = r.json()["hits"]
        hits.extend(j["hits"])
        if not j["hits"] or len(hits) >= int(j.get("total", 0)):
            break
        page += 1
    hits.reverse()  # sort=version lists the newest first
    return [
        {
            "record_id": str(h["id"]),
            "publication_date": h.get("metadata", {}).get("publication_date", ""),
            "files": [
                {"key": f.get("key"), "size": f.get("size"), "checksum": f.get("checksum")}
                for f in h.get("files", [])
            ],
        }
        for h in hits
    ]


def zenodo_release_versions(
    release: ZenodoRelease,
    *,
    cache_dir: str | Path | None = None,
    refresh: bool = False,
) -> list[dict[str, Any]]:
    """Versions of a release, oldest first, cached for ZENODO_VERSIONS_TTL_S.

    When zenodo.org is unreachable a stale cached listing is used; with no
    cache at all the request error is raised.
    """
    cache = Path(cache_dir or zenodo_cache_dir()) / f"versions_{release.record_id}.json"
    cached = None
    if cache.exists():
        try:
            cached = json.loads(cache.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            cached = None
    if cached and not refresh and time.time() - cached.get("fetched", 0) < ZENODO_VERSIONS_TTL_S:
        return cached["versions"]
    try:
        versions = _fetch_versions(release.record_id)
    except (requests.RequestException, KeyError, ValueError) as e:
        if cached:
            print(f"[zenodo] WARN: could not refresh versions of record {release.record_id} "
                  f"({type(e).__name__}); using the cached listing")
            return cached["versions"]
        raise
    cache.parent.mkdir(parents=True, exist_ok=True)
    cache.write_text(json.dumps({"fetched": time.time(), "versions": versions}), encoding="utf-8")
    return versions


def resolve_zenodo_records(
    catalog_key: str,
    version: str | int | None = None,
    *,
    cache_dir: str | Path | None = None,
) -> list[ZenodoRecord]:
    """Zenodo records of `catalog_key` (one per release part) at `version`.

    `version` is 'latest' (default) or vN, numbered from the oldest version (v1).
    If zenodo.org is unreachable and nothing is cached, the latest version falls
    back to the configured record; an explicit version then raises.
    """
    want = parse_zenodo_version(version)
    out: list[ZenodoRecord] = []
    for part in zenodo_release_parts(catalog_key):
        try:
            versions = zenodo_release_versions(part, cache_dir=cache_dir)
        except (requests.RequestException, KeyError, ValueError) as e:
            if want is not None:
                raise ValueError(
                    f"Cannot list the Zenodo versions of {catalog_key} (record {part.record_id}): {e}"
                ) from e
            print(f"[zenodo] WARN: zenodo.org unreachable ({type(e).__name__}); "
                  f"falling back to configured record {part.record_id} for {catalog_key}")
            files = ({"key": part.skymap_filename},) if part.skymap_filename else ()
            out.append(ZenodoRecord(catalog_key, part.record_id, None, None, "", files))
            continue
        n = len(versions)
        idx = n if want is None else want
        if not 1 <= idx <= n:
            listing = ", ".join(f"v{i} (record {v['record_id']}, {v['publication_date']})"
                                for i, v in enumerate(versions, 1))
            raise ValueError(f"{catalog_key} has no Zenodo version v{want}. Available: {listing}")
        v = versions[idx - 1]
        out.append(ZenodoRecord(catalog_key, v["record_id"], idx, n,
                                v["publication_date"], tuple(v["files"])))
    return out


def zenodo_skymap_tarball(
    catalog_key: str, version: str | int | None = None
) -> tuple[ZenodoRecord, str]:
    """(record, filename) of the skymap tarball of `catalog_key` at `version`."""
    for rec in resolve_zenodo_records(catalog_key, version):
        for f in rec.files:
            key = f.get("key") or ""
            if _SKYMAP_TARBALL_RE.search(key):
                return rec, key
    raise ValueError(f"No skymap tarball found in the Zenodo release of {catalog_key} (version {version or 'latest'})")


def zenodo_skymap_url(catalog_key: str, version: str | int | None = None) -> str:
    rec, filename = zenodo_skymap_tarball(catalog_key, version)
    return rec.file_url(filename)


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

