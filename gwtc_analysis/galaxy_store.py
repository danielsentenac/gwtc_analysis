"""The public usegalaxy.org GWTC history seen as an S3 bucket, so that --data-repo galaxy works as --data-repo s3.

The history holds the same files as the `gwtc` bucket, as datasets named after the files: the PE files and skymaps
of the releases (IGWN-GWTC<x>p<y>-...) and the other public files (dataset name: their key in the bucket,
public/zenodo/<record>/<name>, ...). This module lists the datasets (anonymously, cached for a day) and serves
them through the three calls the S3 code uses: list_objects, stat_object, get_object. A release file appears
under the object name <catalog>/<name> (PE files <catalog>/PE/<name>), as in the bucket.
"""
from __future__ import annotations

import json
import os
import threading
import time
from pathlib import Path
from types import SimpleNamespace
from typing import Optional

LISTING_TTL_S = 24 * 3600
_lock = threading.Lock()


def _listing_path() -> Path:
    from .data_repo import zenodo_cache_dir

    return Path(zenodo_cache_dir()).parent / "galaxy" / "history_files.json"


def history_files(refresh: bool = False) -> dict[str, dict]:
    """Dataset name -> {'id', 'ext', 'size'} of the live datasets of the public GWTC history (the last one of a
    name wins: duplicates are copies)."""
    from . import net
    from .data_repo import galaxy_contents_url

    p = _listing_path()
    with _lock:
        if not refresh and p.exists() and time.time() - p.stat().st_mtime < LISTING_TTL_S:
            try:
                return json.loads(p.read_text())
            except ValueError:
                pass
        r = net.get(galaxy_contents_url(), params={"v": "dev", "keys": "id,name,extension,file_size,deleted,"
                                                                          "purged,state,history_content_type"},
                    timeout=(10, 300))
        r.raise_for_status()
        out = {it["name"]: {"id": it["id"], "ext": it.get("extension") or "data", "size": it.get("file_size")}
               for it in r.json()
               if it.get("history_content_type") == "dataset" and not it.get("deleted") and not it.get("purged")
               and it.get("state") == "ok" and it.get("name")}
        p.parent.mkdir(parents=True, exist_ok=True)
        tmp = p.with_name(f"{p.name}.{os.getpid()}.{threading.get_ident()}.part")
        tmp.write_text(json.dumps(out))
        tmp.replace(p)
        return out


def object_name(dataset_name: str) -> Optional[str]:
    """The bucket object name of a history dataset: public/... as is, release files under their catalog."""
    if dataset_name.startswith("public/"):
        return dataset_name
    from .skymap3d import catalog_of_pe_file

    cat = catalog_of_pe_file(dataset_name)
    if cat is None:
        return None
    from .catalog_registry import s3_prefix

    is_pe = "PEDataRelease" in dataset_name and dataset_name.endswith((".h5", ".hdf5"))
    return f"{s3_prefix(cat)}{'PE/' if is_pe else ''}{dataset_name}"


class _Response:
    """A streamed HTTP download with the read/close/release_conn calls of a MinIO response."""

    def __init__(self, url: str):
        from . import net

        self._r = net.get(url, stream=True, timeout=(10, 600))
        self._r.raise_for_status()

    def read(self, n: int = -1) -> bytes:
        if n is None or n < 0:
            return self._r.content
        return self._r.raw.read(n, decode_content=True)

    def close(self) -> None:
        self._r.close()

    def release_conn(self) -> None:
        pass


class GalaxyHistoryStore:
    """list_objects / stat_object / get_object over the history datasets (the bucket argument is ignored)."""

    def __init__(self, refresh: bool = False):
        self._files = history_files(refresh)
        self._by_object = {o: n for n in self._files if (o := object_name(n))}

    def list_objects(self, bucket_name=None, prefix: str = "", recursive: bool = True, **kw):
        prefix = kw.get("prefix", prefix) or ""
        for o, n in sorted(self._by_object.items()):
            if o.startswith(prefix):
                yield SimpleNamespace(object_name=o, size=self._files[n].get("size"))

    def stat_object(self, bucket_name, object_name):
        n = self._by_object.get(object_name)
        if n is None:
            raise KeyError(f"{object_name} is not in the Galaxy history")
        return SimpleNamespace(object_name=object_name, size=self._files[n].get("size"))

    def get_object(self, bucket_name=None, object_name=None, **kw):
        from .data_repo import galaxy_dataset_download_url

        object_name = kw.get("object_name", object_name)
        n = self._by_object.get(object_name)
        if n is None:
            raise KeyError(f"{object_name} is not in the Galaxy history")
        f = self._files[n]
        return _Response(galaxy_dataset_download_url(f["id"], f["ext"]))


def store_client(data_repo: str):
    """The object store of a mirror repository: the MinIO bucket (s3) or the Galaxy history (galaxy)."""
    if data_repo == "galaxy":
        return GalaxyHistoryStore()
    from .gw_stat import _s3_client_from_env

    return _s3_client_from_env()


def repo_label(data_repo: str) -> str:
    return {"s3": "S3 mirror", "galaxy": "Galaxy history"}.get(data_repo, data_repo)
