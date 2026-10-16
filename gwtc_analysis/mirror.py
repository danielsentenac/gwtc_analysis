"""The mirrors (S3 bucket, Galaxy history) of the public files (other than the PE files and skymaps of the releases).

Every public file the modes download (search-sensitivity injections, the area-law release, the PSD and skymap
supplements, the continuous-wave limits, the GW170817 GWTC-1 products) is chosen as without the mirror: same URL,
and for a Zenodo record resolved at run time the same (latest) version. With --data-repo s3 the file of that URL
is downloaded from the `gwtc` bucket when the bucket holds it, under a key derived from the URL (with
--data-repo galaxy, from the public Galaxy history, where the dataset is named after the key):

    https://zenodo.org/records/<record>/files/<name>            -> public/zenodo/<record>/<name>
    https://dcc.ligo.org/public/<nnnn>/<id>/<version>/<name>    -> public/dcc/<id>/<version>/<name>
    https://arxiv.org/e-print/<id>                              -> public/arxiv/<id>.tar.gz

else from its public source. API queries (GWOSC strain and event lists, VizieR, HEASARC) are not mirrored.
"""
from __future__ import annotations

import re
from pathlib import Path
from typing import Callable, Optional
from urllib.parse import urlparse

from . import net

BUCKET = "gwtc"
PREFIX = "public/"

_DATA_REPO = "zenodo"


def set_data_repo(data_repo: Optional[str]) -> None:
    """The repository of the run (--data-repo): s3 makes the downloads below try the mirror first."""
    global _DATA_REPO
    _DATA_REPO = data_repo or "zenodo"


def data_repo() -> str:
    return _DATA_REPO


def mirror_key(url: str) -> Optional[str]:
    """The key of a public file in the bucket, or None for a URL that is not mirrored."""
    u = urlparse(url)
    parts = [p for p in u.path.split("/") if p]
    if u.netloc.endswith("zenodo.org") and len(parts) >= 4 and parts[0] == "records" and parts[2] == "files":
        return f"{PREFIX}zenodo/{parts[1]}/{'/'.join(parts[3:])}"
    if u.netloc == "dcc.ligo.org" and len(parts) >= 5 and parts[0] == "public":
        return f"{PREFIX}dcc/{parts[2]}/{parts[3]}/{'/'.join(parts[4:])}"
    if u.netloc.endswith("arxiv.org") and len(parts) == 2 and parts[0] == "e-print":
        return f"{PREFIX}arxiv/{parts[1]}.tar.gz"
    return None


MIRRORS = ("s3", "galaxy")
BACKUP = "galaxy"     # read when a public source fails, whatever the run's --data-repo


def _from_repo(repo: str, url: str, dest: Path, log: Callable[[str], None], *, otherwise: str) -> bool:
    """Download the copy of `url` held by the mirror `repo` (the S3 bucket, or the Galaxy history, which holds the
    file as a dataset named after its key) to `dest`; False when the mirror lacks it or fails (`otherwise` says
    what happens then, for the log)."""
    key = mirror_key(url)
    if key is None:
        return False
    from .galaxy_store import repo_label, store_client
    from .parameters_estimation import _download_s3_with_progress

    try:
        client = store_client(repo)
        client.stat_object(BUCKET, key)
    except Exception:
        log(f"{key} is not on the {repo_label(repo)} (or it is unreachable): {otherwise}")
        return False
    try:
        dest.parent.mkdir(parents=True, exist_ok=True)
        _download_s3_with_progress(client, bucket=BUCKET, object_name=key, out_path=dest, desc=Path(key).name)
    except Exception as e:                     # a failed transfer
        Path(dest).unlink(missing_ok=True)
        log(f"{repo_label(repo)} failed for {key} ({type(e).__name__}: {e}): {otherwise}")
        return False
    log(f"downloaded from the {repo_label(repo)}: {key}")
    return True


def _from_s3(url: str, dest: Path, log: Callable[[str], None]) -> bool:
    """Download the run's mirror copy of `url` to `dest`; False when the run uses no mirror or the mirror lacks it."""
    if _DATA_REPO not in MIRRORS:
        return False
    return _from_repo(_DATA_REPO, url, dest, log, otherwise=f"downloading it from {urlparse(url).netloc}")


def _from_backup(url: str, dest: Path, log: Callable[[str], None], error: BaseException) -> bool:
    """After the public source of `url` failed with a transient `error` (retries exhausted): the copy in the public
    Galaxy history, whatever --data-repo is (unless the run already tried that history first)."""
    if _DATA_REPO == BACKUP or mirror_key(url) is None or not net.transient(error):
        return False
    from .galaxy_store import repo_label

    log(f"⚠️ {urlparse(url).netloc} failed ({type(error).__name__}): trying the copy on the {repo_label(BACKUP)}")
    return _from_repo(BACKUP, url, dest, log, otherwise="no copy to fall back on")


def on_mirror(url: str) -> bool:
    """Whether the run uses a mirror (s3, galaxy) and it holds the file of `url`."""
    if _DATA_REPO not in MIRRORS or mirror_key(url) is None:
        return False
    from .galaxy_store import store_client

    try:
        store_client(_DATA_REPO).stat_object(BUCKET, mirror_key(url))
        return True
    except Exception:
        return False


def _label() -> str:
    from .galaxy_store import repo_label

    return repo_label(_DATA_REPO)


def _record(url: str, source: str) -> None:
    from .report import record_data_file

    record_data_file(re.sub(r"\?download=1$", "", url), source)


def _source(url: str) -> str:
    """The report's source of a file downloaded from its public source."""
    return urlparse(url).netloc + (f" (not on the {_label()})" if _DATA_REPO in MIRRORS and mirror_key(url) else "")


def _backup_source(url: str) -> str:
    from .galaxy_store import repo_label

    return f"{repo_label(BACKUP)} ({urlparse(url).netloc} failed)"


def fetch_to(url: str, dest: Path, http_download: Callable[[str, Path], None], *,
             log: Callable[[str], None] = print) -> Path:
    """`url` downloaded to `dest` (written whole or not at all): from the run's mirror (s3, galaxy) when it holds
    it, else by `http_download(url, path)`, retried after a transient failure (net.retry); when the source still
    fails, from the copy in the public Galaxy history (BACKUP)."""
    dest = Path(dest)
    tmp = dest.with_name(dest.name + ".part")
    if _from_s3(url, tmp, log):
        _record(url, _label())
    else:
        try:
            net.retry(lambda: http_download(url, tmp), url, log=log)
        except Exception as e:
            tmp.unlink(missing_ok=True)
            if not _from_backup(url, tmp, log, e):
                raise
            _record(url, _backup_source(url))
        else:
            _record(url, _source(url))
    tmp.replace(dest)
    return dest


def fetch_bytes(url: str, http_get: Callable[[str], bytes], *, log: Callable[[str], None] = print,
                tmp_dir: Optional[Path] = None) -> bytes:
    """The content of `url`: from the run's mirror when it holds it, else `http_get(url)` (retried, then the Galaxy
    copy, as in fetch_to)."""
    import tempfile

    if _DATA_REPO in MIRRORS and mirror_key(url):
        with tempfile.TemporaryDirectory(dir=tmp_dir) as d:
            p = Path(d) / "file"
            if _from_s3(url, p, log):
                _record(url, _label())
                return p.read_bytes()
    try:
        data = net.retry(lambda: http_get(url), url, log=log)
    except Exception as e:
        with tempfile.TemporaryDirectory(dir=tmp_dir) as d:
            p = Path(d) / "file"
            if not _from_backup(url, p, log, e):
                raise
            _record(url, _backup_source(url))
            return p.read_bytes()
    _record(url, _source(url))
    return data


def cached(url: str) -> None:
    """Record, for the report, a public file used from a local cache."""
    _record(url, "local cache")
