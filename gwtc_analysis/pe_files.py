"""Which parameter-estimation (PE) file of a release is read, and where it is downloaded from.

The file depends only on the options: it is chosen in the Zenodo listing of the catalog release (its latest
version, or the one of --zenodo-version) by a rule that depends on the analysis, never on the data repository:

- ``SAMPLING`` (population, Hubble-constant, spin, lensing, counterpart, EOS and scalar-cloud modes): the samples
  under the prior they were drawn with, whose distance prior is recorded in the file: ``*_mixed_nocosmo.h5`` for
  GWTC-2.1 and GWTC-3, the single ``*-combined_PEDataRelease.hdf5`` from GWTC-4.0;
- ``CATALOG`` (parameters_estimation): the posteriors quoted in the catalog papers: ``*_mixed_cosmo.h5`` for
  GWTC-2.1 and GWTC-3 (reweighted to a source population uniform in comoving volume and source-frame time),
  the combined file from GWTC-4.0.

--data-repo only says where that file is downloaded from: the S3 mirror or Galaxy serve the file of the same
name when they have it, else it is downloaded from Zenodo. The release files reweighted to the cosmological
prior still record the sampling prior (PowerLaw, alpha=2): `distance_prior_desc` gives the prior their samples
follow.
"""
from __future__ import annotations

import json
import os
import threading
import time
from pathlib import Path
from typing import Callable, Iterable, Optional

SAMPLING = "sampling"
CATALOG = "catalog"

S3_BUCKET = "gwtc"
S3_LISTING_TTL_S = 24 * 3600
_record_ids: dict[str, list] = {}
_rids_lock = threading.Lock()

# the prior of the samples of a cosmologically reweighted release file (GWTC-2.1 and GWTC-3 *_cosmo.h5): checked on
# GW151226, whose nocosmo samples reweighted by this prior over D_L^2 give the distance quantiles of the cosmo file
COSMO_REWEIGHTED_PRIOR = ("UniformSourceFrame(cosmology='Planck15') [cosmologically reweighted release file: "
                          "the prior recorded in it is the sampling prior]")


def is_cosmo_reweighted(filename: str) -> bool:
    """Whether a release file holds samples reweighted to the cosmological prior (``*_cosmo.h5``)."""
    n = Path(str(filename)).name.lower()
    return "cosmo" in n and "nocosmo" not in n


def describe(filename: str) -> tuple[str, str]:
    """(short type, description) of a PE file, for the reports."""
    n = Path(str(filename)).name.lower()
    if "nocosmo" in n:
        return "nocosmo", "nocosmo: the samples under the prior they were drawn with (D_L^2), recorded in the file"
    if is_cosmo_reweighted(n):
        return "cosmo", ("cosmo: reweighted to a population uniform in comoving volume and source-frame time "
                         "(Planck15), the posteriors quoted in the catalog")
    if "combined" in n:
        return "combined", "combined: the single file of the release (GWTC-4.0 on), prior recorded in the file"
    if "unofficial" in n:
        return "unofficial bundle", "unofficial bundle rebuilt from the public GWTC-1 samples (build_unofficial_pe)"
    return "other", "a file given on the command line"


def distance_prior_desc(filename: str, recorded: str) -> str:
    """The luminosity-distance prior the samples of a release file follow: the recorded one, except for the
    cosmologically reweighted files."""
    return COSMO_REWEIGHTED_PRIOR if is_cosmo_reweighted(filename) else recorded


def _score(filename: str, kind: str) -> int:
    f = filename.lower()
    if kind == SAMPLING:
        return 0 if "nocosmo" in f else 1 if "combined" in f else 3 if is_cosmo_reweighted(f) else 2
    if kind == CATALOG:
        return (0 if "mixed_cosmo" in f else 1 if "recombined" in f else 2 if f.endswith(".hdf5")
                else 3 if f.endswith(".h5") else 4)
    raise ValueError(f"unknown PE file kind {kind!r}")


def pick(candidates: list[dict], kind: str = SAMPLING) -> Optional[dict]:
    """The file of an event's candidates ({'filename', ...}) that the rule `kind` reads (first on a tie)."""
    if not candidates:
        return None
    return sorted(candidates, key=lambda c: _score(str(c.get("filename", "")), kind))[0]


def select(names: Iterable[str], cache: Path | str, kind: str = SAMPLING, *, record_ids: Optional[list] = None,
           zenodo_versions: Optional[dict] = None) -> Optional[dict]:
    """The Zenodo entry {'record_id', 'filename', 'url'} of the first of `names` (an event's full name, its
    GPS-rounded variants, its common name) in the PE index of `cache`/index; None if no release has it."""
    from . import parameters_estimation as pe

    if record_ids is None:                 # resolved (and logged) once per process and versions
        key = json.dumps(zenodo_versions or {}, sort_keys=True)
        with _rids_lock:
            if key not in _record_ids:
                _record_ids[key] = pe.zenodo_pe_record_ids(zenodo_versions)
            record_ids = _record_ids[key]
    index = pe.build_zenodo_pe_index(cache_dir=str(Path(cache) / "index"), force_refresh=False,
                                     record_ids=record_ids)
    for nm in names:
        if nm and index.get(nm):
            return pick(index[nm], kind)
    return None


def event_entry(row, cache: Path | str, kind: str = SAMPLING) -> dict:
    """`select` for an event row (event, gps, common_name): its full name, the names rounded by +-2 s (as
    `hubble_constant`), then its common name."""
    from .hubble_constant import _full_name

    names = [row["event"]] + [_full_name(row["gps"] + d) for d in (-2, -1, 1, 2)] + [row.get("common_name") or ""]
    entry = select(names, cache, kind)
    if entry is None:
        raise ValueError(f"no Zenodo PE release found for {row['event']}")
    return entry


def extract_is_current(path: Path, filename: str) -> bool:
    """Whether an extract (samples kept from a PE file) was made from `filename`; extracts that do not record
    their source are not."""
    import h5py

    try:
        with h5py.File(path, "r") as h:
            src = h.attrs.get("source", h.attrs.get("source_file", ""))
    except OSError:
        return False
    src = src.decode() if isinstance(src, bytes) else str(src)
    return src == filename


_s3_lock = threading.Lock()


def s3_objects(cache: Path | str, data_repo: str = "s3") -> dict[str, str]:
    """File name -> object of the PE files on a mirror (the S3 bucket, or the Galaxy history seen as one), listed
    once a day into `cache`/index."""
    from .galaxy_store import store_client

    if data_repo == "galaxy":                 # the history listing has its own daily cache
        return {Path(o.object_name).name: o.object_name for o in store_client("galaxy").list_objects(S3_BUCKET)
                if "PEDataRelease" in o.object_name and o.object_name.endswith((".h5", ".hdf5"))}
    idx = Path(cache) / "index" / "s3_pe_files.json"
    with _s3_lock:
        if idx.exists() and time.time() - idx.stat().st_mtime < S3_LISTING_TTL_S:
            try:
                return json.loads(idx.read_text())
            except ValueError:
                pass
        out = {}
        for o in store_client("s3").list_objects(S3_BUCKET, recursive=True):
            n = o.object_name
            if "PEDataRelease" in n and n.endswith((".h5", ".hdf5")):
                out[Path(n).name] = n
        idx.parent.mkdir(parents=True, exist_ok=True)
        tmp = idx.with_name(f"{idx.name}.{os.getpid()}.{threading.get_ident()}.part")
        tmp.write_text(json.dumps(out))
        tmp.replace(idx)
        return out


def fetch(entry: dict, cache: Path | str, data_repo: Optional[str] = None, *, desc: str = "",
          log: Callable[[str], None] = print) -> tuple[Path, bool]:
    """The file of `entry` in `cache`/files (or in a shared PE cache), else downloaded: from the S3 mirror when
    `data_repo` is s3 and it has the file, else from Zenodo. Returns the path and whether this call downloaded
    it (a file already there belongs to whichever mode fetched it, and is kept)."""
    from . import net, parameters_estimation as pe

    from .mirror import data_repo as run_data_repo
    from .report import note_pe_source

    data_repo = data_repo or run_data_repo()      # the --data-repo of the run
    fname = entry["filename"]
    dest = Path(cache) / "files" / fname
    if dest.exists() and dest.stat().st_size > 0:
        note_pe_source(fname, "PE cache")
        return dest, False
    shared = pe._shared_pe_file(fname)
    if shared is not None:
        note_pe_source(fname, "PE cache of another mode")
        return shared, False
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    if data_repo in ("s3", "galaxy"):
        from .galaxy_store import repo_label, store_client

        label = repo_label(data_repo)
        try:
            obj = (s3_objects(cache) if data_repo == "s3" else s3_objects(cache, "galaxy")).get(fname)
            if obj is not None:
                pe._download_s3_with_progress(store_client(data_repo), bucket=S3_BUCKET, object_name=obj,
                                              out_path=tmp, desc=desc or fname)
                tmp.replace(dest)
                note_pe_source(fname, label)
                return dest, True
            log(f"{fname} is not on the {label}: downloading it from Zenodo")
            note_pe_source(fname, f"Zenodo (not on the {label})")
        except Exception as e:                 # the mirror unreachable or a failed transfer: Zenodo
            tmp.unlink(missing_ok=True)
            log(f"{label} failed for {fname} ({type(e).__name__}: {e}): downloading it from Zenodo")
            note_pe_source(fname, f"Zenodo ({label} failed)")
    else:
        note_pe_source(fname, "Zenodo")
    net.retry(lambda: pe._download_http_with_progress(entry["url"], tmp, desc=desc or fname), entry["url"], log=log)
    tmp.replace(dest)
    return dest, True
