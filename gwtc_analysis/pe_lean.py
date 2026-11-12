"""A PE file read with pesummary one analysis at a time.

pesummary reads the posterior samples of every analysis (label) of a file as Python lists, three to five times
their size on disk, and builds the arrays of all of them again at every access to `samples_dict`: 2 to 4 GB for a
file of ten analyses, of which a run uses one or two. `read_lean` gives an object with the attributes the
parameters_estimation mode uses (labels, approximant, psd, skymap, config, samples_dict), read by pesummary itself:

- the labels and settings from a copy of the file whose sample tables are cut to their first rows and whose sky
  maps, noise spectra and prior samples are left out;
- the samples of a label when they are first asked for, from its table in the file (pesummary returns the table
  as it is; a table it would edit goes through pesummary, from a copy holding that label alone);
- the sky map and the noise spectra of a label when they are first asked for, from a copy holding that label
  alone.

A file that is not laid out as one group per label with a `posterior_samples` table, or that fails on this path,
is read whole by pesummary, as before.
"""
from __future__ import annotations

import tempfile
from collections.abc import Mapping
from pathlib import Path
from typing import Any, Callable, Iterator, Optional

KEEP = 3                                                    # labels whose samples stay in memory
STUB_ROWS = 2                                               # rows of each sample table in the copy read first
HEAVY = ("posterior_samples", "psds", "skymap", "priors", "calibration_envelope", "injection_data")


def _label_groups(h) -> list[str]:
    import h5py

    return [k for k in h if isinstance(h[k], h5py.Group) and isinstance(h[k].get("posterior_samples"), h5py.Dataset)
            and h[k]["posterior_samples"].dtype.names]


def _copy(src, dst, *, label: Optional[str] = None, samples: bool = False, skymap: bool = False,
          psds: bool = False) -> None:
    """Copy the file without what fills the memory. Without `label`: every label, its sample table cut to
    STUB_ROWS rows and without its sky map, noise spectra and prior samples. With `label`: that label alone, with
    its whole sample table (`samples`), its sky map (`skymap`) or its noise spectra (`psds`), and without the
    other heavy groups."""
    labels = set(_label_groups(src))
    for k, v in src.attrs.items():
        dst.attrs[k] = v
    for name in src:
        if name not in labels:
            src.copy(src[name], dst, name=name)
            continue
        if label is not None and name != label:
            continue
        g = dst.create_group(name)
        for k, v in src[name].attrs.items():
            g.attrs[k] = v
        for child in src[name]:
            if child == "posterior_samples":
                ds = src[name][child]
                out = g.create_dataset(child, data=ds[()] if samples else ds[:STUB_ROWS])
                for k, v in ds.attrs.items():
                    out.attrs[k] = v
            elif child in ("skymap", "psds"):
                if (skymap and child == "skymap") or (psds and child == "psds"):
                    src.copy(src[name][child], g, name=child)
            elif child not in HEAVY:
                src.copy(src[name][child], g, name=child)


class _OneLabel:
    """Reads one label of the file with pesummary, from a copy holding that label alone."""

    def __init__(self, path: Path, tmp_dir: Path):
        self._path, self._tmp, self._n = path, tmp_dir, 0

    def read(self, label: str, **what):
        import h5py
        from pesummary.io import read

        self._n += 1
        one = self._tmp / f"label_{self._n}.h5"
        try:
            with h5py.File(self._path, "r") as src, h5py.File(one, "w") as dst:
                _copy(src, dst, label=label, **what)
            return read(str(one))
        finally:
            one.unlink(missing_ok=True)


class _LazySamples(Mapping):
    """label -> samples (pesummary SamplesDict), read from the file when first asked for."""

    def __init__(self, reader: _OneLabel, labels: list[str]):
        self._reader, self._labels = reader, list(labels)
        self._loaded: dict[str, Any] = {}

    def __iter__(self) -> Iterator[str]:
        return iter(self._labels)

    def __len__(self) -> int:
        return len(self._labels)

    def __contains__(self, label) -> bool:
        return label in self._labels

    def _from_table(self, label: str):
        """The samples of `label` as pesummary gives them for a release file: the columns of its table, in the
        file's order (checked against pesummary's own read on GWTC-2.1, GWTC-3 and GWTC-4 files). None for a
        table pesummary would change: a column that is not a number, or a NaN likelihood (it drops those rows)."""
        import h5py
        import numpy as np
        from pesummary.utils.samples_dict import SamplesDict

        with h5py.File(self._reader._path, "r") as h:
            t = h[label]["posterior_samples"][()]
        names = list(t.dtype.names)
        if any(t.dtype[n].kind not in "fiu" for n in names):
            return None
        if "log_likelihood" in names and np.isnan(t["log_likelihood"]).any():
            return None
        return SamplesDict(names, np.array([np.asarray(t[n], float) for n in names]))

    def light(self, label: str):
        """The columns of `label` read one at a time from the file, for a table over all the labels that needs
        a few columns of each; the samples themselves when the label is loaded or pesummary would edit its table."""
        import h5py
        import numpy as np

        if label in self._loaded:
            return self[label]
        with h5py.File(self._reader._path, "r") as h:
            ds = h[label]["posterior_samples"]
            names = list(ds.dtype.names)
            plain = all(ds.dtype[n].kind in "fiu" for n in names) and not (
                "log_likelihood" in names and np.isnan(ds.fields("log_likelihood")[()]).any())
        return _TableColumns(self._reader._path, label, names) if plain else self[label]

    def __getitem__(self, label: str):
        if label not in self._labels:
            raise KeyError(label)
        if label not in self._loaded:
            direct = self._from_table(label)
            if direct is None:                               # a table pesummary would edit: read by pesummary
                sd = self._reader.read(label, samples=True).samples_dict
                direct = sd[label] if label in getattr(sd, "keys", lambda: ())() else sd
            self._loaded[label] = direct
            if len(self._loaded) > KEEP:
                import gc

                while len(self._loaded) > KEEP:
                    self._loaded.pop(next(iter(self._loaded)))
                gc.collect()
        else:
            self._loaded[label] = self._loaded.pop(label)    # most recently used last
        return self._loaded[label]


class _TableColumns(Mapping):
    """The columns of the sample table of one label, each read from the file when asked for."""

    def __init__(self, path: Path, label: str, names: list[str]):
        self._path, self._label, self._names = path, label, names

    def __iter__(self) -> Iterator[str]:
        return iter(self._names)

    def __len__(self) -> int:
        return len(self._names)

    def __getitem__(self, name: str):
        import h5py
        import numpy as np

        if name not in self._names:
            raise KeyError(name)
        with h5py.File(self._path, "r") as h:
            return np.asarray(h[self._label]["posterior_samples"].fields(name)[()], float)


_UNREAD = object()


class _LazySkymaps(dict):
    """label -> sky map of the labels that have one, read from the file when first asked for; a dict, so that a
    supplementary map can be attached to the labels without one."""

    def __init__(self, reader: _OneLabel, labels: list[str]):
        super().__init__({lab: _UNREAD for lab in labels})
        self._reader = reader

    def __getitem__(self, label: str):
        v = super().__getitem__(label)
        if v is _UNREAD:
            v = self._reader.read(label, skymap=True).skymap[label]
            super().__setitem__(label, v)
        return v

    def get(self, label, default=None):
        return self[label] if label in self else default

    def values(self):
        return [self[k] for k in self]

    def items(self):
        return [(k, self[k]) for k in self]


class _LazyPsds(dict):
    """label -> noise spectra of the labels that have a `psds` group, read from the file when first asked for; a
    dict, so that supplementary spectra can be attached. `has` tells whether a label has spectra without
    reading them."""

    def __init__(self, reader: _OneLabel, filled: dict[str, bool]):
        super().__init__({lab: _UNREAD for lab in filled})
        self._reader, self._filled = reader, dict(filled)

    def has(self, label: str) -> Optional[bool]:
        """Whether `label` has a spectrum in the file; None once its entry was read or replaced (look at it)."""
        if label not in self:
            return False
        return self._filled.get(label) if super().__getitem__(label) is _UNREAD else None

    def __getitem__(self, label: str):
        v = super().__getitem__(label)
        if v is _UNREAD:
            v = self._reader.read(label, psds=True).psd[label]
            super().__setitem__(label, v)
        return v

    def get(self, label, default=None):
        return self[label] if label in self else default

    def values(self):
        return [self[k] for k in self]

    def items(self):
        return [(k, self[k]) for k in self]


class LeanPE:
    """The pesummary object of the cut copy, with `samples_dict` and `skymap` read label by label from the file."""

    def __init__(self, base: Any, samples: _LazySamples, skymaps: Optional[_LazySkymaps],
                 psds: Optional[_LazyPsds], tmp: tempfile.TemporaryDirectory):
        self._base, self.samples_dict, self._tmp = base, samples, tmp
        if skymaps:
            self.skymap = skymaps
        if psds:
            self.psd = psds

    def __getattr__(self, name: str):                 # labels, approximant, psd, config...
        return getattr(self._base, name)


def read_lean(path: str | Path, log: Optional[Callable[[str], None]] = None):
    """The PE file `path` as `LeanPE`, or as pesummary's own object when the file cannot be read that way."""
    import h5py
    from pesummary.io import read

    log = log or (lambda m: None)
    path = Path(path)
    tmp = None
    try:
        with h5py.File(path, "r") as src:
            labels = _label_groups(src)
            if not labels:
                return read(str(path))
            try:
                tmp = tempfile.TemporaryDirectory(prefix=".pe_lean_", dir=path.parent)
            except OSError:                           # a read-only folder: the system's temporary folder
                tmp = tempfile.TemporaryDirectory(prefix="pe_lean_")
            stub = Path(tmp.name) / "stub.h5"
            with h5py.File(stub, "w") as dst:
                _copy(src, dst)
            with_map = [lab for lab in labels if "skymap" in src[lab]]
            with_psd = {lab: any(isinstance(v, h5py.Dataset) and v.shape and v.shape[0] > 0
                                 for v in src[lab]["psds"].values())
                        for lab in labels if isinstance(src[lab].get("psds"), h5py.Group)}
        base = read(str(stub))
        stub.unlink(missing_ok=True)
        if sorted(base.labels) != sorted(labels):
            raise ValueError(f"labels {sorted(base.labels)} read for {sorted(labels)}")
        one = _OneLabel(path, Path(tmp.name))
        return LeanPE(base, _LazySamples(one, list(base.labels)),
                      _LazySkymaps(one, [lab for lab in base.labels if lab in with_map]),
                      _LazyPsds(one, {lab: with_psd[lab] for lab in base.labels if lab in with_psd}), tmp)
    except Exception as e:                            # an unexpected layout: the whole file, as before
        log(f"ℹ️ [INFO] {path.name} read whole ({type(e).__name__}: {e})")
        if tmp is not None:
            tmp.cleanup()
        return read(str(path))
