"""Galaxy catalogs for the dark-siren analysis: a standard galaxy file from GLADE+ or from any catalog.

The standard galaxy file is an HDF5 file with a group ``galaxies`` holding one dataset per column,
appended in chunks so that catalogs larger than memory (DES, Rubin) can be converted:

- ``ra``, ``dec``: right ascension and declination [rad];
- ``z``: redshift of the galaxy (CMB frame, corrected for peculiar velocities when available);
- ``sigmaz``: its 1-sigma uncertainty;
- ``m``: apparent magnitude in the band of the analysis.

Its attributes record the band (an icarogw band name, e.g. ``K-glade+``), the source and the selection.
The ``galaxy_catalog`` mode turns this file into the icarogw line-of-sight catalog (see
``dark_catalog_icarogw.py``).

Sources:

- ``fetch_glade_kband``: the GLADE+ galaxies with a Ks magnitude (the selection of the GWTC-3.0, GWTC-4.0
  and GWTC-5.0 cosmology analyses), downloaded from VizieR (VII/291) in declination bands, each band
  cached so that an interrupted download resumes;
- ``convert_catalog``: any catalog in Parquet, FITS, HDF5 or CSV, read in chunks with a column mapping,
  for deeper catalogs (DES Y6, Rubin).
"""
from __future__ import annotations

import io
import json
import os
import time
from pathlib import Path
from typing import Callable, Iterable, Iterator, Optional

import numpy as np
import pandas as pd

GALAXY_COLUMNS = ("ra", "dec", "z", "sigmaz", "m")
GLADE_ASU = "https://vizier.cds.unistra.fr/viz-bin/asu-tsv"
GLADE_TABLE = "VII/291/gladep"
GLADE_KBAND_COLS = ("GLADE+", "Type", "RAJ2000", "DEJ2000", "Kmag", "zhelio", "zcmb", "f_zcmb", "e_z", "e_zhelio")
GLADE_NUMERIC = ("RAJ2000", "DEJ2000", "Kmag", "zhelio", "zcmb", "f_zcmb", "e_z", "e_zhelio")
GLADE_REDSHIFTS = ("zcmb", "zhelio")
GLADE_SIGMAZ = ("quadrature", "measurement", "peculiar")
_GLADE_SIGMAZ_TEXT = {"quadrature": "sqrt(e_z^2 + e_zhelio^2): peculiar-velocity and measurement errors in quadrature",
                      "measurement": "e_zhelio: measurement error", "peculiar": "e_z: peculiar-velocity error"}
# the selection of the GWTC-4.0 analysis (Section 3.2 of arXiv:2509.04348): galaxies with a Ks magnitude and a
# redshift; its threshold map leaves ~5% of the nside-32 pixels empty, as this selection does
GLADE_DEFAULT_SELECTION = dict(types=("G",), redshift="zcmb", sigmaz="quadrature", sigmaz_const=None,
                               sigmaz_relative=False, where=None)
GLADE_MAX_ROWS = 400_000


def _log(msg: str) -> None:
    print(f"[galaxies] {msg}", flush=True)


def cache_dir() -> Path:
    return Path(os.environ.get("GWTC_GALAXY_CACHE", Path.home() / ".cache_gwtc_analysis" / "galaxies"))


# ---------------------------------------------------------------------------
# standard galaxy file
# ---------------------------------------------------------------------------

class GalaxyWriter:
    """Appends galaxies, chunk by chunk, to the standard galaxy file (written as <out>.part, renamed on close)."""

    def __init__(self, out: str | Path, band: str, attrs: Optional[dict] = None):
        import h5py

        self.out = Path(out)
        self.part = self.out.with_name(self.out.name + ".part")
        self.out.parent.mkdir(parents=True, exist_ok=True)
        self.h = h5py.File(self.part, "w")
        g = self.h.create_group("galaxies")
        for c in GALAXY_COLUMNS:
            g.create_dataset(c, shape=(0,), maxshape=(None,), dtype="f8", chunks=(1 << 16,), compression="gzip")
        self.h.attrs["band"] = band
        self.h.attrs["columns"] = json.dumps({"ra": "rad", "dec": "rad", "z": "", "sigmaz": "", "m": "mag"})
        for k, v in (attrs or {}).items():
            self.h.attrs[k] = v if isinstance(v, (int, float, str)) else json.dumps(v)
        self.n = 0

    def append(self, df: pd.DataFrame) -> int:
        """Append the rows of `df` (columns GALAXY_COLUMNS, ra/dec in rad) with all columns finite."""
        d = df[list(GALAXY_COLUMNS)].astype(float)
        d = d[np.isfinite(d.to_numpy()).all(axis=1)]
        k = len(d)
        if k:
            g = self.h["galaxies"]
            for c in GALAXY_COLUMNS:
                g[c].resize((self.n + k,))
                g[c][self.n:] = d[c].to_numpy()
            self.n += k
        return k

    def close(self) -> Path:
        self.h.attrs["n_galaxies"] = self.n
        self.h.close()
        self.part.replace(self.out)
        return self.out


def read_galaxy_file(path: str | Path) -> tuple[dict, dict]:
    """(columns as numpy arrays, attributes) of a standard galaxy file."""
    import h5py

    with h5py.File(path, "r") as h:
        return {c: h["galaxies"][c][:] for c in GALAXY_COLUMNS}, dict(h.attrs)


# ---------------------------------------------------------------------------
# GLADE+ Ks band
# ---------------------------------------------------------------------------

def glade_selection(**kw) -> dict:
    """GLADE+ selection: GLADE_DEFAULT_SELECTION overridden by the arguments that are not None.

    - types: GLADE+ object types kept (G galaxies, Q quasars);
    - redshift: zcmb (CMB frame, peculiar velocities corrected below z = 0.05) or zhelio (heliocentric);
    - sigmaz: quadrature (measurement and peculiar-velocity errors), measurement (e_zhelio), peculiar (e_z); or
      sigmaz_const, a constant (per 1 + z with sigmaz_relative);
    - where: a pandas query on the VizieR columns (RAJ2000, DEJ2000, Kmag, zhelio, zcmb, f_zcmb, e_z, e_zhelio)."""
    sel = dict(GLADE_DEFAULT_SELECTION)
    sel.update({k: v for k, v in kw.items() if v is not None})
    if isinstance(sel["types"], str):
        sel["types"] = tuple(t.strip() for t in sel["types"].split(",") if t.strip())
    sel["types"] = tuple(sel["types"])
    if not sel["types"] or not set(sel["types"]) <= {"G", "Q"}:
        raise ValueError(f"GLADE+ types must be among G, Q, got {sel['types']}")
    if sel["redshift"] not in GLADE_REDSHIFTS:
        raise ValueError(f"GLADE+ redshift must be one of {', '.join(GLADE_REDSHIFTS)}")
    if sel["sigmaz"] not in GLADE_SIGMAZ:
        raise ValueError(f"GLADE+ sigmaz must be one of {', '.join(GLADE_SIGMAZ)}")
    sel["sigmaz_relative"] = bool(sel["sigmaz_relative"])
    return sel


def glade_selection_text(sel: dict) -> str:
    sz = (f"constant {sel['sigmaz_const']:g}" + (" (1 + z)" if sel["sigmaz_relative"] else "")
          if sel.get("sigmaz_const") is not None else sel["sigmaz"])
    return (f"Type {'+'.join(sel['types'])}, Kmag finite, {sel['redshift']} > 0, sigmaz {sz} finite"
            + (f", where {sel['where']}" if sel.get("where") else ""))


def glade_sigmaz(e_z: np.ndarray, e_zhelio: np.ndarray) -> np.ndarray:
    """Redshift uncertainty of a GLADE+ galaxy: the measurement error (spectroscopic, or 2MPZ photometric) and the
    peculiar-velocity error added in quadrature; a missing term counts as zero, both missing gives NaN."""
    a, b = np.asarray(e_z, float), np.asarray(e_zhelio, float)
    both = np.isnan(a) & np.isnan(b)
    s = np.sqrt(np.nan_to_num(a) ** 2 + np.nan_to_num(b) ** 2)
    return np.where(both, np.nan, s)


def glade_kband_frame(raw: pd.DataFrame, galaxies_only: bool = True, selection: Optional[dict] = None) -> pd.DataFrame:
    """Standard columns from a GLADE+ VizieR table: ra/dec in rad, z, sigmaz, m = Kmag (Vega), with the selection
    of `glade_selection` (default: galaxies, zcmb, errors in quadrature; `galaxies_only=False` adds quasars)."""
    sel = selection or glade_selection(types=("G",) if galaxies_only else ("G", "Q"))
    df = raw.copy()
    for c in GLADE_NUMERIC:
        if c in df:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    if "Type" in df:
        df = df[df["Type"].astype(str).str.strip().isin(sel["types"])]
    if sel.get("where"):
        df = df.query(sel["where"])
    z = df[sel["redshift"]].to_numpy(dtype=float)
    if sel.get("sigmaz_const") is not None:
        sz = float(sel["sigmaz_const"]) * ((1 + z) if sel["sigmaz_relative"] else np.ones_like(z))
    elif sel["sigmaz"] == "measurement":
        sz = df["e_zhelio"].to_numpy(dtype=float)
    elif sel["sigmaz"] == "peculiar":
        sz = df["e_z"].to_numpy(dtype=float)
    else:
        sz = glade_sigmaz(df["e_z"], df["e_zhelio"])
    out = pd.DataFrame({"ra": np.deg2rad(df["RAJ2000"].to_numpy()), "dec": np.deg2rad(df["DEJ2000"].to_numpy()),
                        "z": z, "sigmaz": sz, "m": df["Kmag"].to_numpy()})
    return out[(out["z"] > 0) & np.isfinite(out["m"])]


def _vizier_band(dec_lo: float, dec_hi: float, timeout: float, session) -> pd.DataFrame:
    import requests

    params = {"-source": GLADE_TABLE, "-out": ",".join(GLADE_KBAND_COLS), "-out.max": str(GLADE_MAX_ROWS),
              "Kmag": "<99", "DEJ2000": f"{dec_lo:.4f}..{dec_hi:.4f}", "-oc.form": "d"}
    for attempt in range(5):
        try:
            r = session.get(GLADE_ASU, params=params, timeout=timeout)
            r.raise_for_status()
            break
        except requests.RequestException:
            if attempt == 4:
                raise
            time.sleep(10 * (attempt + 1))
    lines = [ln for ln in r.text.splitlines() if ln and not ln.startswith("#")]
    if len(lines) < 3:                          # header, units, dashes: no row
        return pd.DataFrame(columns=list(GLADE_KBAND_COLS))
    df = pd.read_csv(io.StringIO("\n".join([lines[0]] + lines[3:])), sep="\t", dtype=str)
    if len(df) >= GLADE_MAX_ROWS:
        raise ValueError(f"GLADE+ band {dec_lo}..{dec_hi} hit the {GLADE_MAX_ROWS}-row limit: use narrower bands")
    return df


def glade_dec_bands(width: float = 2.0) -> list[tuple[float, float]]:
    """Declination bands covering the sky; VizieR's range is inclusive, so a galaxy exactly on an edge is kept
    once by the dedup on its GLADE+ number."""
    edges = np.arange(-90.0, 90.0 + width / 2, width)
    edges[-1] = 90.0
    return [(float(a), float(b)) for a, b in zip(edges[:-1], edges[1:])]


def fetch_glade_kband(out: str | Path, *, band_width: float = 2.0, galaxies_only: bool = True,
                      selection: Optional[dict] = None, cache: Optional[Path] = None, timeout: float = 300.0,
                      fetch: Optional[Callable[[float, float], pd.DataFrame]] = None) -> Path:
    """Standard galaxy file of the GLADE+ galaxies with a Ks magnitude, from VizieR in declination bands.

    `selection` (see `glade_selection`) chooses the types, redshift, redshift error and an extra cut. Each band is
    cached (raw VizieR columns, CSV) in `cache`; a rerun only downloads the missing bands, and the bands cached
    without a column the selection needs. `fetch(dec_lo, dec_hi)` replaces the VizieR query (tests)."""
    import requests

    sel = selection or glade_selection(types=("G",) if galaxies_only else ("G", "Q"))
    need = {"GLADE+", "Type", "RAJ2000", "DEJ2000", "Kmag", sel["redshift"], "e_z", "e_zhelio"}
    cache = Path(cache or cache_dir() / "glade_kband")
    cache.mkdir(parents=True, exist_ok=True)
    bands = glade_dec_bands(band_width)
    seen: set[str] = set()
    w = GalaxyWriter(out, band="K-glade+", attrs=dict(
        source="GLADE+ (Dalya et al. 2022), VizieR VII/291, galaxies with a Ks magnitude",
        selection=glade_selection_text(sel), glade_selection=json.dumps(sel),
        redshift=("zcmb (CMB frame, peculiar velocities corrected below z = 0.05)" if sel["redshift"] == "zcmb"
                  else "zhelio (heliocentric)"),
        sigmaz=_GLADE_SIGMAZ_TEXT[sel["sigmaz"]] if sel.get("sigmaz_const") is None else
        f"constant {sel['sigmaz_const']:g}" + (" (1 + z)" if sel["sigmaz_relative"] else ""),
        magnitude="Kmag, 2MASS Ks, Vega"))
    try:
        with requests.Session() as s:
            s.headers["User-Agent"] = "gwtc_analysis (https://github.com/danielsentenac/gwtc_analysis)"
            for i, (lo, hi) in enumerate(bands):
                f = cache / f"dec_{lo:+07.2f}_{hi:+07.2f}.csv"
                raw = pd.read_csv(f, dtype=str) if f.exists() else None
                if raw is None or not need <= set(raw.columns):
                    raw = (fetch or (lambda a, b: _vizier_band(a, b, timeout, s)))(lo, hi)
                    raw.to_csv(f.with_suffix(".part"), index=False)
                    f.with_suffix(".part").replace(f)
                if "GLADE+" in raw:
                    ids = raw["GLADE+"].astype(str).str.strip()
                    keep = ~ids.isin(seen)
                    seen.update(ids[keep])
                    raw = raw[keep.to_numpy()]
                k = w.append(glade_kband_frame(raw, selection=sel))
                _log(f"band {i + 1}/{len(bands)} dec {lo:+.1f}..{hi:+.1f}: {len(raw)} rows, {k} kept, total {w.n}")
    except BaseException:
        w.h.close()
        raise
    return w.close()


# ---------------------------------------------------------------------------
# any catalog, in chunks
# ---------------------------------------------------------------------------

def _chunks(path: Path, columns: list[str], chunk_rows: int, fmt: Optional[str] = None,
            hdf5_group: Optional[str] = None) -> Iterator[pd.DataFrame]:
    """Chunks of the `columns` of a catalog file: Parquet (file or directory, e.g. a HATS/LSDB partition tree),
    FITS (first table HDU), HDF5 (datasets of a group) or CSV."""
    fmt = (fmt or "").lower() or (
        "parquet" if path.is_dir() or path.suffix in (".parquet", ".pq") else
        "fits" if path.suffix.lower() in (".fits", ".fit", ".fz") or path.name.lower().endswith(".fits.gz") else
        "hdf5" if path.suffix.lower() in (".h5", ".hdf5", ".hdf") else "csv")
    if fmt == "parquet":
        import pyarrow.dataset as ds

        for b in ds.dataset(str(path), format="parquet").to_batches(columns=columns, batch_size=chunk_rows):
            yield b.to_pandas()
    elif fmt == "fits":
        from astropy.io import fits

        with fits.open(path, memmap=True) as hdul:
            t = next(h for h in hdul if isinstance(h, (fits.BinTableHDU, fits.TableHDU)))
            n = t.header["NAXIS2"]
            for a in range(0, n, chunk_rows):
                rows = t.data[a:a + chunk_rows]
                yield pd.DataFrame({c: np.asarray(rows[c]).astype(float) for c in columns})
    elif fmt == "hdf5":
        import h5py

        with h5py.File(path, "r") as h:
            g = h[hdf5_group] if hdf5_group else h
            n = g[columns[0]].shape[0]
            for a in range(0, n, chunk_rows):
                yield pd.DataFrame({c: g[c][a:a + chunk_rows] for c in columns})
    elif fmt == "csv":
        yield from pd.read_csv(path, usecols=columns, chunksize=chunk_rows)
    else:
        raise ValueError(f"Unknown catalog format {fmt!r} (parquet, fits, hdf5, csv)")


def convert_catalog(path: str | Path, out: str | Path, *, band: str, columns: dict, angle_unit: str = "deg",
                    sigmaz: Optional[float] = None, sigmaz_relative: bool = False, where: Optional[str] = None,
                    chunk_rows: int = 2_000_000, fmt: Optional[str] = None, hdf5_group: Optional[str] = None,
                    source: str = "") -> Path:
    """Standard galaxy file from a catalog read in chunks.

    `columns` maps the standard names to the catalog's: {"ra": ..., "dec": ..., "z": ..., "m": ...,
    "sigmaz": ...} (``sigmaz`` optional: else the constant `sigmaz`, times (1 + z) with `sigmaz_relative`).
    `where` is a pandas query applied to each chunk, on the catalog's column names (quality cuts, star-galaxy
    separation, ...). `band` is the icarogw band of the magnitude."""
    need = ("ra", "dec", "z", "m")
    missing = [c for c in need if c not in columns]
    if missing:
        raise ValueError(f"Column mapping lacks {missing}")
    if "sigmaz" not in columns and sigmaz is None:
        raise ValueError("Give a sigmaz column or a constant sigmaz")
    if angle_unit not in ("deg", "rad"):
        raise ValueError("angle_unit must be deg or rad")
    src_cols = sorted(set(columns.values()) | (_query_names(where) if where else set()))
    w = GalaxyWriter(out, band=band, attrs=dict(source=source or str(path), columns_map=columns,
                                                angle_unit=angle_unit, where=where or "",
                                                sigmaz=columns.get("sigmaz", f"{sigmaz} {'(1+z)' if sigmaz_relative else ''}")))
    try:
        for chunk in _chunks(Path(path), src_cols, chunk_rows, fmt, hdf5_group):
            if where:
                chunk = chunk.query(where)
            to_rad = np.deg2rad if angle_unit == "deg" else (lambda x: x)
            z = chunk[columns["z"]].to_numpy(float)
            sz = (chunk[columns["sigmaz"]].to_numpy(float) if "sigmaz" in columns
                  else np.full(len(chunk), float(sigmaz)) * ((1 + z) if sigmaz_relative else 1.0))
            w.append(pd.DataFrame({"ra": to_rad(chunk[columns["ra"]].to_numpy(float)),
                                   "dec": to_rad(chunk[columns["dec"]].to_numpy(float)),
                                   "z": z, "sigmaz": sz, "m": chunk[columns["m"]].to_numpy(float)}))
            _log(f"{w.n} galaxies written")
    except BaseException:
        w.h.close()
        raise
    return w.close()


def _query_names(where: str) -> set[str]:
    """Column names used in a pandas query (identifiers, backquoted names)."""
    import re
    import keyword

    names = set(re.findall(r"`([^`]+)`", where))
    rest = re.sub(r"`[^`]+`|'[^']*'|\"[^\"]*\"", " ", where)
    names |= {t for t in re.findall(r"[A-Za-z_][A-Za-z0-9_]*", rest)
              if not keyword.iskeyword(t) and t not in ("and", "or", "not", "in", "True", "False", "nan", "inf")}
    return names


def galaxy_summary(path: str | Path) -> dict:
    """Counts and ranges of a standard galaxy file (for the report and the logs)."""
    cols, attrs = read_galaxy_file(path)
    z, m = cols["z"], cols["m"]
    return dict(n=int(len(z)), band=str(attrs.get("band")), z_median=float(np.median(z)) if len(z) else float("nan"),
                z_max=float(np.max(z)) if len(z) else float("nan"),
                sigmaz_median=float(np.median(cols["sigmaz"])) if len(z) else float("nan"),
                m_median=float(np.median(m)) if len(z) else float("nan"), source=str(attrs.get("source", "")),
                selection=selection_text(attrs))


# the selection text of GLADE+ files written before the selection was configurable (the default selection)
_GLADE_OLD_DEFAULT_TEXT = "Type G, Kmag finite, zcmb > 0, sigmaz finite"


def selection_text(attrs: dict) -> str:
    """The galaxy selection of a standard galaxy file, in the wording of `glade_selection_text` for GLADE+."""
    if "glade_selection" in attrs:
        sel = json.loads(attrs["glade_selection"])
        sel["types"] = tuple(sel["types"])
        return glade_selection_text(sel)
    text = str(attrs.get("selection", ""))
    return glade_selection_text(glade_selection()) if text == _GLADE_OLD_DEFAULT_TEXT else text
