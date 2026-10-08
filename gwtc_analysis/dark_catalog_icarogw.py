"""icarogw line-of-sight galaxy catalog for the dark-siren analysis (`galaxy_catalog` mode).

Self-contained on purpose, like h0_icarogw.py: it runs in the icarogw environment (numpy, h5py, healpy,
icarogw), which usually does not have gwtc_analysis installed. It reads the standard galaxy file written by
`gwtc_analysis.galaxies` and drives the pixelated catalog pipeline of icarogw (`icarogw.catalog`, the
functions marked "LVK reviewed"), stage by stage, each stage restartable:

    python dark_catalog_icarogw.py shard     --workdir DIR --galaxies FILE --nside 64 [--nshards 1024]
    python dark_catalog_icarogw.py pixels    --workdir DIR --chunk I --nchunks N
    python dark_catalog_icarogw.py gather    --workdir DIR --nchunks N
    python dark_catalog_icarogw.py prepare   --workdir DIR --chunk I --nchunks N
    python dark_catalog_icarogw.py init      --workdir DIR
    python dark_catalog_icarogw.py interpolate --workdir DIR --chunk I --nchunks N
    python dark_catalog_icarogw.py finish    --workdir DIR
    python dark_catalog_icarogw.py summary   --workdir DIR
    python dark_catalog_icarogw.py status    --workdir DIR

`shard` streams the galaxy file once and spreads the galaxies over `nshards` files by HEALPix pixel range;
`pixels` turns shards into icarogw's per-pixel files (the layout of `create_pixelated_catalogs`) and marks
NaNs; `gather` lists the filled pixels; `prepare` computes the apparent-magnitude threshold (median in the coarser `nside_mthr` pixels) and the
redshift grid of each pixel; `init` builds the common redshift grid and the threshold map; `interpolate`
computes the in-catalog line-of-sight interpolant of each pixel (Schechter function of the band, luminosity
weight epsilon); `finish` assembles the catalog file read by the likelihood (`icarogw_catalog`).

The `pixels`, `prepare` and `interpolate` stages split their work into `nchunks` independent chunks, so they
run in parallel locally or as batch array jobs (Slurm), which is what catalogs deeper than GLADE+ (DES, Rubin)
need. The settings are written once to <workdir>/catalog_settings.json by `shard`.

Each `prepare` and `interpolate` chunk also writes a summary of its pixels (threshold and validity; interpolants)
to <workdir>/done. `init` (with the fixed redshift grid) and `finish` then assemble the catalog from these few
files instead of opening every pixel file in turn, which costs about 70 ms a file on a cluster file system (an hour
for the 46,000 pixels of GLADE+ at nside 64): they reproduce initialize_icarogw_catalog and
icarogw_catalog.build_from_pixelated_files step by step, with icarogw's own objects. Without the summaries (work
directories made before), icarogw's functions are used.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path

import numpy as np

# icarogw copies neighbouring pixel files while other jobs may write them; HDF5 file locking is not needed for
# that, and often fails on cluster file systems
os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")

SETTINGS = "catalog_settings.json"
FIELDS = ("ra", "dec", "z", "sigmaz", "m")


def _log(msg: str) -> None:
    print(f"[catalog] {time.strftime('%H:%M:%S')} {msg}", flush=True)


def _enter(workdir: Path) -> None:
    """icarogw reads `config` from the import path: CPU mode (CUPY=False) from the work directory."""
    workdir.mkdir(parents=True, exist_ok=True)
    os.chdir(workdir)
    if not Path("config.py").exists():
        Path("config.py").write_text("CUPY=False\n")
    sys.path.insert(0, str(workdir))


def settings(workdir: Path) -> dict:
    return json.loads((workdir / SETTINGS).read_text())


def pixel_dir(workdir: Path) -> Path:
    return workdir / "pixels"


def filled_pixels(workdir: Path) -> np.ndarray:
    return np.atleast_1d(np.genfromtxt(pixel_dir(workdir) / "filled_pixels.txt").astype(int))


def chunk_of(items: np.ndarray, chunk: int, nchunks: int) -> np.ndarray:
    """The `chunk`-th of `nchunks` interleaved parts of `items` (interleaved: neighbouring pixels, of similar
    cost, are spread over the chunks)."""
    if not 0 <= chunk < nchunks:
        raise ValueError(f"chunk {chunk} not in [0, {nchunks})")
    return np.asarray(items)[chunk::nchunks]


def cosmo_ref(zmax: float):
    """Reference cosmology of the catalog construction (Planck15, as icarogw's catalog assembly)."""
    import icarogw
    from astropy.cosmology import Planck15

    c = icarogw.cosmology.astropycosmology(zmax=zmax)
    c.build_cosmology(Planck15)
    return c


def _done(workdir: Path, stage: str, chunk: int | None = None) -> Path:
    d = workdir / "done"
    d.mkdir(exist_ok=True)
    return d / (stage if chunk is None else f"{stage}_{chunk}")


# ---------------------------------------------------------------------------
# stages
# ---------------------------------------------------------------------------

def stage_shard(workdir: Path, galaxies: Path, opts: dict, chunk_rows: int = 2_000_000) -> None:
    """Spread the galaxies over `nshards` HDF5 files by pixel range, streaming the galaxy file once."""
    import h5py
    import healpy as hp

    nside, nshards = int(opts["nside"]), int(opts["nshards"])
    npix = hp.nside2npix(nside)
    nshards = min(nshards, npix)
    opts = dict(opts, nshards=nshards, galaxies=str(Path(galaxies).resolve()))
    shards = workdir / "shards"
    shards.mkdir(parents=True, exist_ok=True)
    if _done(workdir, "shard").exists():
        _log("shard: already done")
        return
    for f in shards.glob("shard_*.h5"):
        f.unlink()
    with h5py.File(galaxies, "r") as g:
        band = str(g.attrs["band"])
        n = g["galaxies"]["ra"].shape[0]
        opts.setdefault("band", band)
        handles = {}
        try:
            for a in range(0, n, chunk_rows):
                cols = {c: g["galaxies"][c][a:a + chunk_rows] for c in FIELDS}
                pix = hp.ang2pix(nside, np.pi / 2 - cols["dec"], cols["ra"])     # icarogw.conversions.radec2indeces
                shard = (pix * nshards) // npix
                order = np.argsort(shard, kind="stable")
                bounds = np.searchsorted(shard[order], np.arange(nshards + 1))
                for s in np.nonzero(np.diff(bounds))[0]:
                    sel = order[bounds[s]:bounds[s + 1]]
                    h = handles.get(s)
                    if h is None:
                        h = handles[s] = h5py.File(shards / f"shard_{s:05d}.h5", "w")
                        for c in FIELDS + ("pixel",):
                            h.create_dataset(c, shape=(0,), maxshape=(None,), dtype="i8" if c == "pixel" else "f8",
                                             chunks=(1 << 14,))
                    k = len(sel)
                    for c, v in list(cols.items()) + [("pixel", pix)]:
                        h[c].resize((h[c].shape[0] + k,))
                        h[c][-k:] = v[sel]
                _log(f"shard: {min(a + chunk_rows, n)}/{n} galaxies")
        finally:
            for h in handles.values():
                h.close()
    (workdir / SETTINGS).write_text(json.dumps(opts, indent=1))
    _done(workdir, "shard").touch()


def stage_pixels(workdir: Path, chunk: int, nchunks: int) -> None:
    """Per-pixel files in the layout of icarogw.catalog.create_pixelated_catalogs, NaNs marked."""
    import h5py
    import healpy as hp
    import icarogw

    s = settings(workdir)
    nside, grouping = int(s["nside"]), s["grouping"]
    out = pixel_dir(workdir)
    out.mkdir(exist_ok=True)
    marker = _done(workdir, "pixels", chunk)
    if marker.exists():
        _log(f"pixels {chunk}/{nchunks}: already done")
        return
    shard_ids = chunk_of(np.arange(int(s["nshards"])), chunk, nchunks)
    written = []
    for sid in shard_ids:
        f = workdir / "shards" / f"shard_{sid:05d}.h5"
        if not f.exists():
            continue
        with h5py.File(f, "r") as h:
            pix = h["pixel"][:]
            cols = {c: h[c][:] for c in FIELDS}
        order = np.argsort(pix, kind="stable")
        u, start = np.unique(pix[order], return_index=True)
        stop = np.append(start[1:], len(order))
        for p, a, b in zip(u, start, stop):
            sel = order[a:b]
            path = out / f"pixel_{int(p)}.hdf5"
            with h5py.File(path, "w") as cat:
                cat.attrs["nside"] = nside
                cat.attrs["nest"] = False
                cat.attrs["dOmega_sterad"] = hp.nside2pixarea(nside, degrees=False)
                cat.attrs["dOmega_deg2"] = hp.nside2pixarea(nside, degrees=True)
                cat.attrs["Ntotal_galaxies_original"] = len(sel)
                g = cat.require_group("catalog")
                for c in FIELDS:
                    g.create_dataset(c, data=cols[c][sel], compression="gzip", chunks=True, maxshape=(None,))
            icarogw.catalog.remove_nans_pixelated_files(str(out), int(p), list(FIELDS), grouping)
            written.append(int(p))
    np.savetxt(workdir / "done" / f"pixels_{chunk}.list", np.array(written, dtype=int), fmt="%d")
    marker.touch()
    _log(f"pixels {chunk}/{nchunks}: {len(written)} pixel files from {len(shard_ids)} shards")


def stage_gather(workdir: Path, nchunks: int) -> np.ndarray:
    """filled_pixels.txt from the pixel lists of the `nchunks` chunks of the pixels stage."""
    lists = [workdir / "done" / f"pixels_{c}.list" for c in range(nchunks)]
    missing = [str(p.name) for p in lists if not p.exists()]
    if missing:
        raise RuntimeError(f"pixels stage incomplete: {missing[:5]}{' ...' if len(missing) > 5 else ''}")
    pix = np.sort(np.concatenate([np.atleast_1d(np.loadtxt(p, dtype=int)) for p in lists if p.stat().st_size]))
    np.savetxt(pixel_dir(workdir) / "filled_pixels.txt", pix, fmt="%d")
    _done(workdir, "gather").touch()
    _log(f"gather: {len(pix)} filled pixels")
    return pix


def coarse_pixels(pixels: np.ndarray, nside: int, nside_mthr: int) -> np.ndarray:
    """The nside_mthr pixel of each pixel, as icarogw.catalog.calculate_mthr_pixelated_files finds it (pixel centre,
    RING ordering)."""
    import healpy as hp

    theta, phi = hp.pix2ang(nside, np.asarray(pixels, dtype=int))
    return hp.ang2pix(nside_mthr, theta, phi)


def prepare_chunk(pixels: np.ndarray, nside: int, nside_mthr: int, chunk: int, nchunks: int) -> np.ndarray:
    """The pixels of the `chunk`-th part, whole coarse pixels at a time.

    calculate_mthr_pixelated_files reads (by copying their files) every pixel of the coarse pixel it belongs to:
    with whole coarse pixels per chunk, no chunk reads a file another chunk is writing."""
    big = coarse_pixels(pixels, nside, nside_mthr)
    groups = np.unique(big)
    mine = set(chunk_of(groups, chunk, nchunks).tolist())
    return np.asarray(pixels)[np.isin(big, list(mine))]


def stage_prepare(workdir: Path, chunk: int, nchunks: int) -> None:
    """Apparent-magnitude threshold and redshift grid of each pixel of the chunk (whole coarse pixels)."""
    import icarogw

    s = settings(workdir)
    out = str(pixel_dir(workdir))
    if not _done(workdir, "gather").exists():
        raise RuntimeError("run the gather stage first")
    marker = _done(workdir, "prepare", chunk)
    if marker.exists():
        _log(f"prepare {chunk}/{nchunks}: already done")
        return
    cref = cosmo_ref(float(s["cosmo_zmax"]))
    nint = _nintegration(s)
    todo = prepare_chunk(filled_pixels(workdir), int(s["nside"]), int(s["nside_mthr"]), chunk, nchunks)
    mine = set(int(p) for p in todo)
    for f in pixel_dir(workdir).glob("pixel_*_*.hdf5"):   # temporary copies left by an interrupted run
        if int(f.stem.split("_")[2]) in mine:
            f.unlink()
    for i, p in enumerate(todo):
        icarogw.catalog.calculate_mthr_pixelated_files(out, int(p), "m", s["grouping"], int(s["nside_mthr"]),
                                                       mthr_percentile=float(s["mthr_percentile"]))
        icarogw.catalog.get_redshift_grid_for_files(out, int(p), s["grouping"], cref, Nintegration=nint,
                                                    Numsigma=int(s["numsigma"]), zcut=float(s["zcut"]))
        if (i + 1) % 500 == 0:
            _log(f"prepare {chunk}/{nchunks}: {i + 1}/{len(todo)} pixels")
    _write_prepare_summary(workdir, s["grouping"], todo, chunk)
    marker.touch()
    _log(f"prepare {chunk}/{nchunks}: {len(todo)} pixels")


def _write_prepare_summary(workdir: Path, grouping: str, pixels: np.ndarray, chunk: int) -> None:
    """done/prepare_<chunk>.npz: threshold, validity (a galaxy usable in the interpolant) and galaxy counts (all,
    in the in-catalog term) of each pixel."""
    import h5py

    n = len(pixels)
    mthr, valid = np.empty(n), np.empty(n, dtype=bool)
    n_gal, n_used = np.empty(n, dtype=int), np.empty(n, dtype=int)
    for k, p in enumerate(pixels):
        with h5py.File(pixel_dir(workdir) / f"pixel_{int(p)}.hdf5", "r") as f:
            g = f[grouping]
            mthr[k] = g.attrs["mthr"]
            v = g["valid_galaxies_interpolant"][:]
            valid[k], n_used[k] = bool(np.any(v)), int(np.count_nonzero(v))
            n_gal[k] = int(f.attrs["Ntotal_galaxies_original"])
    np.savez(workdir / "done" / f"prepare_{chunk}.npz", pixels=np.asarray(pixels, dtype=int), mthr=mthr, valid=valid,
             n_gal=n_gal, n_used=n_used)


def _chunk_files(workdir: Path, stage: str, ext: str) -> list[Path] | None:
    """The summary files of every chunk of `stage`, or None when a chunk has none (older work directory)."""
    d = workdir / "done"
    chunks = sorted(int(p.name.split("_")[1]) for p in d.glob(f"{stage}_*") if p.name.split("_")[1].isdigit()
                    and not p.suffix)
    files = [d / f"{stage}_{c}{ext}" for c in chunks]
    return files if files and all(f.exists() for f in files) else None


def _nintegration(s: dict):
    """Integer: icarogw's adaptive redshift grid (refined around every galaxy); "logspace:zmin:N": a fixed grid
    up to zcut, much cheaper to merge for deep catalogs."""
    v = s["nintegration"]
    if isinstance(v, str) and v.startswith("logspace:"):
        _, zmin, n = v.split(":")
        return np.logspace(np.log10(float(zmin)), np.log10(float(s["zcut"])), int(n))
    return int(v)


def stage_init(workdir: Path) -> None:
    """Common redshift grid and apparent-magnitude threshold (MOC) map, in the catalog file."""
    import icarogw

    s = settings(workdir)
    marker = _done(workdir, "init")
    if marker.exists():
        _log("init: already done")
        return
    outfile = workdir / s["outfile"]
    if outfile.exists():
        outfile.unlink()
    summaries = _chunk_files(workdir, "prepare", ".npz")
    if isinstance(_nintegration(s), np.ndarray) and summaries:
        _init_from_summaries(workdir, s, outfile, summaries)
    else:
        icarogw.catalog.initialize_icarogw_catalog(str(pixel_dir(workdir)), str(outfile), s["grouping"])
    marker.touch()
    _log(f"init: {outfile.name}")


def _init_from_summaries(workdir: Path, s: dict, outfile: Path, summaries: list[Path]) -> None:
    """icarogw.catalog.initialize_icarogw_catalog, fixed-grid case, from the prepare summaries."""
    import h5py
    import healpy as hp
    from mhealpy import HealpixMap

    grouping = s["grouping"]
    d = [np.load(f) for f in summaries]
    pix = np.concatenate([x["pixels"] for x in d])
    o = np.argsort(pix)                                    # filled_pixels.txt is sorted
    pix, mthr, valid = pix[o], np.concatenate([x["mthr"] for x in d])[o], np.concatenate([x["valid"] for x in d])[o]
    if not np.array_equal(pix, filled_pixels(workdir)):
        raise RuntimeError("prepare summaries do not cover the filled pixels")
    with h5py.File(pixel_dir(workdir) / f"pixel_{int(pix[0])}.hdf5", "r") as tmpcat:
        nside, nest = tmpcat.attrs["nside"], tmpcat.attrs["nest"]
        dO_sr, dO_deg2 = tmpcat.attrs["dOmega_sterad"], tmpcat.attrs["dOmega_deg2"]
        nint_attr = tmpcat[grouping].attrs["Nintegration"]
        z_grid = tmpcat[grouping]["z_grid"][:]                 # the same fixed grid in every pixel
    keep = valid & np.isfinite(mthr)
    actual_filled_pixels, mth_map = pix[keep], mthr[keep]
    with h5py.File(outfile, "a") as icat:
        icat.attrs["nside"], icat.attrs["nest"] = nside, nest
        icat.attrs["dOmega_sterad"], icat.attrs["dOmega_deg2"] = dO_sr, dO_deg2
        sub = icat.require_group(grouping)
        sub.attrs["zcut"] = nint_attr                          # as icarogw writes it
        sub.create_dataset("z_grid", data=z_grid)
        moc_mthr_map = HealpixMap.moc_from_pixels(nside, actual_filled_pixels, nest=nest, density=False)
        moc_mthr_map[moc_mthr_map.data == 0.] = -np.inf
        mapping_filled_pixels = np.ones_like(actual_filled_pixels)
        for i, pixo in enumerate(actual_filled_pixels):
            theta, phi = hp.pix2ang(nside, pixo, nest=nest)
            p = moc_mthr_map.ang2pix(theta, phi)
            moc_mthr_map[p] = mth_map[i]
            mapping_filled_pixels[i] = p
        sub.create_dataset("mthr_filled_pixels_healpy", data=actual_filled_pixels)
        sub.create_dataset("mthr_filled_pixels_healpy_to_moc_labels", data=mapping_filled_pixels)
        sub.create_dataset("mthr_moc_map", data=moc_mthr_map.data)
        sub.create_dataset("uniq_moc_map", data=moc_mthr_map.uniq)
    np.savetxt(pixel_dir(workdir) / f"{grouping}_common_zgrid.txt", z_grid)
    _log(f"init: {len(actual_filled_pixels)} pixels with galaxies in the in-catalog term (from {len(summaries)} "
         "prepare summaries)")


def stage_interpolate(workdir: Path, chunk: int, nchunks: int) -> None:
    """In-catalog line-of-sight interpolant dN_gal/dz dOmega of each pixel of the chunk."""
    import icarogw

    s = settings(workdir)
    marker = _done(workdir, "interpolate", chunk)
    if marker.exists():
        _log(f"interpolate {chunk}/{nchunks}: already done")
        return
    out = str(pixel_dir(workdir))
    z_grid = np.genfromtxt(pixel_dir(workdir) / f"{s['grouping']}_common_zgrid.txt")
    cref = cosmo_ref(float(s["cosmo_zmax"]))
    todo = chunk_of(filled_pixels(workdir), chunk, nchunks)
    t0 = time.time()
    for i, p in enumerate(todo):
        icarogw.catalog.calculate_interpolant_files(out, z_grid, int(p), s["grouping"], s["subgrouping"], s["band"],
                                                    cref, float(s["epsilon"]), ptype=s["ptype"])
        if (i + 1) % 200 == 0:
            _log(f"interpolate {chunk}/{nchunks}: {i + 1}/{len(todo)} pixels, {time.time() - t0:.0f} s")
    _write_interpolate_summary(workdir, s, todo, chunk)
    marker.touch()
    _log(f"interpolate {chunk}/{nchunks}: {len(todo)} pixels in {time.time() - t0:.0f} s")


def _write_interpolate_summary(workdir: Path, s: dict, pixels: np.ndarray, chunk: int) -> None:
    """done/interpolate_<chunk>.h5: the interpolant of each pixel of the chunk."""
    import h5py

    vals = []
    for p in pixels:
        with h5py.File(pixel_dir(workdir) / f"pixel_{int(p)}.hdf5", "r") as f:
            vals.append(f[s["grouping"]][s["subgrouping"]]["vals_interpolant"][:])
    part = workdir / "done" / f"interpolate_{chunk}.h5.part"
    with h5py.File(part, "w") as h:
        h.create_dataset("pixels", data=np.asarray(pixels, dtype=int))
        h.create_dataset("vals", data=np.array(vals) if vals else np.zeros((0, 0)))
    part.replace(workdir / "done" / f"interpolate_{chunk}.h5")


def _build_from_summaries(cat, workdir: Path, s: dict, summaries: list[Path]) -> None:
    """icarogw_catalog.build_from_pixelated_files, step by step, the interpolants read from the interpolate
    summaries instead of the pixel files."""
    import h5py
    import icarogw
    from astropy.cosmology import Planck15
    from mhealpy import HealpixMap

    with h5py.File(cat.outfile, "r") as icat:
        cat.moc_mthr_map = HealpixMap(data=icat[cat.grouping]["mthr_moc_map"][:], uniq=icat[cat.grouping]["uniq_moc_map"][:])
        cat.z_grid = icat[cat.grouping]["z_grid"][:]
        filled = icat[cat.grouping]["mthr_filled_pixels_healpy"][:]
        moc_pixels = icat[cat.grouping]["mthr_filled_pixels_healpy_to_moc_labels"][:]
    vals = {}
    for f in summaries:
        with h5py.File(f, "r") as h:
            for p, v in zip(h["pixels"][:], h["vals"][:]):
                vals[int(p)] = v
    cat.sky_grid = np.arange(0, len(cat.moc_mthr_map.data), 1).astype(int)
    cat.band, cat.epsilon = s["band"], float(s["epsilon"])
    cat.calc_kcorr = icarogw.catalog.kcorr(cat.band)
    cat.sch_fun = icarogw.cosmology.galaxy_MF(band=cat.band)
    cat.sch_fun.build_effective_number_density_interpolant(cat.epsilon)
    proxy = icarogw.cosmology.astropycosmology(zmax=cat.z_grid[-1] * 2)
    proxy.build_cosmology(Planck15)
    cat.sch_fun.build_MF(proxy)
    dl_proxy = proxy.z2dl(cat.z_grid)
    moc_to_filled = {int(m): int(f) for m, f in zip(moc_pixels, filled)}
    nz, npix = len(cat.z_grid), len(cat.sky_grid)
    cat.dNgal_dzdOm_vals = np.zeros((nz, npix))
    cat.bg_vals = np.empty((nz, npix))
    bg_empty = (cat.sch_fun.background_effective_galaxy_density(-np.inf * np.ones_like(cat.z_grid), cat.z_grid)
                * proxy.dVc_by_dzdOmega_at_z(cat.z_grid))
    dvc = proxy.dVc_by_dzdOmega_at_z(cat.z_grid)
    for pix in cat.sky_grid:
        hp_pix = moc_to_filled.get(int(pix))
        if hp_pix is None:
            cat.bg_vals[:, pix] = bg_empty
        else:
            cat.dNgal_dzdOm_vals[:, pix] = vals[hp_pix]
            mthr = cat.calc_Mthr(cat.z_grid, pix * np.ones_like(cat.z_grid).astype(int), proxy, dl=dl_proxy)
            cat.bg_vals[:, pix] = cat.sch_fun.background_effective_galaxy_density(mthr, cat.z_grid) * dvc


def stage_finish(workdir: Path) -> None:
    """Assemble the catalog file read by the likelihood (icarogw_catalog.load_from_hdf5_file)."""
    import h5py
    import icarogw

    s = settings(workdir)
    marker = _done(workdir, "finish")
    if marker.exists():
        _log("finish: already done")
        return
    outfile = workdir / s["outfile"]
    with h5py.File(outfile, "a") as h:                      # rerun after a crash: drop a partial subgroup
        if s["subgrouping"] in h[s["grouping"]]:
            del h[s["grouping"]][s["subgrouping"]]
    cat = icarogw.catalog.icarogw_catalog(str(outfile), s["grouping"], s["subgrouping"])
    summaries = _chunk_files(workdir, "interpolate", ".h5")
    if summaries:
        _build_from_summaries(cat, workdir, s, summaries)
    else:
        cat.build_from_pixelated_files(str(pixel_dir(workdir)))
    cat.save_to_hdf5_file()
    with h5py.File(outfile, "a") as h:
        h.attrs["gwtc_analysis_settings"] = json.dumps(s)
    marker.touch()
    _log(f"finish: {outfile} ({outfile.stat().st_size / 1e6:.0f} MB), {len(cat.z_grid)} redshifts x "
         f"{cat.dNgal_dzdOm_vals.shape[1]} sky pixels")


def stage_summary(workdir: Path, n_pixels: int = 400) -> dict:
    """summary.json, the apparent-magnitude threshold map and the completeness against redshift (for the report)."""
    import h5py
    import healpy as hp
    import icarogw
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    s = settings(workdir)
    nside = int(s["nside"])
    pix = filled_pixels(workdir)
    mthr = np.full(hp.nside2npix(nside), hp.UNSEEN)
    n_gal = n_used = 0
    summaries = _chunk_files(workdir, "prepare", ".npz")
    d = [np.load(f) for f in summaries] if summaries else []
    if d and all("n_gal" in x for x in d):            # from the prepare summaries: no pass over the pixel files
        for x in d:
            ok = np.isfinite(x["mthr"])
            mthr[x["pixels"][ok]] = x["mthr"][ok]
            n_gal += int(x["n_gal"].sum())
            n_used += int(x["n_used"].sum())
    else:
        for p in pix:
            with h5py.File(pixel_dir(workdir) / f"pixel_{int(p)}.hdf5", "r") as f:
                g = f[s["grouping"]]
                mthr[p] = g.attrs["mthr"] if np.isfinite(g.attrs["mthr"]) else hp.UNSEEN
                n_gal += int(f.attrs["Ntotal_galaxies_original"])
                n_used += int(np.count_nonzero(g["valid_galaxies_interpolant"][:]))
    plots = workdir / "plots"
    plots.mkdir(exist_ok=True)
    hp.mollview(mthr, title=f"Apparent-magnitude threshold ({s['band']}, nside {s['nside_mthr']} medians)",
                unit="m_thr", cmap="viridis")
    hp.graticule()
    plt.savefig(plots / "mthr_map.png", dpi=110)
    plt.close("all")
    cat = icarogw.catalog.icarogw_catalog(str(workdir / s["outfile"]), s["grouping"], s["subgrouping"])
    cat.load_from_hdf5_file()
    cref = cosmo_ref(float(s["cosmo_zmax"]))
    cat.sch_fun.build_MF(cref)
    rng = np.random.default_rng(1)
    sample = rng.choice(pix, size=min(n_pixels, len(pix)), replace=False)
    ra, dec = icarogw.conversions.indices2radec(sample, nside)
    moc = cat.get_NUNIQ_pixel(ra, dec)
    z = np.linspace(1e-3, min(float(s["zcut"]), 0.6), 120)
    _, _, inco, fig, ax = cat.check_differential_effective_galaxies(z, moc, cref)
    ax[0].set_ylabel(r"$dN_{\rm gal}/dz\,d\Omega$ [sr$^{-1}$]")
    ax[1].set_ylabel("completeness")
    ax[1].set_xlabel("z")
    fig.savefig(plots / "completeness.png", dpi=110)
    plt.close("all")
    comp = 1 - inco
    m_ok = mthr[mthr != hp.UNSEEN]
    out = dict(n_galaxies=n_gal, n_galaxies_in_catalog_term=n_used, n_filled_pixels=int(len(pix)),
               sky_fraction=float(len(pix) / hp.nside2npix(nside)), n_redshifts=int(len(cat.z_grid)),
               n_moc_pixels=int(cat.dNgal_dzdOm_vals.shape[1]),
               file_mb=float((workdir / s["outfile"]).stat().st_size / 1e6),
               mthr_median=float(np.median(m_ok)) if len(m_ok) else None,
               mthr_p10_p90=[float(np.percentile(m_ok, 10)), float(np.percentile(m_ok, 90))] if len(m_ok) else None,
               completeness_z=[float(x) for x in z[::10]],
               completeness_median=[float(x) for x in np.median(comp, axis=1)[::10]],
               completeness_p20=[float(x) for x in np.percentile(comp, 20, axis=1)[::10]],
               completeness_p80=[float(x) for x in np.percentile(comp, 80, axis=1)[::10]],
               settings=s)
    (workdir / "summary.json").write_text(json.dumps(out, indent=1))
    _log(f"summary: {n_gal} galaxies ({n_used} in the in-catalog term) in {len(pix)} pixels; median m_thr "
         f"{out['mthr_median']}")
    return out


def status(workdir: Path) -> dict:
    """Which stages and chunks are done."""
    d = workdir / "done"
    names = sorted(p.name for p in d.glob("*")) if d.exists() else []
    out = {}
    for st in ("shard", "gather", "init", "finish"):
        out[st] = st in names
    for st in ("pixels", "prepare", "interpolate"):
        out[st] = sorted(int(n.split("_")[1]) for n in names if n.startswith(st + "_") and not n.endswith(".list"))
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("stage", choices=("shard", "pixels", "gather", "prepare", "init", "interpolate", "finish",
                                      "summary", "status"))
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--galaxies", help="standard galaxy file (shard)")
    ap.add_argument("--settings", help="JSON of the catalog settings (shard)")
    ap.add_argument("--chunk", type=int, default=0)
    ap.add_argument("--nchunks", type=int, default=1)
    a = ap.parse_args(argv)
    workdir = Path(a.workdir).resolve()
    _enter(workdir)
    if a.stage == "shard":
        if not a.galaxies or not a.settings:
            ap.error("shard needs --galaxies and --settings")
        stage_shard(workdir, Path(a.galaxies), json.loads(Path(a.settings).read_text()))
    elif a.stage == "pixels":
        stage_pixels(workdir, a.chunk, a.nchunks)
    elif a.stage == "gather":
        stage_gather(workdir, a.nchunks)
    elif a.stage == "prepare":
        stage_prepare(workdir, a.chunk, a.nchunks)
    elif a.stage == "init":
        stage_init(workdir)
    elif a.stage == "interpolate":
        stage_interpolate(workdir, a.chunk, a.nchunks)
    elif a.stage == "finish":
        stage_finish(workdir)
    elif a.stage == "summary":
        stage_summary(workdir)
    else:
        print(json.dumps(status(workdir)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
