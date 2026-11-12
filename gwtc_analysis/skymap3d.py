"""3D sky maps of the PE releases: credible areas and volumes, distance by direction, host-galaxy candidates.

The PE releases ship, besides the PE files, an archive of FITS sky maps made by ``ligo-skymap-from-samples``
(Singer et al. 2016): a multi-order HEALPix map whose pixels carry the probability density and an ansatz of
the luminosity-distance distribution along that line of sight (DISTMU, DISTSIGMA, DISTNORM). It is the
3D localization (probability per unit volume), and it is small: the GWTC-5.0 map of GW240615_113620 has
21,504 pixels (0.7 MB) where the HEALPix array inside the PE file has 201 million (1.6 GB) and no distance.

Galaxies: GLADE+ (Dalya et al. 2022, VizieR VII/291) is queried through the VizieR ASU service, in cones
covering the 90% region and the distance range of the map; a user catalog (CSV/TSV with ra, dec and
dist [Mpc] or z columns) can be given instead. Each galaxy gets the 3D probability density at its position
and its searched credible volume; the share of host probability assumes the catalog complete, which GLADE+
is not beyond a few hundred Mpc.
"""
from __future__ import annotations

import html
import json
import os
import re
import shutil
import tarfile
from pathlib import Path
from typing import Any, Callable

import numpy as np
import pandas as pd

LogCb = Callable[[str], None] | None

GLADE_ASU = "https://vizier.cds.unistra.fr/viz-bin/asu-tsv"
GLADE_TABLE = "VII/291/gladep"
GLADE_COLS = ("GLADE+", "PGC", "HyperLEDA", "2MASS", "WISExSCOS", "RAJ2000", "DEJ2000", "dL", "zcmb", "f_dL",
              "Bmag", "Kmag", "W1mag", "M*")
GLADE_NAME_COLS = ("GLADE+", "PGC", "HyperLEDA", "2MASS", "WISExSCOS")
GLADE_MAX_ROWS = 200_000          # per cone
GALAXY_MAX_AREA = 100.0          # deg^2: above this, a GLADE+ query is too large to be useful
TOP_GALAXIES = 25


def cache_dir() -> Path:
    return Path(os.environ.get("GWTC_SKYMAP_CACHE", Path.home() / ".cache_gwtc_analysis" / "skymaps"))


def _log(log_cb: LogCb, msg: str) -> None:
    if log_cb:
        log_cb(msg)


# ---------------------------------------------------------------------
# Finding the FITS map of an event
# ---------------------------------------------------------------------

def catalog_of_pe_file(name: str) -> str | None:
    """Catalog key of a PE release file from its IGWN-GWTC<x>p<y> prefix (GWTC-5.0 file → "GWTC-5")."""
    from .catalog_registry import CATALOGS

    m = re.search(r"IGWN-(GWTC\d+p\d+)", name)
    if not m:
        return None
    for key, cat in CATALOGS.items():
        if any(z.skymap_filename and f"IGWN-{m.group(1)}-" in z.skymap_filename for z in cat.zenodo):
            return key
    return None


def catalog_of_event(event: str) -> str | None:
    """The catalog whose observing runs contain the event date (newest catalog first, updates excluded)."""
    from .catalog_registry import CATALOGS, OBSERVING_RUNS

    m = re.search(r"GW(\d{2})(\d{2})(\d{2})(?:_(\d{2})(\d{2})(\d{2}))?", event)
    if not m:
        return None
    from astropy.time import Time

    y, mo, d, hh, mm, ss = (int(x) if x else 0 for x in m.groups())
    gps = Time(f"20{y:02d}-{mo:02d}-{d:02d}T{hh:02d}:{mm:02d}:{ss:02d}", scale="utc").gps
    for key, cat in reversed(list(CATALOGS.items())):
        if cat.update_of:
            continue
        if any(OBSERVING_RUNS[r].start_gps <= gps <= OBSERVING_RUNS[r].end_gps for r in cat.runs if r in OBSERVING_RUNS):
            return cat.products_from or key          # GWTC-1 events: in the GWTC-2.1 release
    return None


def member_waveform(member: str) -> str:
    """Waveform part of a skymap file name (see `gw_stat.skymap_waveform`)."""
    from .gw_stat import skymap_waveform

    return skymap_waveform(member)


def select_member(members: list[str], event: str, label: str | None) -> str | None:
    """The member of `event` for the PE `label` ('C00:IMRPhenomXPHM-SpinTaylor'), else the fallback of
    `gw_stat.choose_skymap` (Mixed, then IMRPhenomXPHM_SpinTaylor, ...)."""
    from .gw_stat import choose_skymap

    mine = [n for n in members if re.search(re.escape(event) + r"(?!\d)", os.path.basename(n))]
    by_wf: dict[str, str] = {}
    for n in sorted(mine):
        by_wf.setdefault(member_waveform(n), n)
    hit = choose_skymap(by_wf, label)
    return hit[1] if hit else None


def _tar_members(tar_path: Path) -> list[str]:
    """FITS members of a skymap tarball, cached as JSON next to it (listing a .tar.gz reads all of it)."""
    idx = tar_path.with_name(tar_path.name + ".members.json")
    if idx.exists() and idx.stat().st_mtime >= tar_path.stat().st_mtime:
        return json.loads(idx.read_text())
    with tarfile.open(tar_path, "r:gz") as tf:
        names = [m.name for m in tf.getmembers() if re.search(r"\.fits(\.gz)?$", m.name, re.I)]
    idx.write_text(json.dumps(names))
    return names


def fetch_fits(event: str, label: str | None, *, pe_file_name: str = "", catalog: str | None = None,
               version: str | None = None, log_cb: LogCb = None) -> Path | None:
    """Local path of the 3D FITS sky map of `event` for the PE `label`, from the Zenodo skymap archive of its
    catalog (downloaded once into the skymap cache); None if the catalog or the event has none."""
    from .gw_stat import download_zenodo_skymaps_tarball

    from .catalog_registry import CATALOGS

    key = catalog_of_pe_file(pe_file_name) or (catalog if catalog in CATALOGS else None) or catalog_of_event(event)
    if key is not None and CATALOGS[key].products_from:
        key = CATALOGS[key].products_from
    if key is None or not CATALOGS[key].has_skymaps:
        _log(log_cb, f"ℹ️ [INFO] No 3D sky-map archive known for {event} (catalog {key or 'unknown'})")
        return None
    root = cache_dir()
    root.mkdir(parents=True, exist_ok=True)
    tar_path = Path(download_zenodo_skymaps_tarball(key, cache_dir=str(root), progress=True, version=version))
    member = select_member(_tar_members(tar_path), event, label)
    if member is None:
        _log(log_cb, f"ℹ️ [INFO] {event} has no FITS sky map in the {key} archive ({tar_path.name})")
        return None
    out = root / "fits" / os.path.basename(member).replace(":", "_")
    if not out.exists():
        out.parent.mkdir(parents=True, exist_ok=True)
        with tarfile.open(tar_path, "r:gz") as tf:
            src = tf.extractfile(member)
            if src is None:
                return None
            part = out.with_name(out.name + ".part")
            with open(part, "wb") as f:
                shutil.copyfileobj(src, f)
            part.replace(out)
    _log(log_cb, f"ℹ️ [INFO] 3D sky map of {event} ({key}, {member_waveform(member)}): {out}")
    return out


# ---------------------------------------------------------------------
# Reading and summarizing
# ---------------------------------------------------------------------

def read_moc(path: str | Path):
    """Multi-order map (astropy Table with UNIQ, PROBDENSITY and, for a 3D map, DISTMU/DISTSIGMA/DISTNORM)."""
    from ligo.skymap.io import read_sky_map

    return read_sky_map(str(path), moc=True)


def is_3d(m) -> bool:
    return all(c in m.colnames for c in ("DISTMU", "DISTSIGMA", "DISTNORM"))


def _pixel_prob(m) -> np.ndarray:
    import astropy_healpix as ah
    import astropy.units as u

    level, _ = ah.uniq_to_level_ipix(m["UNIQ"])
    return np.asarray(m["PROBDENSITY"]) * ah.nside_to_pixel_area(ah.level_to_nside(level)).to_value(u.sr)


def _map_2d(m):
    """The map without its distance columns: ligo.skymap's crossmatch then leaves the volumes out."""
    return m[["UNIQ", "PROBDENSITY"]]


def volume_stats(m, ra_deg=None, dec_deg=None, dist_mpc=None, contours=(), chunk: int = 2_000_000) -> dict:
    """Credible volumes (Mpc^3) of a 3D map and, at positions with a distance, the probability density per unit
    volume (Mpc^-3) and the searched probability in volume: the quantities of ligo.skymap's crossmatch, on the
    same voxels (each pixel times 999 distance shells out to 6 times the mean distance) and with the same
    formulas. crossmatch holds eight arrays of the voxels at once (1.3 GB for 20,000 pixels); here the density is
    the only one kept, with its sort order, and the cumulative sums are taken by chunks of the sorted voxels."""
    import astropy_healpix as ah
    import healpy as hp
    from ligo.skymap import distance, moc

    o = np.argsort(np.asarray(m["PROBDENSITY"], float))[::-1]          # pixels by decreasing density
    uniq = np.asarray(m["UNIQ"])[o]
    dens, mu, sigma, norm = (np.asarray(m[c], float)[o] for c in ("PROBDENSITY", "DISTMU", "DISTSIGMA", "DISTNORM"))
    d_a = moc.uniq2pixarea(uniq)
    have = dist_mpc is not None and np.size(dist_mpc) > 0
    true_dist = np.atleast_1d(np.asarray(dist_mpc, float)) if have else None

    n_r = 1000
    distmean, _ = distance.parameters_to_marginal_moments(dens * d_a, mu, sigma)
    max_r = 6 * distmean
    if have and np.max(true_dist) > max_r:
        max_r = np.max(true_dist)
    d_r = max_r / n_r
    r = d_r * np.arange(1, n_r)
    shell = (np.square(r) + np.square(d_r) / 12) * d_r                 # voxel volume = shell[k] * d_a[p]
    n = len(r)

    with np.errstate(all="ignore"):
        coef = dens * norm / (sigma * np.sqrt(2 * np.pi))
        dp_dv = np.empty((len(uniq), n))
        rows = max(1, chunk // n)
        for a in range(0, len(uniq), rows):
            b = min(a + rows, len(uniq))
            x = (r[None, :] - mu[a:b, None]) / sigma[a:b, None]
            np.square(x, out=x)
            x *= -0.5
            np.exp(x, out=x)
            x *= coef[a:b, None]
            x[np.isnan(x)] = 0                                         # as crossmatch: invalid voxels hold nothing
            dp_dv[a:b] = x
    flat = dp_dv.ravel()
    order = np.argsort(flat)                                           # ascending; read from its end

    targets = pos_of = None
    out = dict(contour_vols=[np.inf] * len(contours))
    if have:
        level, ipix = moc.uniq2nest(uniq)
        top = int(np.max(level))
        max_ipix = ipix << np.int64(2 * (top - level))
        pix = hp.ang2pix(ah.level_to_nside(top), 0.5 * np.pi - np.radians(np.atleast_1d(dec_deg)),
                         np.radians(np.atleast_1d(ra_deg)), nest=True)
        i = np.argsort(max_ipix)
        p_true = i[np.digitize(pix, max_ipix[i]) - 1]
        k_true = np.digitize(true_dist, r) - 1
        k_true[k_true < 0] += n                                        # crossmatch's index -1: the last shell
        targets = p_true.astype(np.int64) * n + k_true
        out["probdensity_vol"] = flat[targets]
        searched = np.full(len(targets), np.nan)
        t_sorted = np.sort(targets)

    todo = {j: float(c) for j, c in enumerate(contours)}
    carry_p = carry_v = 0.0
    total = len(order)
    for a in range(0, total, chunk):
        sel = order[max(total - a - chunk, 0):total - a][::-1]         # the next voxels, densest first
        dv = shell[sel % n] * d_a[sel // n]
        cum_p = np.cumsum(flat[sel] * dv) + carry_p
        cum_v = np.cumsum(dv) + carry_v
        for j, c in list(todo.items()):                                # np.interp of crossmatch, piece by piece
            q = int(np.searchsorted(cum_p, c, side="left"))
            if q < len(cum_p):
                p0, v0 = (cum_p[q - 1], cum_v[q - 1]) if q else (carry_p, carry_v)
                out["contour_vols"][j] = float(v0 + (c - p0) * (cum_v[q] - v0) / (cum_p[q] - p0)
                                               if cum_p[q] > p0 else cum_v[q])
                del todo[j]
        if have:
            hit = np.flatnonzero(np.isin(sel, t_sorted))
            if hit.size:
                where = dict(zip(sel[hit].tolist(), cum_p[hit].tolist()))
                searched = np.array([where.get(t, old) for t, old in zip(targets.tolist(), searched.tolist())])
        carry_p, carry_v = float(cum_p[-1]), float(cum_v[-1])
        if not todo and not have:
            break
    if have:
        out["searched_prob_vol"] = searched
    return out


def summarize(m) -> dict[str, Any]:
    """Credible areas (deg^2) and volumes (Mpc^3), peak direction, distance there and overall, resolution."""
    import astropy_healpix as ah
    from ligo.skymap.postprocess import crossmatch

    r = crossmatch(_map_2d(m), contours=(0.5, 0.9))          # the areas; the volumes below, in less memory
    prob = _pixel_prob(m)
    level, ipix = ah.uniq_to_level_ipix(m["UNIQ"])
    i = int(np.argmax(np.asarray(m["PROBDENSITY"])))
    lon, lat = ah.healpix_to_lonlat(ipix[i], ah.level_to_nside(level[i]), order="nested")
    out = dict(area50=float(r.contour_areas[0]), area90=float(r.contour_areas[1]), ra_peak=float(lon.deg),
               dec_peak=float(lat.deg), npix=len(m), nside_max=int(ah.level_to_nside(level.max())),
               prob_total=float(prob.sum()), is_3d=is_3d(m))
    meta = getattr(m, "meta", {}) or {}
    out["distmean"] = float(meta.get("distmean", np.nan))
    out["diststd"] = float(meta.get("diststd", np.nan))
    if out["is_3d"]:
        out["vol50"], out["vol90"] = volume_stats(m, contours=(0.5, 0.9))["contour_vols"]
        out["distmu_peak"] = float(m["DISTMU"][i])
        out["distsigma_peak"] = float(m["DISTSIGMA"][i])
    return out


def _region_radius(m, peak_ra: float, peak_dec: float, level: float = 0.9) -> float:
    """Largest angular distance (deg) from the peak to a pixel of the `level` credible region."""
    import astropy_healpix as ah
    import astropy.units as u
    from astropy.coordinates import SkyCoord

    dens = np.asarray(m["PROBDENSITY"])
    prob = _pixel_prob(m)
    o = np.argsort(dens)[::-1]
    k = np.searchsorted(np.cumsum(prob[o]), level * prob.sum()) + 1
    lv, ip = ah.uniq_to_level_ipix(m["UNIQ"][o[:k]])
    lon, lat = ah.healpix_to_lonlat(ip, ah.level_to_nside(lv), order="nested")
    sep = SkyCoord(lon, lat).separation(SkyCoord(peak_ra * u.deg, peak_dec * u.deg)).deg
    half_pix = np.degrees(ah.nside_to_pixel_resolution(ah.level_to_nside(lv)).to_value(u.rad))
    return float(np.max(sep + half_pix))


# ---------------------------------------------------------------------
# Galaxies
# ---------------------------------------------------------------------

def region_cones(m, level: float = 0.9, max_cones: int = 40) -> list[tuple[float, float, float]]:
    """(ra, dec, radius) cones in degrees covering the `level` credible region: one per HEALPix cell of the finest
    order that needs at most `max_cones` cells. A region split into distant patches is covered patch by patch."""
    import astropy_healpix as ah
    import astropy.units as u

    dens = np.asarray(m["PROBDENSITY"])
    prob = _pixel_prob(m)
    o = np.argsort(dens)[::-1]
    k = np.searchsorted(np.cumsum(prob[o]), level * prob.sum()) + 1
    lv, ip = ah.uniq_to_level_ipix(m["UNIQ"][o[:k]])
    cells: set[int] = set()
    for order in range(int(lv.max()), -1, -1):
        # a pixel coarser than `order` covers all of its sub-cells at that order
        cells = set()
        for l, i in zip(lv, ip):
            if l >= order:
                cells.add(int(i) >> (2 * (int(l) - order)))
            else:
                n = 4 ** (order - int(l))
                cells.update(range(int(i) * n, (int(i) + 1) * n))
            if len(cells) > max_cones:
                break
        if len(cells) <= max_cones:
            break
    nside = ah.level_to_nside(order)
    lon, lat = ah.healpix_to_lonlat(np.array(sorted(cells)), nside, order="nested")
    # the farthest corner of a HEALPix cell is within ~ its resolution of the centre
    radius = float(ah.nside_to_pixel_resolution(nside).to_value(u.deg)) * 1.05
    return [(float(a), float(d), radius) for a, d in zip(lon.deg, lat.deg)]


def query_glade(cones: list[tuple[float, float, float]], dmin: float, dmax: float, *,
                timeout: float = 120.0) -> pd.DataFrame:
    """GLADE+ galaxies inside the union of cones (ra, dec, radius in deg) with luminosity distance in [dmin, dmax]
    Mpc, one VizieR ASU request per cone (indexed, a second each; an OR of circles in TAP is not). Cached in the
    skymap cache."""
    import hashlib
    import io
    import time

    import requests

    key = json.dumps([[round(c, 5) for c in cone] for cone in cones] + [round(dmin, 2), round(dmax, 2)])
    cache = cache_dir() / "glade" / (hashlib.sha1(key.encode()).hexdigest()[:16] + ".csv")
    if cache.exists():
        return pd.read_csv(cache, dtype={c: str for c in GLADE_NAME_COLS})
    from .galaxies import vizier_get

    frames = []
    queue = list(cones)
    with requests.Session() as s:
        s.headers["User-Agent"] = "gwtc_analysis (https://github.com/danielsentenac/gwtc_analysis)"
        while queue:
            ra, dec, rad = queue.pop(0)
            params = {"-source": GLADE_TABLE, "-c": f"{ra:.5f} {dec:+.5f}", "-c.rd": f"{rad:.5f}",
                      "-out": ",".join(GLADE_COLS), "dL": f"{dmin:.2f}..{dmax:.2f}", "-out.max": str(GLADE_MAX_ROWS),
                      "-oc.form": "d"}
            r = vizier_get(s, params, timeout)          # the CDS, else its CfA mirror
            lines = [ln for ln in r.text.splitlines() if ln and not ln.startswith("#")]
            if len(lines) < 3:                       # header, units, dashes: no row
                continue
            body = "\n".join([lines[0]] + lines[3:])
            if len(lines) - 3 >= GLADE_MAX_ROWS and rad > 0.05:
                # the service stops at GLADE_MAX_ROWS: query the cone again as 7 cones of half the radius (a centre
                # and a hexagonal ring at sqrt(3)/2 of the radius cover it), duplicates dropped below
                d = np.sqrt(3) / 2 * rad
                queue.append((ra, dec, rad / 2))
                for k in range(6):
                    a = np.radians(60 * k)
                    dd = float(np.clip(dec + d * np.sin(a), -89.999, 89.999))
                    queue.append(((ra + d * np.cos(a) / max(np.cos(np.radians(dd)), 1e-3)) % 360, dd, rad / 2))
                continue
            frames.append(pd.read_csv(io.StringIO(body), sep="\t", dtype={c: str for c in GLADE_NAME_COLS}))
    df = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(columns=list(GLADE_COLS))
    for c in df.columns[df.dtypes == object]:
        df[c] = df[c].str.strip().replace({"-": None, "": None})    # GLADE+ writes "-" for no name
    for c in ("RAJ2000", "DEJ2000", "dL", "zcmb", "Bmag", "Kmag", "W1mag", "M*"):
        if c in df:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    df = df.rename(columns={"RAJ2000": "ra", "DEJ2000": "dec", "dL": "dist", "zcmb": "z", "M*": "mstar"})
    if "GLADE+" in df:
        df = df.drop_duplicates("GLADE+").reset_index(drop=True)
    cache.parent.mkdir(parents=True, exist_ok=True)
    tmp = cache.with_suffix(".part")
    df.to_csv(tmp, index=False)
    tmp.replace(cache)
    return df


def load_galaxies(path: str | Path) -> pd.DataFrame:
    """A user galaxy catalog: CSV/TSV with ra and dec (deg) and dist (Mpc) or z columns (any case)."""
    path = Path(path)
    df = pd.read_csv(path, sep=None, engine="python")
    df.columns = [c.strip() for c in df.columns]
    low = {c.lower(): c for c in df.columns}
    ren = {}
    for want, names in (("ra", ("ra", "raj2000", "ra_deg")), ("dec", ("dec", "dej2000", "dec_deg")),
                        ("dist", ("dist", "distance", "dl", "d_l", "lumdist", "luminosity_distance")),
                        ("z", ("z", "redshift", "zcmb"))):
        for n in names:
            if n in low:
                ren[low[n]] = want
                break
    df = df.rename(columns=ren)
    if "ra" not in df or "dec" not in df:
        raise ValueError(f"{path}: needs ra and dec columns")
    if "dist" not in df:
        if "z" not in df:
            raise ValueError(f"{path}: needs a dist [Mpc] or z column")
        from astropy.cosmology import Planck15

        df["dist"] = Planck15.luminosity_distance(df["z"].to_numpy(float)).value
    return df


def rank_galaxies(m, gal: pd.DataFrame) -> pd.DataFrame:
    """Galaxies with the 3D probability density at their position (Mpc^-3), their searched credible area and
    volume, and their share of the host probability among the listed galaxies, sorted by density."""
    import astropy.units as u
    from astropy.coordinates import SkyCoord
    from ligo.skymap.postprocess import crossmatch

    gal = gal.dropna(subset=["ra", "dec", "dist"])
    gal = gal[gal["dist"] > 0].reset_index(drop=True)
    if gal.empty:
        return gal.assign(dp_dv=[], searched_prob=[], searched_prob_vol=[], host_share=[])
    ra, dec, dist = (gal[c].to_numpy(float) for c in ("ra", "dec", "dist"))
    r = crossmatch(_map_2d(m), SkyCoord(ra * u.deg, dec * u.deg))       # on the sky; in volume below
    v = volume_stats(m, ra, dec, dist) if is_3d(m) else dict(probdensity_vol=np.nan, searched_prob_vol=np.nan)
    out = gal.assign(dp_dv=v["probdensity_vol"], searched_prob=np.atleast_1d(r.searched_prob),
                     searched_prob_vol=v["searched_prob_vol"])
    tot = out["dp_dv"].sum()
    out["host_share"] = out["dp_dv"] / tot if tot > 0 else np.nan
    return out.sort_values("dp_dv", ascending=False).reset_index(drop=True)


# ---------------------------------------------------------------------
# Plot and report
# ---------------------------------------------------------------------

def _raster_order(area90: float, max_level: int) -> int:
    """HEALPix order for the zoomed panels: a few hundred pixels or more in the 90% region, at most order 10."""
    order = 6
    while order < min(max_level, 10) and area90 / (41252.96 / (12 * 4 ** order)) < 400:
        order += 1
    return order


def plot(m, s: dict[str, Any], out_png: str | Path, *, title: str = "", galaxies: pd.DataFrame | None = None,
         dist_samples: np.ndarray | None = None) -> Path:
    """Three panels: probability with the 50/90% contours (and the top galaxies), the distance along each line
    of sight (DISTMU) in the 90% region, and the marginal distance distribution of the map against the PE
    samples. Zoomed on the 90% region, or all-sky when that is wider than 40 degrees."""
    import matplotlib

    with matplotlib.rc_context(matplotlib.rcParamsDefault):      # not the style pesummary sets on import
        return _plot(m, s, out_png, title=title, galaxies=galaxies, dist_samples=dist_samples)


def _plot(m, s, out_png, *, title, galaxies, dist_samples) -> Path:
    import astropy_healpix as ah
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import ligo.skymap.plot  # noqa: F401  (registers the astro projections)
    from ligo.skymap import distance, moc
    from ligo.skymap.postprocess import find_greedy_credible_levels

    max_level = int(np.max(ah.uniq_to_level_ipix(m["UNIQ"])[0]))
    order = _raster_order(s["area90"], max_level)
    # one column of the map rasterized at a time (12.6 million pixels each at order 10), kept only while drawn
    nside = ah.level_to_nside(order)
    prob = np.asarray(moc.rasterize(m[["UNIQ", "PROBDENSITY"]], order)["PROBDENSITY"]) \
        * ah.nside_to_pixel_area(nside).value
    cls = (100 * find_greedy_credible_levels(prob)).astype(np.float32)
    radius = _region_radius(m, s["ra_peak"], s["dec_peak"]) * 1.25
    zoom = radius < 40

    fig = plt.figure(figsize=(14, 4.6))
    kw = (dict(projection="astro zoom", center=f"{s['ra_peak']}d {s['dec_peak']}d", radius=f"{max(radius, 0.5)} deg")
          if zoom else dict(projection="astro hours mollweide"))
    ax1 = fig.add_subplot(1, 3, 1, **kw)
    ax1.imshow_hpx((prob, "ICRS"), nested=True, cmap="cylon")
    ax1.contour_hpx((cls, "ICRS"), nested=True, levels=[50, 90], colors="k", linewidths=0.7)
    ax1.grid(alpha=0.4)
    ax1.set_title(f"Probability: 50% {_num(s['area50'])} deg², 90% {_num(s['area90'])} deg²", fontsize=9)
    if galaxies is not None and len(galaxies):
        top = galaxies.head(TOP_GALAXIES)
        ax1.plot(top["ra"], top["dec"], "o", ms=3, mfc="none", mec="tab:blue", mew=0.8,
                 transform=ax1.get_transform("world"), label=f"top {len(top)} galaxies")
        ax1.legend(loc="lower left", fontsize=7)
    panels = [ax1]
    if s["is_3d"]:
        ax2 = fig.add_subplot(1, 3, 2, **kw)
        mu = np.asarray(moc.rasterize(m[["UNIQ", "DISTMU"]], order)["DISTMU"], float)
        mu[(cls > 90) | ~np.isfinite(mu)] = np.nan
        im = ax2.imshow_hpx((mu, "ICRS"), nested=True, cmap="viridis")
        del mu
        ax2.contour_hpx((cls, "ICRS"), nested=True, levels=[90], colors="k", linewidths=0.7)
        ax2.grid(alpha=0.4)
        fig.colorbar(im, ax=ax2, shrink=0.75, pad=0.02).set_label("distance along the line of sight [Mpc]",
                                                                  fontsize=8)
        ax2.set_title(f"Distance by direction (90% volume {s['vol90']:.3g} Mpc³)", fontsize=9)
        panels.append(ax2)

        ax3 = fig.add_subplot(1, 3, 3)
        dmax = s["distmean"] + 5 * s["diststd"] if np.isfinite(s["distmean"]) else 5000.0
        d = np.linspace(dmax / 2000, dmax, 1000)
        # on the multi-order pixels: the same sum over the sky as on the raster, and the same densest pixel
        ok = np.isfinite(np.asarray(m["DISTMU"], float))
        pdf = distance.marginal_pdf(d, _pixel_prob(m)[ok], np.asarray(m["DISTMU"])[ok],
                                    np.asarray(m["DISTSIGMA"])[ok], np.asarray(m["DISTNORM"])[ok])
        ax3.plot(d, pdf, color="k", lw=1.2, label="3D sky map (all directions)")
        i = int(np.argmax(np.asarray(m["PROBDENSITY"])))
        cond = distance.conditional_pdf(d, m["DISTMU"][i], m["DISTSIGMA"][i], m["DISTNORM"][i])
        ax3.plot(d, cond, color="tab:red", lw=1, ls="--", label="along the most probable direction")
        if dist_samples is not None and len(dist_samples):
            ax3.hist(dist_samples, bins=60, range=(0, dmax), density=True, histtype="step", color="tab:blue",
                     label="PE samples")
        ax3.set_xlabel("luminosity distance [Mpc]")
        ax3.set_ylabel("probability density [Mpc⁻¹]")
        ax3.set_xlim(0, dmax)
        ax3.legend(fontsize=7)
        ax3.set_title(f"Distance: {s['distmean']:.0f} ± {s['diststd']:.0f} Mpc", fontsize=9)
    for ax in panels:
        ax.tick_params(labelsize=7)
        ax.coords[0].set_axislabel("RA", fontsize=8)
        ax.coords[1].set_axislabel("Dec", fontsize=8)
    if title:
        fig.suptitle(title, fontsize=10)
    fig.tight_layout()
    out_png = Path(out_png)
    fig.savefig(out_png, dpi=130)
    plt.close(fig)
    return out_png


def report_html(s: dict[str, Any], *, fits_name: str, galaxies: pd.DataFrame | None = None,
                galaxies_tsv: str | None = None, galaxy_note: str = "") -> str:
    """HTML section of the parameters_estimation report."""
    rows = [("90% / 50% credible area", f"{_num(s['area90'])} / {_num(s['area50'])} deg²"),
            ("Most probable direction", f"RA {s['ra_peak']:.3f}°, Dec {s['dec_peak']:+.3f}°")]
    if s["is_3d"]:
        rows += [("90% / 50% credible volume", f"{s['vol90']:.3g} / {s['vol50']:.3g} Mpc³ (luminosity distance)"),
                 ("Distance, all directions", f"{s['distmean']:.0f} ± {s['diststd']:.0f} Mpc"),
                 ("Distance along the most probable direction",
                  f"{s['distmu_peak']:.0f} ± {s['distsigma_peak']:.0f} Mpc (DISTMU ± DISTSIGMA)")]
    rows += [("Map", f"{s['npix']:,} multi-order pixels, finest nside {s['nside_max']} "
                     f"(<code>{html.escape(fits_name)}</code>)")]
    out = ["<h3>3D sky localization</h3>",
           "<p>The FITS sky map of the PE release (made with <code>ligo-skymap-from-samples</code> from the "
           "posterior samples) gives, in each direction, the probability and the distribution of the luminosity "
           "distance along that line of sight (DISTMU, DISTSIGMA, DISTNORM): the probability per unit volume "
           "(Singer et al. 2016).</p>" if s["is_3d"] else
           "<p>The FITS sky map of the PE release has no distance layers: only the 2D localization.</p>",
           "<table style='border-collapse:collapse'>"]
    out += [f"<tr><td style='padding:2px 12px 2px 0'>{k}</td><td>{v}</td></tr>" for k, v in rows]
    out.append("</table>")
    if galaxies is not None:
        inside = int((galaxies["searched_prob_vol"] <= 0.9).sum()) if len(galaxies) else 0
        out.append(f"<h4>Host-galaxy candidates</h4><p>{galaxy_note} {len(galaxies):,} galaxies queried, "
                   f"{inside:,} inside the 90% credible volume. Ranked by the 3D probability density at their "
                   "position; the share assumes the catalog complete, which it is not at large distance, and "
                   "weights all galaxies equally (no luminosity or stellar-mass weighting).</p>")
        if len(galaxies):
            head = galaxies.head(10)
            names = {"GLADE+": "GLADE+", "PGC": "PGC", "WISExSCOS": "WISE×SCOS", "ra": "RA [°]", "dec": "Dec [°]",
                     "dist": "d<sub>L</sub> [Mpc]", "z": "z", "mstar": "M<sub>*</sub> [10¹⁰ M<sub>☉</sub>]"}
            cols = [c for c in names if c in head and head[c].notna().any()]
            out.append("<table style='border-collapse:collapse;font-size:90%'><tr>"
                       + "".join(f"<th style='padding:2px 8px'>{names[c]}</th>" for c in cols)
                       + "<th>dP/dV [Mpc⁻³]</th><th>searched volume</th><th>share</th></tr>")
            for _, g in head.iterrows():
                cells = [html.escape(_fmt(g[c])) for c in cols]
                out.append("<tr>" + "".join(f"<td style='padding:2px 8px'>{c}</td>" for c in cells)
                           + f"<td>{g['dp_dv']:.3g}</td><td>{g['searched_prob_vol']:.3f}</td>"
                           f"<td>{g['host_share']:.2%}</td></tr>")
            out.append("</table>")
        if galaxies_tsv:
            out.append(f"<p>All of them: <code>{html.escape(galaxies_tsv)}</code></p>")
    elif galaxy_note:
        out.append(f"<p>{galaxy_note}</p>")
    return "\n".join(out)


def _num(x: float) -> str:
    """3 significant digits, without exponent below 10^6 (8,720 rather than 8.72e+03)."""
    return f"{x:,.0f}" if 1e3 <= abs(x) < 1e6 else f"{x:.3g}"


def _fmt(v: Any) -> str:
    if isinstance(v, (float, np.floating)):
        return "" if not np.isfinite(v) else f"{v:.4g}" if abs(v) < 1e4 else f"{v:.0f}"
    return "" if v is None else str(v)


# ---------------------------------------------------------------------
# Driver used by parameters_estimation
# ---------------------------------------------------------------------

def run_skymap3d(event: str, label: str | None, outdir: str | Path, *, pe_file_name: str = "",
                 catalog: str | None = None, version: str | None = None, galaxies: str | None = "glade",
                 galaxy_max_area: float = GALAXY_MAX_AREA, dist_samples: np.ndarray | None = None,
                 fits_path: str | Path | None = None, log_cb: LogCb = None) -> dict[str, Any] | None:
    """Fetch the 3D map, write `<event>_<waveform>_skymap3d.png`, a copy of the FITS file, the host-galaxy table
    `<event>_host_galaxies.tsv`; returns dict(files=..., html=..., summary=...) or None without a map.
    `galaxies`: "glade", a catalog file, or None/"none"."""
    path = Path(fits_path) if fits_path else fetch_fits(event, label, pe_file_name=pe_file_name, catalog=catalog,
                                                        version=version, log_cb=log_cb)
    if path is None:
        return None
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    m = read_moc(path)
    s = summarize(m)
    wf = re.sub(r"[^A-Za-z0-9_.+-]", "_", member_waveform(path.name))
    files: list[str] = []
    fits_copy = outdir / path.name
    if fits_copy.resolve() != path.resolve():
        shutil.copyfile(path, fits_copy)
    files.append(str(fits_copy))

    gal_df, gal_tsv, note = None, None, ""
    choice = (galaxies or "none").strip()
    if choice.lower() != "none" and not s["is_3d"]:
        note = "No galaxy cross-match: the map has no distance layers."
    elif choice.lower() == "glade" and s["area90"] > galaxy_max_area:
        note = (f"No GLADE+ cross-match: the 90% area ({s['area90']:.0f} deg²) is above {galaxy_max_area:g} deg² "
                "(raise --galaxy-max-area, or give a catalog file with --galaxies).")
    elif choice.lower() != "none":
        try:
            if choice.lower() == "glade":
                cones = region_cones(m)
                dlo = max(s["distmean"] - 4 * s["diststd"], 0.0)
                dhi = s["distmean"] + 4 * s["diststd"]
                _log(log_cb, f"ℹ️ [INFO] Querying GLADE+ (VizieR) in {len(cones)} cone(s) of {cones[0][2]:.2f}° "
                             f"covering the 90% region, {dlo:.0f}-{dhi:.0f} Mpc")
                raw = query_glade(cones, dlo, dhi)
                note = (f"GLADE+ (Dálya et al. 2022, VizieR VII/291), {len(cones)} cone(s) of {cones[0][2]:.2f}° "
                        f"covering the 90% region, luminosity distance {dlo:.0f}-{dhi:.0f} Mpc:")
            else:
                raw = load_galaxies(choice)
                note = f"Galaxy catalog <code>{html.escape(str(choice))}</code>:"
            gal_df = rank_galaxies(m, raw)
            gal_tsv = str(outdir / f"{event}_host_galaxies.tsv")
            gal_df.to_csv(gal_tsv, sep="\t", index=False)
            files.append(gal_tsv)
            _log(log_cb, f"ℹ️ [OK] {len(gal_df)} galaxies ranked, "
                         f"{int((gal_df['searched_prob_vol'] <= 0.9).sum()) if len(gal_df) else 0} in the 90% volume:"
                         f" {gal_tsv}")
        except Exception as e:                              # network, catalog format: not fatal
            gal_df, note = None, f"Galaxy cross-match failed: {html.escape(str(e))[:300]}"
            _log(log_cb, f"⚠️ [WARN] {note}")

    png = plot(m, s, outdir / f"{event}_{wf}_skymap3d.png", title=f"{event} — 3D localization ({wf})",
               galaxies=gal_df, dist_samples=dist_samples)
    files.insert(0, str(png))
    _log(log_cb, f"ℹ️ [OK] Saved {png}")
    return dict(files=files, summary=s,
                html=report_html(s, fits_name=path.name, galaxies=gal_df,
                                 galaxies_tsv=Path(gal_tsv).name if gal_tsv else None, galaxy_note=note))
