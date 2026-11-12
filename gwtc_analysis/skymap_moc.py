"""Sky areas and credible levels on the multi-order (MOC) form of a skymap, without flattening it.

A multi-order map holds a few 10^4 to 10^5 pixels of mixed sizes; flattened to its finest resolution it becomes
10^8 pixels or more and several GB per map (with the sort and the cumulative sum on top). Working on the
multi-order rows keeps the memory to a few tens of MB and gives the same numbers: a coarse pixel is a block of
equal fine pixels, so its place in the greedy ranking is known without expanding it.
"""
from __future__ import annotations

import gzip
import io
from pathlib import Path
from typing import NamedTuple, Union

import numpy as np

SKY_DEG2 = 41252.96


class MocMap(NamedTuple):
    order: np.ndarray        # HEALPix order of each pixel
    ipix: np.ndarray         # NESTED index of each pixel at its own order
    prob: np.ndarray         # probability in each pixel (not a density)

    @property
    def max_order(self) -> int:
        return int(self.order.max())


def read_moc(source: Union[str, Path, bytes]) -> MocMap:
    """The probability map of a FITS skymap (a path, or the bytes of a .fits or .fits.gz file) as multi-order
    pixels. A flat map on disk is merged into multi-order pixels as it is read."""
    from ligo.skymap import moc
    from ligo.skymap.io import read_sky_map

    if isinstance(source, (bytes, bytearray, memoryview)):
        data = bytes(source)
        if data[:2] == b"\x1f\x8b":
            data = gzip.decompress(data)
        source = io.BytesIO(data)
    else:
        source = str(source)

    t = read_sky_map(source, moc=True)
    uniq = np.asarray(t["UNIQ"], dtype=np.int64)
    order, ipix = moc.uniq2nest(uniq)
    prob = np.asarray(t["PROBDENSITY"], dtype=float) * moc.uniq2pixarea(uniq)
    del t
    return MocMap(np.asarray(order, dtype=np.int64), np.asarray(ipix, dtype=np.int64), prob)


def _ranked(m: MocMap):
    """Pixels by decreasing probability density, with the cumulative probability and the number of finest-order
    pixels each one stands for."""
    n_fine = np.left_shift(np.int64(1), 2 * (m.max_order - m.order))
    # density ∝ prob / n_fine; a stable sort keeps the file order among equal densities
    rank = np.argsort(-(m.prob / n_fine), kind="stable")
    return rank, np.cumsum(m.prob[rank]), n_fine[rank]


def credible_area_deg2(m: MocMap, cred: float = 0.9) -> float:
    """Sky area (deg²) of the greedy `cred` credible region, counted in pixels of the finest order of the map, as
    on the flattened map: the pixel that straddles the level gives the fine pixels that still fit below it."""
    rank, cum, n_fine = _ranked(m)
    fine_area = SKY_DEG2 / (12.0 * 4.0 ** m.max_order)
    k = int(np.searchsorted(cum, cred, side="right"))          # pixels 0..k-1 are entirely inside
    n = int(n_fine[:k].sum())
    if k < len(cum):
        before = float(cum[k - 1]) if k else 0.0
        p_fine = float(m.prob[rank[k]]) / int(n_fine[k])
        if p_fine > 0:
            n += min(int(n_fine[k]), int(np.floor((cred - before) / p_fine)))
    return float(n * fine_area)


def credible_level_at(m: MocMap, ra_deg: float, dec_deg: float) -> float:
    """Greedy credible level, in [0, 1], of the pixel holding (ra, dec): the probability of all the pixels at
    least as dense as this one, itself included."""
    import healpy as hp

    top = m.max_order
    fine = int(hp.ang2pix(2 ** top, np.deg2rad(90.0 - dec_deg), np.deg2rad(ra_deg), nest=True))
    hit = np.flatnonzero(m.ipix == np.right_shift(np.int64(fine), 2 * (top - m.order)))
    if hit.size == 0:
        raise ValueError(f"no pixel of the map holds RA {ra_deg}, Dec {dec_deg}")
    rank, cum, _ = _ranked(m)
    return float(cum[np.flatnonzero(rank == hit[0])[0]])
