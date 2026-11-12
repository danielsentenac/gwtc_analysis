"""The multi-order sky area and credible level against the same numbers on the flattened map."""
import healpy as hp
import numpy as np
import pytest

from gwtc_analysis import skymap_moc as sm

TOP = 5          # finest order of the test map


def _mixed_map():
    """A blob around (RA 30, Dec -30): order-5 pixels inside 25°, order-3 pixels elsewhere."""
    rng = np.random.default_rng(1)
    v0 = hp.ang2vec(np.radians(120.0), np.radians(30.0))
    order, ipix, dens = [], [], []
    for coarse in range(hp.nside2npix(2 ** 3)):
        sep = np.degrees(np.arccos(np.clip(np.dot(hp.pix2vec(8, coarse, nest=True), v0), -1, 1)))
        if sep < 25.0:
            kids = coarse * 4 ** (TOP - 3) + np.arange(4 ** (TOP - 3))
            s = np.degrees(np.arccos(np.clip(np.dot(np.array(hp.pix2vec(2 ** TOP, kids, nest=True)).T, v0), -1, 1)))
            order += [TOP] * len(kids)
            ipix += list(kids)
            dens += list(np.exp(-0.5 * (s / 6.0) ** 2) * (1 + 0.01 * rng.random(len(kids))))
        else:
            order.append(3)
            ipix.append(coarse)
            dens.append(np.exp(-0.5 * (sep / 6.0) ** 2) + 1e-4 * rng.random())
    order, ipix, dens = np.array(order), np.array(ipix), np.array(dens)
    prob = dens / 4.0 ** order
    return sm.MocMap(order, ipix, prob / prob.sum())


def _flat(m):
    """The map at its finest order (NESTED), each coarse pixel shared equally among its fine pixels."""
    out = np.zeros(hp.nside2npix(2 ** TOP))
    for o, i, p in zip(m.order, m.ipix, m.prob):
        n = 4 ** (TOP - o)
        out[i * n:(i + 1) * n] = p / n
    return out


@pytest.mark.parametrize("cred", [0.1, 0.5, 0.9, 0.99, 0.9999])
def test_area_as_flat(cred):
    m = _mixed_map()
    flat = np.sort(_flat(m))[::-1]
    area = np.count_nonzero(np.cumsum(flat) <= cred) * sm.SKY_DEG2 / len(flat)
    assert sm.credible_area_deg2(m, cred) == pytest.approx(area, abs=1.5 * sm.SKY_DEG2 / len(flat))


@pytest.mark.parametrize("ra, dec", [(30.0, -30.0), (33.0, -27.0), (50.0, -30.0), (300.0, 40.0)])
def test_level_at_point_as_flat(ra, dec):
    m = _mixed_map()
    flat = _flat(m)
    pix = hp.ang2pix(2 ** TOP, np.radians(90.0 - dec), np.radians(ra), nest=True)
    level = flat[flat >= flat[pix]].sum()          # ties (the fine pixels of one coarse pixel) all counted
    assert sm.credible_level_at(m, ra, dec) == pytest.approx(level, abs=1e-9)


def test_whole_sky():
    m = _mixed_map()
    assert sm.credible_area_deg2(m, 1.0 + 1e-9) == pytest.approx(sm.SKY_DEG2)
    assert sm.credible_level_at(m, 30.0, -30.0) < 0.05
