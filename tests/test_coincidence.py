"""Multi-messenger coincidences: sky overlap, chance probability, catalog parsing; custom bright-siren counterparts."""
from __future__ import annotations

import healpy as hp
import numpy as np
import pandas as pd
import pytest

from gwtc_analysis import coincidence as co

NSIDE = 64


def _blob(ra, dec, sigma_deg):
    """A Gaussian GW map centred on (ra, dec)."""
    v0 = hp.ang2vec(np.radians(90 - dec), np.radians(ra))
    cos = np.clip(np.dot(np.array(hp.pix2vec(NSIDE, np.arange(hp.nside2npix(NSIDE)))).T, v0), -1, 1)
    p = np.exp((cos - 1) / np.radians(sigma_deg) ** 2)
    return co.GWMap(p / p.sum(), gps=1e9)


def test_overlap_and_chance_fraction():
    """Uniform map: overlap 1 everywhere; a concentrated map: Npix p at the peak for a point, smaller for a broad
    counterpart; q is the sky fraction at least as good, 1 outside the map."""
    flat = co.GWMap(np.full(hp.nside2npix(NSIDE), 1.0 / hp.nside2npix(NSIDE)), gps=1e9)
    assert flat.overlap(10, 20, 0) == pytest.approx(1.0) and flat.overlap(10, 20, 5) == pytest.approx(1.0, rel=1e-3)
    gw = _blob(150.0, -30.0, 3.0)
    peak = gw.overlap(150.0, -30.0, 0)
    assert peak == pytest.approx(hp.nside2npix(NSIDE) * gw.prob.max(), rel=0.05)
    assert gw.overlap(150.0, -30.0, 8.0) < peak                 # a coarse localization dilutes the overlap
    assert gw.q(peak, 0) < 1e-3 and gw.q(gw.overlap(330.0, 30.0, 0), 0) == 1.0
    assert gw.credible_percent(150.0, -30.0) < 3 and gw.credible_percent(330.0, 30.0) > 99
    assert gw.area_deg2(0.9) == pytest.approx(2 * np.pi * np.log(10) * 3.0 ** 2, rel=0.1)   # Gaussian: 2πσ² ln 10


class _FakeCatalogs(co.Catalogs):
    def __init__(self, tables):
        self.sources, self.refresh = tuple(tables), False
        self.tables = tables
        self.rates = {s: 1e-5 for s in tables}


def test_search_event_window_and_probability():
    gw = _blob(150.0, -30.0, 3.0)
    grbs = pd.DataFrame(dict(name=["GRB_in", "GRB_late", "GRB_far", "first", "last"],
                             gps=[1e9 + 2, 1e9 + 30, 1e9 + 1, 9e8, 1.1e9],      # the catalog spans the GW time
                             ra=[150.5, 150.0, 330.0, 0.0, 0.0], dec=[-30.0, -30.0, 30.0, 0.0, 0.0],
                             sigma_deg=[0.0] * 5, redshift=[0.05, np.nan, np.nan, np.nan, np.nan]))
    cands, summary = co.search_event(gw, "GWtest", _FakeCatalogs({"gbm": grbs}), co.WINDOWS, far_gw=1.0)
    names = {c["counterpart"] for c in cands}
    assert names == {"GRB_in", "GRB_far"}                      # +30 s is outside -1/+5 s
    best = next(c for c in cands if c["counterpart"] == "GRB_in")
    assert best["redshift"] == 0.05 and best["credible_percent"] < 50
    rate_dt = 1e-5 * 6.0
    assert best["p_chance"] == pytest.approx(-np.expm1(-rate_dt * best["q"]))
    assert best["raven_far_per_yr"] == pytest.approx(1.0 * rate_dt / best["overlap"])
    assert summary[0]["best"] == "GRB_in" and summary[0]["n_in_window"] == 2


AMON = """<table><tr><th>EVENT</th></tr>
<tr><td>132143_36142394</td><td>19/04/05</td><td>21:21:10.50</td><td>BRONZE</td><td>128.6956</td><td>-53.3144</td>
<td>39.13</td><td>14.09</td><td>8.0e+05</td><td>4.6e-01</td><td>1.87</td><td>early entry</td></tr>
<tr><td>140000_1</td><td>1</td><td>24/01/02</td><td>03:04:05.00</td><td>GOLD</td><td>10.0</td><td>+5.0</td>
<td>60.0</td><td>20.0</td><td>1.0e+05</td><td>5.0e-01</td><td>0.3</td><td>rev 1</td></tr>
<tr><td>140000_1</td><td>0</td><td>24/01/02</td><td>03:04:05.00</td><td>GOLD</td><td>11.0</td><td>+5.5</td>
<td>90.0</td><td>20.0</td><td>1.0e+05</td><td>5.0e-01</td><td>0.3</td><td>rev 0</td></tr>
<tr><td>140000_2</td><td>0</td><td>24/01/02</td><td>03:04:05.00</td><td>GOLD</td><td>10.0</td><td>+5.0</td>
<td>60.0</td><td>20.0</td><td>1.0e+05</td><td>5.0e-01</td><td>0.3</td><td>duplicate</td></tr></table>"""


def test_parse_icecube_alerts():
    """The last revision of each alert, duplicates (same time and position) removed, entries before the start of
    the stream (June 2019) left out."""
    d = co._parse_amon(AMON)
    assert len(d) == 1 and d.loc[0, "ra"] == 10.0 and d.loc[0, "error90_arcmin"] == 60.0 and d.loc[0, "kind"] == "GOLD"


def test_custom_bright_counterpart():
    """An event outside the list needs a position and a redshift; its run comes from its name."""
    from gwtc_analysis.counterpart import _run_of_event, get_counterpart

    assert _run_of_event("GW190425_081805") == "O3a" and _run_of_event("GW230529_181500") == "O4a"
    with pytest.raises(ValueError, match="--ra, --dec, --redshift"):
        get_counterpart("GW190425_081805")
    cp = get_counterpart("GW190425_081805", 251.2, -24.2, redshift=(0.035, 0.003), host="H", population="BNS")
    assert (cp.run, cp.population, cp.z, cp.host, cp.pe_event) == ("O3a", "BNS", 0.035, "H", "GW190425_081805")
    assert get_counterpart("GW170817").host == "NGC 4993"     # the registered ones are unchanged


def test_speed_of_gravity_bound():
    """GW170817-like: 1.74 s from 26 Mpc (5% of the line-of-sight distance), emission 0-10 s after the merger:
    -3e-15 <= dv/c <= +7e-16 (ApJL 848, L13)."""
    gw = _blob(150.0, -30.0, 3.0)
    gw.distmu = np.full(len(gw.prob), 40.0)
    gw.distsigma = np.full(len(gw.prob), 1e-6)            # a sharp distance: every quantile at 40 Mpc
    q = gw.distance_quantiles(150.0, -30.0)
    assert q == pytest.approx([40.0, 40.0, 40.0], rel=1e-3)
    b = co.speed_of_gravity(gw, 150.0, -30.0, 1.74)
    k = co.C_KM_S / (40.0 * co.MPC_KM)
    assert b["gw_speed_dv_over_c_max"] == pytest.approx(1.74 * k) and b["gw_speed_dv_over_c_min"] == pytest.approx(-8.26 * k)
    assert 26.0 / 40.0 * b["gw_speed_dv_over_c_max"] == pytest.approx(6.7e-16, rel=0.02)    # the paper's +7e-16 at 26 Mpc
    assert co.speed_of_gravity(co.GWMap(gw.prob, 1e9), 150.0, -30.0, 1.74) == {}           # no distance layers
