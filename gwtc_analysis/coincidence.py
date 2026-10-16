"""Multi-messenger coincidences of GW events: gamma-ray bursts and neutrinos in time and on the sky.

An offline version of what RAVEN does in low latency for the LVK alerts (Urban 2016; LVK public alerts user guide):
for each GW sky map, the counterparts of external catalogs within a time window of the merger, scored by how well
their position overlaps the GW localization.

Catalogs (public):

- ``gbm``: Fermi-GBM burst catalog (HEASARC ``fermigbrst``, 2008 onwards). Positions with their statistical 1σ
  error and a 3.7° systematic added in quadrature (Connaughton et al. 2015). The catalog replaces the GBM position
  by the counterpart's when one was localized (GRB 170817A carries the position of AT2017gfo, error 0).
- ``swift``: Swift bursts compiled by the burst advocates (HEASARC ``swiftgrbba``, 2004 onwards): the best of the
  UVOT, XRT and BAT 90% positions, the redshift and the host galaxy when measured.
- ``icecube-alerts``: IceCube Gold and Bronze alert tracks (GCN, from the start of the stream on 2019 June 17; the
  last revision of each, duplicates removed), 90% statistical errors. The earlier entries of the GCN table (April
  and May 2019, several hundred) are before the stream started and are left out: O3a events before June 2019 are
  not covered.
- ``icecube-events``: the IceCube point-source event list (HEASARC ``icecubepsc``, 2008 to July 2018: O1, O2), all
  track events with their 1σ errors; queried per window.

Statistic, for a counterpart of localization σ at (ra, dec):

- ``credible_percent``: the greedy credible level of the GW map at its position;
- ``overlap`` I_Ω: Npix Σ p_GW p_ext, the GW probability seen by the counterpart's localization relative to a uniform
  sky (1 for an unrelated direction, Npix p_GW at a point);
- ``q``: the fraction of isotropic positions (same σ) with a larger overlap;
- ``p_chance`` = 1 − exp(−R Δt q): the probability that the external catalog, at its mean rate R, puts at least one
  counterpart this well placed in the window Δt by chance; for the best counterpart of an event and source it is
  the event's chance probability. The summary counts the maps below 1% against the number expected, and corrects
  the smallest probability for the number of maps searched;
- ``raven_far_per_yr`` = FAR_GW R Δt / I_Ω, RAVEN's coincidence false-alarm rate, when the GW FAR is known.
"""
from __future__ import annotations

import hashlib
import html as _html
import io
import json
import os
import re
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Optional

import healpy as hp
import numpy as np
import pandas as pd

NSIDE = 256                            # working resolution of the GW maps (0.23° pixels)
WINDOWS = {"grb": (-1.0, 5.0), "neutrino": (-500.0, 500.0)}     # s around the merger (RAVEN's for GRBs; IceCube's)
SOURCES = {"grb": ("gbm", "swift"), "neutrino": ("icecube-alerts", "icecube-events")}
GBM_SYSTEMATIC_DEG = 3.7
R90_TO_SIGMA = 1.0 / np.sqrt(-2.0 * np.log(0.1))                # 2D Gaussian: r90 = 2.146 σ
TAP_URL = "https://heasarc.gsfc.nasa.gov/xamin/vo/tap/sync"
AMON_URL = "https://gcn.gsfc.nasa.gov/amon_icecube_gold_bronze_events.html"
CACHE_MAX_AGE_DAYS = 7.0               # the catalogs grow: refreshed weekly
P_NOTABLE = 0.01                       # chance probability below which a coincidence is counted in the summary


def _log(msg: str) -> None:
    print(f"[coincidence] {msg}", flush=True)


def cache_dir() -> Path:
    return Path(os.environ.get("GWTC_COINCIDENCE_CACHE", Path.home() / ".cache_gwtc_analysis" / "coincidence"))


def _cached(name: str, fetch, max_age_days: float = CACHE_MAX_AGE_DAYS, refresh: bool = False) -> pd.DataFrame:
    p = cache_dir() / f"{name}.csv"
    if p.exists() and not refresh and (time.time() - p.stat().st_mtime) < max_age_days * 86400:
        return pd.read_csv(p)
    df = fetch()
    p.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(p.with_suffix(".part"), index=False)
    p.with_suffix(".part").replace(p)
    return df


def tap_query(adql: str, timeout: float = 300.0) -> pd.DataFrame:
    """A HEASARC TAP query (VOTable) as a DataFrame."""
    from astropy.io.votable import parse_single_table

    from . import net

    r = net.post(TAP_URL, data=dict(REQUEST="doQuery", LANG="ADQL", QUERY=adql), timeout=timeout)
    r.raise_for_status()
    t = parse_single_table(io.BytesIO(r.content)).to_table()
    df = t.to_pandas()
    for c in df.columns:
        if df[c].dtype == object:
            df[c] = df[c].apply(lambda v: v.decode() if isinstance(v, bytes) else v)
    return df


def _gps_from_mjd(mjd) -> np.ndarray:
    from astropy.time import Time

    return np.asarray(Time(np.asarray(mjd, float), format="mjd", scale="utc").gps, float)


def _gps_from_iso(iso: Iterable[str]) -> np.ndarray:
    from astropy.time import Time

    iso = list(iso)
    out = np.full(len(iso), np.nan)
    ok = [i for i, s in enumerate(iso) if isinstance(s, str) and s]
    if ok:
        out[ok] = Time([iso[i] for i in ok], format="isot", scale="utc").gps
    return out


# ---------------------------------------------------------------------------
# catalogs: name, gps, ra, dec, sigma_deg (+ extra columns); one table per source
# ---------------------------------------------------------------------------
def load_gbm(refresh: bool = False) -> pd.DataFrame:
    def fetch():
        d = tap_query("SELECT name, trigger_time, ra, dec, error_radius, t90 FROM fermigbrst")
        return pd.DataFrame(dict(name=d["name"], gps=_gps_from_mjd(d["trigger_time"]), ra=d["ra"], dec=d["dec"],
                                 error_stat_deg=d["error_radius"], t90_s=d["t90"]))

    d = _cached("fermi_gbm", fetch, refresh=refresh)
    stat = d["error_stat_deg"].fillna(0.0)
    # a zero error is a position taken from a localized counterpart: no GBM systematic then
    d["sigma_deg"] = np.where(stat > 0, np.hypot(stat, GBM_SYSTEMATIC_DEG), 0.0)
    return d


def load_swift(refresh: bool = False) -> pd.DataFrame:
    def fetch():
        d = tap_query("SELECT name, trigger_time, trigger_obs, ra, dec, error_radius, xrt_ra, xrt_dec, "
                      "xrt_error_radius, uvot_ra, uvot_dec, uvot_error_radius, redshift, host_galaxy, bat_t90 "
                      "FROM swiftgrbba")
        # trigger_time is the time of day; the date is in the name (GRB YYMMDD[A-Z])
        iso = []
        for n, t in zip(d["name"], d["trigger_time"]):
            m = re.search(r"(\d{2})(\d{2})(\d{2})", str(n))
            ok = m and isinstance(t, str) and re.match(r"\d{2}:\d{2}:\d{2}", t)
            iso.append(f"20{m.group(1)}-{m.group(2)}-{m.group(3)}T{t.strip()}" if ok else "")
        return d.assign(gps=_gps_from_iso(iso))

    d = _cached("swift", fetch, refresh=refresh)
    # the best 90% position: UVOT (arcsec), XRT (arcsec), BAT (arcmin)
    ra, dec, r90 = d["ra"].astype(float), d["dec"].astype(float), d["error_radius"].astype(float) / 60.0
    for p, unit in (("xrt", 3600.0), ("uvot", 3600.0)):
        ok = d[f"{p}_ra"].notna() & d[f"{p}_error_radius"].notna()
        ra, dec = ra.where(~ok, d[f"{p}_ra"]), dec.where(~ok, d[f"{p}_dec"])
        r90 = r90.where(~ok, d[f"{p}_error_radius"].astype(float) / unit)
    out = pd.DataFrame(dict(name=d["name"], gps=d["gps"], ra=ra, dec=dec, sigma_deg=r90.fillna(0.0) * R90_TO_SIGMA,
                            redshift=pd.to_numeric(d["redshift"], errors="coerce"), host_galaxy=d["host_galaxy"],
                            t90_s=d["bat_t90"]))
    return out.dropna(subset=["gps", "ra", "dec"]).reset_index(drop=True)


def _parse_amon(text: str) -> pd.DataFrame:
    rows = []
    for r in re.findall(r"<tr[^>]*>(.*?)</tr>", text, flags=re.S | re.I):
        c = [_html.unescape(re.sub(r"<[^>]+>", "", x)).strip()
             for x in re.findall(r"<t[dh][^>]*>(.*?)</t[dh]>", r, flags=re.S | re.I)]
        if len(c) < 11 or not re.match(r"\d+_\d+", c[0]):
            continue
        rev = 0
        if not re.match(r"\d{2}/\d{2}/\d{2}", c[1]):           # newer rows carry a revision column
            rev, c = int(c[1]) if c[1].isdigit() else 0, [c[0]] + c[2:]
        date, tod, kind, ra, dec, e90 = c[1], c[2], c[3], c[4], c[5], c[6]
        yy, mm, dd = date.split("/")
        rows.append(dict(name=f"IceCube-{c[0]}", rev=rev, iso=f"20{yy}-{mm}-{dd}T{tod}", kind=kind.upper(),
                         ra=float(ra), dec=float(dec), error90_arcmin=float(e90),
                         energy_tev=float(c[8]) if len(c) > 8 and c[8] else np.nan,
                         signalness=float(c[9]) if len(c) > 9 and c[9] else np.nan))
    d = pd.DataFrame(rows)
    d = d.sort_values("rev").groupby("name", as_index=False).last()
    d = d.assign(gps=_gps_from_iso(d["iso"])).drop(columns=["iso"])
    # the same alert is sometimes listed under two event numbers: one row per time and position
    key = d["gps"].round(0).astype(str) + d["ra"].round(2).astype(str) + d["dec"].round(2).astype(str)
    d = d.loc[~key.duplicated() & (d["gps"] >= ICECUBE_ALERTS_START_GPS)]
    return d.sort_values("gps").reset_index(drop=True)


ICECUBE_ALERTS_START_GPS = 1244764818.0       # 2019-06-17 00:00 UTC, the start of the Gold and Bronze alert stream


def load_icecube_alerts(refresh: bool = False) -> pd.DataFrame:
    def fetch():
        from . import net

        r = net.get(AMON_URL, timeout=120)
        r.raise_for_status()
        return _parse_amon(r.text)

    d = _cached("icecube_alerts", fetch, refresh=refresh)
    return d.assign(sigma_deg=d["error90_arcmin"] / 60.0 * R90_TO_SIGMA)


ICECUBE_EVENTS_END_GPS = 1215043218.0          # 2018-07-08, the end of the point-source event list


def icecube_events(gps0: float, t0: float, t1: float) -> tuple[pd.DataFrame, float]:
    """The IceCube point-source events within [gps0 + t0, gps0 + t1], and their local rate (per second, from the day
    around): their rate (about 1e5 per year) varies with the detector configuration."""
    from astropy.time import Time

    mjd = lambda g: float(Time(g, format="gps").mjd)
    key = f"icecube_events_{gps0:.0f}_{t0:g}_{t1:g}"

    def fetch():
        d = tap_query(f"SELECT time, ra, dec, error_radius, event_energy FROM icecubepsc WHERE time > "
                      f"{mjd(gps0 - 43200):.8f} AND time < {mjd(gps0 + 43200):.8f}")
        return d

    d = _cached(key, fetch, max_age_days=1e6)
    rate = len(d) / 86400.0
    d = d.assign(gps=_gps_from_mjd(d["time"]))
    d = d[(d["gps"] >= gps0 + t0) & (d["gps"] <= gps0 + t1)]
    out = pd.DataFrame(dict(name=[f"IceCube-psc-{g:.2f}" for g in d["gps"]], gps=d["gps"], ra=d["ra"], dec=d["dec"],
                            sigma_deg=d["error_radius"], energy_gev=d["event_energy"]))
    return out.reset_index(drop=True), rate


def catalog_rate(d: pd.DataFrame) -> float:
    """Mean rate of a catalog, per second, over its span."""
    g = d["gps"].dropna()
    return len(g) / max(g.max() - g.min(), 1.0)


# ---------------------------------------------------------------------------
# GW sky map
# ---------------------------------------------------------------------------
@dataclass
class GWMap:
    prob: np.ndarray                       # RING, NSIDE, sums to 1
    gps: Optional[float]
    distmu: Optional[np.ndarray] = None    # RING: the line-of-sight distance ansatz of the map, Mpc
    distsigma: Optional[np.ndarray] = None
    credible: np.ndarray = field(init=False)
    _cache: dict = field(init=False, default_factory=dict)

    def __post_init__(self):
        o = np.argsort(self.prob)[::-1]
        c = np.empty_like(self.prob)
        c[o] = np.cumsum(self.prob[o])
        self.credible = c

    @property
    def nside(self) -> int:
        return hp.npix2nside(len(self.prob))

    def area_deg2(self, level: float = 0.9) -> float:
        return float((self.credible <= level).sum() * hp.nside2pixarea(self.nside, degrees=True))

    def credible_percent(self, ra: float, dec: float) -> float:
        return 100.0 * float(self.credible[hp.ang2pix(self.nside, np.radians(90 - dec), np.radians(ra))])

    def distance_quantiles(self, ra: float, dec: float, q=(0.05, 0.5, 0.95)) -> Optional[np.ndarray]:
        """Quantiles (Mpc) of the distance along the line of sight of (ra, dec): p(r) ∝ r² N(r; μ, σ), r > 0, the
        ansatz of the 3D sky maps; None without distance layers."""
        if self.distmu is None:
            return None
        i = hp.ang2pix(self.nside, np.radians(90 - dec), np.radians(ra))
        mu, sig = float(self.distmu[i]), float(self.distsigma[i])
        if not (np.isfinite(mu) and np.isfinite(sig)) or sig <= 0:
            return None
        r = np.linspace(1e-3, max(mu, 0) + 8 * sig, 20000)
        cdf = np.cumsum(r ** 2 * np.exp(-0.5 * ((r - mu) / sig) ** 2))
        return np.interp(q, cdf / cdf[-1], r)

    def _smoothed(self, sigma_deg: float) -> np.ndarray:
        """The map convolved with a Gaussian of width σ: its value at a position, times Npix, is the overlap of a
        counterpart localized there with that error (the map itself below half a pixel)."""
        key = round(float(sigma_deg or 0.0), 3)
        if key < 0.5 * np.degrees(hp.nside2resol(self.nside)):
            return self.prob
        if key not in self._cache:
            import warnings

            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                s = np.clip(hp.smoothing(self.prob, sigma=np.radians(key)), 0, None)
            self._cache[key] = s / s.sum()
        return self._cache[key]

    def overlap(self, ra: float, dec: float, sigma_deg: float) -> float:
        """I_Ω = Npix Σ p_GW p_ext: 1 for an unrelated direction, Npix p_GW for a point."""
        s = self._smoothed(sigma_deg)
        return len(s) * float(s[hp.ang2pix(self.nside, np.radians(90 - dec), np.radians(ra))])

    def q(self, overlap: float, sigma_deg: float) -> float:
        """Fraction of the sky (isotropic positions, same σ) with an overlap at least `overlap`; exact."""
        s = self._smoothed(sigma_deg)
        if overlap < 1e-6:                     # outside the map (smoothing round-off): every position is as good
            return 1.0
        return float(np.mean(len(s) * s >= overlap * (1 - 1e-9)))


def read_gw_map(path, nside: int = NSIDE) -> GWMap:
    """A GW sky map (FITS, flat or multi-order) at `nside`, with its merger GPS time from the header."""
    import warnings

    from astropy.io import fits
    from ligo.skymap.io import read_sky_map

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        if str(fits.getheader(str(path), 1).get("ORDERING", "")).upper() == "NUNIQ":
            from ligo.skymap.bayestar import rasterize

            s = read_sky_map(str(path), moc=True)
            cols = [c for c in ("UNIQ", "PROBDENSITY", "DISTMU", "DISTSIGMA", "DISTNORM") if c in s.colnames]
            r = rasterize(s[cols], order=int(np.log2(nside)))          # NESTED, probability per pixel
            p = hp.reorder(np.asarray(r["PROB"], float), n2r=True)
            dist = [hp.reorder(np.asarray(r[c], float), n2r=True) for c in ("DISTMU", "DISTSIGMA")] \
                if "DISTMU" in r.colnames else [None, None]
            meta = s.meta
        else:
            try:
                (m, mu, sig, _), meta = read_sky_map(str(path), moc=False, nest=False, distances=True)
                dist = [hp.ud_grade(np.where(np.isfinite(x), x, np.nan), nside) for x in (mu, sig)]
            except (KeyError, ValueError):                             # a map without distance layers
                m, meta = read_sky_map(str(path), moc=False, nest=False)
                dist = [None, None]
            p = hp.ud_grade(np.asarray(m, float), nside, power=-2)
    p = np.clip(p, 0, None)
    return GWMap(p / p.sum(), meta.get("gps_time"), dist[0], dist[1])


# ---------------------------------------------------------------------------
# search
# ---------------------------------------------------------------------------
@dataclass
class Catalogs:
    """The external catalogs of the search, loaded once."""
    sources: tuple[str, ...]
    refresh: bool = False
    tables: dict = field(default_factory=dict)
    rates: dict = field(default_factory=dict)

    def __post_init__(self):
        loaders = {"gbm": load_gbm, "swift": load_swift, "icecube-alerts": load_icecube_alerts}
        for s in self.sources:
            if s in loaders:
                t = loaders[s](refresh=self.refresh)
                self.tables[s], self.rates[s] = t, catalog_rate(t)
                _log(f"{s}: {len(t)} entries, {self.rates[s] * 3.156e7:.0f} per year")

    def in_window(self, source: str, gps0: float, t0: float, t1: float) -> tuple[pd.DataFrame, float]:
        if source == "icecube-events":
            if gps0 > ICECUBE_EVENTS_END_GPS:
                return pd.DataFrame(), np.nan
            return icecube_events(gps0, t0, t1)
        t = self.tables[source]
        return t[(t["gps"] >= gps0 + t0) & (t["gps"] <= gps0 + t1)].reset_index(drop=True), self.rates[source]

    def covers(self, source: str, gps0: float) -> bool:
        if source == "icecube-events":
            return gps0 <= ICECUBE_EVENTS_END_GPS
        g = self.tables[source]["gps"]
        return g.min() <= gps0 <= g.max()


MESSENGER = {s: m for m, ss in SOURCES.items() for s in ss}
C_KM_S, MPC_KM = 299792.458, 3.0856775814913673e19
EMISSION_LAG_S = (0.0, 10.0)           # gamma rays emitted 0-10 s after the merger (the GW170817 paper's assumption)


def speed_of_gravity(gw: GWMap, ra: float, dec: float, dt_s: float, lag_s=EMISSION_LAG_S) -> dict:
    """Bounds on (v_GW - c)/c from a GRB arriving dt_s after the merger, if they are the same source: the gamma
    rays left between lag_s[0] and lag_s[1] s after the merger and travelled D, the 5% quantile of the GW distance
    along their line of sight (conservative), so (v_GW - c)/c lies in [c (dt - lag_max) / D, c (dt - lag_min) / D]
    (LVK, Fermi-GBM and INTEGRAL, ApJL 848, L13: -3e-15 to +7e-16 for GW170817, 1.74 s, D = 26 Mpc)."""
    qd = gw.distance_quantiles(ra, dec)
    if qd is None or qd[0] <= 0:
        return {}
    k = C_KM_S / (qd[0] * MPC_KM)          # per second of delay
    return dict(distance_los_mpc_5=float(qd[0]), distance_los_mpc_50=float(qd[1]),
                gw_speed_dv_over_c_min=float(k * (dt_s - lag_s[1])), gw_speed_dv_over_c_max=float(k * (dt_s - lag_s[0])))


def search_event(gw: GWMap, event: str, catalogs: Catalogs, windows: dict, far_gw: Optional[float] = None) -> tuple[list[dict], list[dict]]:
    """The counterparts of one GW map within the windows: (one row per counterpart, one row per source)."""
    cands, summary = [], []
    for source in catalogs.sources:
        messenger = MESSENGER[source]
        t0, t1 = windows[messenger]
        if gw.gps is None or not catalogs.covers(source, gw.gps):
            summary.append(dict(event=event, source=source, messenger=messenger, covered=False))
            continue
        found, rate = catalogs.in_window(source, gw.gps, t0, t1)
        dt = t1 - t0
        best = None
        for _, c in found.iterrows():
            sig = float(c["sigma_deg"]) if pd.notna(c["sigma_deg"]) else 0.0
            ov = gw.overlap(float(c["ra"]), float(c["dec"]), sig)
            q = gw.q(ov, sig)
            row = dict(event=event, messenger=messenger, source=source, counterpart=c["name"],
                       dt_s=float(c["gps"] - gw.gps), ra=float(c["ra"]), dec=float(c["dec"]), sigma_deg=sig,
                       credible_percent=gw.credible_percent(float(c["ra"]), float(c["dec"])), overlap=ov, q=q,
                       p_chance=float(-np.expm1(-rate * dt * q)),
                       raven_far_per_yr=(far_gw * rate * dt / ov) if far_gw and ov > 0 else np.nan)
            if messenger == "grb":
                row.update(speed_of_gravity(gw, row["ra"], row["dec"], row["dt_s"]))
            for k in ("redshift", "host_galaxy", "energy_tev", "energy_gev", "signalness", "kind", "t90_s"):
                if k in c and pd.notna(c[k]):
                    row[k] = c[k]
            cands.append(row)
            if best is None or ov > best["overlap"]:
                best = row
        summary.append(dict(event=event, source=source, messenger=messenger, covered=True, window_s=f"{t0:g}/{t1:g}",
                            rate_per_window=rate * dt, n_in_window=len(found),
                            best=best["counterpart"] if best else "", best_credible_percent=best["credible_percent"]
                            if best else np.nan, best_overlap=best["overlap"] if best else np.nan,
                            p_chance=best["p_chance"] if best else np.nan))
    return cands, summary


# ---------------------------------------------------------------------------
# driver: the catalog skymaps (Zenodo releases) or given maps
# ---------------------------------------------------------------------------
def _gw_fars(catalogs: Iterable[str]) -> list[tuple[float, float, str]]:
    """(GPS, FAR per year, name) of the GWOSC events of the catalogs, to attach the GW FAR to a map by its time."""
    from .catalog_registry import gwosc_lists, products_catalog
    from .gw_stat import fetch_gwtc_events

    cats = [c for c in catalogs if c != "ALL"]
    out = []
    for lst in gwosc_lists(cats if cats and "ALL" not in catalogs else None):
        try:
            evs = (fetch_gwtc_events(lst) or {}).get("events", {})
        except Exception as e:                       # the FAR is optional
            _log(f"WARN: no FARs from {lst}: {type(e).__name__}: {e}")
            continue
        for name, v in evs.items():
            if v.get("GPS") is not None:
                out.append((float(v["GPS"]), float(v["far"]) if v.get("far") is not None else np.nan,
                            v.get("commonName") or name))
    return out


def _far_for(gps: Optional[float], fars: list) -> float:
    if gps is None or not fars:
        return np.nan
    g = np.array([f[0] for f in fars])
    i = int(np.argmin(np.abs(g - gps)))
    return fars[i][1] if abs(g[i] - gps) < 2.0 else np.nan


def _iter_catalog_maps(catalogs: list[str], skymap_label: str, zenodo_versions: Optional[dict], tmp: Path):
    """(event, catalog, path of a temporary FITS file) for one map per event of the catalogs' skymap releases."""
    import tarfile

    from .gw_stat import build_skymap_index_from_tar, download_zenodo_skymaps_tarball, select_skymap_member
    from .search_skymaps import _expand_catalogs_for_skymaps, _skymap_sources

    cats, only = _skymap_sources(_expand_catalogs_for_skymaps(list(catalogs)))
    prefer = "Mixed" if skymap_label.lower() == "any" else skymap_label
    tmp.mkdir(parents=True, exist_ok=True)
    for catalog in cats:
        tar_path = download_zenodo_skymaps_tarball(catalog, progress=True, verbose=False,
                                                   version=(zenodo_versions or {}).get(products_catalog(catalog)))
        index = build_skymap_index_from_tar(tar_path, verbose=False)
        events = sorted({ev for ev, _ in index})
        if catalog in only:
            from .search_skymaps import _event_selected

            events = [e for e in events if _event_selected(e, only[catalog])]
        with tarfile.open(tar_path, mode="r:gz") as tf:
            for ev in events:
                member = select_skymap_member(index, ev, prefer=prefer)
                if member is None:
                    continue
                f = tf.extractfile(member)
                if f is None:
                    continue
                p = tmp / Path(member).name
                p.write_bytes(f.read())
                try:
                    yield ev, catalog, p
                finally:
                    p.unlink(missing_ok=True)


def _plot_coincidence(map_path: Path, event: str, rows: list[dict], prob: float, plots: Path) -> Optional[Path]:
    """The GW map with the counterparts inside its `prob` credible region."""
    from .read_skymap import plot_skymap_with_ra_dec
    from .search_skymaps import _normalize_plot_output, _safe_filename

    inside = [r for r in rows if r["credible_percent"] <= 100 * prob]
    if not inside:
        return None
    best = min(inside, key=lambda r: r["p_chance"])
    title = _safe_filename(f"{event}_{best['counterpart']}")
    try:
        produced = plot_skymap_with_ra_dec(str(map_path), title, best["ra"], best["dec"], "grey",
                                           contour_levels=(50, 90))
        out = _normalize_plot_output(produced, plots / f"{title}.png")
        return Path(out) if out else None
    except Exception as e:
        _log(f"WARN: plot of {event} failed: {type(e).__name__}: {e}")
        return None


def run_coincidences(*, catalogs: list[str] = (), skymaps: list[str] = (), messengers=("grb", "neutrino"),
                     sources: Optional[list[str]] = None, windows: Optional[dict] = None, prob: float = 0.9,
                     skymap_label: str = "Mixed", zenodo_versions: Optional[dict] = None,
                     out_tsv: Optional[str] = "coincidences.tsv", out_report: Optional[str] = "coincidences.html",
                     plots_dir: str = "coincidence_plots", refresh: bool = False) -> dict:
    """GRB and neutrino coincidences of the catalog skymaps and/or the given skymap files (see the module docstring).
    Writes `out_tsv` (one row per counterpart in a window), <out_tsv>.events.tsv (one row per map and source) and
    the report; returns the summary numbers."""
    from .report import write_simple_html_report

    windows = {**WINDOWS, **(windows or {})}
    sources = tuple(sources or [s for m in messengers for s in SOURCES[m]])
    bad = [s for s in sources if s not in MESSENGER]
    if bad:
        raise ValueError(f"unknown source(s) {bad}; choose among {list(MESSENGER)}")
    cats = Catalogs(sources, refresh=refresh)
    fars = _gw_fars(catalogs)              # by time: the catalog maps, and given maps of catalog events
    plots = Path(plots_dir)
    plots.mkdir(parents=True, exist_ok=True)
    cand_rows, ev_rows, images = [], [], []

    def one(event, catalog, path):
        gw = read_gw_map(path)
        c, s = search_event(gw, event, cats, windows, far_gw=_far_for(gw.gps, fars))
        for r in c + s:
            r["catalog"] = catalog
        for r in s:
            r["gps"], r["area90_deg2"] = gw.gps, gw.area_deg2()
        cand_rows.extend(c)
        ev_rows.extend(s)
        png = _plot_coincidence(Path(path), event, c, prob, plots)
        if png:
            images.append(png)

    for path in skymaps:
        ev = re.search(r"(GW\d{6}(_\d{6})?|S\d{6}[a-z]+)", Path(path).name)
        one(ev.group(1) if ev else Path(path).stem, "file", path)
    if catalogs:
        n = 0
        for ev, catalog, path in _iter_catalog_maps(list(catalogs), skymap_label, zenodo_versions,
                                                    plots / ".tmp_skymaps"):
            try:
                one(ev, catalog, path)
            except Exception as e:
                _log(f"WARN: {ev}: {type(e).__name__}: {e}")
                ev_rows.append(dict(event=ev, catalog=catalog, error=f"{type(e).__name__}: {e}"))
            n += 1
            if n % 25 == 0:
                _log(f"{n} maps searched")

    cand = pd.DataFrame(cand_rows)
    evs = pd.DataFrame(ev_rows)
    if out_tsv:
        cand.to_csv(out_tsv, sep="\t", index=False, float_format="%.6g")
        evs.to_csv(Path(out_tsv).with_suffix(".events.tsv"), sep="\t", index=False, float_format="%.6g")
    result = {}
    paras = []
    if not evs.empty and "covered" in evs:
        cov = evs[evs["covered"] == True]                      # noqa: E712
        for s in sources:
            e = cov[cov["source"] == s]
            if e.empty:
                continue
            p = e["p_chance"].dropna()
            n_in = int((e["n_in_window"] > 0).sum())
            n_inside = int((e["best_credible_percent"] <= 100 * prob).sum())
            n_sig = int((p < P_NOTABLE).sum())
            p_min = float(p.min()) if len(p) else np.nan
            # the smallest chance probability among the maps searched, corrected for their number
            p_global = float(-np.expm1(len(e) * np.log1p(-p_min))) if len(p) else np.nan
            result[s] = dict(maps=len(e), with_counterpart_in_window=n_in, inside_credible=n_inside,
                             below_p_notable=n_sig, expected_below=len(e) * P_NOTABLE, p_min=p_min, p_global=p_global)
            paras.append(f"<b>{s}</b> ({MESSENGER[s]}, window {windows[MESSENGER[s]][0]:g} to "
                         f"{windows[MESSENGER[s]][1]:g} s): {len(e)} GW maps covered, {n_in} with a counterpart in the "
                         f"window, {n_inside} of them inside the {100 * prob:g}% credible region; {n_sig} with a chance "
                         f"probability below {P_NOTABLE:g} ({len(e) * P_NOTABLE:.2f} expected by chance)"
                         + (f"; smallest {p_min:.2g}, {p_global:.2g} after the {len(e)} trials." if len(p) else "."))
    if not cand.empty and "gw_speed_dv_over_c_min" in cand:
        grb = cand[(cand["messenger"] == "grb") & (cand["credible_percent"] <= 100 * prob)
                   & cand["gw_speed_dv_over_c_min"].notna()]
        for r in grb.itertuples():
            result.setdefault("speed_of_gravity", []).append(dict(event=r.event, counterpart=r.counterpart,
                                                                  dv_over_c=(r.gw_speed_dv_over_c_min,
                                                                             r.gw_speed_dv_over_c_max)))
            paras.append(f"<b>Speed of gravity</b>, if {r.counterpart} and {r.event} are the same source: the burst "
                         f"arrived {r.dt_s:.2f} s after the merger (catalog trigger time) from at least "
                         f"{r.distance_los_mpc_5:.0f} Mpc (5% of the GW distance along its line of sight); with the "
                         f"gamma rays emitted {EMISSION_LAG_S[0]:g}–{EMISSION_LAG_S[1]:g} s after the merger, "
                         f"{r.gw_speed_dv_over_c_min:.1e} ≤ (v<sub>GW</sub> − c)/c ≤ {r.gw_speed_dv_over_c_max:.1e} "
                         "(GW170817 paper: −3 × 10⁻¹⁵ to +7 × 10⁻¹⁶ with 1.74 s, the delay to the burst onset, "
                         "and 26 Mpc).")
    if out_report:
        top = cand.sort_values("p_chance").head(30) if not cand.empty else cand
        cols = [c for c in ("event", "source", "counterpart", "dt_s", "credible_percent", "overlap", "p_chance",
                            "raven_far_per_yr", "sigma_deg", "redshift", "host_galaxy") if c in top]
        intro = ("Gamma-ray bursts (Fermi-GBM, Swift) and neutrinos (IceCube) within a time window of each GW merger, "
                 "scored by the overlap of their localization with the GW sky map (an offline version of the LVK's "
                 "RAVEN). p_chance: the probability that the catalog puts a counterpart this well placed in the "
                 "window by chance. A Swift burst with a redshift inside the GW map would be a bright-siren "
                 "candidate.")
        write_simple_html_report(out_report, title="Multi-messenger coincidences", paragraphs=[intro] + paras,
                                 images=[str(i) for i in images],
                                 tables=[("Counterparts in the windows, most significant first",
                                          top[cols].to_html(index=False, float_format=lambda x: f"{x:.3g}",
                                                            na_rep=""))] if not top.empty else [])
        _log(f"report written to {out_report}")
    return result
