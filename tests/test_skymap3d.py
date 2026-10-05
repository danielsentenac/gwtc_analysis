"""Offline tests of the 3D sky maps: file names, catalog routing, a synthetic map and its galaxy ranking."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from gwtc_analysis import skymap3d as s3

pytest.importorskip("ligo.skymap")


def test_member_names():
    g5 = ["parameter_estimation/skymaps/IGWN-GWTC5p0-29ebe06b7_25-GW240615_113620-IMRPhenomXPHM_SpinTaylor_Skymap_PEDataRelease.fits.gz",
          "parameter_estimation/skymaps/IGWN-GWTC5p0-29ebe06b7_25-GW240615_113620-IMRPhenomXPHM_Skymap_PEDataRelease.fits.gz",
          "parameter_estimation/skymaps/IGWN-GWTC5p0-29ebe06b7_25-GW240615_113620-NRSur7dq4_Skymap_PEDataRelease.fits.gz",
          "parameter_estimation/skymaps/IGWN-GWTC5p0-29ebe06b7_25-GW240615_160735-NRSur7dq4_Skymap_PEDataRelease.fits.gz"]
    assert s3.member_waveform(g5[0]) == "IMRPhenomXPHM_SpinTaylor"
    assert s3.select_member(g5, "GW240615_113620", "C00:IMRPhenomXPHM-SpinTaylor") == g5[0]
    assert s3.select_member(g5, "GW240615_113620", "C00:IMRPhenomXPHM") == g5[1]
    assert s3.select_member(g5, "GW240615_160735", "C00:IMRPhenomXPHM-SpinTaylor") == g5[3]
    assert s3.select_member(g5, "GW240616_000000", None) is None
    g21 = ["IGWN-GWTC2p1-v2-GW150914_095045_PEDataRelease_cosmo_reweight_C01:IMRPhenomXPHM.fits",
           "IGWN-GWTC2p1-v2-GW150914_095045_PEDataRelease_cosmo_reweight_C01:Mixed.fits"]
    assert s3.member_waveform(g21[0]) == "IMRPhenomXPHM"
    assert s3.member_waveform(g21[1].replace(":", "_")) == "Mixed"         # cached name
    assert s3.select_member(g21, "GW150914_095045", "C01:Mixed") == g21[1]
    assert s3.select_member(g21, "GW150914_095045", "C01:SEOBNRv4PHM") == g21[1]   # no such map: Mixed


def test_catalog_routing():
    assert s3.catalog_of_pe_file("IGWN-GWTC5p0-29ebe06b7_25-GW240615_113620-combined_PEDataRelease.hdf5") == "GWTC-5"
    assert s3.catalog_of_pe_file("IGWN-GWTC4p1-18965dda8_5-GW230529_181500-combined_PEDataRelease.hdf5") == "GWTC-4.1"
    assert s3.catalog_of_pe_file("IGWN-GWTC3p0-v3-GW200129_065458_PEDataRelease_mixed_cosmo.h5") == "GWTC-3"
    assert s3.catalog_of_pe_file("GW170817_GWTC-1.hdf5") is None
    assert s3.catalog_of_event("GW150914_095045") == "GWTC-2.1"     # GWTC-1 events: GWTC-2.1 release
    assert s3.catalog_of_event("GW200129_065458") == "GWTC-3"
    assert s3.catalog_of_event("GW240615_113620") == "GWTC-5"


def _mock_map(tmp_path, ra0=40.0, dec0=20.0, width=4.0, dmean=400.0, dstd=60.0):
    """A 3D map at nside 32: a Gaussian blob on the sky, the same distance distribution in every direction."""
    import astropy_healpix as ah
    import astropy.units as u
    from astropy.coordinates import SkyCoord
    from astropy.table import Table
    from ligo.skymap.distance import moments_to_parameters
    from ligo.skymap.io import write_sky_map

    level = 5
    nside = ah.level_to_nside(level)
    ipix = np.arange(12 * nside ** 2)
    lon, lat = ah.healpix_to_lonlat(ipix, nside, order="nested")
    sep = SkyCoord(lon, lat).separation(SkyCoord(ra0 * u.deg, dec0 * u.deg)).deg
    prob = np.exp(-0.5 * (sep / width) ** 2)
    prob /= prob.sum()
    mu, sigma, norm = moments_to_parameters(dmean, dstd)
    t = Table(dict(UNIQ=ah.level_ipix_to_uniq(level, ipix),
                   PROBDENSITY=prob / ah.nside_to_pixel_area(nside).to_value(u.sr),
                   DISTMU=np.full(len(ipix), float(mu)), DISTSIGMA=np.full(len(ipix), float(sigma)),
                   DISTNORM=np.full(len(ipix), float(norm))))
    t.meta.update(distmean=dmean, diststd=dstd)
    path = tmp_path / "mock_skymap.fits"
    write_sky_map(str(path), t)
    return path


def test_summary_and_ranking(tmp_path):
    m = s3.read_moc(_mock_map(tmp_path))
    s = s3.summarize(m)
    # a 2D Gaussian of width w: 90% radius w*sqrt(2 ln 10), area pi r^2
    assert s["area90"] == pytest.approx(np.pi * (4.0 * np.sqrt(2 * np.log(10))) ** 2, rel=0.1)
    assert s["is_3d"] and s["vol90"] > 0 and s["ra_peak"] == pytest.approx(40, abs=1) and s["dec_peak"] == pytest.approx(20, abs=1)
    cones = s3.region_cones(m, max_cones=40)
    assert 1 <= len(cones) <= 40
    near = min(cones, key=lambda c: (c[0] - 40) ** 2 + (c[1] - 20) ** 2)
    assert abs(near[0] - 40) < near[2] and abs(near[1] - 20) < near[2]
    gal = pd.DataFrame(dict(ra=[40.0, 42.0, 40.0, 70.0], dec=[20.0, 21.0, 20.0, -10.0], dist=[400.0, 400.0, 900.0, 400.0]))
    r = s3.rank_galaxies(m, gal)
    assert list(r.index) == [0, 1, 2, 3] and r.loc[0, "ra"] == 40 and r.loc[0, "dist"] == 400
    assert r["host_share"].sum() == pytest.approx(1.0) and r.loc[0, "searched_prob_vol"] < 0.1
    assert r.iloc[-1]["searched_prob_vol"] > 0.99


def test_run_with_catalog_file(tmp_path):
    path = _mock_map(tmp_path)
    cat = tmp_path / "gal.csv"
    from astropy.cosmology import Planck15

    z = float(np.interp(400, Planck15.luminosity_distance(np.linspace(0.01, 0.2, 200)).value, np.linspace(0.01, 0.2, 200)))
    pd.DataFrame(dict(RA=[40.0, 45.0], Dec=[20.0, 21.0], redshift=[z, 0.15])).to_csv(cat, index=False)
    out = s3.run_skymap3d("GW000000_000000", "C00:Mock", tmp_path / "out", fits_path=path, galaxies=str(cat),
                          dist_samples=np.random.default_rng(0).normal(400, 60, 2000))
    names = [p.split("/")[-1] for p in out["files"]]
    assert names[0].endswith("_skymap3d.png") and "mock_skymap.fits" in names and "GW000000_000000_host_galaxies.tsv" in names
    ranked = pd.read_csv(tmp_path / "out" / "GW000000_000000_host_galaxies.tsv", sep="\t")
    assert ranked.loc[0, "ra"] == 40 and ranked.loc[0, "dist"] == pytest.approx(400, rel=0.01)
    assert "Host-galaxy candidates" in out["html"] and "credible volume" in out["html"]
    none = s3.run_skymap3d("GW000000_000000", None, tmp_path / "o2", fits_path=path, galaxies="none")
    assert "Host-galaxy" not in none["html"] and len(none["files"]) == 2


def test_cli_flags():
    from gwtc_analysis.cli import build_parser

    a = build_parser().parse_args(["parameters_estimation", "--src-name", "GW240615_113620"])
    assert a.skymap_3d and a.galaxies == "glade" and a.galaxy_max_area == 100
    a = build_parser().parse_args(["parameters_estimation", "--src-name", "X", "--no-skymap-3d", "--galaxies", "none"])
    assert not a.skymap_3d and a.galaxies == "none"
