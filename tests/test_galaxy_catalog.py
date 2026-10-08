"""Tests of the galaxy catalogs of the dark-siren analysis (`galaxies.py`, the `galaxy_catalog` mode and the
galaxy-catalog option of `hubble_constant`).

The icarogw part runs in the icarogw environment: the end-to-end test is run only when GWTC_ICAROGW_PYTHON points
to its Python interpreter.
"""
from __future__ import annotations

import json
import os
import subprocess
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import pytest

from gwtc_analysis import galaxies as gx
from gwtc_analysis import galaxy_catalog as gc
from gwtc_analysis import hubble_constant as hc
from gwtc_analysis.cli import build_parser

from test_hubble_constant import FAKE_LISTS, _write_mixture  # noqa: E402

ICAROGW_PYTHON = os.environ.get("GWTC_ICAROGW_PYTHON")


# ---------------------------------------------------------------------------
# galaxy files
# ---------------------------------------------------------------------------

def test_glade_sigmaz_and_frame():
    """sigmaz adds the peculiar-velocity and measurement errors in quadrature; quasars and z <= 0 are dropped;
    angles go to radians and the magnitude is Kmag."""
    assert gx.glade_sigmaz([0.003, np.nan, np.nan], [0.004, 0.015, np.nan])[:2] == pytest.approx([0.005, 0.015])
    assert np.isnan(gx.glade_sigmaz([np.nan], [np.nan])[0])
    raw = pd.DataFrame({"GLADE+": ["1", "2", "3", "4"], "Type": ["G", "Q", "G", "G"],
                        "RAJ2000": ["180", "10", "90", "45"], "DEJ2000": ["-30", "0", "45", "10"],
                        "Kmag": ["12.5", "13", "11", ""], "zcmb": ["0.05", "0.5", "-0.001", "0.1"],
                        "f_zcmb": ["1"] * 4, "e_z": ["0.001"] * 4, "e_zhelio": ["0.015"] * 4})
    df = gx.glade_kband_frame(raw)
    assert len(df) == 1                                   # quasar, negative redshift and missing Kmag dropped
    r = df.iloc[0]
    assert r["ra"] == pytest.approx(np.pi) and r["dec"] == pytest.approx(-np.pi / 6)
    assert r["m"] == 12.5 and r["sigmaz"] == pytest.approx(np.hypot(0.001, 0.015))


def test_fetch_glade_kband_bands_dedup_and_resume(tmp_path):
    """Declination bands are fetched once (cached), a galaxy on a band edge is kept once."""
    calls = []

    def fake(lo, hi):
        calls.append((lo, hi))
        rows = [] if lo > -86 else [("7", "G", "1.0", f"{hi}", "12", "0.0199", "0.02", "1", "0.001", "0.01")]
        if lo == -90:
            rows.append(("8", "G", "2.0", "-89", "11", "0.0299", "0.03", "1", "0.001", "0.01"))
        return pd.DataFrame(rows, columns=list(gx.GLADE_KBAND_COLS))

    out = gx.fetch_glade_kband(tmp_path / "g.h5", band_width=2.0, cache=tmp_path / "cache", fetch=fake)
    cols, attrs = gx.read_galaxy_file(out)
    assert len(cols["z"]) == 2 and attrs["band"] == "K-glade+"   # galaxy 7 on the -88 edge appears twice in VizieR
    assert len(calls) == 90
    gx.fetch_glade_kband(tmp_path / "g2.h5", band_width=2.0, cache=tmp_path / "cache", fetch=fake)
    assert len(calls) == 90                                       # every band from the cache


@pytest.mark.parametrize("fmt", ["csv", "parquet"])
def test_convert_catalog_in_chunks(tmp_path, fmt):
    """Any catalog: column mapping, degrees to radians, constant relative sigmaz, quality cut, chunked reading."""
    rng = np.random.default_rng(2)
    n = 1000
    df = pd.DataFrame({"RA": rng.uniform(0, 360, n), "DEC": rng.uniform(-60, 0, n), "ZPHOT": rng.uniform(0.05, 1, n),
                       "MAG_R": rng.uniform(18, 24, n), "CLASS": rng.uniform(0, 1, n)})
    src = tmp_path / f"cat.{fmt}"
    df.to_csv(src, index=False) if fmt == "csv" else df.to_parquet(src)
    out = gx.convert_catalog(src, tmp_path / "g.h5", band="r-upglade",
                             columns={"ra": "RA", "dec": "DEC", "z": "ZPHOT", "m": "MAG_R"}, sigmaz=0.05,
                             sigmaz_relative=True, where="CLASS > 0.5 and MAG_R < 23.9", chunk_rows=128)
    cols, attrs = gx.read_galaxy_file(out)
    keep = df.query("CLASS > 0.5 and MAG_R < 23.9")
    assert len(cols["z"]) == len(keep) and attrs["band"] == "r-upglade"
    assert cols["ra"] == pytest.approx(np.deg2rad(keep["RA"].to_numpy()))
    assert cols["sigmaz"] == pytest.approx(0.05 * (1 + keep["ZPHOT"].to_numpy()))
    with pytest.raises(ValueError, match="sigmaz"):
        gx.convert_catalog(src, tmp_path / "x.h5", band="r", columns={"ra": "RA", "dec": "DEC", "z": "ZPHOT", "m": "MAG_R"})


def test_query_names():
    assert gx._query_names("CLASS > 0.5 and `FLAGS GOLD` == 0 and MAG_R < 23.9") == {"CLASS", "FLAGS GOLD", "MAG_R"}


# ---------------------------------------------------------------------------
# settings, Slurm, CLI
# ---------------------------------------------------------------------------

def test_catalog_settings_defaults_are_gwtc4():
    s = gc.catalog_settings("K-glade+")
    assert (s["nside"], s["nside_mthr"], s["mthr_percentile"], s["epsilon"]) == (64, 32, 50.0, 1.0)
    assert s["grouping"] == "K-glade+" and s["subgrouping"] == "eps_1"
    assert s["outfile"] == "catalog_K-glade+_nside64_eps1.hdf5"
    assert gc.catalog_settings("K-glade+", epsilon=0.0)["subgrouping"] == "eps_0"
    with pytest.raises(ValueError):
        gc.catalog_settings("K-glade+", nside=32, nside_mthr=64)


def test_slurm_scripts(tmp_path):
    """Array jobs for the chunked stages, one job otherwise, chained by afterok; the runner copied next to them."""
    sub = gc.write_slurm_scripts(tmp_path, "/env/bin/python", ["shard", "pixels", "gather", "interpolate", "finish"],
                                 16, tmp_path / "g.h5", ["--partition=htc", "mem=4G"], "source /x/conda.sh",
                                 ["--mem=16G"])
    d = tmp_path / "slurm"
    assert (d / "dark_catalog_icarogw.py").exists()
    pixels = (d / "01_pixels.sh").read_text()
    assert "#SBATCH --array=0-15" in pixels and "--chunk $SLURM_ARRAY_TASK_ID --nchunks 16" in pixels
    assert "#SBATCH --partition=htc" in pixels and "#SBATCH --mem=4G" in pixels and "source /x/conda.sh" in pixels
    assert "--mem=16G" not in pixels
    finish = (d / "04_finish.sh").read_text()
    assert finish.index("#SBATCH --mem=4G") < finish.index("#SBATCH --mem=16G")      # the later option wins
    shard = (d / "00_shard.sh").read_text()
    assert "--array" not in shard and "--galaxies" in shard and "settings_input.json" in shard
    gather = (d / "02_gather.sh").read_text()
    assert "--nchunks 16" in gather and "--chunk" not in gather
    s = sub.read_text()
    assert s.count("sbatch --parsable") == 5 and "--dependency=afterok:$dep" in s


def test_cli_galaxy_catalog_arguments():
    a = build_parser().parse_args(["galaxy_catalog", "--input-catalog", "des.parquet", "--columns", "ra=RA", "dec=DEC",
                                   "z=Z", "m=MAG_R", "--band", "r-upglade", "--nintegration", "logspace:0.001:2000",
                                   "--executor", "slurm", "--jobs", "64", "--slurm-option=--partition=htc", "--submit"])
    assert a.columns == ["ra=RA", "dec=DEC", "z=Z", "m=MAG_R"] and a.nintegration == "logspace:0.001:2000"
    assert a.executor == "slurm" and a.jobs == 64 and a.submit and a.slurm_option == ["--partition=htc"]
    assert build_parser().parse_args(["galaxy_catalog"]).source == "glade-kband"


# ---------------------------------------------------------------------------
# hubble_constant with a galaxy catalog
# ---------------------------------------------------------------------------

def _fake_catalog_file(path: Path, **over) -> Path:
    st = gc.catalog_settings("K-glade+", **over)
    with h5py.File(path, "w") as h:
        h.create_group(st["grouping"]).create_group(st["subgrouping"])
        h.attrs["gwtc_analysis_settings"] = json.dumps(st)
    return path


@pytest.fixture
def fake_gwosc(monkeypatch):
    monkeypatch.setattr(hc.gw, "fetch_gwtc_events", lambda c: {"events": FAKE_LISTS.get(c, {})})


def test_prepare_inputs_with_galaxy_catalog(tmp_path, fake_gwosc):
    """With a galaxy catalog, inputs.h5 holds the sky position of the PE samples, isotropic positions for the
    injections and the catalog; extracts without the sky position are not used."""
    _write_mixture(tmp_path / "mix.hdf")
    cat = _fake_catalog_file(tmp_path / "catalog.hdf5")
    cache = tmp_path / "pe"
    (cache / "samples").mkdir(parents=True)
    rng = np.random.default_rng(0)
    for name in ("GW150914_095045", "GW230601_224134"):
        with h5py.File(cache / "samples" / f"{name}.h5", "w") as h:
            g = h.create_group("C00:IMRPhenomXPHM-SpinTaylor")
            for k, v in (("mass_1", rng.uniform(30, 40, 500)), ("mass_2", rng.uniform(20, 30, 500)),
                         ("luminosity_distance", rng.uniform(300, 2000, 500)), ("ra", rng.uniform(0, 6.28, 500)),
                         ("dec", rng.uniform(-1.5, 1.5, 500))):
                g.create_dataset(k, data=v)
            g.attrs["prior:luminosity_distance"] = "PowerLaw(alpha=2, minimum=10, maximum=10000)"
    hc.prepare_inputs(tmp_path / "work", "gwtc4", tmp_path / "mix.hdf", 0.25, 10.0, 3.0, hc.H0_DEFAULT_EXCLUDE,
                      cache, keep_pe_files=False, galaxy_catalog=cat)
    with h5py.File(tmp_path / "work" / "inputs.h5") as h:
        assert h.attrs["galaxy_catalog"] == str(cat.resolve())
        assert h.attrs["catalog_grouping"] == "K-glade+" and h.attrs["catalog_subgrouping"] == "eps_1"
        g = h["GW150914_095045"]
        assert len(g["ra"]) == len(g["dl"]) and len(g["dec"]) == len(g["dl"])
        gi = h["_injections"]
        assert len(gi["ra"]) == len(gi["prior"])
        assert abs(np.mean(np.sin(gi["dec"][:]))) < 0.1                 # isotropic
    sel = json.loads((tmp_path / "work" / "selection.json").read_text())
    assert sel["galaxy_catalog"]["band"] == "K-glade+"
    # an extract without the sky position does not count as cached when the sky is needed
    with h5py.File(cache / "samples" / "GW150914_095045.h5", "a") as h:
        del h["C00:IMRPhenomXPHM-SpinTaylor"]["ra"]
    assert not hc._extract_has_sky(cache / "samples" / "GW150914_095045.h5")


def test_galaxy_catalog_info_checks(tmp_path):
    with pytest.raises(ValueError, match="not found"):
        hc.galaxy_catalog_info(tmp_path / "nope.hdf5")
    with h5py.File(tmp_path / "raw.hdf5", "w"):
        pass
    with pytest.raises(ValueError, match="not a finished"):
        hc.galaxy_catalog_info(tmp_path / "raw.hdf5")
    info = hc.galaxy_catalog_info(_fake_catalog_file(tmp_path / "c.hdf5"))
    assert info["settings"]["nside"] == 64


def test_published_dark_siren_needs_the_same_catalog():
    pub = hc.H0_SENSITIVITY_RELEASES["gwtc4"]["published"]["dark_plp"]
    assert hc._matches_published_catalog(gc.catalog_settings("K-glade+"), pub)
    assert not hc._matches_published_catalog(gc.catalog_settings("K-glade+", epsilon=0.0), pub)
    assert not hc._matches_published_catalog(gc.catalog_settings("K-glade+", nside=128), pub)


# ---------------------------------------------------------------------------
# end to end in the icarogw environment
# ---------------------------------------------------------------------------

@pytest.mark.skipif(not ICAROGW_PYTHON, reason="set GWTC_ICAROGW_PYTHON to the icarogw interpreter")
def test_build_catalog_and_dark_likelihood(tmp_path):
    """A mock catalog through every stage (2 chunks), then the dark-siren likelihood of icarogw evaluated on it."""
    rng = np.random.default_rng(5)
    n = 1500
    w = gx.GalaxyWriter(tmp_path / "gal.h5", band="K-glade+", attrs=dict(source="mock"))
    z = rng.uniform(0.005, 0.2, n)
    w.append(pd.DataFrame({"ra": rng.uniform(0, 2 * np.pi, n), "dec": np.arcsin(rng.uniform(-1, 1, n)), "z": z,
                           "sigmaz": 0.01 * (1 + z), "m": rng.uniform(8, 14, n)}))
    w.close()
    out = gc.run_galaxy_catalog(workdir=tmp_path / "cat", galaxies=tmp_path / "gal.h5", nside=8, nside_mthr=4,
                                nshards=8, jobs=2, icarogw_python=ICAROGW_PYTHON,
                                out_report_html=tmp_path / "cat.html")
    assert out is not None and out.exists() and (tmp_path / "cat.html").exists()
    summ = json.loads((tmp_path / "cat" / "summary.json").read_text())
    assert summ["n_galaxies"] == n and summ["n_filled_pixels"] > 0
    # inputs.h5 for two events and a few injections, then the likelihood at one point
    info = hc.galaxy_catalog_info(out)
    with h5py.File(tmp_path / "inputs.h5", "w") as h:
        for e in range(2):
            g = h.create_group(f"GW00000{e}_000000")
            k = 400
            dl = rng.uniform(200, 800, k)
            for name, v in (("m1", rng.uniform(30, 40, k)), ("m2", rng.uniform(20, 30, k)), ("dl", dl),
                            ("prior", dl ** 2), ("ra", rng.normal(1.0, 0.05, k)), ("dec", rng.normal(0.3, 0.05, k))):
                g.create_dataset(name, data=v)
        gi = h.create_group("_injections")
        m = 3000
        dl = rng.uniform(100, 3000, m)
        for name, v in (("mass_1", rng.uniform(10, 80, m)), ("mass_2", rng.uniform(5, 40, m)),
                        ("luminosity_distance", dl), ("prior", np.full(m, 1e-6)),
                        ("ra", rng.uniform(0, 2 * np.pi, m)), ("dec", np.arcsin(rng.uniform(-1, 1, m)))):
            gi.create_dataset(name, data=v)
        gi.attrs.update(ntotal=1e5, Tobs=1.0)
        h.attrs.update(galaxy_catalog=str(info["path"]), catalog_grouping=info["grouping"],
                       catalog_subgrouping=info["subgrouping"])
    code = (
        "import sys, numpy as np\n"
        f"sys.path.insert(0, {str(Path(hc.__file__).parent)!r})\n"
        "import h0_icarogw as r\n"
        f"like, rate, cat, inj = r.build_likelihood(300, 1.0, 'plp', {str(tmp_path / 'inputs.h5')!r})\n"
        "p = dict(H0=70., Om0=0.3065, alpha=3.4, beta=1.1, mmin=5., mmax=87., delta_m=4.8, mu_g=34., sigma_g=3.6,\n"
        "         lambda_peak=0.04, gamma=2.7, kappa=3., zp=2.)\n"
        "like.parameters.update(p)\n"
        "print('LNL', like.log_likelihood(), type(rate).__name__)\n")
    env = dict(os.environ, HDF5_USE_FILE_LOCKING="FALSE")
    r = subprocess.run([ICAROGW_PYTHON, "-c", code], capture_output=True, text=True, cwd=tmp_path, env=env)
    line = [ln for ln in r.stdout.splitlines() if ln.startswith("LNL")]
    assert r.returncode == 0 and line, r.stderr[-2000:]
    _, lnl, name = line[0].split()
    assert name == "CBC_catalog_vanilla_rate" and (np.isfinite(float(lnl)) or float(lnl) == -np.inf)


def test_prepare_chunks_hold_whole_coarse_pixels():
    """The threshold of a pixel reads every pixel of its coarse (nside_mthr) pixel: the prepare chunks must hold
    whole coarse pixels, so that no chunk reads a file another chunk is writing (they raced before)."""
    healpy = pytest.importorskip("healpy")
    import importlib.util

    spec = importlib.util.spec_from_file_location("dark_catalog_icarogw", gc.RUNNER)
    runner = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(runner)
    pix = np.sort(np.random.default_rng(3).choice(healpy.nside2npix(64), 30000, replace=False))
    parts = [runner.prepare_chunk(pix, 64, 32, c, 4) for c in range(4)]
    assert np.array_equal(np.sort(np.concatenate(parts)), pix)
    coarse = [set(runner.coarse_pixels(p, 64, 32).tolist()) for p in parts]
    assert all(not (coarse[i] & coarse[j]) for i in range(4) for j in range(i + 1, 4))
    # the coarse pixel is that of the pixel centre, as icarogw finds it
    theta, phi = healpy.pix2ang(64, pix[:100])
    assert np.array_equal(runner.coarse_pixels(pix[:100], 64, 32), healpy.ang2pix(32, theta, phi))


@pytest.mark.skipif(not ICAROGW_PYTHON, reason="set GWTC_ICAROGW_PYTHON to the icarogw interpreter")
def test_assembly_from_summaries_matches_icarogw(tmp_path):
    """init and finish from the chunk summaries give exactly the catalog of icarogw's own functions, which read
    every pixel file in turn."""
    import shutil

    rng = np.random.default_rng(7)
    n = 2500
    w = gx.GalaxyWriter(tmp_path / "gal.h5", band="K-glade+", attrs=dict(source="mock"))
    z = rng.uniform(0.002, 0.3, n)
    w.append(pd.DataFrame({"ra": rng.uniform(0, 2 * np.pi, n), "dec": np.arcsin(rng.uniform(-1, 1, n)), "z": z,
                           "sigmaz": np.where(rng.random(n) < 0.3, 2e-4, 0.015), "m": rng.uniform(8, 14, n)}))
    w.close()
    kw = dict(galaxies=tmp_path / "gal.h5", nside=8, nside_mthr=4, nshards=8, jobs=2, nintegration="logspace:0.0001:600",
              icarogw_python=ICAROGW_PYTHON, out_report_html=None)
    gc.run_galaxy_catalog(stages=["shard", "pixels", "gather", "prepare"], workdir=tmp_path / "fast", **kw)
    shutil.copytree(tmp_path / "fast", tmp_path / "ref")
    for f in (tmp_path / "ref" / "done").glob("prepare_*.npz"):
        f.unlink()
    rest = ["init", "interpolate", "finish"]
    fast = gc.run_galaxy_catalog(stages=rest, workdir=tmp_path / "fast", **kw)
    gc.run_galaxy_catalog(stages=["init", "interpolate"], workdir=tmp_path / "ref", **kw)
    for f in (tmp_path / "ref" / "done").glob("interpolate_*.h5"):
        f.unlink()
    ref = gc.run_galaxy_catalog(stages=["finish"], workdir=tmp_path / "ref", **kw)

    def datasets(path):
        out = {}
        with h5py.File(path, "r") as h:
            h.visititems(lambda k, v: out.__setitem__(k, v[()]) if isinstance(v, h5py.Dataset) else None)
        return out

    a, b = datasets(fast), datasets(ref)
    assert sorted(a) == sorted(b) and len(b) == 7
    for k in b:
        assert np.allclose(a[k], b[k], rtol=1e-12, atol=0, equal_nan=True), k


# ---------------------------------------------------------------------------
# configurable selection and processing
# ---------------------------------------------------------------------------

def _raw_glade():
    return pd.DataFrame({"GLADE+": ["1", "2", "3"], "Type": ["G", "Q", "G"], "RAJ2000": ["10", "20", "30"],
                         "DEJ2000": ["0", "5", "10"], "Kmag": ["12", "15", "13"], "zhelio": ["0.049", "0.6", "0.2"],
                         "zcmb": ["0.05", "0.6", "0.2"], "f_zcmb": ["1", "0", "0"], "e_z": ["0.002", "", ""],
                         "e_zhelio": ["0.001", "0.01", "0.015"]})


def test_glade_selection_options():
    """Types, redshift column, redshift error and an extra cut on the VizieR columns."""
    raw = _raw_glade()
    assert len(gx.glade_kband_frame(raw, selection=gx.glade_selection())) == 2          # quasar dropped
    assert len(gx.glade_kband_frame(raw, selection=gx.glade_selection(types="G,Q"))) == 3
    df = gx.glade_kband_frame(raw, selection=gx.glade_selection(redshift="zhelio"))
    assert df["z"].tolist() == pytest.approx([0.049, 0.2])
    df = gx.glade_kband_frame(raw, selection=gx.glade_selection(sigmaz="measurement"))
    assert df["sigmaz"].tolist() == pytest.approx([0.001, 0.015])
    df = gx.glade_kband_frame(raw, selection=gx.glade_selection(sigmaz="peculiar"))
    assert df["sigmaz"].iloc[0] == pytest.approx(0.002) and np.isnan(df["sigmaz"].iloc[1])   # NaN: dropped later
    df = gx.glade_kband_frame(raw, selection=gx.glade_selection(sigmaz_const=0.01, sigmaz_relative=True))
    assert df["sigmaz"].tolist() == pytest.approx([0.0105, 0.012])
    assert len(gx.glade_kband_frame(raw, selection=gx.glade_selection(where="f_zcmb == 1"))) == 1
    with pytest.raises(ValueError):
        gx.glade_selection(types="G,X")
    with pytest.raises(ValueError):
        gx.glade_selection(sigmaz="other")
    assert gx.glade_selection_text(gx.glade_selection()) == "Type G, Kmag finite, zcmb > 0, sigmaz quadrature finite"


def test_fetch_glade_refetches_bands_missing_a_column(tmp_path):
    """A band cached without a column the selection needs (zhelio) is downloaded again; the selection is recorded."""
    cache = tmp_path / "cache"
    cache.mkdir()
    old = _raw_glade().drop(columns=["zhelio"])
    for lo, hi in gx.glade_dec_bands(90.0):
        old.to_csv(cache / f"dec_{lo:+07.2f}_{hi:+07.2f}.csv", index=False)
    calls = []

    def fake(lo, hi):
        calls.append(lo)
        return _raw_glade() if lo == -90 else _raw_glade().iloc[:0]

    gx.fetch_glade_kband(tmp_path / "a.h5", band_width=90.0, cache=cache, fetch=fake)
    assert calls == []                                    # zcmb: the cached bands suffice
    out = gx.fetch_glade_kband(tmp_path / "b.h5", band_width=90.0, cache=cache, fetch=fake,
                               selection=gx.glade_selection(redshift="zhelio"))
    assert calls == [-90.0, 0.0]
    cols, attrs = gx.read_galaxy_file(out)
    assert cols["z"].tolist() == pytest.approx([0.049, 0.2])
    assert json.loads(attrs["glade_selection"])["redshift"] == "zhelio"
    assert gx.selection_text(attrs).startswith("Type G, Kmag finite, zhelio > 0")


def test_existing_galaxy_file_with_another_selection_is_refused(tmp_path):
    def fake(lo, hi):
        return _raw_glade() if lo == -90 else _raw_glade().iloc[:0]

    gx.fetch_glade_kband(tmp_path / "galaxies.h5", band_width=90.0, cache=tmp_path / "c", fetch=fake)
    gc.run_galaxy_catalog(stages=["galaxies"], workdir=tmp_path, out_report_html=None)        # same selection: kept
    with pytest.raises(ValueError, match="another GLADE"):
        gc.run_galaxy_catalog(stages=["galaxies"], workdir=tmp_path, glade_types="G,Q", out_report_html=None)


def test_ptype_and_zmin_settings():
    assert gc.catalog_settings("K-glade+", ptype="gaussian_nocom")["ptype"] == "gaussian_nocom"
    with pytest.raises(ValueError, match="ptype"):
        gc.catalog_settings("K-glade+", ptype="lorentzian")
    assert gc.catalog_settings("K-glade+", zmin=0.05, zcut=0.35)["zmin"] == 0.05
    with pytest.raises(ValueError, match="zmin"):
        gc.catalog_settings("K-glade+", zmin=0.5, zcut=0.35)
    with pytest.raises(ValueError, match="logarithmic"):
        gc.catalog_settings("K-glade+", zmin=0.05, nintegration=10)


def test_runner_grid_starts_at_zmin():
    import importlib.util

    spec = importlib.util.spec_from_file_location("dark_catalog_icarogw", gc.RUNNER)
    runner = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(runner)
    s = gc.catalog_settings("K-glade+", zmin=0.05, zcut=0.35, nintegration="logspace:0.0001:100")
    g = runner._nintegration(s)
    assert g[0] == pytest.approx(0.05) and g[-1] == pytest.approx(0.35) and len(g) == 100
    g = runner._nintegration(gc.catalog_settings("K-glade+"))
    assert g[0] == pytest.approx(1e-4)


def test_published_match_uses_all_given_settings():
    """The comparison needs the settings the paper gives (threshold map, redshift range, selection), not those it
    leaves open (ptype); catalogs from older versions, without the new keys, still match."""
    pub = hc.H0_SENSITIVITY_RELEASES["gwtc4"]["published"]["dark_plp"]
    base = gc.catalog_settings("K-glade+")
    old = {k: v for k, v in base.items() if k not in ("zmin",)}
    assert hc._matches_published_catalog(old, pub)
    assert hc._matches_published_catalog(gc.catalog_settings("K-glade+", ptype="gaussian_nocom"), pub)
    assert not hc._matches_published_catalog(gc.catalog_settings("K-glade+", zmin=0.01), pub)
    assert not hc._matches_published_catalog(gc.catalog_settings("K-glade+", mthr_percentile=90), pub)
    assert hc._matches_published_catalog(
        dict(base, galaxy_selection="Type G, Kmag finite, zcmb > 0, sigmaz quadrature finite"), pub)
    assert not hc._matches_published_catalog(
        dict(base, galaxy_selection="Type G+Q, Kmag finite, zcmb > 0, sigmaz quadrature finite"), pub)


def test_settings_file(tmp_path):
    """--settings: the file's values are defaults, the command line overrides them, keys and values are checked."""
    from gwtc_analysis.cli import _apply_settings_file

    f = tmp_path / "s.json"
    f.write_text(json.dumps({"zcut": 0.35, "zmin": "0.05", "glade-types": "G,Q", "ptype": "gaussian_nocom",
                             "stages": ["galaxies", "shard"], "slurm_option": ["--mem=4G"]}))
    p = build_parser()
    a = _apply_settings_file(p, ["galaxy_catalog", "--settings", str(f), "--zcut", "0.3"])
    assert a.zcut == 0.3 and a.zmin == 0.05 and a.glade_types == "G,Q" and a.ptype == "gaussian_nocom"
    assert a.stages == ["galaxies", "shard"] and a.slurm_option == ["--mem=4G"]
    f.write_text(json.dumps({"zcutt": 0.3}))
    with pytest.raises(ValueError, match="unknown option"):
        _apply_settings_file(build_parser(), ["galaxy_catalog", "--settings", str(f)])
    f.write_text(json.dumps({"ptype": "lorentzian"}))
    with pytest.raises(ValueError, match="not among"):
        _apply_settings_file(build_parser(), ["galaxy_catalog", "--settings", str(f)])
    f.write_text(json.dumps({"neff-pe": 20, "pe-samples": 3000}))
    a = _apply_settings_file(build_parser(), ["hubble_constant", "--settings", str(f)])
    assert a.neff_pe == 20.0 and a.pe_samples == 3000


def test_likelihood_settings(tmp_path):
    """--neff-pe and --neff-inj are recorded in likelihood.json, read by the runner, and cannot change once runs exist."""
    import importlib.util

    spec = importlib.util.spec_from_file_location("h0_icarogw", Path(hc.__file__).with_name("h0_icarogw.py"))
    runner = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(runner)
    assert runner.likelihood_settings(tmp_path) == {"neff_pe": 10, "neff_inj": None}
    assert hc.set_likelihood_settings(tmp_path, neff_pe=20) == {"neff_pe": 20.0, "neff_inj": None}
    assert hc.set_likelihood_settings(tmp_path, neff_inj=600) == {"neff_pe": 20.0, "neff_inj": 600}
    assert runner.likelihood_settings(tmp_path) == {"neff_pe": 20.0, "neff_inj": 600}
    (tmp_path / "result").mkdir()
    (tmp_path / "result" / "plp_seed1_result.json").write_text("{}")
    assert hc.set_likelihood_settings(tmp_path) == {"neff_pe": 20.0, "neff_inj": 600}        # unchanged: fine
    with pytest.raises(ValueError, match="another work directory"):
        hc.set_likelihood_settings(tmp_path, neff_pe=50)


@pytest.mark.skipif(not ICAROGW_PYTHON, reason="set GWTC_ICAROGW_PYTHON to the icarogw interpreter")
def test_catalog_with_zmin_and_gaussian_nocom(tmp_path):
    """zmin: the grid starts there and below it icarogw sees no catalog (completeness correction alone);
    the Gaussian-as-posterior redshift probability builds too."""
    rng = np.random.default_rng(8)
    n = 1200
    w = gx.GalaxyWriter(tmp_path / "gal.h5", band="K-glade+", attrs=dict(source="mock"))
    z = rng.uniform(0.005, 0.2, n)
    w.append(pd.DataFrame({"ra": rng.uniform(0, 2 * np.pi, n), "dec": np.arcsin(rng.uniform(-1, 1, n)), "z": z,
                           "sigmaz": np.full(n, 0.005), "m": rng.uniform(8, 13, n)}))
    w.close()
    out = gc.run_galaxy_catalog(workdir=tmp_path / "cat", galaxies=tmp_path / "gal.h5", nside=8, nside_mthr=4,
                                nshards=8, jobs=2, zmin=0.05, zcut=0.3, nintegration="logspace:0.001:400",
                                ptype="gaussian_nocom", icarogw_python=ICAROGW_PYTHON, out_report_html=None)
    st = json.loads((tmp_path / "cat" / "catalog_settings.json").read_text())
    assert st["zmin"] == 0.05 and st["ptype"] == "gaussian_nocom"
    code = (
        "import sys, numpy as np, json\n"
        f"sys.path.insert(0, {str(Path(hc.__file__).parent)!r})\n"
        "import h0_icarogw as r, icarogw\n"
        f"cat = r.load_galaxy_catalog({str(out)!r}, 'K-glade+', 'eps_1')\n"
        "cw = icarogw.wrappers.FlatLambdaCDM_wrap(zmax=20.0); cw.update(H0=70., Om0=0.3065)\n"
        "cat.sch_fun.build_MF(cw.cosmology)\n"
        "zz = np.array([0.02, 0.1]); pix = np.zeros(2, dtype=int)\n"
        "gc_, bg = cat.effective_galaxy_number_interpolant(zz, pix, cw.cosmology, average=False)\n"
        "full = cat.sch_fun.background_effective_galaxy_density(-np.inf * np.ones(2), zz) * "
        "cw.cosmology.dVc_by_dzdOmega_at_z(zz)\n"
        "print('OUT', json.dumps(dict(z0=float(cat.z_grid[0]), gc=gc_.tolist(), ratio=(bg / full).tolist())))\n")
    env = dict(os.environ, HDF5_USE_FILE_LOCKING="FALSE")
    r = subprocess.run([ICAROGW_PYTHON, "-c", code], capture_output=True, text=True, cwd=tmp_path, env=env)
    line = [ln for ln in r.stdout.splitlines() if ln.startswith("OUT")]
    assert r.returncode == 0 and line, r.stderr[-2000:]
    d = json.loads(line[0][4:])
    assert d["z0"] == pytest.approx(0.05)
    assert d["gc"][0] == 0.0 and d["ratio"][0] == pytest.approx(1.0)     # below zmin: no catalog, full correction
