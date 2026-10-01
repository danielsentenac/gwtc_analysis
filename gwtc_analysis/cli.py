from __future__ import annotations

import argparse
from typing import List, Optional

from .catalogs import RATES_DEFAULT_RELEASE, RATES_SENSITIVITY_RELEASES, run_catalog_statistics, run_merger_rates
from .hubble_constant import (H0_DEFAULT_EXCLUDE, H0_DEFAULT_MASS_MODEL, H0_DEFAULT_RELEASE, H0_SENSITIVITY_RELEASES,
                              MASS_MODELS as H0_MASS_MODELS, STAGES as H0_STAGES, run_hubble_constant)
from .bright_siren import COUNTERPARTS as BRIGHT_SIREN_COUNTERPARTS, run_bright_siren
from .event_selection import run_event_selection
from .source_classes import (MASS_GAP as SC_MASS_GAP, NS_MAX_MASS as SC_NS_MAX, PISN_GAP_MIN as SC_PISN_MIN,
                             PRESETS as SOURCE_PRESETS)
from .search_skymaps import run_search_skymaps
from .parameters_estimation import run_parameters_estimation
from .unofficial_pe import build_unofficial_pe_bundle, get_unofficial_pe_spec, list_unofficial_pe_specs
from .gw_stat import ALLOWED_CATALOGS as ALLOWED_CATALOGS
from .catalog_registry import DEFAULT_H0_RELEASE, catalog_help, catalog_runs_help, release_runs_help, update_help

DEFAULT_H0_RELEASE_TEXT = f"{DEFAULT_H0_RELEASE}: the published analysis"
from .data_repo import parse_zenodo_version, zenodo_catalogs
import sys


def _format_allowed_catalogs() -> str:
    return ", ".join(ALLOWED_CATALOGS)


def _validate_catalogs(catalogs: list[str]) -> None:
    bad = [c for c in catalogs if c != "ALL" and c not in ALLOWED_CATALOGS]
    if bad:
        raise ValueError(
            "Unknown catalog(s): "
            + ", ".join(bad)
            + ". Allowed catalogs are: "
            + _format_allowed_catalogs()
        )


def _split_csv(s: str) -> List[str]:
    return [x.strip() for x in s.split(",") if x.strip()]


def _parse_catalogs(items: Optional[List[str]]) -> List[str]:
    """Parse catalogs passed as space-separated items, each item optionally comma-separated."""
    if not items:
        return []
    out: List[str] = []
    for it in items:
        out.extend(_split_csv(it))
    return out


def _parse_zenodo_versions(items: Optional[List[str]], data_repo: str) -> Optional[dict[str, str]]:
    """Parse --zenodo-version CAT=VER items into {catalog: version}."""
    if not items:
        return None
    if data_repo != "zenodo":
        raise ValueError("--zenodo-version only applies with --data-repo zenodo")
    out: dict[str, str] = {}
    for it in _parse_catalogs(items):
        cat, sep, ver = it.partition("=")
        cat, ver = cat.strip(), ver.strip()
        if not sep or not cat or not ver:
            raise ValueError(f"Invalid --zenodo-version {it!r}: expected CATALOG=VERSION, e.g. GWTC-3=v2")
        if cat not in zenodo_catalogs():
            raise ValueError(
                f"No Zenodo release for catalog {cat!r} in --zenodo-version. "
                f"Catalogs with Zenodo releases: {', '.join(zenodo_catalogs())}"
            )
        parse_zenodo_version(ver)  # validate the format early
        out[cat] = ver
    return out


def _add_zenodo_version_arg(p: argparse.ArgumentParser) -> None:
    p.add_argument(
        "--zenodo-version",
        nargs="+",
        default=None,
        metavar="CATALOG=VERSION",
        help=(
            "With --data-repo zenodo, read an older Zenodo release version of a catalog instead of the latest "
            "(e.g. --zenodo-version GWTC-3=v2 GWTC-4=v1). Versions are numbered from the oldest (v1); "
            "list them with the zenodo_releases mode."
        ),
    )


def _fraction_or_auto(x: str):
    """'auto' or a fraction in (0, 1]."""
    if str(x).lower() == "auto":
        return "auto"
    try:
        f = float(x)
    except ValueError:
        raise argparse.ArgumentTypeError(f"expected 'auto' or a number in (0, 1], got {x!r}") from None
    if not 0 < f <= 1:
        raise argparse.ArgumentTypeError(f"expected a number in (0, 1], got {f}")
    return f


def _none_if_empty(x):
    """Argparse with nargs can yield [] instead of None."""
    if x is None:
        return None
    if isinstance(x, list) and len(x) == 0:
        return None
    return x


def build_parser() -> argparse.ArgumentParser:
    fmt = argparse.ArgumentDefaultsHelpFormatter

    p = argparse.ArgumentParser(
        prog="gwtc_analysis",
        description=(
            "GWTC analysis tool.\n\n"
            "Use one of the MODE subcommands below. Each mode has its own detailed help:\n"
            "  gwtc_analysis catalog_statistics -h\n"
            "  gwtc_analysis rates -h\n"
            "  gwtc_analysis hubble_constant -h\n"
            "  gwtc_analysis bright_siren -h\n"
            "  gwtc_analysis event_selection -h\n"
            "  gwtc_analysis search_skymaps -h\n"
            "  gwtc_analysis parameters_estimation -h\n"
            "  gwtc_analysis build_unofficial_pe -h\n"
            "  gwtc_analysis zenodo_releases -h\n"
            "  gwtc_analysis check_catalogs -h\n"
        ),
        formatter_class=argparse.RawTextHelpFormatter,
    )

    from . import __version__
    from . import catalog_registry as _reg_v
    p.add_argument("--version", action="version",
                   version=f"gwtc_analysis {__version__}\n" + _reg_v.coverage_text(__version__).replace("`", ""),
                   help="Show the version and the catalogs it covers, then exit.")

    sub = p.add_subparsers(dest="mode", required=True, metavar="MODE")

    # ---------------------------------------------------------------------
    # catalog_statistics
    # ---------------------------------------------------------------------
    p_cat = sub.add_parser(
        "catalog_statistics",
        help="Build per-event TSV table and HTML summary report (with PIE plots) from GW catalogs.",
        description=(
            "Fetch events statistics for one or more catalogs.\n\n"
            "Outputs:\n"
            "  --out-events : TSV table of events\n"
            "  --out-report : HTML report (tables + plots)\n\n"
            "Optional additions:\n"
            "  --include-detectors : detector network via GWOSC calls\n"
            "  --include-area      : sky localization area Axx\n"
        ),
        formatter_class=argparse.RawTextHelpFormatter,
    )
    p_cat.add_argument(
        "--catalogs",
        required=True,
        nargs="+",
        help=f"Catalog keys, space-separated (e.g. {catalog_help()}). ALL takes them all except the updates "
             f"({update_help()}), which are used only when named.",
    )
    p_cat.add_argument("--out-events", default="catalogs_statistics.tsv", help="Output TSV path (per-event table).")
    p_cat.add_argument("--out-report", default="catalogs_statistics.html", help="Output HTML report path.")
    p_cat.add_argument("--include-detectors", action="store_true", help="Include detector network via GWOSC v2 calls.")
    p_cat.add_argument("--include-area", action="store_true", help="Compute sky localization area Axx if skymaps are available.")
    p_cat.add_argument("--area-cred", type=float, default=0.9, help="Credible level for sky area: 0.9→A90, 0.5→A50, 0.95→A95.")
    p_cat.add_argument("--plots-dir", default="cat_plots", help="Directory for plots (default: cat_plots).")
    p_cat.add_argument("--data-repo", choices=["galaxy", "zenodo", "s3"], default="zenodo", help="Where to read data from: galaxy | zenodo | s3.")
    _add_zenodo_version_arg(p_cat)

    # ---------------------------------------------------------------------
    # rates
    # ---------------------------------------------------------------------
    p_rate = sub.add_parser(
        "rates",
        help="Estimate BNS / NSBH / BBH merger rates (per Gpc^3 per year) from the catalogs.",
        description=(
            "Merger rates R = N / <VT> per population.\n\n"
            "N: GWOSC candidates (confident and marginal lists) inside the time span of the LVK\n"
            "search-sensitivity injections with FAR below --far-threshold, classified by median\n"
            "source-frame masses. <VT>: sensitive volume-time from the injections, reweighted to\n"
            "fixed population models. The injection file of --sensitivity-release is retrieved\n"
            "automatically from Zenodo (latest version of the record) and cached.\n\n"
            "Outputs:\n"
            "  --out-rates  : TSV of rates per population\n"
            "  --out-events : TSV of the events counted\n"
            "  --out-report : HTML report with observed and selection-corrected mass distributions\n"
        ),
        formatter_class=argparse.RawTextHelpFormatter,
    )
    p_rate.add_argument("--out-rates", default="merger_rates.tsv", help="Output TSV of rates per population.")
    p_rate.add_argument("--out-events", default="merger_rates_events.tsv", help="Output TSV of the events counted.")
    p_rate.add_argument("--out-report", default="merger_rates.html", help="Output HTML report path.")
    p_rate.add_argument("--plots-dir", default="rates_plots", help="Directory for plots (default: rates_plots).")
    p_rate.add_argument(
        "--sensitivity-release",
        choices=list(RATES_SENSITIVITY_RELEASES),
        default=RATES_DEFAULT_RELEASE,
        help="LVK search-sensitivity release retrieved automatically from Zenodo: "
        + "; ".join(f"{k} = {v[1]}" for k, v in RATES_SENSITIVITY_RELEASES.items()) + ".",
    )
    p_rate.add_argument("--sensitivity-file", default=None, help="Local LVK injection HDF file to use instead of --sensitivity-release.")
    p_rate.add_argument("--far-threshold", type=float, default=1.0, help="FAR threshold [1/yr] for both injections and events.")
    p_rate.add_argument("--ns-max-mass", type=float, default=2.5, help="Maximum neutron-star mass [Msun] separating NS from BH.")
    p_rate.add_argument("--bbh-kappa", type=float, default=2.9, help="BBH rate evolution R ∝ (1+z)^kappa.")
    p_rate.add_argument("--bbh-z-ref", type=float, default=0.2, help="Redshift at which the evolving BBH rate is reported.")
    p_rate.add_argument("--catalogs", nargs="+", default=None,
                        help=f"Catalog keys ({catalog_help()}, or ALL): events and injections are restricted to "
                             f"their observing runs ({catalog_runs_help()}). Default: the runs of the real-injection "
                             "mixture (O3 onward).")
    p_rate.add_argument("--snr-threshold", type=float, default=10.0,
                        help="Network SNR threshold for the semi-analytic O1+O2 injections (with GWTC-1).")

    # ---------------------------------------------------------------------
    # hubble_constant
    # ---------------------------------------------------------------------
    p_h0 = sub.add_parser(
        "hubble_constant",
        help="Estimate the Hubble constant from the BBH mass spectrum (spectral siren, icarogw).",
        description=(
            "Spectral-siren H0: the BBH mass distribution (--mass-model: Power Law + Peak or Multi Peak) and\n"
            "the Madau-Dickinson rate evolution fitted together with H0 (flat LCDM, Om0 = 0.3065), with icarogw\n"
            "and bilby/dynesty. The default setup reproduces the GWTC-4.0 cosmology paper (arXiv:2509.04348, v3):\n"
            "H0 = 105.5 (+46.4 / -35.8) km/s/Mpc (plp), 72.3 (+42.5 / -25.6) km/s/Mpc (mltp).\n\n"
            "Stages (--stages, default all, in this order):\n"
            "  prepare : select the events, download their PE samples from Zenodo (tens of GB, cached in\n"
            "            --pe-cache; only the extracted samples are kept unless --keep-pe-files), prepare\n"
            "            the found injections of --sensitivity-release -> <workdir>/inputs.h5\n"
            "  sample  : one dynesty run per --seeds value (hours each; resumable), --parallel at a time;\n"
            "            with --inj-fraction auto a probe first chooses the fastest reliable injection subset\n"
            "  combine : merge the runs, posterior + corner plot + stability diagnostics\n"
            "  reweight: when the runs used a subset of the injections, reweight their posterior to all of them\n"
            "  report  : HTML report and TSV summary\n\n"
            "icarogw needs its own environment: pass its interpreter with --icarogw-python. The sampler\n"
            "gwtc_analysis/h0_icarogw.py is standalone, so runs can also be started by hand on other\n"
            "machines sharing <workdir>:  python h0_icarogw.py run --workdir DIR --seed N\n"
        ),
        formatter_class=argparse.RawTextHelpFormatter,
    )
    p_h0.add_argument("--stages", nargs="+", choices=list(H0_STAGES), default=list(H0_STAGES),
                      help="Stages to run (default: all).")
    p_h0.add_argument("--workdir", default="hubble_constant_run", help="Work directory (inputs, runs, posterior).")
    p_h0.add_argument("--out-report", default="hubble_constant.html", help="Output HTML report path.")
    p_h0.add_argument("--out-summary", default="hubble_constant.tsv", help="Output TSV of the posterior quantiles.")
    p_h0.add_argument(
        "--sensitivity-release",
        choices=list(H0_SENSITIVITY_RELEASES),
        default=H0_DEFAULT_RELEASE,
        help="LVK search-sensitivity release (and matching catalogs and runs): "
        + "; ".join(f"{k} = {v['label']}" for k, v in H0_SENSITIVITY_RELEASES.items()) + ".",
    )
    p_h0.add_argument("--catalogs", nargs="+", default=None,
                      help=f"Catalog keys ({catalog_help()}, or ALL): events and injections are restricted to "
                           "their observing runs. Default: all the runs of --sensitivity-release "
                           f"({release_runs_help()}; {DEFAULT_H0_RELEASE_TEXT}).")
    p_h0.add_argument("--sensitivity-file", default=None,
                      help="Local LVK injection mixture file (semi-analytic O1+O2 + real) instead of the release's.")
    p_h0.add_argument("--far-threshold", type=float, default=0.25,
                      help="FAR threshold [1/yr] for the events and the real injections.")
    p_h0.add_argument("--snr-threshold", type=float, default=10.0,
                      help="Network SNR threshold for the semi-analytic O1+O2 injections.")
    p_h0.add_argument("--min-mass", type=float, default=3.0,
                      help="Minimum source-frame mass [Msun] of both components (potential neutron stars excluded).")
    p_h0.add_argument("--exclude", nargs="*", default=list(H0_DEFAULT_EXCLUDE), help="Events left out.")
    p_h0.add_argument("--pe-cache", default=None,
                      help="PE cache directory (files/, samples/, index/); default ~/.cache_gwtc_analysis/pe_catalog "
                           "or $GWTC_PE_CACHE.")
    p_h0.add_argument("--keep-pe-files", action="store_true", help="Keep the full PE files after extraction.")
    p_h0.add_argument("--mass-model", choices=list(H0_MASS_MODELS), default=H0_DEFAULT_MASS_MODEL,
                      help="BBH primary-mass model: " + "; ".join(f"{k} = {v}" for k, v in H0_MASS_MODELS.items())
                      + ". Use one work directory per model.")
    p_h0.add_argument("--seeds", nargs="+", type=int, default=[1], help="One sampler run per seed.")
    p_h0.add_argument("--parallel", type=int, default=1,
                      help="Seeds run at the same time on this machine (each with --npool processes; "
                           "logs in <workdir>/logs).")
    p_h0.add_argument("--nlive", type=int, default=100, help="dynesty live points per run.")
    p_h0.add_argument("--npool", type=int, default=4,
                      help="Worker processes per run: random walks of one seed run at the same time.")
    p_h0.add_argument("--naccept", type=int, default=60, help="dynesty accepted steps per MCMC walk.")
    p_h0.add_argument("--pe-samples", type=int, default=1500, help="PE samples per event.")
    p_h0.add_argument("--inj-fraction", type=_fraction_or_auto, default="auto",
                      help="Fraction of the found injections used by the sampler runs: 'auto' (a probe chooses the "
                           "fastest reliable subset, the posterior being then reweighted to all the injections), or a "
                           "number in (0, 1], 1 = all the injections, as in the paper.")
    p_h0.add_argument("--min-ess-fraction", type=float, default=0.5,
                      help="With --inj-fraction auto: smallest predicted effective-sample-size fraction accepted for "
                           "the reweighting to all the injections.")
    p_h0.add_argument("--probe-points", type=int, default=30,
                      help="With --inj-fraction auto: finite-likelihood prior points used by the probe.")
    p_h0.add_argument("--reweight-pe-samples", type=int, default=None,
                      help="PE samples per event of the reweighting target (default: those of the runs).")
    p_h0.add_argument("--icarogw-python", default=None,
                      help="Python interpreter of the icarogw environment (default: the current one).")

    # ---------------------------------------------------------------------
    # bright_siren
    # ---------------------------------------------------------------------
    p_bs = sub.add_parser(
        "bright_siren",
        help="Estimate the Hubble constant from an event with an identified host galaxy (bright siren).",
        description=(
            "Bright-siren H0: the luminosity distance of the GW signal, at the sky position of the\n"
            "counterpart, against the Hubble-flow velocity of the host galaxy (v_H = v_r - <v_p>).\n"
            "Flat H0 prior; for sources uniform in volume and a GW-limited detection, the selection term\n"
            "cancels the volume prior on the distance. The default setup follows LVK 2017\n"
            "(arXiv:1710.05835): H0 = 70.0 (+12.0 / -8.0) km/s/Mpc (maximum a posteriori, 68%).\n\n"
            "GW170817 uses the bundle built from public GWTC-1 products (build_unofficial_pe).\n"
            "--spectral-posterior combines the result with a spectral-siren posterior (hubble_constant).\n"
        ),
        formatter_class=argparse.RawTextHelpFormatter,
    )
    p_bs.add_argument("--src-name", default="GW170817", choices=list(BRIGHT_SIREN_COUNTERPARTS),
                      help="Event with an identified host galaxy.")
    p_bs.add_argument("--pe-label", nargs="+", default=None,
                      help="PE label(s) to use (default: all the labels of the PE file, LowSpin first).")
    p_bs.add_argument("--pe-file", default=None, help="PE file to read instead of the event's bundle.")
    p_bs.add_argument("--cache-dir", default=".cache_gwosc",
                      help="Cache root of the unofficial PE bundle (as in build_unofficial_pe).")
    p_bs.add_argument("--v-recession", nargs=2, type=float, metavar=("V", "SIGMA"), default=None,
                      help="Recession velocity of the host and its uncertainty, km/s (default for GW170817: "
                           "3327 72, the NGC 4993 group in the CMB frame).")
    p_bs.add_argument("--v-peculiar", nargs=2, type=float, metavar=("V", "SIGMA"), default=None,
                      help="Peculiar velocity of the host and its uncertainty, km/s (default for GW170817: 310 150).")
    p_bs.add_argument("--redshift", nargs=2, type=float, metavar=("Z", "SIGMA"), default=None,
                      help="Hubble-flow redshift of the host and its uncertainty, instead of the velocities "
                           "(default for GW190521: 0.438 0.0015).")
    p_bs.add_argument("--selection", choices=("auto", "euclidean", "injections"), default="auto",
                      help="Selection term: euclidean (GW-limited, nearby sources: beta ∝ H0^3), injections (LVK "
                           "sensitivity injections of the event's run), auto (euclidean below z = 0.05).")
    p_bs.add_argument("--sensitivity-release", choices=list(H0_SENSITIVITY_RELEASES), default=None,
                      help=f"Injections of the selection term (default: {DEFAULT_H0_RELEASE}).")
    p_bs.add_argument("--sensitivity-file", default=None, help="Local sensitivity file instead of the release.")
    p_bs.add_argument("--far-threshold", type=float, default=0.25, help="Found injections: FAR below this, per year.")
    p_bs.add_argument("--snr-threshold", type=float, default=10.0,
                      help="Found semi-analytic O1+O2 injections: network SNR above this.")
    p_bs.add_argument("--pe-cache", default=None,
                      help="PE cache of the events read from Zenodo (default: that of hubble_constant).")
    p_bs.add_argument("--sky-radius", type=float, default=3.0,
                      help="For samples not fixed to the counterpart's position: keep those within this angle (deg).")
    p_bs.add_argument("--spectral-posterior", default=None,
                      help="Spectral-siren H0 posterior to combine with: a hubble_constant work directory or a "
                           "posterior TSV with an H0 column.")
    p_bs.add_argument("--h0-range", nargs=2, type=float, metavar=("MIN", "MAX"), default=[10.0, 200.0],
                      help="Flat H0 prior range, km/s/Mpc (that of the spectral siren by default).")
    p_bs.add_argument("--out-report", default="bright_siren.html", help="Output HTML report path.")
    p_bs.add_argument("--out-summary", default="bright_siren.tsv",
                      help="Output TSV of the H0 summary (the posterior grid goes to <name>.posterior.tsv).")
    p_bs.add_argument("--plots-dir", default="bright_siren_plots", help="Directory for the plots.")

    # ---------------------------------------------------------------------
    # event_selection
    # ---------------------------------------------------------------------
    p_sel = sub.add_parser(
        "event_selection",
        help="Select GW events based on physical criteria (mass, distance, spin, source class) and write a TSV.",
        description=(
            "Select events by simple cuts on source-frame masses, luminosity distance and chi_eff (GWOSC medians),\n"
            "and/or a --preset class of sources. Cuts are optional; if a cut is not provided, it is not applied.\n"
        ),
        formatter_class=argparse.RawTextHelpFormatter,
    )
    p_sel.add_argument("--catalogs", required=True, nargs="+", help=f"Catalog keys, space-separated (e.g. {catalog_help()}). ALL takes them all except the updates "
             f"({update_help()}), which are used only when named.")
    p_sel.add_argument("--out-selection", default="event_selection.tsv", help="Output TSV path for the selected events.")
    p_sel.add_argument("--m1-min", type=float, default=None, help="Minimum primary mass (source frame).")
    p_sel.add_argument("--m1-max", type=float, default=None, help="Maximum primary mass (source frame).")
    p_sel.add_argument("--m2-min", type=float, default=None, help="Minimum secondary mass (source frame).")
    p_sel.add_argument("--m2-max", type=float, default=None, help="Maximum secondary mass (source frame).")
    p_sel.add_argument("--dl-min", type=float, default=None, help="Minimum luminosity distance (Mpc).")
    p_sel.add_argument("--dl-max", type=float, default=None, help="Maximum luminosity distance (Mpc).")
    p_sel.add_argument("--chi-eff-min", type=float, default=None, help="Minimum effective spin chi_eff.")
    p_sel.add_argument("--chi-eff-max", type=float, default=None, help="Maximum effective spin chi_eff.")
    p_sel.add_argument("--preset", choices=list(SOURCE_PRESETS), default=None,
                       help="Class of sources (the cuts apply on top): neutron-stars (a component below "
                            "--ns-max-mass), mass-gap (a component in --mass-gap), hierarchical (primary above "
                            "--pisn-gap-min, or chi_eff < 0 at 90%%: earlier-generation black holes).")
    p_sel.add_argument("--ns-max-mass", type=float, default=SC_NS_MAX, help="Maximum neutron-star mass (M_sun).")
    p_sel.add_argument("--mass-gap", nargs=2, type=float, metavar=("LO", "HI"), default=list(SC_MASS_GAP),
                       help="Lower mass gap between neutron stars and black holes (M_sun).")
    p_sel.add_argument("--pisn-gap-min", type=float, default=SC_PISN_MIN,
                       help="Lower edge of the pair-instability mass gap (M_sun; ~45-65 in the literature).")
    p_sel.add_argument("--out-plot", default=None,
                       help="Optional PNG of the selected events among all the events of the catalogs (m2 and D_L "
                            "against m1).")

    # ---------------------------------------------------------------------
    # search_skymaps
    # ---------------------------------------------------------------------
    p_sky = sub.add_parser(
        "search_skymaps",
        help="Search GW sky localizations for a given sky position (RA/Dec).",
        description=(
            "Given a sky position (RA/Dec in degrees) and catalogs,\n"
            "report which events contain that position above the requested credible level.\n"
            "Plotting: Hit skymaps are produced.\n"
        ),
        formatter_class=argparse.RawTextHelpFormatter,
    )
    p_sky.add_argument("--catalogs", required=True, nargs="+", help=f"Catalog keys, space-separated (e.g. {catalog_help()}). ALL takes them all except the updates "
             f"({update_help()}), which are used only when named.")
    p_sky.add_argument("--ra-deg", type=float, required=True, help="Right ascension (deg).")
    p_sky.add_argument("--dec-deg", type=float, required=True, help="Declination (deg).")
    p_sky.add_argument("--prob", type=float, default=0.9, help="Credible-level threshold (0–1). Common values: 0.9, 0.5, 0.95.")
    p_sky.add_argument("--skymap-label", default="Mixed", help="Label selector used to filter skymap (default: Mixed).")
    p_sky.add_argument("--out-events", default="search_skymaps.tsv", help="Output TSV file (default: search_skymaps.tsv).")
    p_sky.add_argument("--out-report", default="search_skymaps.html", help="Optional output HTML report path for hits.")
    p_sky.add_argument("--plots-dir", default="sky_plots", help="Directory for hit plots (default: sky_plots).")
    p_sky.add_argument("--data-repo", choices=["galaxy", "zenodo", "s3"], default="zenodo", help="Where to read data from: galaxy | zenodo | s3.")
    _add_zenodo_version_arg(p_sky)

    # ---------------------------------------------------------------------
    # parameters_estimation
    # ---------------------------------------------------------------------
    p_pe = sub.add_parser(
        "parameters_estimation",
        help="Generate parameter-estimation plots (posteriors, strain, waveforms) for one event.",
        description=(
            "Generate PE plots (posteriors, skymap, strain overlays, PSD) for a single event.\n"
        ),
        formatter_class=argparse.RawTextHelpFormatter,
    )
    p_pe.add_argument("--out-report", default="parameters_estimation.html", help="Output HTML report path.")
    p_pe.add_argument("--src-name", dest="src_name", required=True, help="Source event name (e.g. GW231223_032836).")
    p_pe.add_argument("--data-repo", choices=["galaxy", "zenodo", "s3"], default="zenodo", help="Where to read data from: galaxy | zenodo | s3.")
    _add_zenodo_version_arg(p_pe)
    p_pe.add_argument("--pe-vars", nargs="+", default=None, help=("Extra posterior sample variables to plot (space-separated). Example: --pe-vars chi_eff chi_p luminosity_distance."))
    p_pe.add_argument("--pe-pairs", nargs="+", default=None, help=("Extra 2D posterior pairs to plot as 'x:y' tokens. Example: --pe-pairs mass_1_source:mass_2_source chi_eff:chi_p."))
    p_pe.add_argument("--plots-dir", default="pe_plots", help="Directory for output PE plots (default: pe_plots).")
    p_pe.add_argument("--start", type=float, default=0.2, help="Default seconds before GPS time for overlay and q-transform windows.")
    p_pe.add_argument("--stop", type=float, default=0.1, help="Default seconds after GPS time for overlay and q-transform windows.")
    p_pe.add_argument("--fmin", type=float, default=20.0, help="Default low frequency bound (Hz) used for overlay filtering and q-transform range.")
    p_pe.add_argument("--fmax", type=float, default=300.0, help="Default high frequency bound (Hz) used for overlay filtering and q-transform range.")
    p_pe.add_argument("--fs-low", dest="fmin", type=float, help=argparse.SUPPRESS)
    p_pe.add_argument("--fs-high", dest="fmax", type=float, help=argparse.SUPPRESS)
    p_pe.add_argument("--overlay-start", type=float, default=None, help="Override seconds before GPS time for the whitened overlay window.")
    p_pe.add_argument("--overlay-stop", type=float, default=None, help="Override seconds after GPS time for the whitened overlay window.")
    p_pe.add_argument("--overlay-fmin", type=float, default=None, help="Override low frequency bound (Hz) for overlay whitening/bandpass.")
    p_pe.add_argument("--overlay-fmax", type=float, default=None, help="Override high frequency bound (Hz) for overlay whitening/bandpass.")
    p_pe.add_argument("--q-start", type=float, default=None, help="Override seconds before GPS time for the q-transform window.")
    p_pe.add_argument("--q-stop", type=float, default=None, help="Override seconds after GPS time for the q-transform window.")
    p_pe.add_argument("--q-fmin", type=float, default=None, help="Override low frequency bound (Hz) for the q-transform.")
    p_pe.add_argument("--q-fmax", type=float, default=None, help="Override high frequency bound (Hz) for the q-transform.")
    p_pe.add_argument("--q-fscale", choices=["linear", "log"], default="log", help="Frequency axis scaling for q-transform plots (default: log).")

    # Renamed options (no legacy names)
    p_pe.add_argument(
        "--pe-label",
        default=None,
        help=(
            "PE label used to select posterior samples and metadata. "
            "If omitted and --waveform-engine is provided, the tool selects the closest PE label "
            "by substring match in the PE label. If both are omitted, defaults to Mixed."
        ),
    )
    p_pe.add_argument(
        "--waveform-engine",
        default=None,
        help="Waveform engine used to generate a time-domain waveform for strain overlay. If omitted, a sensible default engine is used for overlays.",
    )

    # ---------------------------------------------------------------------
    # build_unofficial_pe
    # ---------------------------------------------------------------------
    supported_unofficial = ", ".join(list_unofficial_pe_specs()) or "(none)"
    p_unoff = sub.add_parser(
        "build_unofficial_pe",
        help="Build an unofficial PESummary-compatible PE bundle from locally cached source files.",
        description=(
            "Build an unofficial PEDataRelease-style HDF5 bundle for a supported special-case event.\n\n"
            f"Supported events: {supported_unofficial}\n"
        ),
        formatter_class=argparse.RawTextHelpFormatter,
    )
    p_unoff.add_argument("--src-name", dest="src_name", required=True, help="Source event name (e.g. GW170817).")
    p_unoff.add_argument("--cache-dir", default=".cache_gwosc", help="Cache root where unofficial_pe/<bundle>.h5 will be written.")
    p_unoff.add_argument("--force", action="store_true", help="Force rebuilding the unofficial bundle even if a cached copy already exists and is up to date.")

    # ---------------------------------------------------------------------
    # zenodo_releases
    # ---------------------------------------------------------------------
    p_chk = sub.add_parser(
        "check_catalogs",
        help="Compare the catalogs published by GWOSC and Zenodo with the catalog registry.",
        description=(
            "Reports the GWOSC event lists and observing runs that gwtc_analysis/catalog_registry.py does not\n"
            "describe, and proposes for each new list a draft registry entry (its runs, and its Zenodo\n"
            "records found from the PE links of its events, with concept IDs, latest versions and skymap\n"
            "tarballs). Also reports registry records with a newer Zenodo version (used at run time anyway).\n"
        ),
        formatter_class=argparse.RawTextHelpFormatter,
    )
    p_chk.add_argument("--out-json", default=None, help="Optional JSON file with the full report.")
    p_chk.add_argument("--sample-events", type=int, default=3,
                       help="Events of each new list whose PE links are used to find its Zenodo records.")

    p_zen = sub.add_parser(
        "zenodo_releases",
        help="List the Zenodo release versions of each catalog (for --zenodo-version).",
        description=(
            "List the versions of the Zenodo PE/skymap releases, numbered from the oldest (v1).\n"
            "With --data-repo zenodo the latest version is used, unless --zenodo-version selects another.\n"
        ),
        formatter_class=argparse.RawTextHelpFormatter,
    )
    p_zen.add_argument("--catalogs", nargs="+", default=["ALL"], help="Catalog keys, space-separated (e.g. GWTC-3 GWTC-4). ALL key takes them all.")

    return p


def _print_zenodo_releases(catalogs: list[str]) -> None:
    from .data_repo import zenodo_release_parts, zenodo_release_versions

    if "ALL" in catalogs:
        catalogs = zenodo_catalogs()
    for cat in catalogs:
        parts = zenodo_release_parts(cat)
        for n, part in enumerate(parts, 1):
            title = cat if len(parts) == 1 else f"{cat} (part {n} of {len(parts)})"
            print(title)
            versions = zenodo_release_versions(part)
            for i, v in enumerate(versions, 1):
                latest = "  (latest, default)" if i == len(versions) else ""
                print(f"  v{i}  record {v['record_id']:>9}  {v['publication_date']}{latest}")


def main(argv=None) -> int:
    try:
        p = build_parser()
        args = p.parse_args(argv)

        if args.mode == "catalog_statistics":
            catalogs = _parse_catalogs(args.catalogs)
            _validate_catalogs(catalogs)
            run_catalog_statistics(
                catalogs=catalogs,
                out_events_tsv=args.out_events,
                out_report_html=args.out_report,
                include_detectors=args.include_detectors,
                include_area=args.include_area,
                area_cred=args.area_cred,
                data_repo=args.data_repo,
                plots_dir=args.plots_dir,
                zenodo_versions=_parse_zenodo_versions(args.zenodo_version, args.data_repo),
            )
            return 0

        if args.mode == "rates":
            if args.far_threshold <= 0 or args.ns_max_mass <= 1:
                raise ValueError("--far-threshold must be > 0 and --ns-max-mass > 1")
            run_merger_rates(
                out_rates_tsv=args.out_rates,
                out_events_tsv=args.out_events,
                out_report_html=args.out_report,
                plots_dir=args.plots_dir,
                sensitivity_file=args.sensitivity_file,
                sensitivity_release=args.sensitivity_release,
                far_threshold=args.far_threshold,
                ns_max_mass=args.ns_max_mass,
                bbh_kappa=args.bbh_kappa,
                bbh_z_ref=args.bbh_z_ref,
                catalogs=_parse_catalogs(args.catalogs) if args.catalogs else None,
                snr_threshold=args.snr_threshold,
            )
            return 0

        if args.mode == "hubble_constant":
            if args.far_threshold <= 0 or args.pe_samples < 10 or args.parallel < 1 or not 0 < args.min_ess_fraction <= 1:
                raise ValueError("--far-threshold must be > 0, --pe-samples >= 10, --parallel >= 1 and "
                                 "--min-ess-fraction in (0, 1]")
            run_hubble_constant(
                stages=args.stages,
                workdir=args.workdir,
                out_report_html=args.out_report,
                out_summary_tsv=args.out_summary,
                sensitivity_release=args.sensitivity_release,
                sensitivity_file=args.sensitivity_file,
                catalogs=_parse_catalogs(args.catalogs) if args.catalogs else None,
                far_threshold=args.far_threshold,
                snr_threshold=args.snr_threshold,
                min_mass=args.min_mass,
                exclude=args.exclude,
                pe_cache=args.pe_cache,
                keep_pe_files=args.keep_pe_files,
                seeds=args.seeds,
                parallel=args.parallel,
                mass_model=args.mass_model,
                nlive=args.nlive,
                npool=args.npool,
                naccept=args.naccept,
                pe_samples=args.pe_samples,
                inj_fraction=args.inj_fraction,
                min_ess_fraction=args.min_ess_fraction,
                probe_points=args.probe_points,
                reweight_pe_samples=args.reweight_pe_samples,
                icarogw_python=args.icarogw_python,
            )
            return 0

        if args.mode == "bright_siren":
            if args.h0_range[0] <= 0 or args.h0_range[1] <= args.h0_range[0] or args.sky_radius <= 0:
                raise ValueError("--h0-range must be 0 < MIN < MAX and --sky-radius > 0")
            run_bright_siren(
                src_name=args.src_name,
                pe_labels=args.pe_label,
                pe_file=args.pe_file,
                cache_dir=args.cache_dir,
                v_recession=tuple(args.v_recession) if args.v_recession else None,
                v_peculiar=tuple(args.v_peculiar) if args.v_peculiar else None,
                redshift=tuple(args.redshift) if args.redshift else None,
                selection=args.selection,
                sensitivity_release=args.sensitivity_release,
                sensitivity_file=args.sensitivity_file,
                far_threshold=args.far_threshold,
                snr_threshold=args.snr_threshold,
                pe_cache=args.pe_cache,
                sky_radius_deg=args.sky_radius,
                spectral_posterior=args.spectral_posterior,
                h0_range=tuple(args.h0_range),
                out_report_html=args.out_report,
                out_summary_tsv=args.out_summary,
                plots_dir=args.plots_dir,
            )
            return 0

        if args.mode == "event_selection":
            catalogs = _parse_catalogs(args.catalogs)
            _validate_catalogs(catalogs)
            run_event_selection(
                catalogs=catalogs,
                out_tsv=args.out_selection,
                m1_min=args.m1_min,
                m1_max=args.m1_max,
                m2_min=args.m2_min,
                m2_max=args.m2_max,
                dl_min=args.dl_min,
                dl_max=args.dl_max,
                chi_eff_min=args.chi_eff_min,
                chi_eff_max=args.chi_eff_max,
                preset=args.preset,
                ns_max_mass=args.ns_max_mass,
                mass_gap=tuple(args.mass_gap),
                pisn_gap_min=args.pisn_gap_min,
                out_plot=args.out_plot,
            )
            return 0

        if args.mode == "search_skymaps":
            catalogs = _parse_catalogs(args.catalogs)
            _validate_catalogs(catalogs)
            run_search_skymaps(
                catalogs=catalogs,
                out_events_tsv=args.out_events,
                out_report_html=args.out_report,
                ra_deg=args.ra_deg,
                dec_deg=args.dec_deg,
                prob=args.prob,
                plots_dir=args.plots_dir,
                data_repo=args.data_repo,
                skymap_label=args.skymap_label,
                zenodo_versions=_parse_zenodo_versions(args.zenodo_version, args.data_repo),
            )
            return 0

        if args.mode == "parameters_estimation":
            out = run_parameters_estimation(
                src_name=args.src_name,
                plots_dir=args.plots_dir,
                start=args.start,
                stop=args.stop,
                fmin=args.fmin,
                fmax=args.fmax,
                overlay_start=args.overlay_start,
                overlay_stop=args.overlay_stop,
                overlay_fmin=args.overlay_fmin,
                overlay_fmax=args.overlay_fmax,
                q_start=args.q_start,
                q_stop=args.q_stop,
                q_fmin=args.q_fmin,
                q_fmax=args.q_fmax,
                q_fscale=args.q_fscale,
                pe_label=_none_if_empty(args.pe_label),
                waveform_engine=_none_if_empty(args.waveform_engine),
                out_report_html=args.out_report,
                data_repo=args.data_repo,
                pe_vars=args.pe_vars,
                pe_pairs=args.pe_pairs,
                zenodo_versions=_parse_zenodo_versions(args.zenodo_version, args.data_repo),
            )
            # Small manifest (like your previous behavior)
            for k, v in out.items():
                if isinstance(v, list):
                    print(f"[pe] {k}: {len(v)} file(s)")
                else:
                    print(f"[pe] {k}: {v}")
            return 0

        if args.mode == "build_unofficial_pe":
            out = build_unofficial_pe_bundle(
                args.src_name,
                cache_dir=args.cache_dir,
                log_cb=print,
                force_rebuild=bool(args.force),
            )
            if out is None:
                supported = ", ".join(list_unofficial_pe_specs()) or "(none)"
                if get_unofficial_pe_spec(args.src_name) is None:
                    raise ValueError(
                        f"No unofficial PE bundle recipe is registered for {args.src_name}. Supported events: {supported}"
                    )
                raise ValueError(
                    f"Unofficial PE bundle recipe for {args.src_name} exists, but required source files are missing. "
                    "See warnings above for the expected paths."
                )
            print(out)
            return 0

        if args.mode == "check_catalogs":
            from .catalog_check import run_check_catalogs

            run_check_catalogs(out_json=args.out_json, sample_events=args.sample_events)
            return 0

        if args.mode == "zenodo_releases":
            catalogs = _parse_catalogs(args.catalogs)
            bad = [c for c in catalogs if c != "ALL" and c not in zenodo_catalogs()]
            if bad:
                raise ValueError(
                    f"No Zenodo release for: {', '.join(bad)}. "
                    f"Catalogs with Zenodo releases: {', '.join(zenodo_catalogs())}"
                )
            _print_zenodo_releases(catalogs)
            return 0

        raise SystemExit(f"Unsupported mode {args.mode}")

    except ValueError as e:
        # User error → clean message, no traceback
        print(f"Error: {e}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
