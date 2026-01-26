from __future__ import annotations

import argparse
from typing import List, Optional

from .catalogs import run_catalog_statistics
from .event_selection import run_event_selection
from .search_skymaps import run_search_skymaps
from .parameters_estimation import run_parameters_estimation
from .gw_stat import ALLOWED_CATALOGS as ALLOWED_CATALOGS
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
            "  gwtc_analysis event_selection -h\n"
            "  gwtc_analysis search_skymaps -h\n"
            "  gwtc_analysis parameters_estimation -h\n"
        ),
        formatter_class=argparse.RawTextHelpFormatter,
    )

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
        help="Catalog keys, space-separated (e.g. GWTC-1 GWTC-2.1 GWTC-3 GWTC-4). ALL key takes them all.",
    )
    p_cat.add_argument("--out-events", default="catalogs_statistics.tsv", help="Output TSV path (per-event table).")
    p_cat.add_argument("--out-report", default="catalogs_statistics.html", help="Output HTML report path.")
    p_cat.add_argument("--include-detectors", action="store_true", help="Include detector network via GWOSC v2 calls.")
    p_cat.add_argument("--include-area", action="store_true", help="Compute sky localization area Axx if skymaps are available.")
    p_cat.add_argument("--area-cred", type=float, default=0.9, help="Credible level for sky area: 0.9→A90, 0.5→A50, 0.95→A95.")
    p_cat.add_argument("--plots-dir", default="cat_plots", help="Directory for plots (default: cat_plots).")
    p_cat.add_argument("--data-repo", choices=["galaxy", "zenodo", "s3"], default="zenodo", help="Where to read data from: galaxy | zenodo | s3.")

    # ---------------------------------------------------------------------
    # event_selection
    # ---------------------------------------------------------------------
    p_sel = sub.add_parser(
        "event_selection",
        help="Select GW events based on physical criteria (mass, distance) and write a TSV.",
        description=(
            "Select events by simple cuts on source-frame masses and luminosity distance.\n"
            "Cuts are optional; if a cut is not provided, it is not applied.\n"
        ),
        formatter_class=argparse.RawTextHelpFormatter,
    )
    p_sel.add_argument("--catalogs", required=True, nargs="+", help="Catalog keys, space-separated (e.g. GWTC-1 GWTC-2.1 GWTC-3 GWTC-4). ALL key takes them all.")
    p_sel.add_argument("--out-selection", default="event_selection.tsv", help="Output TSV path for the selected events.")
    p_sel.add_argument("--m1-min", type=float, default=None, help="Minimum primary mass (source frame).")
    p_sel.add_argument("--m1-max", type=float, default=None, help="Maximum primary mass (source frame).")
    p_sel.add_argument("--m2-min", type=float, default=None, help="Minimum secondary mass (source frame).")
    p_sel.add_argument("--m2-max", type=float, default=None, help="Maximum secondary mass (source frame).")
    p_sel.add_argument("--dl-min", type=float, default=None, help="Minimum luminosity distance (Mpc).")
    p_sel.add_argument("--dl-max", type=float, default=None, help="Maximum luminosity distance (Mpc).")

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
    p_sky.add_argument("--catalogs", required=True, nargs="+", help="Catalog keys, space-separated (e.g. GWTC-1 GWTC-2.1 GWTC-3 GWTC-4). ALL key takes them all.")
    p_sky.add_argument("--ra-deg", type=float, required=True, help="Right ascension (deg).")
    p_sky.add_argument("--dec-deg", type=float, required=True, help="Declination (deg).")
    p_sky.add_argument("--prob", type=float, default=0.9, help="Credible-level threshold (0–1). Common values: 0.9, 0.5, 0.95.")
    p_sky.add_argument("--waveform", default="Mixed", help="Waveform/approximant selector used to filter skymap (Mixed | IMRPhenomXPHM | SEOBNRv4PHM | SEOBNRv5PHM | IMRPhenomPv2 | IMRPhenomD | IMRPhenomM | IMRPhenomXPHM-SpinTaylor | NRSur7dq4 | NRSur7dq4HM | TaylorF2 | TaylorF2Ecc).")
    p_sky.add_argument("--out-events", default="search_skymaps.tsv", help="Output TSV file (default: search_skymaps.tsv).")
    p_sky.add_argument("--out-report", default="search_skymaps.html", help="Optional output HTML report path for hits.")
    p_sky.add_argument("--plots-dir", default="sky_plots", help="Directory for hit plots (default: sky_plots).")
    p_sky.add_argument("--data-repo", choices=["galaxy", "zenodo", "s3"], default="zenodo", help="Where to read data from: galaxy | zenodo | s3.")

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
    p_pe.add_argument("--pe-vars", nargs="+", default=None, help=("Extra posterior sample variables to plot (space-separated). Example: --pe-vars chi_eff chi_p luminosity_distance."))
    p_pe.add_argument("--pe-pairs", nargs="+", default=None, help=("Extra 2D posterior pairs to plot as 'x:y' tokens. Example: --pe-pairs mass_1_source:mass_2_source chi_eff:chi_p."))
    p_pe.add_argument("--plots-dir", default="pe_plots", help="Directory for output PE plots (default: pe_plots).")
    p_pe.add_argument("--start", type=float, default=0.2, help="Seconds before GPS time for strain window.")
    p_pe.add_argument("--stop", type=float, default=0.1, help="Seconds after GPS time for strain window.")
    p_pe.add_argument("--fs-low", type=float, default=20.0, help="Bandpass low frequency (Hz).")
    p_pe.add_argument("--fs-high", type=float, default=300.0, help="Bandpass high frequency (Hz).")

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

    return p


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
                waveform=args.waveform,
            )
            return 0

        if args.mode == "parameters_estimation":
            out = run_parameters_estimation(
                src_name=args.src_name,
                plots_dir=args.plots_dir,
                start=args.start,
                stop=args.stop,
                fs_low=args.fs_low,
                fs_high=args.fs_high,
                pe_label=_none_if_empty(args.pe_label),
                waveform_engine=_none_if_empty(args.waveform_engine),
                out_report_html=args.out_report,
                data_repo=args.data_repo,
                pe_vars=args.pe_vars,
                pe_pairs=args.pe_pairs,
            )
            # Small manifest (like your previous behavior)
            for k, v in out.items():
                if isinstance(v, list):
                    print(f"[pe] {k}: {len(v)} file(s)")
                else:
                    print(f"[pe] {k}: {v}")
            return 0

        raise SystemExit(f"Unsupported mode {args.mode}")

    except ValueError as e:
        # User error → clean message, no traceback
        print(f"Error: {e}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())

