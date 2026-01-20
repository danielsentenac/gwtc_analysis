from __future__ import annotations

import argparse
from typing import List, Optional

from .catalogs import run_catalog_statistics
from .event_selection import run_event_selection
from .search_skymaps import run_search_skymaps
from .parameters_estimation import run_parameters_estimation


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
            "Fetch events for one or more catalogs and compute derived columns.\n"
            "Outputs:\n"
            "  --out-events : TSV table of events\n"
            "  --out-report : HTML report (tables + plots)\n\n"
            "Optional additions:\n"
            "  --include-detectors : detector network via GWOSC calls\n"
            "  --include-area      : sky localization area Axx (requires skymaps)\n"
        ),
        formatter_class=fmt,
    )
    p_cat.add_argument(
        "--catalogs",
        required=True,
        nargs="+",
        help="Catalog keys. Space-separated; commas also accepted (e.g. GWTC-4.0,GWTC-3-confident).",
    )
    p_cat.add_argument("--out-events", required=True, help="Output TSV path (per-event table).")
    p_cat.add_argument("--out-report", required=True, help="Output HTML report path.")

    p_cat.add_argument("--include-detectors", action="store_true", help="Include detector network via GWOSC v2 calls.")
    p_cat.add_argument("--include-area", action="store_true", help="Compute sky localization area Axx if skymaps are available.")
    p_cat.add_argument("--area-cred", type=float, default=0.9, help="Credible level for sky area: 0.9→A90, 0.5→A50, 0.95→A95.")

    p_cat.add_argument("--skymaps-gwtc21", nargs="+", default=None, help="GWTC-2.1 skymaps collection/files (Galaxy may pass a list of file paths).")
    p_cat.add_argument("--skymaps-gwtc3", nargs="+", default=None, help="GWTC-3 skymaps collection/files (Galaxy may pass a list of file paths).")
    p_cat.add_argument("--skymaps-gwtc4", nargs="+", default=None, help="GWTC-4.0 skymaps collection/files (Galaxy may pass a list of file paths).")

    p_cat.add_argument("--events-json", default=None, help="Offline mode: path to local jsonfull-like events JSON (skip GWOSC network calls).")

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
        formatter_class=fmt,
    )
    p_sel.add_argument("--catalogs", required=True, nargs="+", help="Catalog keys. Space-separated; commas also accepted.")
    p_sel.add_argument("--out-selection", required=True, help="Output TSV path for the selected events.")

    p_sel.add_argument("--m1-min", type=float, default=None, help="Minimum primary mass (source frame).")
    p_sel.add_argument("--m1-max", type=float, default=None, help="Maximum primary mass (source frame).")
    p_sel.add_argument("--m2-min", type=float, default=None, help="Minimum secondary mass (source frame).")
    p_sel.add_argument("--m2-max", type=float, default=None, help="Maximum secondary mass (source frame).")
    p_sel.add_argument("--dl-min", type=float, default=None, help="Minimum luminosity distance (Mpc).")
    p_sel.add_argument("--dl-max", type=float, default=None, help="Maximum luminosity distance (Mpc).")

    p_sel.add_argument("--events-json", default=None, help="Offline mode: path to jsonfull-like events JSON (skip network calls).")

    # ---------------------------------------------------------------------
    # search_skymaps
    # ---------------------------------------------------------------------
    p_sky = sub.add_parser(
        "search_skymaps",
        help="Search GW sky localizations for a given sky position (RA/Dec).",
        description=(
            "Given a sky position (RA/Dec in degrees) and a list of skymap FITS files,\n"
            "report which events contain that position above the requested credible level.\n\n"
            "Plotting:\n"
            "  Hit plots are produced. If --plots-dir is not given, plots are written to the current directory.\n"
        ),
        formatter_class=fmt,
    )
    p_sky.add_argument("--catalogs", required=True, nargs="+", help="Catalog keys. Space-separated; commas also accepted.")
    p_sky.add_argument("--ra-deg", type=float, required=True, help="Right ascension (deg).")
    p_sky.add_argument("--dec-deg", type=float, required=True, help="Declination (deg).")
    p_sky.add_argument("--prob", type=float, default=0.9, help="Credible-level threshold (0–1). Common values: 0.9, 0.5, 0.95.")
    p_sky.add_argument("--skymaps", nargs="+", required=True, help="List of skymap files (FITS or FITS.gz). Galaxy collections may expand to a list of files.")
    p_sky.add_argument("--out-events", default=None, help="Optional output TSV path for hits.")
    p_sky.add_argument("--out-report", default=None, help="Optional output HTML report path for hits.")
    p_sky.add_argument("--plots-dir", default=None, help="Directory for hit plots (default: current directory).")
    p_sky.add_argument("--events-json", default=None, help="Offline mode: path to jsonfull-like events JSON (skip GWOSC network calls).")

    # ---------------------------------------------------------------------
    # parameters_estimation
    # ---------------------------------------------------------------------
    p_pe = sub.add_parser(
        "parameters_estimation",
        help="Generate parameter-estimation plots (posteriors, strain, waveforms) for one event.",
        description=(
            "Generate PE plots (posteriors, skymap, strain overlays, PSD) for a single event.\n"
            "Plots are always produced. If --plots-dir is not given, plots are written to the current directory.\n"
        ),
        formatter_class=fmt,
    )
    p_pe.add_argument("--src-name", dest="src_name", required=True, help="Source event name (e.g. GW231223_032836).")
    p_pe.add_argument("--plots-dir", default=None, help="Directory for output PE plots (default: current directory).")
    p_pe.add_argument("--start", type=float, default=0.2, help="Seconds before GPS time for strain window.")
    p_pe.add_argument("--stop", type=float, default=0.1, help="Seconds after GPS time for strain window.")
    p_pe.add_argument("--fs-low", type=float, default=20.0, help="Bandpass low frequency (Hz).")
    p_pe.add_argument("--fs-high", type=float, default=300.0, help="Bandpass high frequency (Hz).")
    p_pe.add_argument("--sample-method", default="Mixed", help="Posterior sample label/model selector (e.g. Mixed, IMRPhenomXPHM, SEOBNRv4PHM).")
    p_pe.add_argument("--strain-approximant", default="IMRPhenomXPHM", help="Waveform model used to generate time-domain waveform for strain overlay.")

    return p


def main(argv=None) -> int:
    p = build_parser()
    args = p.parse_args(argv)

    if args.mode == "catalog_statistics":
        catalogs = _parse_catalogs(args.catalogs)
        run_catalog_statistics(
            catalogs=catalogs,
            out_events_tsv=args.out_events,
            out_report_html=args.out_report,
            include_detectors=args.include_detectors,
            include_area=args.include_area,
            area_cred=args.area_cred,
            skymaps_dirs={
                "GWTC-2.1-confident": _none_if_empty(args.skymaps_gwtc21),
                "GWTC-3-confident": _none_if_empty(args.skymaps_gwtc3),
                "GWTC-4.0": _none_if_empty(args.skymaps_gwtc4),
            },
            events_json=args.events_json,
        )
        return 0

    if args.mode == "event_selection":
        catalogs = _parse_catalogs(args.catalogs)
        run_event_selection(
            catalogs=catalogs,
            out_tsv=args.out_selection,
            events_json=args.events_json,
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
        run_search_skymaps(
            catalogs=catalogs,
            out_events_tsv=args.out_events,
            out_report_html=args.out_report,
            skymaps_dirs={},  # not used when --skymaps list is given
            events_json=args.events_json,
            ra_deg=args.ra_deg,
            dec_deg=args.dec_deg,
            prob=args.prob,
            plots_dir=args.plots_dir,
            skymaps=args.skymaps,
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
            sample_method=args.sample_method,
            strain_approximant=args.strain_approximant,
        )
        # Small manifest (like your previous behavior)
        for k, v in out.items():
            if isinstance(v, list):
                print(f"[pe] {k}: {len(v)} file(s)")
            else:
                print(f"[pe] {k}: {v}")
        return 0

    raise SystemExit(f"Unsupported mode {args.mode}")


if __name__ == "__main__":
    raise SystemExit(main())

