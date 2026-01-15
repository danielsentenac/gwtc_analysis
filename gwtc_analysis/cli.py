from __future__ import annotations

import argparse
from typing import List

from .catalogs import run_catalog_statistics
from .event_selection import run_event_selection
from .search_skymaps import run_search_skymaps


def _split_csv(s: str) -> List[str]:
    return [x.strip() for x in s.split(",") if x.strip()]


def _none_if_empty(x):
    """Argparse with nargs can yield [] instead of None."""
    if x is None:
        return None
    if isinstance(x, list) and len(x) == 0:
        return None
    return x


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="gwtc_analysis", description="GWTC analysis tool (Galaxy-friendly).")
    p.add_argument(
        "--mode",
        required=True,
        choices=["catalog_statistics", "event_selection","search_skymaps"],
        help="Which analysis to run.",
    )
    p.add_argument(
        "--catalogs",
        required=True,
        nargs="+",
        help="Catalog keys (space-separated). Galaxy may pass multiple. Commas are also accepted.",
    )

    # Outputs
    p.add_argument("--out-events", default=None, help="Output TSV (catalog_statistics)")
    p.add_argument("--out-report", default=None, help="Output HTML report (catalog_statistics)")
    p.add_argument("--out-selection", default=None, help="Output TSV for event_selection")

    # catalog_statistics options
    p.add_argument("--include-detectors", action="store_true", help="Include detector network via GWOSC v2 calls")
    p.add_argument("--include-a90", action="store_true", help="Compute A90 sky area if skymaps are available")
    p.add_argument("--a90-cred", type=float, default=0.9, help="Credible level for sky area (default 0.9)")

    # Skymaps inputs for A90 (Galaxy collections typically map to a list of files)
    p.add_argument("--skymaps-gwtc21", nargs="+", default=None, help="GWTC-2.1 skymaps collection/files")
    p.add_argument("--skymaps-gwtc3", nargs="+", default=None, help="GWTC-3 skymaps collection/files")
    p.add_argument("--skymaps-gwtc4", nargs="+", default=None, help="GWTC-4 skymaps collection/files")

    # event_selection filters
    p.add_argument("--m1-min", type=float, default=None)
    p.add_argument("--m1-max", type=float, default=None)
    p.add_argument("--m2-min", type=float, default=None)
    p.add_argument("--m2-max", type=float, default=None)
    p.add_argument("--dl-min", type=float, default=None)
    p.add_argument("--dl-max", type=float, default=None)
    
    # -----  new mode: search_skymaps  -----
    p.add_argument("--ra-deg",  type=float, help="Right ascension (deg)")
    p.add_argument("--dec-deg", type=float, help="Declination (deg)")
    p.add_argument("--prob",    type=float, default=0.9, help="Credible-level threshold (0–1)")
    p.add_argument("--make-plots", choices=["none", "hits", "all"], default="hits")
    p.add_argument("--plots-dir", default=None, help="Directory for skymap plots (default: no plots)")
    # allow Galaxy collections to be passed as a plain list
    p.add_argument("--skymaps", nargs="*", default=None, help="List of skymap files (FITS or FITS.gz)")

    # Offline mode
    p.add_argument("--events-json", default=None, help="Offline mode: path to jsonfull-like events JSON")

    return p


def main(argv=None) -> int:
    p = build_parser()
    args = p.parse_args(argv)

    raw_cats: List[str] = []
    for item in args.catalogs:
        raw_cats.extend(_split_csv(item))
    catalogs = raw_cats

    if args.mode == "catalog_statistics":
        if not args.out_events or not args.out_report:
            raise SystemExit("--out-events and --out-report are required for mode=catalog_statistics")

        run_catalog_statistics(
            catalogs=catalogs,
            out_events_tsv=args.out_events,
            out_report_html=args.out_report,
            include_detectors=args.include_detectors,
            include_a90=args.include_a90,
            a90_cred=args.a90_cred,
            skymaps_dirs={
                "GWTC-2.1-confident": _none_if_empty(args.skymaps_gwtc21),
                "GWTC-3-confident": _none_if_empty(args.skymaps_gwtc3),
                "GWTC-4.0": _none_if_empty(args.skymaps_gwtc4),
            },
            events_json=args.events_json,
        )
        return 0

    if args.mode == "event_selection":
        if not args.out_selection:
            raise SystemExit("--out-selection is required for mode=event_selection")

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
        if args.ra_deg is None or args.dec_deg is None:
            raise SystemExit("--ra-deg and --dec-deg are required for search_skymaps")
        run_search_skymaps(
            catalogs=catalogs,
            out_events_tsv=args.out_events,
            out_report_html=args.out_report,
            skymaps_dirs={},          # not used when --skymaps list is given
            events_json=args.events_json,
            ra_deg=args.ra_deg,
            dec_deg=args.dec_deg,
            prob=args.prob,
            make_plots=args.make_plots,
            plots_dir=args.plots_dir,
            skymaps=args.skymaps,     # single list from Galaxy
        )
        return 0

    raise SystemExit(f"Unsupported mode {args.mode}")


if __name__ == "__main__":
    raise SystemExit(main())
