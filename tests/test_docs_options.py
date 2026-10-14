"""Every option of every mode is in an options table of its documentation page, and the reports list them."""
from __future__ import annotations

import argparse
import re
from pathlib import Path

import pytest

from gwtc_analysis import cli, report

DOCS = Path(__file__).resolve().parent.parent / "docs"
PAGES = {
    "catalog_statistics": ["modes/catalog-statistics.md"],
    "rates": ["modes/rates.md"],
    "hubble_constant": ["modes/hubble-constant.md", "modes/hubble-constant-bright.md"],
    "galaxy_catalog": ["modes/galaxy-catalog.md"],
    "counterpart": ["modes/counterpart.md"],
    "continuous_waves": ["modes/continuous-waves.md"],
    "area_law": ["modes/area-law.md"],
    "scalar_cloud": ["modes/scalar-cloud.md"],
    "stochastic": ["modes/stochastic.md"],
    "neutron_star_eos": ["modes/neutron-star-eos.md"],
    "spin_population": ["modes/spin-population.md"],
    "lens_pairs": ["modes/lens-pairs.md"],
    "event_selection": ["modes/event-selection.md"],
    "search_skymaps": ["modes/search-skymaps.md"],
    "parameters_estimation": ["modes/parameters-estimation.md"],
    "build_unofficial_pe": ["modes/unofficial-pe.md"],
    "check_catalogs": ["data-sources.md"],
    "zenodo_releases": ["data-sources.md"],
}
PARSER = cli.build_parser()


def test_every_mode_has_a_page():
    assert set(PARSER._mode_parsers) == set(PAGES)


@pytest.mark.parametrize("mode", sorted(PAGES))
def test_options_in_a_docs_table(mode):
    sp = PARSER._mode_parsers[mode]
    opts = {a.option_strings[-1] for a in sp._actions
            if a.option_strings and a.dest != "help" and a.help != argparse.SUPPRESS}
    rows = "\n".join(l for page in PAGES[mode] for l in (DOCS / page).read_text().splitlines()
                     if l.lstrip().startswith("|"))
    documented = set(re.findall(r"`(--[a-z0-9-]+)", rows))
    if mode == "zenodo_releases":           # its single option is given in the text
        documented |= {"--catalogs"}
    assert sorted(opts - documented) == []


def test_report_lists_all_options(tmp_path):
    p = cli.build_parser()
    args = p.parse_args(["hubble_constant", "--method", "bright", "--sky-radius", "4"])
    given = dict(vars(args))
    cli._resolve_h0_method(args, p._mode_parsers["hubble_constant"])
    cli._set_report_options(args, p, given)
    try:
        out = tmp_path / "r.html"
        report.write_simple_html_report(out, title="t", paragraphs=["x"])
        html = out.read_text()
        assert "Options of this run" in html
        for opt in ("--event", "--population", "--selection", "--h0-range", "--settings"):
            assert f"<code>{opt}</code>" in html
        assert "<code>--nlive</code>" not in html                    # an option of the spectral siren only
        assert "<b>4.0</b>" in html and "<b>hubble_constant_bright</b>" not in html   # set, and resolved defaults
    finally:
        report._RUN_OPTIONS = None
