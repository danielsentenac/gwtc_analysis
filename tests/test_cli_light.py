"""The parser loads no analysis module: `gwtc_analysis -h` stays fast. Its constants match the modules'."""
from __future__ import annotations

import subprocess
import sys


import pytest

def test_cli_import_loads_no_heavy_package():
    code = ("import sys, gwtc_analysis.cli as c; c.build_parser(); "
            "print(','.join(m for m in ('numpy', 'pandas', 'healpy', 'ligo.skymap', 'lalsimulation', 'matplotlib') "
            "if m in sys.modules))")
    out = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, check=True).stdout.strip()
    assert out == "", f"loaded by the parser: {out}"


def test_cli_defaults_match_the_modules():
    from gwtc_analysis import cli_defaults as d
    from gwtc_analysis.counterpart import COUNTERPARTS
    from gwtc_analysis.ns_eos import EVENTS

    assert tuple(COUNTERPARTS) == d.COUNTERPART_EVENTS
    assert tuple(EVENTS) == d.NS_EOS_EVENTS


@pytest.mark.parametrize("argv", [
    ["search_skymaps", "--catalogs", "GWTC-3", "--ra-deg", "1", "--dec-deg", "1", "--prob", "1.5"],
    ["search_skymaps", "--catalogs", "GWTC-3", "--ra-deg", "1", "--dec-deg", "1", "--prob", "0"],
    ["catalog_statistics", "--catalogs", "GWTC-3", "--include-area", "--area-cred", "90"],
    ["event_selection", "--catalogs", "GWTC-3", "--area-max", "100", "--area-cred", "-0.5"],
])
def test_credible_levels_outside_0_1_are_refused(argv, capsys):
    from gwtc_analysis import cli

    with pytest.raises(SystemExit) as e:
        cli.main(argv)
    assert e.value.code == 2
    assert "above 0 and at most 1" in capsys.readouterr().err


def test_counterpart_help_leaves_out_the_option_it_refuses(capsys):
    from gwtc_analysis import cli

    with pytest.raises(SystemExit):
        cli.main(["counterpart", "-h"])
    assert "--pe-distance-prior" not in capsys.readouterr().out
    with pytest.raises(SystemExit):
        cli.main(["hubble_constant", "-h"])
    assert "--pe-distance-prior" in capsys.readouterr().out
