"""The parser loads no analysis module: `gwtc_analysis -h` stays fast. Its constants match the modules'."""
from __future__ import annotations

import subprocess
import sys


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
