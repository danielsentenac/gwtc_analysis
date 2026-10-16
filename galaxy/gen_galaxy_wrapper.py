"""Generate the Galaxy tool wrapper (gwtc_analysis.xml) from the command line of the package.

    python galaxy/gen_galaxy_wrapper.py            # writes galaxy/gwtc_analysis.xml

Every mode of cli.py becomes a branch of the tool's mode selector, every option a parameter (same name, type,
default, choices and help), so that the wrapper follows the package. The rules below only say how the options that
are not plain values are handled in Galaxy: output paths (Galaxy datasets), input files (Galaxy data inputs), the
working directories and caches (the job directory), and the cluster settings (left out). Modes left out of the
tool: galaxy_catalog and the spectral and dark methods of hubble_constant, which run icarogw (not on conda) for
hours on many CPUs.
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path
from xml.sax.saxutils import escape, quoteattr

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from gwtc_analysis.cli import build_parser, _h0_method_dests  # noqa: E402
from gwtc_analysis.catalog_registry import allowed_catalogs  # noqa: E402

DOCS = "https://danielsentenac.github.io/gwtc_analysis/modes"

# (branch value, CLI mode, extra fixed arguments, label, docs page)
BRANCHES = [
    ("catalog_statistics", "catalog_statistics", "", "Catalog statistics", "catalog-statistics"),
    ("event_selection", "event_selection", "", "Event selection", "event-selection"),
    ("search_skymaps", "search_skymaps", "", "Search skymaps (sky position, GRB and neutrino coincidences)",
     "search-skymaps"),
    ("parameters_estimation", "parameters_estimation", "", "Parameter estimation plots of one event",
     "parameters-estimation"),
    ("rates", "rates", "", "Merger rates (BNS, NSBH, BBH)", "rates"),
    ("hubble_constant_bright", "hubble_constant", "--method bright", "Hubble constant: bright siren",
     "hubble-constant-bright"),
    ("hubble_constant_joint", "hubble_constant", "--method joint", "Hubble constant: joint posterior",
     "hubble-constant"),
    ("counterpart", "counterpart", "", "Electromagnetic counterpart in the GW posterior", "counterpart"),
    ("neutron_star_eos", "neutron_star_eos", "", "Neutron-star equation of state", "neutron-star-eos"),
    ("spin_population", "spin_population", "", "BBH effective-spin population", "spin-population"),
    ("lens_pairs", "lens_pairs", "", "Strongly lensed pairs", "lens-pairs"),
    ("stochastic", "stochastic", "", "Predicted compact-binary background", "stochastic"),
    ("area_law", "area_law", "", "Hawking's area law with GW250114", "area-law"),
    ("continuous_waves", "continuous_waves", "", "Continuous-wave upper limits", "continuous-waves"),
    ("scalar_cloud", "scalar_cloud", "--stages prepare fit report", "Black hole in a scalar-field cloud (best fit)",
     "scalar-cloud"),
    ("build_unofficial_pe", "build_unofficial_pe", "", "GW170817 PE bundle from the GWTC-1 products",
     "unofficial-pe"),
    ("check_catalogs", "check_catalogs", "", "Check for new catalogs", "../data-sources"),
    ("zenodo_releases", "zenodo_releases", "", "Zenodo release versions", "../data-sources"),
]

# options fixed by the wrapper or left out (paths inside the job, caches, cluster and environment settings)
FIXED = {"plots_dir": "--plots-dir plots", "workdir": "--workdir work", "cache_dir": "--cache-dir cache",
         "pe_cache": "--pe-cache pe_cache"}
SKIP = {"help", "method", "stages", "settings", "keep_pe_files", "icarogw_python", "executor", "slurm_option",
        "slurm_env_setup", "slurm_assembly_option", "submit", "python", "fs_low", "fs_high", "jobs", "workers",
        "npool", "out_plot", "refresh", "refresh_catalogs"}
SLOTS = {("scalar_cloud", "jobs"): "--jobs \\${GALAXY_SLOTS:-1}", ("lens_pairs", "workers"): "--workers \\${GALAXY_SLOTS:-4}"}
# options naming input files: Galaxy data inputs (format, multiple)
FILES = {"sensitivity_file": ("h5", False), "population_posterior": ("tabular", False), "pe_file": ("h5", False),
         "sky_map": ("fits", False), "skymap": ("fits", True), "spectral_posterior": ("tabular", False),
         "rates": ("tabular", False), "inputs": ("tabular", True)}
# outputs: option -> (Galaxy output, local file name, format)
OUTPUTS = {"out_report": ("report", "report.html", "html"), "out_summary": ("table", "table.tsv", "tabular"),
           "out_events": ("table", "table.tsv", "tabular"), "out_rates": ("table", "table.tsv", "tabular"),
           "out_selection": ("table", "table.tsv", "tabular"), "out_json": ("json", "report.json", "json")}
LABELS = {"report": "Report", "table": "Table", "table2": "Events", "json": "Report (JSON)", "text": "Versions",
          "bundle": "GW170817 PE bundle"}
SECOND = {("rates", "out_events"): ("table2", "table2.tsv", "tabular")}     # a second table of a mode
NUMBERS = re.compile(r"^[-+0-9.eE ]+$")
SEL, MULTI, OPT = ' selected="true"', ' multiple="true"', ' optional="true"'


def _actions(mode: str, extra: str, sp: argparse.ArgumentParser):
    acts = [a for a in sp._actions if a.option_strings and a.dest not in SKIP and _pname(a) not in SKIP
            or (mode, a.dest) in SLOTS]
    if mode == "hubble_constant":
        allowed = _h0_method_dests(extra.split()[-1], sp)
        acts = [a for a in acts if a.dest in allowed]
    return acts


def _pname(a) -> str:
    """The Galaxy name of an option's parameter: Galaxy names a param after its argument (--no-skymap-3d ->
    no_skymap_3d), which differs from the argparse dest for negated flags and aliases."""
    return a.option_strings[-1].lstrip("-").replace("-", "_")


def _label(a) -> str:
    return a.option_strings[-1].lstrip("-").replace("-", " ").capitalize()


def _help(a) -> str:
    h = (a.help or "").replace("%(default)s", str(a.default)).replace("%%", "%")
    return h if h and h != argparse.SUPPRESS else ""


def param_xml(mode: str, a, ind: str) -> str:
    """The <param> of an option."""
    name, opt, h = _pname(a), a.option_strings[-1], _help(a)
    attrs = f'argument="{opt}" label={quoteattr(_label(a))} help={quoteattr(h)}'
    if a.dest in FILES:
        fmt, multi = FILES[a.dest]
        attrs_d = attrs.replace(f'argument="{opt}" ', "")
        return (f'{ind}<param name="{name}" {attrs_d} type="data" format="{fmt}" '
                f'optional="{"false" if a.required else "true"}"{MULTI if multi else ""}/>')
    if a.dest == "catalogs":
        default = a.default if isinstance(a.default, list) else [a.default]
        opts = "".join(f'\n{ind}  <option value="{c}"{SEL if c in default else ""}>{c}</option>'
                       for c in allowed_catalogs())
        return f'{ind}<param {attrs} type="select" multiple="true" display="checkboxes">{opts}\n{ind}</param>'
    if isinstance(a, (argparse._StoreTrueAction, argparse._StoreFalseAction)):
        return f'{ind}<param {attrs} type="boolean" truevalue="{opt}" falsevalue="" checked="false"/>'
    if a.choices:
        multi = a.nargs in ("+", "*")
        default = a.default if isinstance(a.default, list) else [a.default]
        opts = "".join(f'\n{ind}  <option value="{c}"{SEL if c in default else ""}>{c}</option>'
                       for c in a.choices)
        extra = MULTI + ("" if a.nargs == "+" and a.default else OPT) if multi else ""
        if not multi and a.default is None:
            extra = OPT
        return f'{ind}<param {attrs} type="select"{extra}>{opts}\n{ind}</param>'
    if a.type in (int, float) and a.nargs is None:
        t = "integer" if a.type is int else "float"
        val = "" if a.default is None else a.default
        return (f'{ind}<param {attrs} type="{t}" value="{val}"'
                f'{OPT if a.default is None else ""}/>')
    # text: single values, lists (space-separated) and fixed-length number lists
    val = a.default if not isinstance(a.default, (list, tuple)) else " ".join(map(str, a.default))
    val = "" if val is None else str(val)
    optional = OPT if a.default is None or a.nargs in ("*", "+") or a.nargs in (2, 3) else ""
    if isinstance(a.nargs, int):
        h2 = f"{a.nargs} numbers separated by a space. " + h
        return (f'{ind}<param {attrs.replace(quoteattr(h), quoteattr(h2))} type="text" value="{escape(val)}"'
                f'{optional}>\n{ind}  <validator type="regex" message="numbers separated by spaces">'
                f'^[-+0-9.eE ]*$</validator>\n{ind}</param>')
    return f'{ind}<param {attrs} type="text" value="{escape(val)}"{optional}/>'


def command_xml(mode: str, a, branch: str) -> str:
    """The Cheetah lines passing an option."""
    v, opt = f"$mode.{_pname(a)}", a.option_strings[-1]
    if (mode, a.dest) in SLOTS:
        return SLOTS[(mode, a.dest)]
    if a.dest in FILES:
        if FILES[a.dest][1]:
            return (f"#if {v}:\n    {opt}\n    #for $f in {v}:\n        '$f'\n    #end for\n#end if")
        return f"#if {v}:\n    {opt} '{v}'\n#end if"
    if isinstance(a, (argparse._StoreTrueAction, argparse._StoreFalseAction)):
        return f"{v}"
    if a.dest == "catalogs" or (a.choices and a.nargs in ("+", "*")):
        return (f"#if {v}:\n    {opt}\n    #for $c in str({v}).split(','):\n        '$c'\n    #end for\n#end if")
    if isinstance(a.nargs, int):
        return f"#if str({v}).strip():\n    {opt} {v}\n#end if"
    if a.nargs in ("+", "*"):
        return (f"#if str({v}).strip():\n    {opt}\n    #for $t in str({v}).split():\n        '$t'\n    #end for\n"
                f"#end if")
    if a.type in (int, float) or a.choices:
        return f"#if str({v}) not in ('', 'None'):\n    {opt} '{v}'\n#end if"
    return f"#if str({v}).strip():\n    {opt} '{v}'\n#end if"


def build() -> str:
    p = build_parser()
    sub = next(a for a in p._actions if isinstance(a, argparse._SubParsersAction))
    helps = {ca.dest: ca.help for ca in sub._choices_actions}
    cmd, whens, outputs_of = [], [], {}
    for branch, mode, extra, label, page in BRANCHES:
        sp = p._mode_parsers[mode]
        acts = _actions(mode, extra, sp)
        dests = {a.dest for a in sp._actions}
        if mode == "hubble_constant":                     # the fixed options and outputs of this method only
            dests &= _h0_method_dests(extra.split()[-1], sp)
        lines = [f"gwtc_analysis {mode}" + (f" {extra}" if extra else "")]
        for a in acts:
            if a.dest in FIXED or a.dest in OUTPUTS or (mode, a.dest) in SECOND:
                continue
            lines.append(command_xml(mode, a, branch))
        lines += [FIXED[d] for d in FIXED if d in dests]
        outs = []
        for d in sorted(dests):
            o = SECOND.get((mode, d)) or OUTPUTS.get(d)
            if o:
                lines.append(f"--{d.replace('_', '-')} '{o[1]}'")
                outs.append(o)
        if mode == "zenodo_releases":
            lines.append("> report.txt")
            outs.append(("text", "report.txt", "txt"))
        if mode == "build_unofficial_pe":
            lines.append("&& cp cache/unofficial_pe/*.h5 bundle.h5")
            outs.append(("bundle", "bundle.h5", "h5"))
        outputs_of[branch] = outs
        lines.append("&& " + " && ".join(f"mv '{o[1]}' '${o[0]}'" for o in dict.fromkeys(outs)))
        body = "\n".join(lines)
        cmd.append(f'#elif str($mode.mode) == "{branch}":\n' + "\n".join("    " + l for l in body.splitlines()))
        params = "\n".join(param_xml(mode, a, "        ") for a in acts if a.dest not in FIXED
                           and a.dest not in OUTPUTS and (mode, a.dest) not in SECOND and (mode, a.dest) not in SLOTS)
        whens.append(f'      <when value="{branch}">\n{params}\n      </when>')
    command = "\n".join(cmd).replace("#elif", "#if", 1) + "\n#end if"
    names = sorted({o[0] for outs in outputs_of.values() for o in outs})
    # outputs
    out_xml = []
    for n in names:
        modes = [b for b, outs in outputs_of.items() if any(o[0] == n for o in outs)]
        fmt = next(o[2] for outs in outputs_of.values() for o in outs if o[0] == n)
        out_xml.append(f'    <data name="{n}" format="{fmt}" label="${{tool.name}} on ${{on_string}}: '
                       f'{LABELS[n]}">\n'
                       f"      <filter>mode['mode'] in {tuple(modes)!r}</filter>\n    </data>")
    plot_modes = [b for b, m, *_ in BRANCHES if "plots_dir" in {a.dest for a in p._mode_parsers[m]._actions}]
    work_modes = [b for b, m, *_ in BRANCHES if "workdir" in {a.dest for a in p._mode_parsers[m]._actions}]
    table_modes = [b for b, outs in outputs_of.items() if any(o[0] == "table" for o in outs)]
    out_xml.append(
        '    <collection name="tables" type="list" label="${tool.name} on ${on_string}: Extra tables">\n'
        '      <discover_datasets pattern="table\\.(?P&lt;designation&gt;.+)\\.tsv" format="tabular" directory="."/>\n'
        f"      <filter>mode['mode'] in {tuple(table_modes)!r}</filter>\n    </collection>")
    out_xml.append(
        '    <collection name="plots" type="list" label="${tool.name} on ${on_string}: Plots">\n'
        '      <discover_datasets pattern="(?P&lt;designation&gt;.+)\\.png" format="png" directory="plots"/>\n'
        '      <discover_datasets pattern="(?P&lt;designation&gt;.+)\\.png" format="png" directory="work/plots"/>\n'
        f"      <filter>mode['mode'] in {tuple(sorted(set(plot_modes) | set(work_modes)))!r}</filter>\n"
        "    </collection>")
    selector = "\n".join(f'        <option value="{b}"{SEL if i == 0 else ""}>{escape(l)}</option>'
                         for i, (b, m, e, l, pg) in enumerate(BRANCHES))
    helprst = "\n".join(f"* **{escape(l)}** (``gwtc_analysis {m}{(' ' + e) if e else ''}``): "
                        f"{escape(helps.get(m) or '')} `Documentation <{DOCS}/{pg}/>`__"
                        for b, m, e, l, pg in BRANCHES)
    return TEMPLATE.format(command=command, selector=selector,
                           whens="\n".join(whens), outputs="\n".join(out_xml), help=helprst)


TEMPLATE = '''<tool id="gwtc_analysis" name="GWTC" version="@TOOL_VERSION@+galaxy@VERSION_SUFFIX@" profile="@PROFILE@">
  <description>Analysis of Gravitational-Wave Transient Catalogs</description>
  <!-- Generated by galaxy/gen_galaxy_wrapper.py from the command line of gwtc_analysis: do not edit by hand. -->
  <macros>
    <import>macros.xml</import>
  </macros>
  <expand macro="requirements"/>
  <command detect_errors="exit_code"><![CDATA[
## pesummary calls getpass.getuser(), which fails in containers whose uid has no passwd entry
export USER="\\${{USER:-galaxy}}" &&
mkdir -p plots work &&
{command}
]]></command>
  <inputs>
    <conditional name="mode">
      <param name="mode" type="select" label="Analysis">
{selector}
      </param>
{whens}
    </conditional>
  </inputs>
  <outputs>
{outputs}
  </outputs>
  <expand macro="tests"/>
  <help><![CDATA[
GWTC analysis
=============

This tool runs **gwtc_analysis** on the public data of the LIGO-Virgo-KAGRA Gravitational-Wave Transient
Catalogs (GWTC-1 to GWTC-5.0): event lists from GWOSC, parameter-estimation files and skymaps from the Zenodo
releases, and the LVK search-sensitivity injections.

Analyses
--------

{help}

Data repository
---------------

Every analysis that downloads files has a **Data repository** option: ``zenodo`` (the public sources), ``s3``
(a MinIO mirror) or ``galaxy`` (the public usegalaxy.org GWTC history). The files read are the same whatever the
repository (the latest Zenodo release versions, or the one selected); a mirror that lacks a file falls back to
Zenodo. Each HTML report lists the files read and the options of the run.

Outputs
-------

An HTML report (plots embedded), the main table of the analysis, extra tables when the analysis writes several,
and the plots as a collection.
]]></help>
  <expand macro="citations"/>
</tool>
'''


if __name__ == "__main__":
    out = Path(__file__).with_name("gwtc_analysis.xml")
    out.write_text(build())
    print(f"wrote {out}")
