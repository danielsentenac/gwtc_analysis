"""Generate the Galaxy tool wrapper (gwtc_analysis.xml) from the command line of the package.

    python galaxy/gen_galaxy_wrapper.py            # writes galaxy/gwtc_analysis.xml

Every mode of cli.py becomes a branch of the tool's mode selector, every option a parameter (same name, type,
default, choices and help), so that the wrapper follows the package. The rules below only say how the options that
are not plain values are handled in Galaxy: output paths (Galaxy datasets), input files (Galaxy data inputs), the
working directories and caches (the job directory), and the cluster settings (left out). Modes left out of the
tool: galaxy_catalog and the spectral and dark methods of hubble_constant, which run icarogw (not on conda) for
hours on many CPUs, and check_catalogs, a maintenance check of the package. The tests are written below (TESTS),
with the reason why some modes have none.
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
    ("primordial_black_holes", "primordial_black_holes", "", "Primordial black holes: limits on the dark-matter fraction",
     "primordial-black-holes"),
    ("area_law", "area_law", "", "Hawking's area law with GW250114", "area-law"),
    ("continuous_waves", "continuous_waves", "", "Continuous-wave upper limits", "continuous-waves"),
    ("scalar_cloud", "scalar_cloud", "--stages prepare fit report", "Black hole in a scalar-field cloud (best fit)",
     "scalar-cloud"),
    ("build_unofficial_pe", "build_unofficial_pe", "", "GW170817 PE bundle from the GWTC-1 products",
     "unofficial-pe"),
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
           "out_selection": ("table", "table.tsv", "tabular"), "out_json": ("json", "report.json", "json"),
           "out_limits": ("table", "table.tsv", "tabular")}
LABELS = {"report": "Report", "table": "Table", "table2": "Events", "json": "Report (JSON)", "text": "Versions",
          "bundle": "GW170817 PE bundle"}
SECOND = {("rates", "out_events"): ("table2", "table2.tsv", "tabular")}     # a second table of a mode
# the description of a branch, when the help of its mode covers more than the branch
DESCRIPTIONS = {"hubble_constant_bright": "Hubble constant from an event with an identified host galaxy (bright siren).",
                "hubble_constant_joint": "Joint Hubble-constant posterior of independent results."}
# Galaxy labels of the options (the flag itself is shown by argument=); an option not listed gets its flag
PARAM_LABELS = {
    "area_cred": "Credible level of the sky area", "bbh_kappa": "BBH rate evolution exponent kappa",
    "bbh_z_ref": "Redshift of the BBH rate", "boson_label": "Label of the boson-mass band",
    "boson_mass": "Boson-mass band (eV)", "catalogs": "Catalogs", "chi_eff_max": "Maximum effective spin",
    "chi_eff_min": "Minimum effective spin", "counterpart_class": "Class of the binary",
    "counterparts": "Messengers", "data_repo": "Data repository", "dec": "Declination of the counterpart (deg)",
    "dec_deg": "Declination (deg)", "dl_max": "Maximum luminosity distance (Mpc)",
    "dl_min": "Minimum luminosity distance (Mpc)", "event": "Event", "events": "Events",
    "exclude": "Events left out", "far_threshold": "False-alarm-rate threshold (per year)",
    "fmax": "Maximum frequency (Hz)", "fmin": "Minimum frequency (Hz)", "force": "Rebuild the bundle",
    "fpp_threshold": "False-positive-probability threshold", "galaxies": "Galaxy catalog",
    "galaxy_max_area": "Largest sky area for the galaxy query (deg²)",
    "grb_window": "GRB time window (s)", "h0_range": "H0 prior range (km/s/Mpc)",
    "high_z": "BBH rate beyond the detected events", "host": "Host galaxy", "include_area": "Include the sky areas",
    "include_detectors": "Include the detector networks", "inputs": "Results to combine",
    "lambda_max": "Upper bound of the tidal-deformability prior", "m1_max": "Maximum primary mass (M☉)",
    "m1_min": "Minimum primary mass (M☉)", "m2_max": "Maximum secondary mass (M☉)",
    "m2_min": "Minimum secondary mass (M☉)", "mass_gap": "Lower mass gap (M☉)", "mass_model": "Mass model",
    "max_injections": "Number of injections used", "min_mass": "Minimum component mass (M☉)",
    "n_background": "Number of background injections", "n_draws": "Number of posterior draws",
    "neutrino_window": "Neutrino time window (s)", "no_correlation": "No correlation with the mass ratio",
    "skymap_3d": "Skip the 3D sky map", "ns_max_mass": "Maximum neutron-star mass (M☉)",
    "overlay_fmax": "Overlay: maximum frequency (Hz)", "overlay_fmin": "Overlay: minimum frequency (Hz)",
    "overlay_start": "Overlay: seconds before the merger", "overlay_stop": "Overlay: seconds after the merger",
    "pe_distance_prior": "Distance prior of the PE samples", "pe_file": "PE file",
    "pe_label": "PE analysis (label)", "pe_pairs": "Extra 2D posterior pairs", "pe_samples": "PE samples per event",
    "pe_vars": "Extra posterior variables", "pisn_gap_min": "Lower edge of the pair-instability gap (M☉)",
    "population": "Population of the sources", "population_posterior": "Population posterior (TSV)",
    "preset": "Class of sources", "prob": "Credible level", "q_fmax": "Q-transform: maximum frequency (Hz)",
    "q_fmin": "Q-transform: minimum frequency (Hz)", "q_fscale": "Q-transform: frequency axis",
    "q_start": "Q-transform: seconds before the merger", "q_stop": "Q-transform: seconds after the merger",
    "ra": "Right ascension of the counterpart (deg)", "ra_deg": "Right ascension (deg)",
    "rates": "Merger rates (TSV of the rates mode)", "redshift": "Host redshift and uncertainty",
    "sample_events": "Sample events", "sample_models": "Models to sample", "sampling_frequency": "Sampling frequency (Hz)",
    "seed": "Random seed", "selection": "Selection term", "sensitivity_file": "Injection file",
    "sensitivity_release": "Search-sensitivity release", "sky_map": "Sky map (FITS)",
    "sky_radius": "Sky radius (deg)", "skymap": "Sky maps (FITS)", "skymap_label": "Sky-map analysis (label)",
    "snr_threshold": "SNR threshold of the O1–O2 injections", "sources": "External catalogs",
    "spectral_posterior": "Spectral-siren posterior (TSV)", "spin_prior": "Spin prior of the PE analyses",
    "src_name": "Event name", "start": "Seconds before the merger", "steps": "Sampler steps",
    "stop": "Seconds after the merger", "tau_d": "Time from black-hole formation to merger (yr)",
    "v_peculiar": "Peculiar velocity of the host (km/s)", "v_recession": "Recession velocity of the host (km/s)",
    "viewing_angle": "Viewing-angle constraint (deg)", "waveform_engine": "Waveform model of the overlay",
    "with_imr": "Also show the full-signal (IMR) analysis", "workers": "Parallel downloads",
    "z_horizon": "Redshift of the farthest events", "zenodo_version": "Zenodo release versions",
    "masses": "PBH masses (M☉)", "mass_window": "Half-width of the mass window (in ln m)",
    "attribute": "Events that may be primordial", "suppression": "Suppression of the early-binary rate",
    "sub_solar": "Leave out the sub-solar searches",
}
# help texts written for Galaxy, where the package's mention work directories or values a data input cannot take
HELP_OVERRIDES = {
    ("hubble_constant", "inputs"): "At least two independent H0 results to combine: posterior TSV files of H0 "
                                   "samples (a column H0), or an H0 grid with a p column.",
    ("rates", "population_posterior"): "Posterior TSV of a hubble_constant run with the same mass model: the BBH "
                                       "rate (with fullpop4, all the rates) is then computed over its samples.",
    ("stochastic", "spectral_posterior"): "Posterior TSV of a hubble_constant spectral siren (mass model plp or mltp).",
    ("counterpart", "sky_map"): "FITS sky map for the searched probability of the position. Leave empty for the "
                                "event's LVK map (when its PE file comes from Zenodo), else a kernel estimate on "
                                "the PE samples.",
    ("parameters_estimation", "galaxies"): "Galaxy catalog cross-matched with the 3D sky map: GLADE+ (via VizieR), "
                                           "none, or a CSV/TSV dataset with ra, dec and dist [Mpc] or z columns.",
}
# event names: a select when the mode knows a single event, a regex otherwise
EVENT_CHOICES = {("build_unofficial_pe", "src_name"): ["GW170817"]}
VALIDATORS = {("parameters_estimation", "src_name"): ("^GW[0-9]{6}(_[0-9]{6})?$",
                                                     "an event name such as GW170817 or GW231223_032836")}
# a nested conditional choosing between groups of options (mode -> name, selector, label, [(value, label, dests)])
GROUPS = {"search_skymaps": ("search", "by", "Search by", [
    ("position", "Sky position (RA, Dec)", ["ra_deg", "dec_deg"]),
    ("counterparts", "GRB and neutrino coincidences", ["counterparts", "sources", "grb_window", "neutrino_window",
                                                       "skymap"])])}
GROUP_REQUIRED = {"ra_deg", "dec_deg", "counterparts"}
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
    return PARAM_LABELS.get(_pname(a)) or PARAM_LABELS.get(a.dest) or a.option_strings[-1].lstrip("-").replace(
        "-", " ").capitalize()


def _help(a, mode: str = "") -> str:
    if (mode, a.dest) in HELP_OVERRIDES:
        return HELP_OVERRIDES[(mode, a.dest)]
    h = (a.help or "").replace("%(default)s", str(a.default)).replace("%%", "%")
    h = re.sub(r"^\[(bright|joint)\] ", "", h)
    if mode in GROUPS:                    # the nested conditional already separates the groups
        h = re.sub(r"; not with --counterparts", "", re.sub(r"^With --counterparts: ", "", h))
        h = h[:1].upper() + h[1:]
    return h if h and h != argparse.SUPPRESS else ""


def param_xml(mode: str, a, ind: str) -> str:
    """The <param> of an option."""
    name, opt, h = _pname(a), a.option_strings[-1], _help(a, mode)
    attrs = f'argument="{opt}" label={quoteattr(_label(a))} help={quoteattr(h)}'
    required = a.required or a.dest in GROUP_REQUIRED or (mode, a.dest) == ("hubble_constant", "inputs")
    if (mode, a.dest) == ("parameters_estimation", "galaxies"):
        return (f'{ind}<conditional name="galaxies">\n'
                f'{ind}  <param name="source" type="select" label="Galaxy catalog" help={quoteattr(h)}>\n'
                f'{ind}    <option value="glade" selected="true">GLADE+ (VizieR)</option>\n'
                f'{ind}    <option value="none">None</option>\n'
                f'{ind}    <option value="history">From the history</option>\n'
                f'{ind}  </param>\n{ind}  <when value="glade"/>\n{ind}  <when value="none"/>\n'
                f'{ind}  <when value="history">\n'
                f'{ind}    <param name="file" type="data" format="csv,tabular" label="Galaxy catalog file" '
                f'help="CSV or TSV with ra, dec and dist [Mpc] or z columns."/>\n'
                f'{ind}  </when>\n{ind}</conditional>')
    if (mode, a.dest) in EVENT_CHOICES:
        opts = "".join(f'\n{ind}  <option value="{c}"{SEL if i == 0 else ""}>{c}</option>'
                       for i, c in enumerate(EVENT_CHOICES[(mode, a.dest)]))
        return f'{ind}<param {attrs} type="select">{opts}\n{ind}</param>'
    if a.dest in FILES:
        fmt, multi = FILES[a.dest]
        attrs_d = attrs.replace(f'argument="{opt}" ', "")
        extra = ' min="2"' if (mode, a.dest) == ("hubble_constant", "inputs") else ""
        return (f'{ind}<param name="{name}" {attrs_d} type="data" format="{fmt}" '
                f'optional="{"false" if required else "true"}"{MULTI if multi else ""}{extra}/>')
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
        extra = MULTI + ("" if (a.nargs == "+" and a.default) or required else OPT) if multi else ""
        if not multi and a.default is None and not required:
            extra = OPT
        return f'{ind}<param {attrs} type="select"{extra}>{opts}\n{ind}</param>'
    if a.type in (int, float) and a.nargs is None:
        t = "integer" if a.type is int else "float"
        val = "" if a.default is None else a.default
        return (f'{ind}<param {attrs} type="{t}" value="{val}"'
                f'{OPT if a.default is None and not required else ""}/>')
    # text: single values, lists (space-separated) and fixed-length number lists
    val = a.default if not isinstance(a.default, (list, tuple)) else " ".join(map(str, a.default))
    val = "" if val is None else str(val)
    optional = OPT if (a.default is None or a.nargs in ("*", "+") or a.nargs in (2, 3)) and not required else ""
    if (mode, a.dest) in VALIDATORS:
        rx, msg = VALIDATORS[(mode, a.dest)]
        return (f'{ind}<param {attrs} type="text" value="{escape(val)}"{optional}>\n'
                f'{ind}  <validator type="regex" message={quoteattr("Expected " + msg)}>{escape(rx)}</validator>\n'
                f'{ind}</param>')
    if isinstance(a.nargs, int):
        h2 = f"{a.nargs} numbers separated by a space. " + h
        return (f'{ind}<param {attrs.replace(quoteattr(h), quoteattr(h2))} type="text" value="{escape(val)}"'
                f'{optional}>\n{ind}  <validator type="regex" message="numbers separated by spaces">'
                f'^[-+0-9.eE ]*$</validator>\n{ind}</param>')
    return f'{ind}<param {attrs} type="text" value="{escape(val)}"{optional}/>'


def command_xml(mode: str, a, branch: str, prefix: str = "$mode") -> str:
    """The Cheetah lines passing an option (`prefix`: the conditional holding its parameter)."""
    v, opt = f"{prefix}.{_pname(a)}", a.option_strings[-1]
    if (mode, a.dest) in SLOTS:
        return SLOTS[(mode, a.dest)]
    if (mode, a.dest) == ("parameters_estimation", "galaxies"):
        return (f"#if str({prefix}.galaxies.source) == 'history':\n    {opt} '{prefix}.galaxies.file'\n#else:\n"
                f"    {opt} '{prefix}.galaxies.source'\n#end if")
    if (mode, a.dest) in EVENT_CHOICES:
        return f"{opt} '{v}'"
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
        grouped = {d for _, _, ds in GROUPS.get(mode, ("", "", "", []))[3] for d in ds}
        for a in acts:
            if a.dest in FIXED or a.dest in OUTPUTS or (mode, a.dest) in SECOND or a.dest in grouped:
                continue
            lines.append(command_xml(mode, a, branch))
        if mode in GROUPS:
            gname, gsel, _, groups = GROUPS[mode]
            by = {a.dest: a for a in acts}
            for i, (val, _, ds) in enumerate(groups):
                lines.append(("#if" if i == 0 else "#elif") + f" str($mode.{gname}.{gsel}) == '{val}':")
                for d in ds:
                    lines += ["    " + l for l in command_xml(mode, by[d], branch, f"$mode.{gname}").splitlines()]
            lines.append("#end if")
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
            lines.append("&& cp cache/unofficial_pe/*.h5 '$bundle'")
        outputs_of[branch] = outs + ([("bundle", "bundle.h5", "h5")] if mode == "build_unofficial_pe" else [])
        if outs:
            lines.append("&& " + " && ".join(f"mv '{o[1]}' '${o[0]}'" for o in dict.fromkeys(outs)))
        body = "\n".join(lines)
        cmd.append(f'#elif str($mode.mode) == "{branch}":\n' + "\n".join("    " + l for l in body.splitlines()))
        params = "\n".join(param_xml(mode, a, "        ") for a in acts if a.dest not in FIXED and a.dest not in grouped
                           and a.dest not in OUTPUTS and (mode, a.dest) not in SECOND and (mode, a.dest) not in SLOTS)
        if mode in GROUPS:
            gname, gsel, glabel, groups = GROUPS[mode]
            by = {a.dest: a for a in acts}
            opts = "\n".join(f'            <option value="{v}"{SEL if i == 0 else ""}>{escape(l)}</option>'
                             for i, (v, l, _) in enumerate(groups))
            whens_g = "\n".join(f'          <when value="{v}">\n'
                                + "\n".join(param_xml(mode, by[d], "            ") for d in ds)
                                + "\n          </when>" for v, _, ds in groups)
            params = (f'        <conditional name="{gname}">\n          <param name="{gsel}" type="select" '
                      f'label="{escape(glabel)}">\n{opts}\n          </param>\n{whens_g}\n        </conditional>\n'
                      + params)
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
                        f"{escape(DESCRIPTIONS.get(b) or helps.get(m) or '')} `Documentation <{DOCS}/{pg}/>`__"
                        for b, m, e, l, pg in BRANCHES)
    return TEMPLATE.format(command=command, selector=selector,
                           whens="\n".join(whens), outputs="\n".join(out_xml), help=helprst, tests=TESTS)


# the tool tests (the comment at the top says why some modes have none)
TESTS = """  <tests>
    <!-- Modes without a test, and why: rates and primordial_black_holes download the LVK search-sensitivity
         injections (~900 MB); neutron_star_eos the PE files of two events and a joint fit; spin_population the PE
         files of ~140 events and an MCMC; lens_pairs the PE files of every O4a event; scalar_cloud a PE file and a
         fit of several minutes. -->
    <test expect_num_outputs="4">
      <conditional name="mode">
        <param name="mode" value="catalog_statistics"/>
        <param name="catalogs" value="GWTC-2.1"/>
        <param name="data_repo" value="zenodo"/>
      </conditional>
      <output name="table" ftype="tabular">
        <assert_contents>
          <has_text text="event_id"/>
          <has_n_lines min="40"/>
        </assert_contents>
      </output>
      <output name="report" ftype="html">
        <assert_contents>
          <has_text text="Options of this run"/>
        </assert_contents>
      </output>
      <output_collection name="plots" type="list" count="4">
        <element name="GWTC-2.1_histograms" ftype="png">
          <assert_contents>
            <has_size min="10000"/>
          </assert_contents>
        </element>
        <element name="GWTC-2.1_m1_m2_snr" ftype="png">
          <assert_contents>
            <has_size min="10000"/>
          </assert_contents>
        </element>
        <element name="GWTC-2.1_remnants" ftype="png">
          <assert_contents>
            <has_size min="10000"/>
          </assert_contents>
        </element>
        <element name="source_types_pie" ftype="png">
          <assert_contents>
            <has_size min="10000"/>
          </assert_contents>
        </element>
      </output_collection>
    </test>
    <test expect_num_outputs="2">
      <conditional name="mode">
        <param name="mode" value="event_selection"/>
        <param name="catalogs" value="GWTC-2.1"/>
        <param name="preset" value="hierarchical"/>
      </conditional>
      <output name="table" ftype="tabular">
        <assert_contents>
          <has_text text="event_id"/>
          <has_n_lines min="2"/>
        </assert_contents>
      </output>
    </test>
    <test expect_num_outputs="4">
      <conditional name="mode">
        <param name="mode" value="search_skymaps"/>
        <param name="catalogs" value="GWTC-2.1"/>
        <conditional name="search">
          <param name="by" value="position"/>
          <param name="ra_deg" value="180.0"/>
          <param name="dec_deg" value="0.0"/>
        </conditional>
        <param name="data_repo" value="galaxy"/>
      </conditional>
      <output name="table" ftype="tabular">
        <assert_contents>
          <has_text text="event_id"/>
          <has_n_lines min="2"/>
        </assert_contents>
      </output>
    </test>
    <test expect_num_outputs="2">
      <conditional name="mode">
        <param name="mode" value="parameters_estimation"/>
        <param name="src_name" value="GW190828_063405"/>
        <param name="data_repo" value="galaxy"/>
      </conditional>
      <output name="report" ftype="html">
        <assert_contents>
          <has_text text="PE files of this run"/>
          <has_text text="IGWN-GWTC2p1-v2-GW190828_063405_PEDataRelease_mixed_cosmo.h5"/>
        </assert_contents>
      </output>
    </test>
    <test expect_num_outputs="4">
      <conditional name="mode">
        <param name="mode" value="hubble_constant_joint"/>
        <param name="inputs" value="h0_fullpop4_posterior.tsv,h0_plp_posterior.tsv"/>
      </conditional>
      <output name="table" ftype="tabular">
        <assert_contents>
          <has_text text="joint"/>
          <has_n_lines n="4"/>
        </assert_contents>
      </output>
    </test>
    <test expect_num_outputs="4">
      <conditional name="mode">
        <param name="mode" value="stochastic"/>
        <param name="spectral_posterior" value="h0_plp_posterior.tsv"/>
        <param name="rates" value="merger_rates.tsv"/>
      </conditional>
      <output name="table" ftype="tabular">
        <assert_contents>
          <has_text text="omega_25Hz_median"/>
        </assert_contents>
      </output>
      <output_collection name="plots" type="list" count="1">
        <element name="omega_gw" ftype="png">
          <assert_contents>
            <has_size min="10000"/>
          </assert_contents>
        </element>
      </output_collection>
    </test>
    <test expect_num_outputs="4">
      <conditional name="mode">
        <param name="mode" value="hubble_constant_bright"/>
      </conditional>
      <output name="report" ftype="html">
        <assert_contents>
          <has_text text="Unofficial-GWTC1-GW170817_PEDataRelease.h5"/>
        </assert_contents>
      </output>
    </test>
    <test expect_num_outputs="4">
      <conditional name="mode">
        <param name="mode" value="counterpart"/>
        <param name="event" value="GW170817"/>
      </conditional>
      <output name="table" ftype="tabular">
        <assert_contents>
          <has_text text="C02:IMRPhenomPv2_NRTidal-LowSpin"/>
          <has_text text="viewing_angle_median"/>
          <has_n_lines n="3"/>
        </assert_contents>
      </output>
      <output_collection name="plots" type="list" count="2">
        <element name="distance_GW170817" ftype="png">
          <assert_contents>
            <has_size min="10000"/>
          </assert_contents>
        </element>
        <element name="distance_inclination_GW170817" ftype="png">
          <assert_contents>
            <has_size min="10000"/>
          </assert_contents>
        </element>
      </output_collection>
    </test>
    <test expect_num_outputs="1">
      <conditional name="mode">
        <param name="mode" value="build_unofficial_pe"/>
        <param name="src_name" value="GW170817"/>
      </conditional>
      <output name="bundle" ftype="h5">
        <assert_contents>
          <has_size min="100000000"/>
        </assert_contents>
      </output>
    </test>
    <test expect_num_outputs="4">
      <conditional name="mode">
        <param name="mode" value="continuous_waves"/>
      </conditional>
      <output name="report" ftype="html">
        <assert_contents>
          <has_text text="Public data files of this run"/>
        </assert_contents>
      </output>
    </test>
    <test expect_num_outputs="4">
      <conditional name="mode">
        <param name="mode" value="area_law"/>
      </conditional>
      <output name="table" ftype="tabular">
        <assert_contents>
          <has_n_lines min="2"/>
        </assert_contents>
      </output>
    </test>
    <test expect_num_outputs="1">
      <conditional name="mode">
        <param name="mode" value="zenodo_releases"/>
        <param name="catalogs" value="GWTC-3"/>
      </conditional>
      <output name="text">
        <assert_contents>
          <has_text text="22685054"/>
        </assert_contents>
      </output>
    </test>
  </tests>"""


TEMPLATE = '''<tool id="gwtc_analysis" name="GWTC" version="@TOOL_VERSION@+galaxy@VERSION_SUFFIX@" profile="@PROFILE@">
  <description>Analysis of Gravitational-Wave Transient Catalogs</description>
  <!-- Generated by galaxy/gen_galaxy_wrapper.py (gwtc_analysis repository) from the command line of the package:
       regenerate it there when the options change, so that the wrapper follows them. -->
  <macros>
    <import>macros.xml</import>
  </macros>
  <expand macro="requirements"/>
  <command detect_errors="exit_code"><![CDATA[
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
{tests}
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
