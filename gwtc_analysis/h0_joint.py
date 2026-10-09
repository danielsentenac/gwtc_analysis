"""Joint Hubble constant (``hubble_constant --method joint``): the product of independent H0 posteriors.

The inputs are work directories of the other methods (or posterior TSV files):

- ``spectral`` or ``dark``: an icarogw run (``posterior_reweighted.tsv`` or ``posterior.tsv``), its H0 samples turned
  into a density by a Gaussian kernel estimate reflected at the prior bounds;
- ``bright``: the posterior grid of a bright siren (``posterior_grid.tsv``, ``bright.json``).

All of them have the same flat prior on H0 (10-200 km/s/Mpc), so the joint posterior of independent data is the
normalized product of the posteriors. Independence is checked: a dark siren already contains the spectral siren of
its events, so at most one spectral or dark input is accepted, and no event may enter two inputs (e.g. GW190521 as a
bright siren and in a spectral siren: leave it out of the latter with ``--exclude``).
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Iterable, Optional

import numpy as np
import pandas as pd

from .h0_bright import BRIGHT_FILE, GRID_FILE, H0_PRIOR, _fmt, _normalize, _plot, summarize
from .report import write_simple_html_report

ICAROGW_POSTERIORS = ("posterior_reweighted.tsv", "posterior.tsv")


def _log(msg: str) -> None:
    print(f"[hubble_constant] {msg}", flush=True)


def siren_kind(path: str | Path) -> str:
    """'dark siren (<band>, ...)' when the hubble_constant work directory of `path` used a galaxy catalog, else
    'spectral siren'."""
    path = Path(path).expanduser()
    sel = (path if path.is_dir() else path.parent) / "selection.json"
    if sel.exists():
        cat = json.loads(sel.read_text()).get("galaxy_catalog")
        if cat:
            return f"dark siren ({cat.get('band', 'galaxy catalog')}, mass spectrum + galaxies)"
    return "spectral siren"


def sample_density(h0: np.ndarray, samples: np.ndarray, bounds: tuple[float, float] = H0_PRIOR) -> np.ndarray:
    """Density of H0 samples on the grid (Gaussian KDE, reflected at the prior bounds)."""
    from scipy.stats import gaussian_kde

    lo, hi = bounds
    kde = gaussian_kde(samples)
    p = kde(h0) + kde(2 * lo - h0) + kde(2 * hi - h0)
    return np.where((h0 >= lo) & (h0 <= hi), p, 0.0)


def _grid_density(h0: np.ndarray, x: np.ndarray, p: np.ndarray) -> np.ndarray:
    return np.interp(h0, x, p, left=0.0, right=0.0)


def read_input(path: str | Path, h0: np.ndarray) -> dict:
    """One input of the joint posterior: dict(name, method, kind, events, density on `h0`)."""
    path = Path(path).expanduser()
    if path.is_dir() and (path / BRIGHT_FILE).exists():
        info = json.loads((path / BRIGHT_FILE).read_text())
        if tuple(info.get("h0_range", H0_PRIOR)) != H0_PRIOR:
            raise ValueError(f"{path}: the bright siren used the H0 range {info['h0_range']}, not the common prior "
                             f"{H0_PRIOR[0]:g}-{H0_PRIOR[1]:g}")
        g = pd.read_csv(path / GRID_FILE, sep="\t")
        va = info.get("viewing_angle")
        kind = f"bright siren ({info['event']}, {info['label']}" + (f", viewing angle {va[0]:g} ± {va[1]:g}°" if va
                                                                    else "") + ")"
        return dict(name=path.name, method="bright", kind=kind, events={info["event"]} | set(info.get("events", [])),
                    density=_grid_density(h0, g["H0"].to_numpy(float), g["p"].to_numpy(float)))
    if path.is_dir():
        post = next((path / n for n in ICAROGW_POSTERIORS if (path / n).exists()), None)
        if post is None:
            raise ValueError(f"{path}: not a hubble_constant work directory (no {BRIGHT_FILE}, "
                             f"{' or '.join(ICAROGW_POSTERIORS)})")
        kind = siren_kind(path)
        summ = path / "summary.json"
        model = json.loads(summ.read_text()).get("mass_model") if summ.exists() else None
        events = set()
        if (path / "events.tsv").exists():
            ev = pd.read_csv(path / "events.tsv", sep="\t")
            events = set(ev["event"]) | set(ev.get("common_name", pd.Series(dtype=str)).dropna())
        return dict(name=path.name, method="dark" if kind.startswith("dark") else "spectral",
                    kind=kind + (f", {model}" if model else "") + ("" if post.name == ICAROGW_POSTERIORS[0] else
                                                                  ", not reweighted"),
                    events=events, density=sample_density(h0, pd.read_csv(post, sep="\t")["H0"].to_numpy(float)))
    df = pd.read_csv(path, sep="\t")
    if "H0" not in df:
        raise ValueError(f"{path} has no H0 column")
    if "p" in df:
        return dict(name=path.name, method="file", kind=f"posterior grid {path.name}", events=set(),
                    density=_grid_density(h0, df["H0"].to_numpy(float), df["p"].to_numpy(float)))
    return dict(name=path.name, method="file", kind=f"posterior samples {path.name}", events=set(),
                density=sample_density(h0, df["H0"].to_numpy(float)))


def check_independent(inputs: list[dict]) -> None:
    """At most one spectral or dark input, and no event in two inputs."""
    icarogw = [i["name"] for i in inputs if i["method"] in ("spectral", "dark")]
    if len(icarogw) > 1:
        raise ValueError(f"{', '.join(icarogw)}: at most one spectral or dark siren in a joint posterior (the dark "
                         "siren already contains the spectral siren of the same events)")
    for a in range(len(inputs)):
        for b in range(a + 1, len(inputs)):
            common = inputs[a]["events"] & inputs[b]["events"]
            if common:
                raise ValueError(f"{', '.join(sorted(common))} in both {inputs[a]['name']} and {inputs[b]['name']}: "
                                 "the inputs are not independent (leave the event out of the spectral or dark siren "
                                 "with --exclude)")


def run_h0_joint(
    inputs: Iterable[str | Path],
    workdir: str | Path = "hubble_constant_joint",
    out_report_html: Optional[str | Path] = "hubble_constant_joint.html",
    out_summary_tsv: Optional[str | Path] = "hubble_constant_joint.tsv",
) -> pd.DataFrame:
    """Joint H0 posterior of independent hubble_constant results, with the table of each one and of their product."""
    inputs = list(inputs)
    if len(inputs) < 2:
        raise ValueError("the joint posterior needs at least two inputs (--inputs DIR DIR ...)")
    lo, hi = H0_PRIOR
    h0 = np.linspace(lo, hi, int(round((hi - lo) / 0.05)) + 1)
    parts = [read_input(p, h0) for p in inputs]
    check_independent(parts)
    joint = np.prod([_normalize(h0, i["density"]) for i in parts], axis=0)
    if not np.any(joint > 0):
        raise ValueError("the posteriors do not overlap: their product is zero on the H0 grid")
    joint = _normalize(h0, joint)
    rows = [dict(analysis=i["kind"], input=i["name"], n_events=len(i["events"]) or None,
                 **summarize(h0, i["density"])) for i in parts]
    rows.append(dict(analysis="joint", input=" × ".join(i["name"] for i in parts),
                     n_events=sum(len(i["events"]) for i in parts) or None, **summarize(h0, joint)))
    for r in rows:
        _log(f"{r['analysis']}: H0 = {_fmt(r)}")
    table = pd.DataFrame(rows)

    workdir = Path(workdir).expanduser()
    workdir.mkdir(parents=True, exist_ok=True)
    grid = pd.DataFrame({"H0": h0, "p": joint, **{f"p_{i['name']}": _normalize(h0, i["density"]) for i in parts}})
    grid.to_csv(workdir / "posterior_joint.tsv", sep="\t", index=False, float_format="%.6g")
    (workdir / "joint.json").write_text(json.dumps(dict(method="joint", inputs=[str(Path(p).expanduser().resolve())
                                                                                 for p in inputs]), indent=1))
    if out_summary_tsv:
        Path(out_summary_tsv).parent.mkdir(parents=True, exist_ok=True)
        table.to_csv(out_summary_tsv, sep="\t", index=False, float_format="%.6g")
    if out_report_html:
        colors = ("#2a78d6", "#1a9e77", "#7a5195", "#52514e")
        curves = [(f"{r['analysis'].split(' (')[0]} ({i['name']}): {r['map']:.0f}", i["density"], colors[k % 4],
                   ("--", ":", "-.")[k % 3]) for k, (r, i) in enumerate(zip(rows, parts))]
        curves.append((f"joint: {rows[-1]['map']:.0f}, 68% {rows[-1]['hpd68_low']:.0f}–{rows[-1]['hpd68_high']:.0f}",
                       joint, "#eb6834", "-"))
        images = [_plot(h0, curves, None, "", "Joint Hubble constant", workdir / "plots" / "h0_joint.png")]
        paras = [f"H<sub>0</sub> = <b>{_fmt(rows[-1])}</b>, the product of {len(parts)} independent posteriors with "
                 f"the same flat prior ({lo:g}–{hi:g} km/s/Mpc): "
                 + "; ".join(f"{i['kind']} ({i['name']})" for i in parts) + ".",
                 "Spectral- and dark-siren samples are turned into a density by a kernel estimate reflected at the "
                 "prior bounds; bright-siren posteriors are already on a grid. The inputs are checked to share no "
                 "event, and at most one of them is a spectral or dark siren (the dark siren contains the spectral "
                 "siren of its events)."]
        Path(out_report_html).parent.mkdir(parents=True, exist_ok=True)
        write_simple_html_report(out_report_html, title="Hubble constant (joint)", paragraphs=paras, images=images,
                                 tables=[("H0 by method and joint", table.to_html(
                                     index=False, float_format=lambda x: f"{x:.3g}", na_rep=""))])
        _log(f"report written to {out_report_html}")
    return table
