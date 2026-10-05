"""Remnant and energetics of one event, from the PE samples: final mass and spin, radiated energy, peak luminosity.

The PE files of GWTC-2.1 and later store, for each sample, the source-frame final mass `final_mass_source`,
the final spin `final_spin`, the energy radiated in gravitational waves `radiated_energy` (M_sun c^2, source
frame) and the peak luminosity `peak_luminosity` (10^56 erg/s). They are computed from the component
masses and spins with fits to numerical-relativity simulations, not measured separately: the test of the
area law with an independent remnant is the `area_law` mode (GW250114).
"""
from __future__ import annotations

from pathlib import Path
from typing import Any, Optional

import numpy as np
import pandas as pd

MSUN_C2_ERG = 1.7877e54
FIELDS = {                                   # field: (name, unit)
    "final_mass_source": ("final mass", "M☉"),
    "final_spin": ("final spin", ""),
    "radiated_energy": ("radiated energy", "M☉c²"),
    "peak_luminosity": ("peak luminosity", "10⁵⁶ erg/s"),
}


def _col(s: Any, k: str) -> Optional[np.ndarray]:
    if k not in set(s.keys()):
        return None
    x = np.asarray(s[k], dtype=float)
    return x[np.isfinite(x)]


def summarize(samples_dict: Any) -> pd.DataFrame:
    """One row per (label, quantity) with the remnant fields, plus the radiated fraction of the total mass."""
    rows = []
    for label in samples_dict.keys():
        s = samples_dict[label]
        if not hasattr(s, "keys"):
            continue
        qs = {f: _col(s, f) for f in FIELDS}
        m1, m2, e = _col(s, "mass_1_source"), _col(s, "mass_2_source"), qs["radiated_energy"]
        if e is not None and m1 is not None and m2 is not None and len(e) == len(m1) == len(m2):
            qs["radiated_fraction"] = e / (m1 + m2)
        for f, x in qs.items():
            if x is None or not len(x):
                continue
            name, unit = FIELDS.get(f, ("radiated fraction", ""))
            rows.append(dict(label=label, quantity=name, unit=unit, median=float(np.median(x)),
                             low_90=float(np.percentile(x, 5)), high_90=float(np.percentile(x, 95))))
    return pd.DataFrame(rows, columns=["label", "quantity", "unit", "median", "low_90", "high_90"])


def plot(samples: Any, label: str, src_name: str, out_png: Path) -> Optional[Path]:
    present = [f for f in FIELDS if _col(samples, f) is not None]
    if not present:
        return None
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    surface, ink, ink2, grid, blue = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df", "#2a78d6"
    fig, axes = plt.subplots(1, len(present), figsize=(3.2 * len(present), 3.2), dpi=150, squeeze=False)
    fig.patch.set_facecolor(surface)
    for ax, f in zip(axes[0], present):
        x = _col(samples, f)
        name, unit = FIELDS[f]
        ax.set_facecolor(surface)
        ax.hist(x, bins=40, density=True, color=blue, alpha=0.85, lw=0, zorder=2)
        q = np.percentile(x, [5, 50, 95])
        for v, ls in zip(q, (":", "-", ":")):
            ax.axvline(v, color=ink2, lw=0.9, ls=ls, zorder=3)
        if f == "final_spin":
            ax.axvline(0.686, color=ink2, lw=0.8, ls=(0, (3, 2)), zorder=3)
            ax.annotate("0.686", xy=(0.686, 1), xycoords=("data", "axes fraction"), xytext=(2, -10),
                        textcoords="offset points", fontsize=7, color=ink2)
        ax.set_title(f"{name}: {q[1]:.3g} (+{q[2] - q[1]:.2g} / −{q[1] - q[0]:.2g})", color=ink, fontsize=9.5, loc="left")
        tex = {"M☉": "M$_\\odot$", "M☉c²": "M$_\\odot$c$^2$", "10⁵⁶ erg/s": "10$^{56}$ erg/s"}.get(unit, unit)
        ax.set_xlabel(f"{name}" + (f" ({tex})" if tex else ""), color=ink2)
        ax.set_yticks([])
        ax.grid(axis="x", color=grid, lw=0.7, zorder=0)
        for sp in ("top", "right", "left"):
            ax.spines[sp].set_visible(False)
        ax.spines["bottom"].set_color(grid)
        ax.tick_params(colors=ink2, labelsize=8)
    fig.suptitle(f"Remnant and energetics — {src_name} ({label}); median and 90% interval", color=ink, fontsize=10,
                 x=0.01, ha="left")
    fig.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, facecolor=surface)
    plt.close(fig)
    return out_png


def report_html(table: pd.DataFrame, label: str) -> str:
    import html

    if table.empty:
        return ("<h2>Remnant and energetics</h2><p>The PE file stores no remnant quantities (final mass and spin, "
                "radiated energy, peak luminosity).</p>")
    main = table[table["label"] == label]
    lines = ["<h2>Remnant and energetics</h2>",
             "<p>Final mass and spin, energy radiated in gravitational waves and peak luminosity, from the PE samples "
             "(source frame). They come from fits to numerical-relativity simulations applied to the component masses "
             "and spins, not from a separate measurement of the remnant; the area-law test with an independent "
             "ringdown remnant is the <code>area_law</code> mode.</p>"]
    if not main.empty:
        get = lambda q: main[main["quantity"] == q]
        parts = []
        for q, fmt in (("final mass", "{:.1f} M☉"), ("final spin", "{:.2f}"), ("radiated energy", "{:.2f} M☉c²"),
                       ("peak luminosity", "{:.2f} × 10⁵⁶ erg/s")):
            r = get(q)
            if not r.empty:
                parts.append(f"{q} {fmt.format(r['median'].iloc[0])}")
        e = get("radiated energy")
        if not e.empty:
            parts.append(f"i.e. {e['median'].iloc[0] * MSUN_C2_ERG:.2g} erg")
        lines.append(f"<p><b>{html.escape(label)}: {html.escape('; '.join(parts))}.</b></p>")
    lines.append(table.to_html(index=False, escape=True, float_format=lambda x: f"{x:.3g}"))
    return "\n".join(lines)


def analyse(samples_dict: Any, label: str, src_name: str, outdir: Path) -> tuple[pd.DataFrame, list[Path], str]:
    table = summarize(samples_dict)
    files: list[Path] = []
    if not table.empty:
        tsv = Path(outdir) / f"{src_name}_remnant.tsv"
        table.to_csv(tsv, sep="\t", index=False, float_format="%.4g")
        files.append(tsv)
        png = plot(samples_dict[label], label, src_name, Path(outdir) / f"remnant_{src_name}.png")
        if png:
            files.append(png)
    return table, files, report_html(table, label)
