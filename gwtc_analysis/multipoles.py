"""Higher multipoles and precession of one event, from the SNRs stored in its PE samples.

The GWTC-4.0 and later PE files store, for each sample, the network SNR of the (3,3), (4,4) and (2,1)
multipoles orthogonal to the (2,2) one (Mills & Fairhurst 2021, arXiv:2007.04313), and the precession SNR
rho_p, the SNR of the weaker of the two harmonics of a precessing signal (Fairhurst et al. 2020,
arXiv:1908.05707). Without the multipole (or precession), rho^2 is chi^2 distributed with 2 degrees of
freedom: rho follows a Rayleigh distribution, P(rho > x) = exp(-x^2 / 2), 11% at 2.1 and 1% at 3.

The posterior medians are compared with these thresholds: >= 3 clear evidence, 2.1-3 a hint. This is the
noise-only scale; the comparison with the prior-only distribution of rho made in the LVK papers needs prior
samples that the files do not contain. Earlier catalogs (GWTC-1 to GWTC-3) do not store these SNRs.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any, Optional

import numpy as np
import pandas as pd

FIELDS = {
    "network_33_multipole_snr": "ρ₃₃",
    "network_44_multipole_snr": "ρ₄₄",
    "network_21_multipole_snr": "ρ₂₁",
    "network_precessing_snr": "ρ_p",
}
_PLOT_LABELS = {"ρ₃₃": "ρ$_{33}$", "ρ₄₄": "ρ$_{44}$", "ρ₂₁": "ρ$_{21}$", "ρ_p": "ρ$_{\\rm p}$"}
HINT, CLEAR = 2.1, 3.0


def evidence(median: float) -> str:
    if not np.isfinite(median):
        return ""
    return "clear" if median >= CLEAR else "hint" if median >= HINT else "none"


def summarize(samples_dict: Any) -> pd.DataFrame:
    """One row per (label, quantity) with the SNR fields: median, 90% interval, P(ρ > 3), evidence."""
    rows = []
    for label in samples_dict.keys():
        s = samples_dict[label]
        keys = set(s.keys()) if hasattr(s, "keys") else set()
        for field, name in FIELDS.items():
            if field not in keys:
                continue
            x = np.asarray(s[field], dtype=float)
            x = x[np.isfinite(x)]
            if not len(x):
                continue
            med = float(np.median(x))
            rows.append(dict(label=label, quantity=name, field=field, median=med,
                             low_90=float(np.percentile(x, 5)), high_90=float(np.percentile(x, 95)),
                             p_noise_at_median=float(np.exp(-med ** 2 / 2)), frac_above_3=float(np.mean(x > CLEAR)),
                             evidence=evidence(med)))
    return pd.DataFrame(rows, columns=["label", "quantity", "field", "median", "low_90", "high_90",
                                       "p_noise_at_median", "frac_above_3", "evidence"])


def plot(samples: Any, label: str, src_name: str, out_png: Path) -> Optional[Path]:
    """Posterior of each SNR of `label`, with the noise-only (Rayleigh) distribution."""
    present = [f for f in FIELDS if f in set(samples.keys())]
    if not present:
        return None
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    surface, ink, ink2, grid, blue = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df", "#2a78d6"
    fig, axes = plt.subplots(1, len(present), figsize=(3.2 * len(present), 3.2), dpi=150, squeeze=False)
    fig.patch.set_facecolor(surface)
    for ax, field in zip(axes[0], present):
        x = np.asarray(samples[field], dtype=float)
        x = x[np.isfinite(x)]
        hi = max(5.0, float(np.percentile(x, 99.5)) * 1.1)
        bins = np.linspace(0, hi, 41)
        ax.set_facecolor(surface)
        ax.hist(x, bins=bins, density=True, color=blue, alpha=0.85, lw=0, zorder=2, label="posterior")
        r = np.linspace(0, hi, 300)
        ax.plot(r, r * np.exp(-r ** 2 / 2), color=ink2, lw=1.2, ls=(0, (3, 2)), zorder=3, label="noise only")
        for v in (HINT, CLEAR):
            ax.axvline(v, color=ink2, lw=0.7, ls=":", zorder=1)
        med = float(np.median(x))
        lab = _PLOT_LABELS[FIELDS[field]]
        ax.set_title(f"{lab}: median {med:.1f} ({evidence(med)})", color=ink, fontsize=9.5, loc="left")
        ax.set_xlabel(lab, color=ink2)
        ax.set_yticks([])
        ax.grid(axis="x", color=grid, lw=0.7, zorder=0)
        for sp in ("top", "right", "left"):
            ax.spines[sp].set_visible(False)
        ax.spines["bottom"].set_color(grid)
        ax.tick_params(colors=ink2, labelsize=8)
    axes[0][0].legend(frameon=False, fontsize=7.5, labelcolor=ink2)
    fig.suptitle(f"Higher multipoles and precession — {src_name} ({label}); dotted: 2.1 and 3", color=ink,
                 fontsize=10, x=0.01, ha="left")
    fig.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, facecolor=surface)
    plt.close(fig)
    return out_png


def report_html(table: pd.DataFrame, label: str) -> str:
    """HTML section: the selected label's SNRs, and the other labels for comparison."""
    import html

    if table.empty:
        return ("<h2>Higher multipoles and precession</h2><p>The PE file stores no multipole or precession SNRs "
                "(they are stored from GWTC-4.0 on).</p>")
    fmt = lambda x: f"{x:.2g}" if abs(x) < 0.01 else f"{x:.2f}"
    main = table[table["label"] == label]
    lines = ["<h2>Higher multipoles and precession</h2>",
             "<p>SNR of the (3,3), (4,4) and (2,1) multipoles beyond the (2,2) one, and precession SNR ρ<sub>p</sub>, "
             "from the PE samples. Without the effect, ρ follows a Rayleigh distribution: P(ρ &gt; 2.1) = 11%, "
             "P(ρ &gt; 3) = 1%. Medians ≥ 3: clear evidence; 2.1–3: a hint. This is the noise-only scale, not the "
             "comparison with the prior made in the LVK papers.</p>"]
    if not main.empty:
        clear = main[main["evidence"] == "clear"]["quantity"].tolist()
        hint = main[main["evidence"] == "hint"]["quantity"].tolist()
        verdict = (f"clear evidence for {', '.join(clear)}" if clear else "no clear evidence") + \
                  (f"; a hint of {', '.join(hint)}" if hint else "")
        lines.append(f"<p><b>{html.escape(label)}: {html.escape(verdict)}.</b></p>")
    cols = ["label", "quantity", "median", "low_90", "high_90", "p_noise_at_median", "frac_above_3", "evidence"]
    lines.append(table[cols].to_html(index=False, escape=True, float_format=fmt))
    if table["label"].nunique() > 1:
        spread = table.groupby("quantity")["median"].agg(["min", "max"])
        wide = spread[(spread["max"] - spread["min"]) > 1.0]
        if not wide.empty:
            lines.append("<p>The waveform models disagree (medians differing by more than 1): "
                         + html.escape(", ".join(f"{q} {r['min']:.1f}–{r['max']:.1f}" for q, r in wide.iterrows()))
                         + ".</p>")
    return "\n".join(lines)


def analyse(samples_dict: Any, label: str, src_name: str, outdir: Path) -> tuple[pd.DataFrame, list[Path], str]:
    """Summary table (also written as TSV), plot of `label`, and the report section."""
    table = summarize(samples_dict)
    files: list[Path] = []
    if not table.empty:
        tsv = Path(outdir) / f"{src_name}_multipoles_precession.tsv"
        table.to_csv(tsv, sep="\t", index=False, float_format="%.4g")
        files.append(tsv)
        png = plot(samples_dict[label], label, src_name, Path(outdir) / f"multipoles_precession_{src_name}.png")
        if png:
            files.append(png)
    return table, files, report_html(table, label)
