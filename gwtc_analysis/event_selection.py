from __future__ import annotations

from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

from . import gw_stat as gw
from . import source_classes as sc

from .catalog_registry import gwosc_aliases

# Catalog key -> GWOSC jsonfull list of its confident events (catalog_registry)
CATALOG_ALIASES: dict[str, str] = gwosc_aliases()


def _as_float_or_nan(x) -> float:
    try:
        return float(x)
    except Exception:
        return float("nan")


def _plot_selection(df: pd.DataFrame, mask: pd.Series, bounds: dict, catalogs: list[str], out_png: Path,
                    preset: Optional[str] = None, regions: Optional[dict] = None) -> Path:
    """The selected events among all the events of the catalogs: m2 against m1, and D_L (chi_eff for the
    hierarchical preset) against m1; the preset's regions are shaded."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    surface, ink, ink2, grid, grey, orange = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df", "#b8b6b0", "#eb6834"
    regions = regions or {}
    second = (("chi_eff", "Effective spin $\\chi_{\\rm eff}$", (bounds.get("chi_eff_min"), bounds.get("chi_eff_max")))
              if preset == "hierarchical" else
              ("luminosity_distance", "Luminosity distance $D_L$ (Mpc)", (bounds.get("dl_min"), bounds.get("dl_max"))))
    ok = df["mass_1_source"].notna() & df["mass_2_source"].notna() & df[second[0]].notna()
    d, sel = df[ok], mask[ok]
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.4), dpi=150)
    fig.patch.set_facecolor(surface)
    for ax, y, ylabel, (lo, hi) in ((axes[0], "mass_2_source", "Secondary mass $m_2$ (source frame, M$_\\odot$)",
                                     (bounds.get("m2_min"), bounds.get("m2_max"))),
                                    (axes[1],) + second):
        ax.set_facecolor(surface)
        ax.scatter(d.loc[~sel, "mass_1_source"], d.loc[~sel, y], s=16, color=grey, lw=0, zorder=2,
                   label=f"not selected ({int((~sel).sum())})")
        ax.scatter(d.loc[sel, "mass_1_source"], d.loc[sel, y], s=34, color=orange, edgecolor=surface, lw=1.2,
                   zorder=3, label=f"selected ({int(sel.sum())})")
        ax.set_xscale("log")
        if y != "chi_eff":
            ax.set_yscale("log")
        for v in (bounds.get("m1_min"), bounds.get("m1_max")):
            if v is not None:
                ax.axvline(v, color=ink2, lw=0.9, ls=(0, (4, 3)), zorder=1)
        for v in (lo, hi):
            if v is not None:
                ax.axhline(v, color=ink2, lw=0.9, ls=(0, (4, 3)), zorder=1)
        ax.set_xlabel("Primary mass $m_1$ (source frame, M$_\\odot$)", color=ink2)
        ax.set_ylabel(ylabel, color=ink2)
        ax.grid(color=grid, lw=0.7, zorder=0)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        for sp in ("left", "bottom"):
            ax.spines[sp].set_color(grid)
        ax.tick_params(colors=ink2, labelsize=9)
    shade = dict(color=orange, alpha=0.10, lw=0, zorder=0)
    lims = [(ax.get_xlim(), ax.get_ylim()) for ax in axes]
    if "gap" in regions:
        axes[0].axvspan(*regions["gap"], **shade); axes[0].axhspan(*regions["gap"], **shade)
        axes[1].axvspan(*regions["gap"], **shade)
    if "ns_max" in regions:
        axes[0].axhspan(lims[0][1][0], regions["ns_max"], **shade)
    if "pisn" in regions:
        for ax, (xl, _) in zip(axes, lims):
            ax.axvspan(regions["pisn"], xl[1], **shade)
        axes[1].axhspan(lims[1][1][0], 0, color=orange, alpha=0.05, lw=0, zorder=0)
    for ax, (xl, yl) in zip(axes, lims):
        ax.set_xlim(xl); ax.set_ylim(yl)
    axes[0].legend(frameon=False, fontsize=9, labelcolor=ink2, loc="upper left")
    crit = ", ".join(([f"preset {preset}"] if preset else [])
                     + [f"{k.replace('_min', ' ≥ ').replace('_max', ' ≤ ')}{v:g}" for k, v in bounds.items() if v is not None])
    fig.suptitle(f"Event selection in {', '.join(catalogs)}" + (f": {crit}" if crit else ""), color=ink, fontsize=11,
                 x=0.01, ha="left")
    fig.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, facecolor=surface)
    plt.close(fig)
    return out_png


def run_event_selection(
    *,
    catalogs: list[str],
    out_tsv: str | Path,
    m1_min: Optional[float] = None,
    m1_max: Optional[float] = None,
    m2_min: Optional[float] = None,
    m2_max: Optional[float] = None,
    dl_min: Optional[float] = None,
    dl_max: Optional[float] = None,
    chi_eff_min: Optional[float] = None,
    chi_eff_max: Optional[float] = None,
    preset: Optional[str] = None,
    ns_max_mass: float = sc.NS_MAX_MASS,
    mass_gap: tuple[float, float] = sc.MASS_GAP,
    pisn_gap_min: float = sc.PISN_GAP_MIN,
    out_plot: Optional[str | Path] = None,
) -> None:
    """Select GWTC events based on source-frame component masses and luminosity distance.

    Uses:
      - mass_1_source
      - mass_2_source
      - luminosity_distance

    `preset` selects a class of sources (`source_classes.apply_preset`): neutron-stars, mass-gap or
    hierarchical; its columns are added to the TSV. The bounds apply on top of it.

    Writes TSV with selected events (at least event_id), and with `out_plot` a PNG of the selected events
    among all the events of the catalogs. The redshift column is GWOSC's: inferred from the luminosity
    distance for the Planck 2015 cosmology, it is the one the source-frame masses were derived with,
    m_src = m_det / (1 + z).
    """

    out_tsv = Path(out_tsv)
    out_tsv.parent.mkdir(parents=True, exist_ok=True)

    # Expand ALL catalog selector
    requested = list(catalogs)
    if "ALL" in catalogs:
        from .catalog_registry import expand_all

        catalogs = expand_all(catalogs)          # the default catalogs, without the updates (GWTC-4.1)

    # Fetch per-catalog event tables
    dfs: list[pd.DataFrame] = []
    for cat in catalogs:
        resolved = CATALOG_ALIASES.get(cat, cat)
        if resolved != cat:
            print(f"[event_selection] Catalog alias applied: {cat} → {resolved}")

        raw = gw.fetch_gwtc_events(catalog=resolved)
        df_cat = gw.events_to_dataframe(raw["events"])
        df_cat["catalog_key"] = cat  # keep user-facing key stable
        dfs.append(df_cat)

    if not dfs:
        pd.DataFrame(columns=["event_id", "catalog_key"]).to_csv(out_tsv, sep="\t", index=False)
        return

    # ---- Combine without pd.concat (future-proof for pandas dtype warnings) ----
    kept: list[pd.DataFrame] = []
    for d in dfs:
        if d is None or d.empty:
            continue
        if not d.notna().to_numpy().any():
            continue
        kept.append(d)

    if not kept:
        pd.DataFrame(columns=["event_id", "catalog_key"]).to_csv(out_tsv, sep="\t", index=False)
        return

    all_cols = sorted(set().union(*(d.columns for d in kept)))

    records: list[dict] = []
    for d in kept:
        d2 = d.reindex(columns=all_cols)
        records.extend(d2.to_dict(orient="records"))

    df_all = pd.DataFrame.from_records(records, columns=all_cols)

    # Ensure numeric
    if preset is not None and preset not in sc.PRESETS:
        raise ValueError(f"unknown preset {preset!r}; choose from {', '.join(sc.PRESETS)}")
    num_cols = ["mass_1_source", "mass_2_source", "luminosity_distance", "redshift", "chi_eff"]
    num_cols += [c for c in df_all.columns if c.endswith(("_lo90", "_hi90"))]
    for col in num_cols:
        if col in df_all.columns:
            df_all[col] = df_all[col].apply(_as_float_or_nan)
        else:
            df_all[col] = np.nan

    # the primary is the heavier component (GWOSC lists a few events the other way round)
    swap = df_all["mass_2_source"] > df_all["mass_1_source"]
    for a, b in (("mass_1_source", "mass_2_source"), ("mass_1_source_lo90", "mass_2_source_lo90"),
                 ("mass_1_source_hi90", "mass_2_source_hi90")):
        if a in df_all and b in df_all and swap.any():
            df_all.loc[swap, [a, b]] = df_all.loc[swap, [b, a]].to_numpy()

    # Build mask
    mask = pd.Series(True, index=df_all.index)

    if m1_min is not None:
        mask &= df_all["mass_1_source"] >= float(m1_min)
    if m1_max is not None:
        mask &= df_all["mass_1_source"] <= float(m1_max)

    if m2_min is not None:
        mask &= df_all["mass_2_source"] >= float(m2_min)
    if m2_max is not None:
        mask &= df_all["mass_2_source"] <= float(m2_max)

    if dl_min is not None:
        mask &= df_all["luminosity_distance"] >= float(dl_min)
    if dl_max is not None:
        mask &= df_all["luminosity_distance"] <= float(dl_max)
    if chi_eff_min is not None:
        mask &= df_all["chi_eff"] >= float(chi_eff_min)
    if chi_eff_max is not None:
        mask &= df_all["chi_eff"] <= float(chi_eff_max)

    extra, regions = pd.DataFrame(index=df_all.index), {}
    if preset:
        pmask, extra, desc = sc.apply_preset(df_all, preset, ns_max_mass=ns_max_mass, mass_gap=mass_gap,
                                             pisn_gap_min=pisn_gap_min)
        mask &= pmask.fillna(False)
        regions = {"neutron-stars": {"ns_max": ns_max_mass}, "mass-gap": {"gap": mass_gap},
                   "hierarchical": {"pisn": pisn_gap_min}}[preset]
        print(f"[event_selection] preset {preset}: {desc}")

    cols = ["event_id", "catalog_key", "mass_1_source", "mass_2_source", "luminosity_distance", "redshift"]
    if preset == "hierarchical" or chi_eff_min is not None or chi_eff_max is not None:
        cols.append("chi_eff")
    out = pd.concat([df_all.loc[mask, cols], extra.loc[mask, [c for c in extra.columns if c not in cols]]], axis=1)

    # Stable order for tests/users
    out = out.sort_values(["catalog_key", "event_id"]).reset_index(drop=True)

    out.to_csv(out_tsv, sep="\t", index=False)

    if out_plot:
        bounds = dict(m1_min=m1_min, m1_max=m1_max, m2_min=m2_min, m2_max=m2_max, dl_min=dl_min, dl_max=dl_max,
                      chi_eff_min=chi_eff_min, chi_eff_max=chi_eff_max)
        png = _plot_selection(df_all, mask, bounds, requested, Path(out_plot), preset=preset, regions=regions)
        print(f"[event_selection] {int(mask.sum())} of {len(df_all)} events selected; plot written to {png}")
