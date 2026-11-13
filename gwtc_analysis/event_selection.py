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


DETECTOR_NAMES = {"H": "H1", "L": "L1", "V": "V1", "K": "K1", "G": "G1"}
DETECTOR_MATCH = ("all", "exact", "any")


def parse_detectors(tokens) -> list[str]:
    """Detector names from what the user wrote: `H L V`, `H1 L1`, `HLV`, `H,L` -> ["H1", "L1", "V1"] (in the order
    H1, L1, V1, K1, G1, without repeats)."""
    found = set()
    for tok in ([tokens] if isinstance(tokens, str) else tokens):
        letters = str(tok).upper().replace("1", "").replace(",", "").replace("-", "").replace(" ", "")
        if not letters or any(c not in DETECTOR_NAMES for c in letters):
            raise ValueError(f"unknown detector {tok!r}: use H, L, V (also K, G), e.g. H L V or HLV")
        found.update(DETECTOR_NAMES[c] for c in letters)
    return [d for d in DETECTOR_NAMES.values() if d in found]


def detectors_mask(networks: pd.Series, detectors: list[str], match: str = "all") -> pd.Series:
    """Which networks (lists of detector names, or None) satisfy the selection: `all` of the detectors are in the
    network, the network is `exact`ly these detectors, or `any` of them is in it. No network: False."""
    if match not in DETECTOR_MATCH:
        raise ValueError(f"detectors_match must be one of {', '.join(DETECTOR_MATCH)}, got {match!r}")
    want = set(detectors)
    test = {"all": want.issubset, "exact": want.__eq__, "any": lambda s: bool(want & s)}[match]
    return networks.apply(lambda x: isinstance(x, (list, tuple)) and len(x) > 0 and bool(test(set(x))))


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
                     + [v if isinstance(v, str) else f"{k.replace('_min', ' ≥ ').replace('_max', ' ≤ ')}{v:g}"
                        for k, v in bounds.items() if v is not None])
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
    area_max: Optional[float] = None,
    area_cred: float = 0.9,
    data_repo: str = "zenodo",
    zenodo_versions: Optional[dict[str, str]] = None,
    detectors: Optional[list[str]] = None,
    detectors_match: str = "all",
    q_min: Optional[float] = None,
    q_max: Optional[float] = None,
    snr_min: Optional[float] = None,
    p_astro_min: Optional[float] = None,
    far_max: Optional[float] = None,
    runs: Optional[list[str]] = None,
) -> None:
    """Select GWTC events based on source-frame component masses and luminosity distance.

    Uses:
      - mass_1_source
      - mass_2_source
      - luminosity_distance

    `preset` selects a class of sources (`source_classes.apply_preset`): neutron-stars, mass-gap or
    hierarchical; its columns are added to the TSV. The bounds apply on top of it.

    `area_max` keeps the events whose sky area at the `area_cred` credible level (deg², column A90_deg2 for 0.9)
    is at most this, computed from the skymaps (`data_repo`, `zenodo_versions` as catalog_statistics) of the events
    passing the other cuts only; events without a skymap are dropped.

    `q_min`, `q_max` bound the mass ratio q = m2 / m1 of the GWOSC medians (after the m1 >= m2 ordering; column
    mass_ratio); `snr_min` the network matched-filter SNR (column network_snr); `p_astro_min` the probability of
    astrophysical origin and `far_max` the false-alarm rate (per year) of GWOSC, which lists for each event the
    values of the pipeline that ranks it best (columns p_astro, far_per_yr); `runs` keeps the events of these
    observing runs (O1 ... O4b, from the GPS time; column run). An event without the value a cut needs is dropped.

    `detectors` (names as `parse_detectors` gives them) keeps the events by their detector network, the one GWOSC
    lists for the event (one query per event passing the other cuts): with `detectors_match` "all" the network
    holds all of them, "exact" it is these detectors and no other, "any" it holds at least one. The network is
    added to the TSV (column detectors); events without one are dropped.

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
    num_cols = ["mass_1_source", "mass_2_source", "luminosity_distance", "redshift", "chi_eff", "snr", "far", "p_astro",
                "gps"]
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
            df_all.loc[swap, [a, b]] = df_all.loc[swap, [b, a]].to_numpy(copy=True)   # a copy: pandas 3 can give a view

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

    df_all["mass_ratio"] = df_all["mass_2_source"] / df_all["mass_1_source"]
    df_all["network_snr"], df_all["far_per_yr"] = df_all["snr"], df_all["far"]
    for value, col, keep in ((q_min, "mass_ratio", "ge"), (q_max, "mass_ratio", "le"), (snr_min, "network_snr", "ge"),
                             (p_astro_min, "p_astro", "ge"), (far_max, "far_per_yr", "le")):
        if value is not None:
            mask &= (df_all[col] >= float(value)) if keep == "ge" else (df_all[col] <= float(value))
    if runs:
        from .catalog_registry import observing_runs, run_of

        known = list(observing_runs())
        unknown = [r for r in runs if r not in known]
        if unknown:
            raise ValueError(f"unknown observing run(s) {', '.join(unknown)}; choose from {', '.join(known)}")
        df_all["run"] = df_all["gps"].apply(lambda g: run_of(g) if np.isfinite(g) else None)
        mask &= df_all["run"].isin(list(runs))

    extra, regions = pd.DataFrame(index=df_all.index), {}
    if preset:
        pmask, extra, desc = sc.apply_preset(df_all, preset, ns_max_mass=ns_max_mass, mass_gap=mass_gap,
                                             pisn_gap_min=pisn_gap_min)
        mask &= pmask.fillna(False)
        regions = {"neutron-stars": {"ns_max": ns_max_mass}, "mass-gap": {"gap": mass_gap},
                   "hierarchical": {"pisn": pisn_gap_min}}[preset]
        print(f"[event_selection] preset {preset}: {desc}")

    if detectors:
        detectors = parse_detectors(detectors)
        df_all["detectors"] = None
        nets = pd.Series([None] * len(df_all), index=df_all.index, dtype=object)
        if mask.any():
            print(f"[event_selection] detector networks of {int(mask.sum())} event(s) from GWOSC")
            sub, _ = gw.add_detectors_and_virgo_flag(df_all.loc[mask], progress=True, plot_network_pie=False)
            for i, v in zip(sub.index, sub["detectors"]):
                nets[i] = v
            missing = int(sum(not isinstance(v, list) for v in nets[mask]))
            if missing:
                print(f"[event_selection] {missing} event(s) without a detector network at GWOSC: not selected")
        df_all["detectors"] = nets.apply(lambda x: ",".join(x) if isinstance(x, list) else "")
        mask &= detectors_mask(nets, detectors, detectors_match)

    area_col = f"A{int(round(100 * area_cred))}_deg2"
    if area_max is not None:
        if not 0 < area_cred < 1:
            raise ValueError(f"area_cred must be in (0, 1), got {area_cred}")
        from .catalogs import add_credible_areas

        df_all[area_col] = np.nan
        if mask.any():
            sub = df_all.loc[mask, ["event_id", "catalog_key"]].copy()
            print(f"[event_selection] sky areas of {len(sub)} event(s) from {data_repo}")
            add_credible_areas(sub, sorted(sub["catalog_key"].unique()), data_repo=data_repo, cred=area_cred,
                               column=area_col, zenodo_versions=zenodo_versions, tag="event_selection")
            df_all.loc[mask, area_col] = sub[area_col].to_numpy()
            missing = int(df_all.loc[mask, area_col].isna().sum())
            if missing:
                print(f"[event_selection] {missing} event(s) without a skymap: not selected")
        mask &= df_all[area_col] <= float(area_max)

    cols = ["event_id", "catalog_key", "mass_1_source", "mass_2_source", "luminosity_distance", "redshift"]
    if preset == "hierarchical" or chi_eff_min is not None or chi_eff_max is not None:
        cols.append("chi_eff")
    cols += [c for c, used in (("mass_ratio", q_min is not None or q_max is not None),
                               ("network_snr", snr_min is not None), ("p_astro", p_astro_min is not None),
                               ("far_per_yr", far_max is not None), ("run", bool(runs))) if used]
    if detectors:
        cols.append("detectors")
    if area_max is not None:
        cols.append(area_col)
    out = pd.concat([df_all.loc[mask, cols], extra.loc[mask, [c for c in extra.columns if c not in cols]]], axis=1)

    # Stable order for tests/users
    out = out.sort_values(["catalog_key", "event_id"]).reset_index(drop=True)

    out.to_csv(out_tsv, sep="\t", index=False)

    if out_plot:
        bounds = dict(m1_min=m1_min, m1_max=m1_max, m2_min=m2_min, m2_max=m2_max, dl_min=dl_min, dl_max=dl_max,
                      chi_eff_min=chi_eff_min, chi_eff_max=chi_eff_max)
        if area_max is not None:
            bounds[area_col.replace("_deg2", "_max")] = area_max   # in the title only: "A90 ≤ 100"
        bounds.update(q_min=q_min, q_max=q_max, SNR_min=snr_min, p_astro_min=p_astro_min, FAR_max=far_max)
        if runs:
            bounds["runs"] = f"runs {' '.join(runs)}"
        if detectors:                                              # in the title only: "detectors H1 L1 (all)"
            bounds["detectors"] = f"detectors {' '.join(detectors)} ({detectors_match})"
        png = _plot_selection(df_all, mask, bounds, requested, Path(out_plot), preset=preset, regions=regions)
        print(f"[event_selection] {int(mask.sum())} of {len(df_all)} events selected; plot written to {png}")
