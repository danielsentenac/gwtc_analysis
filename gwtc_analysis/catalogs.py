from __future__ import annotations

import os
from pathlib import Path
from typing import List, Optional

import numpy as np
import pandas as pd
import difflib
import textwrap


from . import gw_stat as gw
from .report import write_simple_html_report


ALLOWED_CATALOGS = [
    "GWTC-1",
    "GWTC-2.1",
    "GWTC-3",
    "GWTC-4",
    "ALL",
]

CATALOG_STATISTICS_ALIASES = {
    "GWTC-4": "GWTC-4.0",
    "GWTC-3": "GWTC-3-confident",
    "GWTC-2.1": "GWTC-2.1-confident",
    "GWTC-1": "GWTC-1-confident",
}

def _pick_first_existing_col(df: pd.DataFrame, candidates: list[str]) -> str | None:
    for c in candidates:
        if c in df.columns:
            return c
    return None


def _plot_m1_m2_snr_scatter(
    df: pd.DataFrame,
    out_png: Path,
    *,
    catalogs_label: str = "",
) -> Path | None:
    """
    Scatter plot of m1 vs m2 colored by network SNR.
    Returns out_png if created, else None.
    """
    plt = _ensure_matplotlib()

    if "mass_1_source" not in df.columns or "mass_2_source" not in df.columns:
        return None

    m1 = pd.to_numeric(df["mass_1_source"], errors="coerce")
    m2 = pd.to_numeric(df["mass_2_source"], errors="coerce")

    snr_col = _pick_first_existing_col(
        df,
        ["network_snr", "snr_network", "network_matched_filter_snr", "snr"],
    )
    if snr_col is None:
        return None

    snr = pd.to_numeric(df[snr_col], errors="coerce")

    mask = np.isfinite(m1) & np.isfinite(m2) & np.isfinite(snr)
    if not bool(mask.any()):
        return None

    m1v = m1[mask].to_numpy()
    m2v = m2[mask].to_numpy()
    snrv = snr[mask].to_numpy()

    fig, ax = plt.subplots(figsize=(12, 6))
    sc = ax.scatter(m1v, m2v, c=snrv, s=45, alpha=0.9)

    ax.set_xlabel(r"Mass 1 (M$_\odot$)", fontsize=14)
    ax.set_ylabel(r"Mass 2 (M$_\odot$)", fontsize=14)
    ax.set_title(f"{catalogs_label}\nsource masses distribution", fontsize=16)

    ax.grid(True, alpha=0.25)
    ax.axhline(0, linewidth=2)
    ax.axvline(0, linewidth=2)

    cb = fig.colorbar(sc, ax=ax)
    cb.ax.set_title("Network SNR", fontsize=12, pad=10)

    fig.tight_layout()
    fig.savefig(out_png, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out_png


def _plot_histograms_panel(
    df: pd.DataFrame,
    out_png: Path,
    *,
    catalogs_label: str = "",
) -> Path | None:
    """
    Three-panel histogram figure:
      - total mass (m1+m2)
      - luminosity distance
      - network SNR
    Returns out_png if created, else None.
    """
    plt = _ensure_matplotlib()

    if "mass_1_source" not in df.columns or "mass_2_source" not in df.columns:
        return None

    m1 = pd.to_numeric(df["mass_1_source"], errors="coerce")
    m2 = pd.to_numeric(df["mass_2_source"], errors="coerce")
    mtot = (m1 + m2)

    dist_col = _pick_first_existing_col(
        df,
        ["luminosity_distance", "luminosity_distance_mpc", "distance", "dist_mpc"],
    )
    if dist_col is None:
        return None
    dl = pd.to_numeric(df[dist_col], errors="coerce")

    snr_col = _pick_first_existing_col(
        df,
        ["network_snr", "snr_network", "network_matched_filter_snr", "snr"],
    )
    if snr_col is None:
        return None
    snr = pd.to_numeric(df[snr_col], errors="coerce")

    mtot = mtot[np.isfinite(mtot)]
    dl = dl[np.isfinite(dl)]
    snr = snr[np.isfinite(snr)]

    # If everything is empty, skip
    if len(mtot) == 0 and len(dl) == 0 and len(snr) == 0:
        return None

    # White background / black foreground for readability
    plt.rcParams.update(
        {
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.facecolor": "white",
            "text.color": "black",
            "axes.labelcolor": "black",
            "axes.titlecolor": "black",
            "xtick.color": "black",
            "ytick.color": "black",
            "axes.edgecolor": "black",
        }
    )

    fig = plt.figure(figsize=(17, 8))
    gs = fig.add_gridspec(2, 3, height_ratios=[1.0, 1.0], wspace=0.55, hspace=0.35)

    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[0, 2])

    tx1 = fig.add_subplot(gs[1, 0]); tx1.axis("off")
    tx2 = fig.add_subplot(gs[1, 1]); tx2.axis("off")
    tx3 = fig.add_subplot(gs[1, 2]); tx3.axis("off")

    hist_kw = dict(edgecolor="black", linewidth=1.2)

    if len(mtot) > 0:
        ax1.hist(mtot, bins=7, **hist_kw)
    ax1.set_title("Total Mass Histogram", fontsize=18)
    ax1.set_xlabel(r"Mass (M$_\odot$)", fontsize=14)
    ax1.set_ylabel("Count", fontsize=14)
    ax1.tick_params(axis="both", labelsize=12)
    ax1.grid(True, axis="y", alpha=0.25)

    if len(dl) > 0:
        ax2.hist(dl, bins=8, **hist_kw)
    ax2.set_title("Luminosity Distance Histogram", fontsize=18)
    ax2.set_xlabel("Distance (Mpc)", fontsize=14)
    ax2.set_ylabel("Count", fontsize=14)
    ax2.tick_params(axis="both", labelsize=12)
    ax2.grid(True, axis="y", alpha=0.25)

    if len(snr) > 0:
        ax3.hist(snr, bins=10, **hist_kw)
    ax3.set_title("Network SNR Histogram", fontsize=18)
    ax3.set_xlabel("SNR", fontsize=14)
    ax3.set_ylabel("Count", fontsize=14)
    ax3.tick_params(axis="both", labelsize=12)
    ax3.grid(True, axis="y", alpha=0.25)

    tx1.text(
        -0.15, 0.85,
        textwrap.fill(f"Distribution of total mass\nfor events contained in:\n{catalogs_label}", width=38),
        va="top", fontsize=12, color="black"
    )
    tx2.text(
        -0.15, 0.85,
        textwrap.fill(f"Distribution of luminosity distance\nfor events contained in:\n{catalogs_label}", width=38),
        va="top", fontsize=12, color="black"
    )
    tx3.text(
        -0.15, 0.85,
        textwrap.fill(f"Distribution of network SNR\nfor events contained in:\n{catalogs_label}", width=38),
        va="top", fontsize=12, color="black"
    )
    fig.savefig(out_png, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out_png

def _plot_source_type_pie(df: pd.DataFrame, out_png: Path, column: str = "binary_type") -> Path | None:
    s = df.get(column)
    if s is None:
        return None
    s = s.dropna()
    if s.empty:
        return None

    counts = s.value_counts()
    if counts.sum() > 0:
        colors = {
            "BBH": "#9467bd",
            "BH-NS": "#ff7f0e",
            "NS-NS": "#2ca02c",
        }
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(5, 4))
    wedges, _, autotexts = ax.pie(
        counts.values,
        labels=None,
        colors=[colors[k] for k in counts.index],
        startangle=90,
        counterclock=False,
        wedgeprops=dict(width=0.42),
        autopct=lambda p: f"{p:.1f}%" if p >= 3 else "",
        pctdistance=0.78,
    )

    # Center label
    total = int(counts.sum())
    ax.text(
        0, 0,
        f"N = {total}",
        ha="center", va="center",
        fontsize=13,
        fontweight="bold",
    )

    # Legend
    ax.legend(
        wedges,
        [f"{k} ({counts[k]})" for k in counts.index],
        title="Source type",
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        frameon=False,
    )
    ax.set_title("Source type distribution")
    fig.savefig(out_png, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out_png


def _format_allowed_catalogs() -> str:
    return ", ".join(ALLOWED_CATALOGS)
    
def _validate_catalogs(catalogs: list[str]) -> None:
    bad = [c for c in catalogs if c not in ALLOWED_CATALOGS]
    if not bad:
        return

    suggestions = []
    for b in bad:
        m = difflib.get_close_matches(b, ALLOWED_CATALOGS, n=1, cutoff=0.6)
        if m:
            suggestions.append(f"{b} → did you mean {m[0]}?")

    msg = (
        "Unknown catalog(s): " + ", ".join(bad)
        + ". Allowed catalogs are: " + ", ".join(ALLOWED_CATALOGS)
    )
    if suggestions:
        msg += ". Suggestions: " + "; ".join(suggestions)

    raise ValueError(msg)

def _safe_mkdir(p: str | Path) -> None:
    Path(p).mkdir(parents=True, exist_ok=True)

def _ensure_matplotlib():
    import matplotlib
    matplotlib.use("Agg", force=True)  # headless-safe
    import matplotlib.pyplot as plt
    return plt
def _plot_network_counts(df: pd.DataFrame, out_png: str | Path) -> Optional[str]:
    plt = _ensure_matplotlib()
    if "n_det" not in df.columns:
        return None
    vc = df["n_det"].value_counts().sort_index()
    if vc.empty:
        return None
    fig, ax = plt.subplots(figsize=(6,4))
    ax.bar(vc.index.astype(int).astype(str), vc.values)
    ax.set_xlabel("Number of detectors")
    ax.set_ylabel("Count")
    ax.set_title("Detector count distribution")
    fig.tight_layout()
    fig.savefig(out_png, dpi=150)
    plt.close(fig)
    return str(out_png)

def _plot_area_cdf(
    df: pd.DataFrame,
    out_png: Path,
    *,
    column: str,
    catalog_label: str | None = None,
    source_type: str | None = None,   # e.g. "BBH"
    from_zenodo: bool = False,
) -> Path | None:
    import numpy as np
    import pandas as pd
    plt = _ensure_matplotlib()
    
    if column not in df.columns:
        return None

    d = df.copy()

    # Optional filter by source type (MMODA example does BBH only)
    if source_type is not None and "binary_type" in d.columns:
        d = d[d["binary_type"] == source_type]

    # Need detector count for MMODA-style grouping
    if "n_det" not in d.columns:
        # fall back to single curve
        s = pd.to_numeric(d[column], errors="coerce").dropna()
        if s.empty:
            return None
        xs = np.sort(s.values)
        ys = np.arange(1, len(xs) + 1) / len(xs)
        fig, ax = plt.subplots()
        ax.step(xs, ys, where="post")
        ax.set_xscale("log")
        ax.set_xlabel(f"{column.replace('_', ' ')} [deg$^2$]")
        ax.set_ylabel("Cumulative fraction")
        ax.set_title("Sky localization (CDF)")
        fig.savefig(out_png, dpi=150, bbox_inches="tight")
        plt.close(fig)
        return out_png

    fig, ax = plt.subplots()

    # Exclusive groups: by number of detectors ONLY (matches MMODA)
    det_counts = sorted(pd.to_numeric(d["n_det"], errors="coerce").dropna().unique())
    plotted = 0

    for n in det_counts:
        n = int(n)
        g = d[d["n_det"] == n]
        s = pd.to_numeric(g[column], errors="coerce").dropna()
        if len(s) < 2:
            continue

        xs = np.sort(s.values)
        ys = np.arange(1, len(xs) + 1) / len(xs)

        q05, q50, q95 = np.percentile(xs, [5, 50, 95])
        label = f"{n} detector{'s' if n != 1 else ''} (N={len(xs)}; med={q50:.0f}, 5–95%={q05:.0f}–{q95:.0f})"

        ax.step(xs, ys, where="post", label=label)
        plotted += 1

    if plotted == 0:
        return None

    ax.set_xscale("log")
    ax.set_xlabel(rf"$A_{{{int(round(100*0.9))}}}$  [deg$^2$]" if "A90" in column else f"{column} [deg$^2$]")
    ax.set_ylabel("Cumulative fraction")

    # MMODA-like multi-line title
    title_lines = []
    if source_type is not None:
        title_lines.append(f"{source_type} sky localization")
    else:
        title_lines.append("Sky localization")

    if catalog_label:
        title_lines.append(f"[{catalog_label!r}]")

    if from_zenodo:
        title_lines.append(f"{column.split('_')[0]} computed from Zenodo PE skymaps")

    ax.set_title("\n".join(title_lines))
    ax.legend(loc="upper left", frameon=True, fontsize=8)
    ax.grid(True, which="both", linestyle="--", alpha=0.5)
    fig.savefig(out_png, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out_png



def run_catalog_statistics(
    catalogs: List[str],
    out_events_tsv: str | Path,
    out_report_html: Optional[str | Path] = None,
    include_detectors: bool = True,
    include_area: bool = False,
    area_cred: float = 0.9,
    area_column: str | None = None,
    ns_threshold: float = 3.0,
    data_repo: str = "s3",
    plots_dir: Optional[str] = None,
) -> None:
    """
    Fetch events from GWOSC jsonfull for one or more catalogs and compute basic derived columns.

    Notes
    -----
    - Catalog names passed on the CLI are "user-facing" (GWTC-1, GWTC-2.1, GWTC-3, GWTC-4, ALL).
      Internally, some of them are aliased to the GWOSC/GWTC identifiers needed by APIs.
    - For localization areas (--include-area), supported data repos are:
        * s3
        * zenodo
        * galaxy   (expects Galaxy to have staged per-catalog collections under the working directory)

    Outputs
    -------
    out_events_tsv : TSV
        Per-event table with masses, distance, detector network (optional), and credible area (optional).
    out_report_html : HTML (optional)
        Single-file HTML report with summary tables and embedded plots.
    """

    if plots_dir is None:
        plots_dir = "cat_plots"

    plot_dir = Path(plots_dir)
    plot_dir.mkdir(parents=True, exist_ok=True)

    # Expand ALL selector
    catalogs = list(catalogs or [])
    if "ALL" in catalogs:
        catalogs = [c for c in ALLOWED_CATALOGS if c != "ALL"]

    if not catalogs:
        raise ValueError("No catalogs selected")

    _validate_catalogs(catalogs)

    if data_repo not in {"s3", "zenodo", "galaxy"}:
        raise ValueError(f"Unsupported data_repo={data_repo!r}. Use one of: s3, zenodo, galaxy")

    if area_column is None:
        level = int(round(100 * area_cred))
        area_column = f"A{level}_deg2"

    dfs: list[pd.DataFrame] = []

    for cat in catalogs:
        resolved_cat = CATALOG_STATISTICS_ALIASES.get(cat, cat)
        if resolved_cat != cat:
            print(f"Catalog alias applied for statistics: {cat} → {resolved_cat}")

        raw = gw.fetch_gwtc_events(catalog=resolved_cat)
        df_cat = gw.events_to_dataframe(raw["events"])
        df_cat["catalog_key"] = cat   # keep user-facing key stable
        dfs.append(df_cat)

    # Avoid pd.concat warning (future behavior change): build records with a union schema.
    kept: list[pd.DataFrame] = []
    for d in dfs:
        if d is None or d.empty:
            continue
        if not d.notna().to_numpy().any():
            continue
        kept.append(d)

    if not kept:
        raise RuntimeError("No valid catalog dataframes to combine")

    all_cols = sorted(set().union(*(d.columns for d in kept)))

    records: list[dict] = []
    for d in kept:
        d2 = d.reindex(columns=all_cols)
        records.extend(d2.to_dict(orient="records"))

    df0 = pd.DataFrame.from_records(records, columns=all_cols)
    df = gw.prepare_catalog_df(df0, ns_threshold=ns_threshold)

    # ------------------------------------------------------------------
    # Detectors network (requires GWOSC v2 calls)
    # ------------------------------------------------------------------
    fig_network = None
    if include_detectors:
        df, fig_network = gw.add_detectors_and_virgo_flag(
            df, progress=True, verbose=False, plot_network_pie=True
        )
        df["n_det"] = df["detectors"].apply(lambda x: len(x) if isinstance(x, list) else np.nan)
    else:
        df["detectors"] = np.nan
        df["n_det"] = np.nan
        df["has_V1"] = np.nan

    # ------------------------------------------------------------------
    # Credible area (optional)
    # ------------------------------------------------------------------
    if include_area:
        df[area_column] = np.nan

        zenodo_cache_dir = ".cache_gwosc"

        for cat in catalogs:
            m = df["catalog_key"].eq(cat)
            if not m.any():
                continue

            if data_repo == "zenodo":
                tmp = gw.add_localization_area_from_zenodo(
                    df.loc[m],
                    catalog_key=cat,
                    cred=area_cred,
                    column=area_column,
                    cache_dir=zenodo_cache_dir,
                    progress=True,
                    verbose=False,
                )
                df.loc[m, area_column] = tmp[area_column].values
                continue

            if data_repo == "s3":
                tmp = gw.add_localization_area_from_s3(
                    df.loc[m],
                    catalog_key=cat,
                    cred=area_cred,
                    column=area_column,
                    bucket="gwtc",
                    base_prefix="",
                    progress=True,
                    verbose=False,
                )
                df.loc[m, area_column] = tmp[area_column].values
                continue

            if data_repo == "galaxy":
                # Convention: Galaxy collections are staged under the working directory as:
                #   GWTC-2.1-SKYMAPS, GWTC-3-SKYMAPS, GWTC-4-SKYMAPS
                base = Path(f"{cat}-SKYMAPS")
                if not base.exists():
                    # GWTC-1 is covered by GWTC-2.1 on Galaxy
                    if cat == "GWTC-1":
                        base = Path("GWTC-2.1-SKYMAPS")
                    if not base.exists():
                        raise RuntimeError(
                            f"Galaxy mode: cannot find skymaps directory for {cat}. "
                            f"Expected {cat}-SKYMAPS (or GWTC-2.1-SKYMAPS for GWTC-1). "
                            "Make sure the Galaxy collection is staged in the working directory."
                        )

                if hasattr(gw, "add_localization_area_from_galaxy"):
                    tmp = gw.add_localization_area_from_galaxy(
                        df.loc[m],
                        catalog_key=cat,
                        skymap_dir=str(base),
                        cred=area_cred,
                        column=area_column,
                        progress=True,
                        verbose=False,
                    )
                else:
                    # Backward-compatible fallback: the directory-based helper already supports Galaxy collections.
                    tmp = gw.add_localization_area_from_directory(
                        df.loc[m],
                        skymap_dir=str(base),
                        cred=area_cred,
                        column=area_column,
                        progress=True,
                        verbose=False,
                    )

                df.loc[m, area_column] = tmp[area_column].values
                continue

            raise ValueError(
                f"include_area requires a supported data_repo. Got {data_repo!r} (expected s3/zenodo/galaxy)."
            )

    # ------------------------------------------------------------------
    # Write per-event table
    # ------------------------------------------------------------------
    out_events_tsv = Path(out_events_tsv)
    _safe_mkdir(out_events_tsv.parent if out_events_tsv.parent != Path("") else ".")
    df_out = df.copy()

    if "detectors" in df_out.columns:
        df_out["detectors"] = df_out["detectors"].apply(lambda x: ",".join(x) if isinstance(x, list) else "")

    keep = [
        "event_id", "catalog_key", "version",
        "mass_1_source", "mass_2_source", "chirp_mass_source", "total_mass_source", "final_mass_source",
        "luminosity_distance", "redshift", "chi_eff", "chi_p", "snr", "far", "p_astro",
        "binary_type", "detectors", "n_det", "has_V1",
    ]
    if include_area:
        keep.append(area_column)
    keep = [c for c in keep if c in df_out.columns]
    df_out[keep].to_csv(out_events_tsv, sep="\t", index=False)

    # ------------------------------------------------------------------
    # Report (optional)
    # ------------------------------------------------------------------
    if not out_report_html:
        return

    out_report_html = Path(out_report_html)
    _safe_mkdir(out_report_html.parent if out_report_html.parent != Path("") else ".")

    def _rel_to_html(p: Path) -> Path:
        rel = os.path.relpath(str(p), start=str(out_report_html.parent))
        return Path(rel)

    n_total = len(df_out)
    per_cat = df_out["catalog_key"].value_counts().rename_axis("catalog").reset_index(name="N")
    per_type = df_out["binary_type"].value_counts().rename_axis("binary_type").reset_index(name="N")
    per_net = df_out["n_det"].value_counts(dropna=True).sort_index().rename_axis("n_det").reset_index(name="N")

    tables = []
    tables.append(("Counts by catalog", per_cat.to_html(index=False, escape=True)))
    tables.append(("Counts by binary type", per_type.to_html(index=False, escape=True)))
    if include_detectors:
        tables.append(("Counts by detector number", per_net.to_html(index=False, escape=True)))

    img_paths: list[Path] = []
    cat_label = ", ".join(catalogs)
    cat_tag = "_".join(catalogs)

    if include_detectors:
        if fig_network is None:
            raise RuntimeError(
                "include_detectors=True but fig_network is None. "
                "Call gw.add_detectors_and_virgo_flag(..., plot_network_pie=True)."
            )
        p_net = plot_dir / "network_pie.png"
        fig_network.savefig(p_net, dpi=150, bbox_inches="tight")
        img_paths.append(_rel_to_html(p_net))

    p_type = _plot_source_type_pie(
        df_out,
        plot_dir / "source_types_pie.png",
        column="binary_type",
    )
    if p_type:
        img_paths.append(_rel_to_html(p_type))

    p_m1m2 = _plot_m1_m2_snr_scatter(
        df_out,
        plot_dir / f"{cat_tag}_m1_m2_snr.png",
        catalogs_label=cat_label,
    )
    if p_m1m2:
        img_paths.append(_rel_to_html(p_m1m2))

    p_hists = _plot_histograms_panel(
        df_out,
        plot_dir / f"{cat_tag}_histograms.png",
        catalogs_label=cat_label,
    )
    if p_hists:
        img_paths.append(_rel_to_html(p_hists))

    if include_area and area_column in df_out.columns:
        p_area_all = _plot_area_cdf(
            df_out,
            plot_dir / f"{area_column}_cdf.png",
            column=area_column,
            catalog_label=cat_label,
            source_type=None,
            from_zenodo=(data_repo == "zenodo"),
        )
        if p_area_all:
            img_paths.append(_rel_to_html(p_area_all))

        p_area_bbh = _plot_area_cdf(
            df_out,
            plot_dir / f"{area_column}_cdf_bbh.png",
            column=area_column,
            catalog_label=cat_label,
            source_type="BBH",
            from_zenodo=(data_repo == "zenodo"),
        )
        if p_area_bbh:
            img_paths.append(_rel_to_html(p_area_bbh))

    paragraphs = [
        f"Catalogs: {', '.join(catalogs)}",
        f"Total events after basic cleaning: {n_total}",
        f"Per-event table written to: {out_events_tsv.name}",
    ]

    if include_area and area_column in df_out.columns:
        got = int(pd.notna(df_out[area_column]).sum())
        paragraphs.append(f"Credible area computed at cred={area_cred}: {got}/{n_total} events.")

        vals = pd.to_numeric(df_out[area_column], errors="coerce").dropna()
        vals = vals[(vals > 0) & np.isfinite(vals)]
        if not vals.empty:
            p10, p50, p90 = np.percentile(vals.values, [10, 50, 90])
            paragraphs.append(
                f"{area_column}: N={len(vals)}, p10={p10:.1f} deg², median={p50:.1f} deg², p90={p90:.1f} deg²"
            )

    write_simple_html_report(
        out_report_html,
        title="GWTC catalog statistics",
        paragraphs=paragraphs,
        images=img_paths,
        tables=tables,
    )
