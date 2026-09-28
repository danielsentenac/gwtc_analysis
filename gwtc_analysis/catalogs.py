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
    "GWTC-5",
    "ALL",
]

CATALOG_STATISTICS_ALIASES = {
    "GWTC-5": "GWTC-5.0",
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
    zenodo_versions: Optional[dict[str, str]] = None,
) -> None:
    """
    Fetch events from GWOSC jsonfull for one or more catalogs and compute basic derived columns.

    `zenodo_versions` maps catalog keys to the Zenodo release version read with
    data_repo="zenodo" (e.g. {"GWTC-3": "v2"}); other catalogs use the latest.

    Notes
    -----
    - Catalog names passed on the CLI are "user-facing" (GWTC-1, GWTC-2.1, GWTC-3, GWTC-4, GWTC-5, ALL).
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

    # prepare_catalog_df drops events without source-frame component masses
    # (they cannot be placed on mass-based statistics). Report how many were
    # dropped per catalog so the kept count is not surprising.
    mass_drop_rows: list[dict] = []
    if "catalog_key" in df0.columns:
        before = df0["catalog_key"].value_counts()
        after = df["catalog_key"].value_counts() if "catalog_key" in df.columns else {}
        for cat in catalogs:
            n_before = int(before.get(cat, 0))
            n_after = int(after.get(cat, 0))
            n_drop = n_before - n_after
            mass_drop_rows.append(
                {"catalog": cat, "total": n_before, "kept": n_after, "dropped_no_masses": n_drop}
            )
            if n_drop > 0:
                print(
                    f"[catalogs] {cat}: {n_after}/{n_before} events kept "
                    f"({n_drop} dropped: no source-frame masses in GWOSC metadata)"
                )
            else:
                print(f"[catalogs] {cat}: {n_after}/{n_before} events kept")

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
                    zenodo_version=(zenodo_versions or {}).get(cat),
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
                #   GWTC-2.1-SKYMAPS, GWTC-3-SKYMAPS, GWTC-4-SKYMAPS, GWTC-5-SKYMAPS
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
        # Images are embedded as data URIs: the report needs a path readable from here,
        # not one relative to the HTML file (that broke reports written outside the cwd).
        return Path(p)

    n_total = len(df_out)
    per_cat = df_out["catalog_key"].value_counts().rename_axis("catalog").reset_index(name="N")
    per_type = df_out["binary_type"].value_counts().rename_axis("binary_type").reset_index(name="N")
    per_net = df_out["n_det"].value_counts(dropna=True).sort_index().rename_axis("n_det").reset_index(name="N")

    tables = []
    tables.append(("Counts by catalog", per_cat.to_html(index=False, escape=True)))
    if mass_drop_rows:
        drop_df = pd.DataFrame(mass_drop_rows, columns=["catalog", "total", "kept", "dropped_no_masses"])
        tables.append((
            "Events kept vs. dropped (events without source-frame masses are excluded "
            "from mass-based statistics)",
            drop_df.to_html(index=False, escape=True),
        ))
    tables.append(("Counts by binary type", per_type.to_html(index=False, escape=True)))
    if include_detectors:
        tables.append(("Counts by detector number", per_net.to_html(index=False, escape=True)))

    # Loudest events: top 10% by network SNR — the best candidates to inspect
    # the parameter estimation for. List their names so they can be fed to
    # `parameters_estimation --src-name <event>`.
    top_snr_names: list[str] = []
    top_snr_threshold: float | None = None
    snr_col_rep = _pick_first_existing_col(
        df_out, ["network_snr", "snr_network", "network_matched_filter_snr", "snr"]
    )
    name_col = "common_name" if "common_name" in df_out.columns else "event_id"
    if snr_col_rep is not None:
        snr_num = pd.to_numeric(df_out[snr_col_rep], errors="coerce")
        valid = df_out.assign(_snr=snr_num).dropna(subset=["_snr"])
        if not valid.empty:
            top_snr_threshold = float(valid["_snr"].quantile(0.90))
            top = valid[valid["_snr"] >= top_snr_threshold].sort_values("_snr", ascending=False)
            top_snr_names = [str(x) for x in top[name_col].tolist()]
            cols = [c for c in (name_col, "_snr", "mass_1_source", "mass_2_source", "catalog_key") if c in top.columns]
            top_tbl = top[cols].rename(columns={name_col: "event", "_snr": "network_snr"})
            tables.append((
                f"Loudest events — top 10% by network SNR (≥ {top_snr_threshold:.1f}), "
                "best candidates for PE inspection",
                top_tbl.to_html(index=False, escape=True),
            ))

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

    total_drop = sum(r["dropped_no_masses"] for r in mass_drop_rows)
    if total_drop > 0:
        per_cat_drop = ", ".join(
            f"{r['catalog']} {r['kept']}/{r['total']}" for r in mass_drop_rows
        )
        paragraphs.append(
            f"{total_drop} catalog event(s) were dropped from the statistics because the "
            f"GWOSC metadata has no source-frame component masses ({per_cat_drop} kept). "
            "See the 'Events kept vs. dropped' table for the per-catalog breakdown."
        )

    if top_snr_names:
        paragraphs.append(
            f"Loudest events (top 10% by network SNR ≥ {top_snr_threshold:.1f}, "
            f"n={len(top_snr_names)}) — best candidates for PE inspection: "
            + ", ".join(top_snr_names)
        )

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


# ---------------------------------------------------------------------------
# Merger rates: R = N / <VT> from the LVK search-sensitivity injections
# ---------------------------------------------------------------------------
#
# The sensitive volume-time <VT> of a population is estimated from an LVK
# injection campaign (simulated signals added to the real data and searched by
# the real pipelines), reweighted from the injected distribution to the
# population by importance sampling:
#
#   <VT> = T * sum_found[ w * p_pop(m1, m2, spins, z) / p_draw ] / N_generated
#   p_pop(z) ∝ (1+z)^(kappa-1) dVc/dz
#
# N counts the GWOSC candidates inside the injections' time span that pass the
# same FAR threshold (confident and marginal lists), classified by their median
# source-frame masses. The rate posterior uses a Jeffreys prior: Gamma(N+1/2)/<VT>.

# GWTC-4.0 cumulative search sensitivity, real O3 + O4a injections (zenodo 16740128).
RATES_SENSITIVITY_RECORD = "16740128"
RATES_SENSITIVITY_FILE = "mixture-real_o3_o4a-cartesian_spins_20250503134659UTC.hdf"
RATES_EVENT_LISTS = (
    "GWTC-2.1-confident", "GWTC-2.1-marginal",
    "GWTC-3-confident", "GWTC-3-marginal",
    "GWTC-4.0", "GWTC-5.0",
)
_YEAR_S = 3.15576e7


def _rates_sensitivity_path(sensitivity_file: str | Path | None) -> Path:
    """Local injection file: the given path, or the default release (downloaded once)."""
    if sensitivity_file:
        p = Path(sensitivity_file).expanduser()
        if not p.exists():
            raise ValueError(f"Sensitivity file not found: {p}")
        return p
    from .data_repo import zenodo_cache_dir

    dest = zenodo_cache_dir() / f"zenodo_{RATES_SENSITIVITY_RECORD}_{RATES_SENSITIVITY_FILE}"
    if not (dest.exists() and dest.stat().st_size > 0):
        url = f"https://zenodo.org/records/{RATES_SENSITIVITY_RECORD}/files/{RATES_SENSITIVITY_FILE}?download=1"
        print(f"[rates] downloading LVK search-sensitivity injections (~400 MB): {url}")
        gw._download_with_byte_progress(url, dest)
    return dest


def _injection_segments(times: np.ndarray, max_gap_days: float = 7.0) -> list[tuple[float, float]]:
    """Observing periods covered by the injections: their times split at gaps longer than max_gap_days."""
    t = np.sort(np.asarray(times, dtype=float))
    cut = np.where(np.diff(t) > max_gap_days * 86400.0)[0]
    starts = np.concatenate([[t[0]], t[cut + 1]])
    ends = np.concatenate([t[cut], [t[-1]]])
    return list(zip(starts.tolist(), ends.tolist()))


def _load_found_injections(path: Path, far_threshold: float) -> dict:
    """Injections found by any search with FAR below `far_threshold` [1/yr]."""
    import h5py

    lnpdraw_key = "lnpdraw_mass1_source_mass2_source_redshift_spin1x_spin1y_spin1z_spin2x_spin2y_spin2z"
    with h5py.File(path, "r") as h:
        searches = [s.decode() if isinstance(s, bytes) else str(s) for s in h.attrs["searches"]]
        ev = h["events"]
        far = np.min([ev[f"{s}_far"][:] for s in searches], axis=0)
        found = far < far_threshold
        times = ev["time_geocenter"][:]
        out = {k: ev[k][:][found] for k in ("mass1_source", "mass2_source", "redshift", "weights",
                                             "spin1x", "spin1y", "spin1z", "spin2x", "spin2y", "spin2z")}
        out["lnpdraw"] = ev[lnpdraw_key][:][found]
        out["T_yr"] = float(h.attrs["total_analysis_time"]) / _YEAR_S
        out["N_gen"] = float(h.attrs["total_generated"])
        out["segments"] = _injection_segments(times)
    out["a1"] = np.sqrt(out["spin1x"] ** 2 + out["spin1y"] ** 2 + out["spin1z"] ** 2)
    out["a2"] = np.sqrt(out["spin2x"] ** 2 + out["spin2y"] ** 2 + out["spin2z"] ** 2)
    from astropy.cosmology import Planck15

    zg = np.linspace(0.0, float(out["redshift"].max()) * 1.01 + 1e-3, 4000)
    dvdz = 4 * np.pi * Planck15.differential_comoving_volume(zg).to("Gpc3/sr").value
    out["dvdz"] = np.interp(out["redshift"], zg, dvdz)
    return out


def _ln_iso_spin(a: np.ndarray, amax: float) -> np.ndarray:
    """Isotropic spin with magnitude uniform in [0, amax], as a density in Cartesian components."""
    with np.errstate(divide="ignore"):
        return np.where(a < amax, -np.log(4 * np.pi * np.maximum(a, 1e-12) ** 2 * amax), -np.inf)


def _ln_uniform(x: np.ndarray, lo: float, hi: float) -> np.ndarray:
    return np.where((x >= lo) & (x <= hi), -np.log(hi - lo), -np.inf)


def _ln_powerlaw(x: np.ndarray, alpha: float, lo: float, hi: float) -> np.ndarray:
    """p(x) ∝ x^-alpha on [lo, hi]."""
    norm = (hi ** (1 - alpha) - lo ** (1 - alpha)) / (1 - alpha)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where((x >= lo) & (x <= hi), -alpha * np.log(np.maximum(x, 1e-12)) - np.log(norm), -np.inf)


def _ln_power_law_peak(m1: np.ndarray, m2: np.ndarray, *, alpha=3.4, mmin=5.1, mmax=87.0, lam=0.04,
                       mu=34.0, sig=3.6, beta=1.1, delta=4.8) -> np.ndarray:
    """GWTC-3 'Power Law + Peak' BBH mass model with low-mass smoothing, normalized numerically."""
    trapz = getattr(np, "trapezoid", None) or np.trapz

    def smooth(m):
        m = np.asarray(m, dtype=float)
        out = np.zeros_like(m)
        out[m >= mmin + delta] = 1.0
        mid = (m > mmin) & (m < mmin + delta)
        x = m[mid] - mmin
        with np.errstate(over="ignore"):
            out[mid] = 1.0 / (np.exp(delta / x + delta / (x - delta)) + 1.0)
        return out

    g = np.linspace(mmin, 300.0, 60000)
    pl_norm = trapz(np.where(g <= mmax, g ** -alpha, 0.0), g)

    def shape1(m):
        pl = np.where((m >= mmin) & (m <= mmax), np.maximum(m, 1e-12) ** -alpha, 0.0) / pl_norm
        peak = np.exp(-0.5 * ((m - mu) / sig) ** 2) / (sig * np.sqrt(2 * np.pi))
        return ((1 - lam) * pl + lam * peak) * smooth(m)

    norm1 = trapz(shape1(g), g)
    q = g ** beta * smooth(g)
    cq = np.concatenate([[0.0], np.cumsum(0.5 * (q[1:] + q[:-1]) * np.diff(g))])   # ∫_mmin^m m2^beta S dm2
    with np.errstate(divide="ignore", invalid="ignore"):
        lp = (np.log(shape1(m1) / norm1)
              + np.log(np.where(m2 <= m1, np.maximum(m2, 1e-12) ** beta * smooth(m2), 0.0) / np.interp(m1, g, cq)))
    return np.where(np.isfinite(lp), lp, -np.inf)


def _sensitive_vt(inj: dict, ln_mass: np.ndarray, amax1: float, amax2: float, kappa: float = 0.0) -> tuple[float, float]:
    """(<VT> [Gpc^3 yr], effective sample size) of a population over the found injections."""
    lnp = (ln_mass + _ln_iso_spin(inj["a1"], amax1) + _ln_iso_spin(inj["a2"], amax2)
           + np.log(inj["dvdz"]) + (kappa - 1.0) * np.log1p(inj["redshift"]))
    with np.errstate(over="ignore", invalid="ignore"):
        x = inj["weights"] * np.exp(lnp - inj["lnpdraw"])
    x = np.where(np.isfinite(x), x, 0.0)
    s = float(x.sum())
    return inj["T_yr"] * s / inj["N_gen"], (s * s / float((x * x).sum()) if s > 0 else 0.0)


def _rate_quantiles(n: int, vt: float, levels=(0.05, 0.5, 0.95)) -> np.ndarray:
    """Rate posterior quantiles for n detections and <VT>, Jeffreys prior: Gamma(n + 1/2) / VT."""
    from scipy.stats import gamma

    return gamma.ppf(np.asarray(levels), n + 0.5) / vt


def _rates_events(segments: list[tuple[float, float]], far_threshold: float, ns_max_mass: float) -> pd.DataFrame:
    """GWOSC candidates inside the injection `segments` with FAR < threshold, one row per event (lowest FAR kept).

    Only periods covered by the injections count: e.g. GW230518 (engineering run ER15, before O4a) is excluded.
    """
    best: dict[str, dict] = {}
    for cat in RATES_EVENT_LISTS:
        try:
            raw = gw.fetch_gwtc_events(cat)
        except Exception as e:  # a list may not exist on GWOSC yet
            print(f"[rates] WARN: could not fetch {cat}: {type(e).__name__}: {e}")
            continue
        for v in (raw.get("events") or {}).values():
            gps, far = v.get("GPS"), v.get("far")
            if gps is None or far is None:
                continue
            gps, far = float(gps), float(far)
            if not any(a <= gps <= b for a, b in segments) or far >= far_threshold:
                continue
            name = v.get("commonName")
            if name in best and best[name]["far_per_yr"] <= far:
                continue
            best[name] = dict(event=name, catalog=cat, gps=gps, far_per_yr=far, p_astro=v.get("p_astro"),
                              mass_1_source=v.get("mass_1_source"), mass_2_source=v.get("mass_2_source"))
    df = pd.DataFrame(list(best.values()), columns=["event", "catalog", "gps", "far_per_yr", "p_astro",
                                                    "mass_1_source", "mass_2_source"])

    def cls(r):
        m1, m2 = r["mass_1_source"], r["mass_2_source"]
        if m1 is None or m2 is None or pd.isna(m1) or pd.isna(m2):
            return "unknown"
        m1, m2 = max(float(m1), float(m2)), min(float(m1), float(m2))
        return "BNS" if m1 < ns_max_mass else ("NSBH" if m2 < ns_max_mass else "BBH")

    df["class"] = df.apply(cls, axis=1) if len(df) else []
    return df.sort_values("gps").reset_index(drop=True)


def _plot_rates_mass_distribution(inj: dict, events: pd.DataFrame, out_png: Path) -> Path:
    """Observed m1 counts and selection-corrected dR/dln(m1) per mass bin (two panels)."""
    import matplotlib.pyplot as plt

    surface, ink, ink2, grid = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df"
    blue, orange = "#2a78d6", "#eb6834"
    m1, m2, z = inj["mass1_source"], inj["mass2_source"], inj["redshift"]
    with np.errstate(divide="ignore", invalid="ignore"):
        ln_pair = np.where((m2 >= 1) & (m2 <= m1), -np.log(np.maximum(m1 - 1, 1e-9)), -np.inf)  # m2 | m1 ~ U[1, m1]
    edges = np.array([1.0, 2.5, 5, 7, 9, 11, 13.5, 17, 21, 26, 32, 39, 47, 57, 70, 100])
    dln = np.diff(np.log(edges))
    vt = np.empty(len(dln))
    for i, (lo, hi) in enumerate(zip(edges[:-1], edges[1:])):
        with np.errstate(divide="ignore", invalid="ignore"):
            ln_m1 = np.where((m1 >= lo) & (m1 < hi), -np.log(m1) - np.log(np.log(hi / lo)), -np.inf)
        vt[i] = _sensitive_vt(inj, ln_m1 + ln_pair, 0.99, 0.99, kappa=0.0)[0]
    known = events[events["class"] != "unknown"]
    em1 = np.maximum(known["mass_1_source"].astype(float), known["mass_2_source"].astype(float)).to_numpy()
    counts, _ = np.histogram(em1, edges)
    with np.errstate(divide="ignore", invalid="ignore"):
        lo_q, med, hi_q = (np.array([_rate_quantiles(int(n), v)[k] for n, v in zip(counts, vt)]) / dln for k in range(3))

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8.2, 7.2), dpi=150, sharex=True,
                                   gridspec_kw=dict(height_ratios=[1, 1.35]))
    fig.patch.set_facecolor(surface)
    for ax in (ax1, ax2):
        ax.set_facecolor(surface)
        ax.set_xscale("log")
        ax.grid(False)   # the gwpy/pesummary styles switch on a full grid
        ax.grid(axis="y", which="major", color=grid, lw=0.8, zorder=0)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        for s in ("left", "bottom"):
            ax.spines[s].set_color(grid)
        ax.tick_params(colors=ink2, labelsize=9)
    ax1.bar(edges[:-1], counts, width=np.diff(edges) * 0.94, align="edge", color=blue, linewidth=0, zorder=2)
    ax1.set_ylabel("Detected events", color=ink2, fontsize=10)
    ax1.set_title(f"Observed: {len(em1)} detections (FAR < threshold)", color=ink, fontsize=10.5, loc="left")
    det = counts > 0
    ax2.errorbar(np.sqrt(edges[1:] * edges[:-1])[det], med[det], yerr=[med[det] - lo_q[det], hi_q[det] - med[det]],
                 fmt="o", ms=5, color=orange, ecolor=orange, elinewidth=1.5, capsize=0, zorder=3)
    ax2.set_yscale("log")
    ax2.set_ylabel("Merger rate dR/d ln m$_1$ (Gpc$^{-3}$ yr$^{-1}$)", color=ink2, fontsize=10)
    ax2.set_title("Selection-corrected: detections ÷ sensitive volume-time of each mass bin",
                  color=ink, fontsize=10.5, loc="left")
    ax2.set_xlabel("Primary mass $m_1$ (source frame, M$_\\odot$)", color=ink2, fontsize=10)
    ax2.set_xticks([1, 2, 5, 10, 20, 35, 50, 100])
    ax2.set_xticklabels(["1", "2", "5", "10", "20", "35", "50", "100"])
    ax2.set_xlim(1, 100)
    fig.text(0.01, 0.008, "Median masses; 90% Poisson intervals; m$_2$ uniform in [1, m$_1$]; no redshift evolution.",
             fontsize=7.5, color=ink2)
    fig.tight_layout(rect=(0, 0.025, 1, 1))
    fig.savefig(out_png, facecolor=surface)
    plt.close(fig)
    return out_png


def run_merger_rates(
    out_rates_tsv: str | Path = "merger_rates.tsv",
    out_events_tsv: Optional[str | Path] = "merger_rates_events.tsv",
    out_report_html: Optional[str | Path] = "merger_rates.html",
    plots_dir: Optional[str | Path] = "rates_plots",
    sensitivity_file: Optional[str | Path] = None,
    far_threshold: float = 1.0,
    ns_max_mass: float = 2.5,
    bbh_kappa: float = 2.9,
    bbh_z_ref: float = 0.2,
) -> pd.DataFrame:
    """Merger rates per population, R = N / <VT>, from catalog events and LVK injections.

    Populations (fixed shapes): BNS with both masses uniform in [1, ns_max_mass];
    NSBH with the BH mass ~ m^-2.35 on [ns_max_mass, 40] and the NS uniform in
    [1, ns_max_mass]; BBH with the GWTC-3 Power Law + Peak model, reported
    without redshift evolution and with R ∝ (1+z)^bbh_kappa at z = bbh_z_ref.
    Returns the rates table (also written to `out_rates_tsv`).
    """
    inj = _load_found_injections(_rates_sensitivity_path(sensitivity_file), far_threshold)
    events = _rates_events(inj["segments"], far_threshold, ns_max_mass)
    counts = events["class"].value_counts().to_dict()
    n_unknown = int(counts.get("unknown", 0))
    m1, m2 = inj["mass1_source"], inj["mass2_source"]
    ns = ns_max_mass

    with np.errstate(divide="ignore", invalid="ignore"):
        bns = np.where(m1 >= m2, np.log(2.0) + _ln_uniform(m1, 1.0, ns) + _ln_uniform(m2, 1.0, ns), -np.inf)
        nsbh = _ln_powerlaw(m1, 2.35, ns, 40.0) + _ln_uniform(m2, 1.0, ns)
        bbh = _ln_power_law_peak(m1, m2)

    rows = []
    for pop, ln_mass, amax1, amax2, kappa, zref, model in (
        ("BNS", bns, 0.4, 0.4, 0.0, 0.0, f"m1, m2 uniform in [1, {ns:g}] Msun; |spin| < 0.4; constant rate"),
        ("NSBH", nsbh, 0.99, 0.4, 0.0, 0.0, f"m_BH ~ m^-2.35 on [{ns:g}, 40], m_NS uniform in [1, {ns:g}]; constant rate"),
        ("BBH", bbh, 0.99, 0.99, bbh_kappa, bbh_z_ref, f"Power Law + Peak (GWTC-3); R ∝ (1+z)^{bbh_kappa:g}"),
        ("BBH", bbh, 0.99, 0.99, 0.0, 0.0, "Power Law + Peak (GWTC-3); constant rate"),
    ):
        n = int(counts.get(pop, 0))
        vt, neff = _sensitive_vt(inj, ln_mass, amax1, amax2, kappa)
        lo, med, hi = _rate_quantiles(n, vt) * (1.0 + zref) ** kappa
        rows.append(dict(population=pop, model=model, z_ref=zref, n_detected=n, vt_gpc3_yr=vt, n_eff=neff,
                         rate_median=med, rate_05=lo, rate_95=hi))
    rates = pd.DataFrame(rows)

    Path(out_rates_tsv).parent.mkdir(parents=True, exist_ok=True)
    rates.to_csv(out_rates_tsv, sep="\t", index=False, float_format="%.6g")
    if out_events_tsv:
        events.to_csv(out_events_tsv, sep="\t", index=False)
    periods = ", ".join(f"{a:.0f}-{b:.0f}" for a, b in inj["segments"])
    print(f"[rates] {len(events)} events with FAR < {far_threshold:g}/yr in the injection periods (GPS {periods}; "
          f"{inj['T_yr']:.2f} yr analysed); {n_unknown} without masses left out")
    for r in rows:
        print(f"[rates] {r['population']:4s} N={r['n_detected']:3d} <VT>={r['vt_gpc3_yr']:.4g} Gpc3 yr "
              f"R={r['rate_median']:.3g} [{r['rate_05']:.3g}, {r['rate_95']:.3g}] Gpc^-3 yr^-1"
              + (f" at z={r['z_ref']:g}" if r["z_ref"] else "") + f"  ({r['model']})")

    if out_report_html:
        out_report_html = Path(out_report_html)
        plot_dir = Path(plots_dir or "rates_plots")
        plot_dir.mkdir(parents=True, exist_ok=True)
        _ensure_matplotlib()
        png = _plot_rates_mass_distribution(inj, events, plot_dir / "rates_mass_distribution.png")
        show = rates.copy()
        show["rate [90%] (Gpc^-3 yr^-1)"] = [f"{m:.3g} [{a:.3g}, {b:.3g}]" for m, a, b in
                                            zip(show.rate_median, show.rate_05, show.rate_95)]
        show = show[["population", "model", "z_ref", "n_detected", "vt_gpc3_yr", "rate [90%] (Gpc^-3 yr^-1)"]]
        paragraphs = [
            f"Merger rates R = N / &lt;VT&gt; from {len(events)} GWOSC candidates with FAR &lt; {far_threshold:g}/yr "
            f"inside the span of the LVK search-sensitivity injections ({inj['T_yr']:.2f} yr of analysis time). "
            f"Events are classified by their median source-frame masses (NS below {ns:g} Msun); "
            f"{n_unknown} candidate(s) without masses are left out.",
            "Intervals are 90% Poisson (Jeffreys prior) for fixed population shapes; the full LVK analyses fit "
            "the shapes together with the rates, so their intervals are wider and model-dependent.",
        ]
        tables = [("Merger rates", show.to_html(index=False, escape=True, float_format=lambda x: f"{x:.4g}")),
                  ("Events used", events.to_html(index=False, escape=True))]
        write_simple_html_report(
            out_report_html, title="GWTC merger rates", paragraphs=paragraphs,
            images=[png], tables=tables,
        )
    return rates
