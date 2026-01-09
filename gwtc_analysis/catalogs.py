from __future__ import annotations

import os
from pathlib import Path
from typing import List, Optional, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from . import gw_stat as gw
from .report import write_simple_html_report

def _safe_mkdir(p: str | Path) -> None:
    Path(p).mkdir(parents=True, exist_ok=True)

def _plot_network_counts(df: pd.DataFrame, out_png: str | Path) -> Optional[str]:
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

def _plot_a90_cdf(df: pd.DataFrame, out_png: str | Path, column: str = "A90_deg2") -> Optional[str]:
    if column not in df.columns:
        return None
    x = df[column].dropna().astype(float).to_numpy()
    if x.size == 0:
        return None
    x = np.sort(x)
    y = np.arange(1, len(x)+1)/len(x)
    fig, ax = plt.subplots(figsize=(6,4))
    ax.plot(x, y, lw=1.6)
    ax.set_xscale("log")
    ax.set_xlabel(r"$A_{90}$ [deg$^2$]")
    ax.set_ylabel("Cumulative fraction")
    ax.set_title("Sky localization (CDF)")
    ax.grid(True, which="both", alpha=0.35)
    fig.tight_layout()
    fig.savefig(out_png, dpi=150)
    plt.close(fig)
    return str(out_png)

def run_catalog_statistics(
    catalogs: List[str],
    out_events_tsv: str | Path,
    out_report_html: Optional[str | Path] = None,
    include_detectors: bool = True,
    include_a90: bool = False,
    a90_cred: float = 0.9,
    a90_column: str = "A90_deg2",
    skymaps_dirs: Optional[dict] = None,
    ns_threshold: float = 3.0,
    events_json: Optional[str | Path] = None,
) -> None:
    """
    Fetch events from GWOSC jsonfull for one or more catalogs and compute basic derived columns.

    Outputs
    -------
    out_events_tsv : TSV
        Per-event table with masses, distance, detector network (optional), and A90 (optional).
    out_report_html : HTML (optional)
        Single-file HTML report with summary tables and embedded plots.
    """
    if not catalogs:
        raise ValueError("No catalogs selected")

    import json

    dfs = []

    if events_json is not None:
        raw = json.loads(Path(events_json).read_text(encoding="utf-8"))
        df_cat = gw.events_to_dataframe(raw["events"])
        df_cat["catalog_key"] = "OFFLINE"
        dfs.append(df_cat)

    else:
        for cat in catalogs:
            raw = gw.fetch_gwtc_events(catalog=cat)
            df_cat = gw.events_to_dataframe(raw["events"])
            df_cat["catalog_key"] = cat
            dfs.append(df_cat)

    df0 = pd.concat(dfs, ignore_index=True)
    df = gw.prepare_catalog_df(df0, ns_threshold=ns_threshold)
    # In offline mode, network calls are not available.
    if events_json is not None:
        include_detectors = False
        include_a90 = False

    # Detectors network (requires GWOSC v2 calls)
    fig_network = None
    if include_detectors:
        df, fig_network = gw.add_detectors_and_virgo_flag(df, progress=True, verbose=False, plot_network_pie=False)
        df["n_det"] = df["detectors"].apply(lambda x: len(x) if isinstance(x, list) else np.nan)
    else:
        df["detectors"] = np.nan
        df["n_det"] = np.nan
        df["has_V1"] = np.nan

    # A90 (optional): use explicit Galaxy collection directories if provided.
    # For pure Galaxy tools, we expect skymaps collections to be passed as directory paths.
    if include_a90:
        df[a90_column] = np.nan

    skymaps_dirs = skymaps_dirs or {}

    for cat in catalogs:
        m = df["catalog_key"].eq(cat)
        if not m.any():
            continue

        skydir = skymaps_dirs.get(cat)

        if skydir is None:
            # fallback: try mapping from a base galaxy_inputs directory (interactive-style),
            # if gw_stat provides that helper.
            if hasattr(gw, "add_localization_area_from_galaxy_inputs"):
                tmp = gw.add_localization_area_from_galaxy_inputs(
                    df.loc[m],
                    catalog_key=cat,
                    cred=a90_cred,
                    column=a90_column,
                    base_dir="galaxy_inputs",
                    progress=True,
                    verbose=False,
                )
                df.loc[m, a90_column] = tmp[a90_column].values
                continue

            raise RuntimeError(
                f"A90 requested but no skymaps directory was provided for {cat}. "
                "Provide the corresponding skymaps collection input."
            )

        # Use the directory directly (collection directory)
        if hasattr(gw, "add_localization_area_from_directory"):
            tmp = gw.add_localization_area_from_directory(
                df.loc[m],
                skymap_dir=skydir,
                cred=a90_cred,
                column=a90_column,
                progress=True,
                verbose=False,
            )
            df.loc[m, a90_column] = tmp[a90_column].values
        else:
            raise RuntimeError("gw_stat missing add_localization_area_from_directory()")

    # Write per-event table
    out_events_tsv = Path(out_events_tsv)
    _safe_mkdir(out_events_tsv.parent if out_events_tsv.parent != Path('') else ".")
    df_out = df.copy()

    # Normalize detectors list to comma-joined for TSV
    if "detectors" in df_out.columns:
        df_out["detectors"] = df_out["detectors"].apply(lambda x: ",".join(x) if isinstance(x, list) else "")

    # Keep a practical column subset (still fairly rich)
    keep = [
        "event_id","catalog_key","version",
        "mass_1_source","mass_2_source","chirp_mass_source","total_mass_source","final_mass_source",
        "luminosity_distance","redshift","chi_eff","chi_p","snr","far","p_astro",
        "binary_type","detectors","n_det","has_V1",
    ]
    if include_a90:
        keep.append(a90_column)
    keep = [c for c in keep if c in df_out.columns]
    df_out[keep].to_csv(out_events_tsv, sep="\t", index=False)

    # Report
    out_report_html = Path(out_report_html)
    _safe_mkdir(out_report_html.parent if out_report_html.parent != Path('') else ".")

    # Summary tables
    n_total = len(df_out)
    per_cat = df_out["catalog_key"].value_counts().rename_axis("catalog").reset_index(name="N")
    per_type = df_out["binary_type"].value_counts().rename_axis("binary_type").reset_index(name="N")
    per_net = df_out["n_det"].value_counts(dropna=True).sort_index().rename_axis("n_det").reset_index(name="N")

    tables = []
    tables.append(("Counts by catalog", per_cat.to_html(index=False, escape=True)))
    tables.append(("Counts by binary type", per_type.to_html(index=False, escape=True)))
    if include_detectors:
        tables.append(("Counts by detector number", per_net.to_html(index=False, escape=True)))

    # Plots
    img_paths = []
    workdir = out_report_html.parent
    p1 = _plot_network_counts(df_out, workdir/"network_counts.png")
    if p1: img_paths.append(p1)
    if include_a90:
        p2 = _plot_a90_cdf(df_out, workdir/"a90_cdf.png", column=a90_column)
        if p2: img_paths.append(p2)

    paragraphs = [
        f"Catalogs: {', '.join(catalogs)}",
        f"Total events after basic cleaning: {n_total}",
        f"Per-event table written to: {out_events_tsv.name}",
    ]
    if include_a90:
        got = int(pd.notna(df_out[a90_column]).sum())
        paragraphs.append(f"A90 computed at cred={a90_cred}: {got}/{n_total} events.")

    write_simple_html_report(
        out_report_html,
        title="GWTC catalog statistics",
        paragraphs=paragraphs,
        images=img_paths,
        tables=tables,
    )
