from __future__ import annotations

import os
from pathlib import Path
from typing import List, Optional, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import difflib

from . import gw_stat as gw
from .report import write_simple_html_report


ALLOWED_CATALOGS = [
    "GWTC-2.1",
    "GWTC-3",
    "GWTC-4",
]

def _plot_source_type_pie(df: pd.DataFrame, out_png: Path, column: str = "binary_type") -> Path | None:
    s = df.get(column)
    if s is None:
        return None
    s = s.dropna()
    if s.empty:
        return None

    counts = s.value_counts()

    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    ax.pie(counts.values, labels=counts.index, autopct="%1.0f%%", startangle=90)
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
    import matplotlib.pyplot as plt

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
    skymaps_dirs: Optional[dict] = None,  # values can be str (dir) OR list[str] (files)
    ns_threshold: float = 3.0,
    events_json: Optional[str | Path] = None,
    data_repo="local",
) -> None:
    """
    Fetch events from GWOSC jsonfull for one or more catalogs and compute basic derived columns.

    Outputs
    -------
    out_events_tsv : TSV
        Per-event table with masses, distance, detector network (optional), and credible area (optional).
    out_report_html : HTML (optional)
        Single-file HTML report with summary tables and embedded plots.
    """
    if not catalogs:
        raise ValueError("No catalogs selected")

    # If reading from GWOSC, validate catalog names early
    if events_json is None:
        _validate_catalogs(catalogs)

    if area_column is None:
        level = int(round(100 * area_cred))
        area_column = f"A{level}_deg2"

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
        include_area = False

    # Detectors network (requires GWOSC v2 calls)
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

    # Credible area (optional)
    if include_area:
        df[area_column] = np.nan

    skymaps_dirs = skymaps_dirs or {}

    # Cache dir for Zenodo downloads (kept hidden; can be centralized later)
    zenodo_cache_dir = ".cache_gwosc"

    # --- Per-catalog credible area fill ---
    for cat in catalogs:
        m = df["catalog_key"].eq(cat)
        if not m.any():
            continue

        # If area not requested, nothing to do here
        if not include_area:
            continue

        # --- ZENODO: do not require skymaps_dirs ---
        if data_repo == "zenodo":
            if not hasattr(gw, "add_localization_area_from_zenodo"):
                raise RuntimeError("gw_stat missing add_localization_area_from_zenodo()")

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


        # --- local data-repo (directory/collection required) ---
        skydir = skymaps_dirs.get(cat)

        if skydir is None:
            # fallback: try mapping from a base galaxy_inputs directory (interactive-style),
            # if gw_stat provides that helper.
            if hasattr(gw, "add_localization_area_from_galaxy_inputs"):
                tmp = gw.add_localization_area_from_galaxy_inputs(
                    df.loc[m],
                    catalog_key=cat,
                    cred=area_cred,
                    column=area_column,
                    base_dir="galaxy_inputs",
                    progress=True,
                    verbose=False,
                )
                df.loc[m, area_column] = tmp[area_column].values
                continue

            raise RuntimeError(
                f"Credible area requested but no skymaps directory was provided for {cat}. "
                "Provide the corresponding skymaps collection input."
            )

        # --- Normalize skymaps input: directory path OR list of files (Galaxy collection) ---
        skydir_path: Path
        if isinstance(skydir, (str, Path)):
            skydir_path = Path(skydir)
        elif isinstance(skydir, (list, tuple)):
            # Galaxy may pass a 1-element list containing a directory
            if len(skydir) == 1 and Path(skydir[0]).is_dir():
                skydir_path = Path(skydir[0])
            else:
                # Materialize the collection into a local directory of symlinks/copies
                workdir_for_links = (
                    Path(out_events_tsv).parent if Path(out_events_tsv).parent != Path("") else Path(".")
                )
                skydir_path = _materialize_skymap_collection_dir(list(skydir), workdir_for_links, tag=cat)
        else:
            raise RuntimeError(f"Unsupported skymaps_dirs entry for {cat}: {type(skydir)}")

        if not skydir_path.exists():
            raise RuntimeError(f"Skymaps directory does not exist for {cat}: {skydir_path}")

        skydir = str(skydir_path)

        # Use the directory directly (collection directory)
        if hasattr(gw, "add_localization_area_from_directory"):
            tmp = gw.add_localization_area_from_directory(
                df.loc[m],
                skymap_dir=skydir,
                cred=area_cred,
                column=area_column,
                progress=True,
                verbose=False,
            )
            df.loc[m, area_column] = tmp[area_column].values
        else:
            raise RuntimeError("gw_stat missing add_localization_area_from_directory()")

    # Write per-event table
    out_events_tsv = Path(out_events_tsv)
    _safe_mkdir(out_events_tsv.parent if out_events_tsv.parent != Path("") else ".")
    df_out = df.copy()

    # Normalize detectors list to comma-joined for TSV
    if "detectors" in df_out.columns:
        df_out["detectors"] = df_out["detectors"].apply(lambda x: ",".join(x) if isinstance(x, list) else "")

    # Keep a practical column subset (still fairly rich)
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

    # Report (optional)
    if not out_report_html:
        return

    out_report_html = Path(out_report_html)
    _safe_mkdir(out_report_html.parent if out_report_html.parent != Path("") else ".")

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

    # Detector network PIE only (no histogram)
    if include_detectors:
        if fig_network is None:
            raise RuntimeError(
                "include_detectors=True but fig_network is None. "
                "Call gw.add_detectors_and_virgo_flag(..., plot_network_pie=True)."
            )
        p_net = workdir / "network_pie.png"
        fig_network.savefig(p_net, dpi=150, bbox_inches="tight")
        img_paths.append(p_net)

    # Source type PIE (BBH/BNS/NSBH/...)
    p_type = _plot_source_type_pie(
        df_out,
        workdir / "source_types_pie.png",
        column="binary_type",
    )
    if p_type:
        img_paths.append(p_type)

    # Area CDF
    if include_area:
        p_area = _plot_area_cdf(
            df_out,
            workdir / "area_cdf_bbh.png",
            column=area_column,
            catalog_label="GWTC-4",
            source_type="BBH",
            from_zenodo=(data_repo == "zenodo"),
        )
        if p_area:
            img_paths.append(p_area)


    paragraphs = [
        f"Catalogs: {', '.join(catalogs)}",
        f"Total events after basic cleaning: {n_total}",
        f"Per-event table written to: {out_events_tsv.name}",
    ]
    if include_area:
        got = int(pd.notna(df_out[area_column]).sum())
        paragraphs.append(
            f"Credible area computed at cred={area_cred}: {got}/{n_total} events."
        )

    write_simple_html_report(
        out_report_html,
        title="GWTC catalog statistics",
        paragraphs=paragraphs,
        images=img_paths,
        tables=tables,
    )




    
def _materialize_skymap_collection_dir(files: list[str], workdir: Path, tag: str) -> Path:
    """
    Create a directory under workdir and symlink all skymap files into it.
    This lets us reuse add_localization_area_from_directory().
    """
    out = workdir / f"skymaps_{tag}"
    out.mkdir(parents=True, exist_ok=True)

    for f in files:
        src = Path(f)
        if not src.exists():
            continue
        # keep original filename to preserve GW event id parsing
        dst = out / src.name
        if dst.exists():
            continue
        try:
            dst.symlink_to(src)
        except OSError:
            # fallback: hardlink/copy if symlink not allowed
            import shutil
            shutil.copy2(src, dst)

    return out

