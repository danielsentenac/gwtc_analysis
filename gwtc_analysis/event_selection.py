from __future__ import annotations

from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

from . import gw_stat as gw

# Catalog aliases for GWOSC jsonfull endpoints
# (GWTC-4 is currently served as GWTC-4.0; GWTC-3/2.1 have "confident" endpoints for jsonfull)
CATALOG_ALIASES: dict[str, str] = {
    "GWTC-4": "GWTC-4.0",
    "GWTC-3": "GWTC-3-confident",
    "GWTC-2.1": "GWTC-2.1-confident",
    "GWTC-1": "GWTC-1-confident",
}


def _as_float_or_nan(x) -> float:
    try:
        return float(x)
    except Exception:
        return float("nan")


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
) -> None:
    """Select GWTC events based on source-frame component masses and luminosity distance.

    Uses:
      - mass_1_source
      - mass_2_source
      - luminosity_distance

    Writes TSV with selected events (at least event_id).
    """

    out_tsv = Path(out_tsv)
    out_tsv.parent.mkdir(parents=True, exist_ok=True)

    # Expand ALL catalog selector
    if "ALL" in catalogs:
        catalogs = [c for c in gw.ALLOWED_CATALOGS if c != "ALL"]

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
    for col in ["mass_1_source", "mass_2_source", "luminosity_distance"]:
        if col in df_all.columns:
            df_all[col] = df_all[col].apply(_as_float_or_nan)
        else:
            df_all[col] = np.nan

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

    out = df_all.loc[
        mask,
        ["event_id", "catalog_key", "mass_1_source", "mass_2_source", "luminosity_distance"],
    ].copy()

    # Stable order for tests/users
    out = out.sort_values(["catalog_key", "event_id"]).reset_index(drop=True)

    out.to_csv(out_tsv, sep="\t", index=False)
