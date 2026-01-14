from __future__ import annotations

import json
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

from . import gw_stat as gw


def _as_float_or_nan(x):
    try:
        return float(x)
    except Exception:
        return np.nan


def run_event_selection(
    *,
    catalogs: list[str],
    out_tsv: str | Path,
    events_json: Optional[str | Path] = None,
    m1_min: Optional[float] = None,
    m1_max: Optional[float] = None,
    m2_min: Optional[float] = None,
    m2_max: Optional[float] = None,
    dl_min: Optional[float] = None,
    dl_max: Optional[float] = None,
) -> None:
    """
    Select GWTC events based on source-frame component masses and luminosity distance.

    Uses:
      - mass_1_source
      - mass_2_source
      - luminosity_distance

    Writes TSV with selected events (at least event_id).
    """

    out_tsv = Path(out_tsv)
    out_tsv.parent.mkdir(parents=True, exist_ok=True)

    dfs: list[pd.DataFrame] = []

    if events_json is not None:
        raw = json.loads(Path(events_json).read_text(encoding="utf-8"))
        df = gw.events_to_dataframe(raw["events"])
        df["catalog_key"] = "OFFLINE"
        dfs.append(df)
    else:
        # Online mode: fetch each catalog through your existing gw_stat helper
        for cat in catalogs:
            raw = gw.fetch_gwtc_events(catalog=cat)
            df = gw.events_to_dataframe(raw["events"])
            df["catalog_key"] = cat
            dfs.append(df)

    if not dfs:
        out = pd.DataFrame(columns=["event_id", "catalog_key"])
        out.to_csv(out_tsv, sep="\t", index=False)
        return

    df_all = pd.concat(dfs, ignore_index=True)

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

    out = df_all.loc[mask, ["event_id", "catalog_key", "mass_1_source", "mass_2_source", "luminosity_distance"]].copy()

    # Stable order for tests/users
    out = out.sort_values(["catalog_key", "event_id"]).reset_index(drop=True)

    out.to_csv(out_tsv, sep="\t", index=False)

