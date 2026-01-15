# gwtc_analysis/search_skymaps.py
from __future__ import annotations

import csv
import json
import re
from pathlib import Path
from typing import Optional, Iterable, Any

import numpy as np
import healpy as hp
from ligo.skymap.io.fits import read_sky_map

from .read_skymap import plot_skymap_with_ra_dec

def run_search_skymaps(
    *,
    catalogs: list[str],
    out_events_tsv: str,
    out_report_html: Optional[str],
    skymaps_dirs: dict[str, Optional[str]],
    events_json: dict | str,
    ra_deg: float,
    dec_deg: float,
    prob: float,
    make_plots: str = "hits",
    plots_dir: Optional[str] = None,
    skymaps: Optional[list[str]] = None,   # NEW: Galaxy collection support
) -> None:
    ...
    selected_event_ids = _extract_event_ids(events_json)
    requested_percent = 100.0 * prob
    ...
    rows: list[dict[str, Any]] = []

    # NEW: pick input iterator
    if skymaps:
        # Galaxy: explicit list of files (collection elements)
        skymap_paths = []
        for s in skymaps:
            p = Path(s)
            if p.is_dir():
                skymap_paths.extend(list(_iter_skymaps(p)))
            else:
                skymap_paths.append(p)
        # In this mode, "catalog" is whatever user passed; we’ll record the first catalog (or empty)
        catalog_label = catalogs[0] if catalogs else ""
        # normalize plots_dir -> plots_path
        plots_path = Path(plots_dir) if plots_dir else Path("plots")
        plots_path.mkdir(parents=True, exist_ok=True)
        for skymap_path in skymap_paths:
            _process_one_skymap(
                rows=rows,
                catalog=catalog_label,
                skymap_path=skymap_path,
                selected_event_ids=selected_event_ids,
                ra_deg=ra_deg,
                dec_deg=dec_deg,
                prob=prob,
                requested_percent=requested_percent,
                make_plots=make_plots,
                plots_path=plots_path,
            )
    else:
        # Existing behavior: directory per catalog
        for catalog in catalogs:
            base = skymaps_dirs.get(catalog)
            if not base:
                continue
            base_path = Path(base)
            if not base_path.exists():
                rows.append(
                    dict(
                        catalog=catalog,
                        event_id="",
                        skymap_path=str(base_path),
                        status="error",
                        error="skymaps directory not found",
                    )
                )
                continue

            for skymap_path in _iter_skymaps(base_path):
                _process_one_skymap(
                    rows=rows,
                    catalog=catalog,
                    skymap_path=skymap_path,
                    selected_event_ids=selected_event_ids,
                    ra_deg=ra_deg,
                    dec_deg=dec_deg,
                    prob=prob,
                    requested_percent=requested_percent,
                    make_plots=make_plots,
                    plots_path=plots_path,
                )

 
    out_tsv_path = Path(out_events_tsv)
    _write_tsv(out_tsv_path, rows)

    if out_report_html:
        Path(out_report_html).write_text(
            "<html><body><p>search_skymaps: TSV written; plots in ./plots (if enabled).</p></body></html>\n"
        )

def _process_one_skymap(
    *,
    rows: list[dict[str, Any]],
    catalog: str,
    skymap_path: Path,
    selected_event_ids: set[str],
    ra_deg: float,
    dec_deg: float,
    prob: float,
    requested_percent: float,
    make_plots: str,
    plots_path: Path,
) -> None:
    event_id = _event_id_from_filename(skymap_path)

    if selected_event_ids and event_id not in selected_event_ids:
        return

    plot_png = ""
    try:
        cls_at_point = credible_level_at_radec_percent(skymap_path, ra_deg, dec_deg)
        inside = cls_at_point <= requested_percent

        do_plot = (make_plots == "all") or (make_plots == "hits" and inside)
        if do_plot:
            safe_title = _safe_filename(event_id)
            produced = plot_skymap_with_ra_dec(
                str(skymap_path),
                safe_title,
                ra_deg,
                dec_deg,
                "grey",
                contour_levels=(50, 90, prob * 100.0),
            )
            plot_png = _normalize_plot_output(produced, plots_path / f"{safe_title}.png")

        rows.append(
            dict(
                catalog=catalog,
                event_id=event_id,
                skymap_path=str(skymap_path),
                requested_credible_percent=requested_percent,
                credible_at_point_percent=cls_at_point,
                inside_requested_credible=bool(inside),
                plot_png=str(plot_png) if plot_png else "",
                status="ok",
                error="",
            )
        )
    except Exception as e:
        rows.append(
            dict(
                catalog=catalog,
                event_id=event_id,
                skymap_path=str(skymap_path),
                plot_png=str(plot_png) if plot_png else "",
                status="error",
                error=str(e),
            )
        )


def credible_level_at_radec_percent(skymap_path: Path, ra_deg: float, dec_deg: float) -> float:
    """
    Exact greedy credible level at the HEALPix pixel containing (ra,dec).
    Returns percent in [0,100]. Smaller => more inside.
    """
    m = read_sky_map(str(skymap_path), moc=False)
    if isinstance(m, tuple):
        m = m[0]
    m = np.asarray(m, dtype=float)

    nside = hp.get_nside(m)
    theta = np.deg2rad(90.0 - dec_deg)
    phi = np.deg2rad(ra_deg)
    ipix = hp.ang2pix(nside, theta, phi, nest=False)

    order = np.argsort(m)[::-1]
    cdf = np.cumsum(m[order])

    cls = np.empty_like(m)
    cls[order] = 100.0 * cdf
    return float(cls[ipix])


def _iter_skymaps(base: Path) -> Iterable[Path]:
    for p in base.rglob("*"):
        if not p.is_file():
            continue
        n = p.name.lower()
        if n.endswith(".fits") or n.endswith(".fits.gz"):
            yield p


def _event_id_from_filename(f: Path) -> str:
    m = re.search(r"(GW\d{6,}[A-Za-z0-9_-]*)", f.name)
    if m:
        return m.group(1)
    name = f.name
    if name.endswith(".fits.gz"):
        return name[:-8]
    if name.endswith(".fits"):
        return name[:-5]
    return f.stem


def _safe_filename(s: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", s).strip("_") or "event"


def _normalize_plot_output(produced: str, final_path: Path) -> str:
    produced_path = Path(produced)
    candidates = [
        produced_path,
        produced_path.with_suffix(".png"),
    ]
    found = None
    for c in candidates:
        if c.exists() and c.is_file():
            found = c
            break
    if found is None:
        return ""

    final_path.parent.mkdir(parents=True, exist_ok=True)
    if final_path.exists():
        final_path.unlink()
    found.rename(final_path)
    return str(final_path)


def _extract_event_ids(events_json: dict | str) -> set[str]:
    data: Any = events_json
    if isinstance(events_json, str):
        p = Path(events_json)
        if p.exists():
            data = json.loads(p.read_text())
        else:
            return set()

    ids: set[str] = set()

    def add_id(x: Any) -> None:
        if isinstance(x, str) and x.strip():
            ids.add(x.strip())

    if isinstance(data, dict):
        if "id" in data:
            add_id(data["id"])
        if "event_id" in data:
            add_id(data["event_id"])
        if "events" in data:
            ev = data["events"]
            if isinstance(ev, list):
                for item in ev:
                    if isinstance(item, str):
                        add_id(item)
                    elif isinstance(item, dict):
                        add_id(item.get("id") or item.get("event_id") or item.get("name"))
    elif isinstance(data, list):
        for item in data:
            if isinstance(item, str):
                add_id(item)
            elif isinstance(item, dict):
                add_id(item.get("id") or item.get("event_id") or item.get("name"))

    return ids


def _write_tsv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("")
        return

    keys: list[str] = []
    seen = set()
    for r in rows:
        for k in r.keys():
            if k not in seen:
                keys.append(k)
                seen.add(k)

    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow(r)

