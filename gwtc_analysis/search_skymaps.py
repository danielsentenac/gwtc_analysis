# gwtc_analysis/search_skymaps.py
from __future__ import annotations

import csv
import json
import re
import warnings
from pathlib import Path
from typing import Optional, Iterable, Any

import healpy as hp
import numpy as np
from ligo.skymap.io.fits import read_sky_map

from .read_skymap import plot_skymap_with_ra_dec


def _tqdm_or_none():
    try:
        from tqdm.auto import tqdm  # type: ignore
        return tqdm
    except Exception:
        return None


def run_search_skymaps(
    *,
    catalogs: list[str],
    out_events_tsv: str | None,
    out_report_html: Optional[str],
    skymaps_dirs: dict[str, Optional[str]],
    events_json: dict | str,
    ra_deg: float,
    dec_deg: float,
    prob: float,
    plots_dir: Optional[str] = None,
    skymaps: Optional[list[str]] = None,  # Galaxy collection support
    data_repo: str = "local",
    waveform: str = "Mixed",
) -> None:
    selected_event_ids = _extract_event_ids(events_json)
    requested_percent = 100.0 * prob

    rows: list[dict[str, Any]] = []

    if out_events_tsv is None:
        out_events_tsv = "search_skymaps.tsv"

    if plots_dir is None:
        plots_dir = "sky_plots"

    plots_path = Path(plots_dir)
    plots_path.mkdir(parents=True, exist_ok=True)

    wf = (waveform or "").strip()
    wf_filter_on = bool(wf) and wf.lower() != "any"

    # -------------------------
    # Zenodo repo
    # -------------------------
    if data_repo == "zenodo":
        from .gw_stat import (
            download_zenodo_skymaps_tarball,
            build_skymap_index_from_tar,
            select_skymap_member,
        )
        import re
        import tarfile

        # waveform selection:
        # - if waveform="any": prefer Mixed, else prefer requested waveform
        prefer = "Mixed" if waveform.lower() == "any" else waveform

        tqdm = _tqdm_or_none()
        processed = 0
        matched = 0
        errors = 0
        miss = 0

        # We iterate by event ids (fast) rather than by tar members (slow).
        # Normalize event ids like "GW200105_162426"
        event_ids = sorted(selected_event_ids) if selected_event_ids else []

        pbar = None
        if tqdm is not None:
            pbar = tqdm(total=len(event_ids) if event_ids else None, unit=" event", desc="Zenodo skymaps")

        for catalog in catalogs:
            tar_path = download_zenodo_skymaps_tarball(catalog, progress=True, verbose=False)
            index = build_skymap_index_from_tar(tar_path, verbose=False)

            # If no events_json filtering is provided, fall back to scanning all events in index (heavier)
            if not event_ids:
                # build a unique set of event keys from the index
                event_keys = sorted({ev for (ev, _approx) in index.keys()})
            else:
                # map selected_event_ids -> "GW\d{6}_\d{6}" keys used in index
                event_keys = []
                for ev in event_ids:
                    m = re.search(r"(GW\d{6}_\d{6})", str(ev))
                    if m:
                        event_keys.append(m.group(1))

            with tarfile.open(tar_path, mode="r:gz") as tf:
                for ev_key in event_keys:
                    member_name = select_skymap_member(index, ev_key, prefer=prefer)
                    if member_name is None:
                        miss += 1
                        if pbar is not None:
                            pbar.update(1)
                            try:
                                pbar.set_postfix(matched=matched, miss=miss, errors=errors)
                            except Exception:
                                pass
                        continue

                    try:
                        f = tf.extractfile(member_name)
                        if f is None:
                            raise OSError(f"could not extract {member_name}")

                        data = f.read()

                        before_len = len(rows)
                        _process_one_skymap_zenodo_member(
                            rows=rows,
                            catalog=catalog,
                            member_name=member_name,
                            member_bytes=data,
                            selected_event_ids=selected_event_ids,
                            ra_deg=ra_deg,
                            dec_deg=dec_deg,
                            prob=prob,
                            requested_percent=requested_percent,
                            plots_path=plots_path,
                        )

                        # infer hit/error from appended row
                        if len(rows) > before_len:
                            last = rows[-1]
                            if last.get("status") == "ok" and last.get("inside_requested_credible") is True:
                                matched += 1
                            if last.get("status") == "error":
                                errors += 1

                    except Exception:
                        errors += 1
                        rows.append(
                            dict(
                                catalog=catalog,
                                event_id=ev_key,
                                skymap_path=f"zenodo:{catalog}:{member_name}",
                                status="error",
                                error="extract_or_process_failed",
                            )
                        )

                    processed += 1
                    if pbar is not None:
                        pbar.update(1)
                        try:
                            pbar.set_postfix(matched=matched, miss=miss, errors=errors)
                        except Exception:
                            pass
                        try:
                            pbar.set_description(f"Zenodo skymaps ({catalog})")
                        except Exception:
                            pass

        if pbar is not None:
            pbar.close()

        out_tsv_path = Path(out_events_tsv)
        _write_tsv(out_tsv_path, rows)

        if out_report_html:
            _write_search_skymaps_report(
                out_path=Path(out_report_html),
                rows=rows,
                tsv_path=Path(out_events_tsv),
            )
        return
    # -------------------------
    # S3 repo
    # -------------------------
    if data_repo == "s3":
        from .gw_stat import _s3_client_from_env, _catalog_to_s3_prefix

        client = _s3_client_from_env()
        bucket = "gwtc"

        tqdm = _tqdm_or_none()
        matched = 0
        errors = 0

        pbar = tqdm(total=None, unit=" skymap", desc="S3 skymaps") if tqdm is not None else None

        try:
            for catalog in catalogs:
                prefix = _catalog_to_s3_prefix(catalog)

                # Safety: never allow empty / root prefix scans
                if not prefix or prefix == "/":
                    raise RuntimeError(
                        f"Refusing to scan S3 with unsafe prefix={prefix!r} for catalog={catalog!r}"
                    )

                for obj in client.list_objects(bucket_name=bucket, prefix=prefix, recursive=True):
                    key = obj.object_name
                    if not (key.endswith(".fits") or key.endswith(".fits.gz")):
                        continue

                    if wf_filter_on and wf not in key:
                        continue

                    try:
                        resp = client.get_object(bucket_name=bucket, object_name=key)
                        try:
                            data = resp.read()
                        finally:
                            resp.close()
                            resp.release_conn()

                        before_len = len(rows)
                        _process_one_skymap_s3_object(
                            rows=rows,
                            catalog=catalog,
                            s3_key=key,
                            s3_bytes=data,
                            selected_event_ids=selected_event_ids,
                            ra_deg=ra_deg,
                            dec_deg=dec_deg,
                            prob=prob,
                            requested_percent=requested_percent,
                            plots_path=plots_path,
                        )

                        if len(rows) > before_len:
                            last = rows[-1]
                            if (
                                last.get("status") == "ok"
                                and last.get("inside_requested_credible") is True
                            ):
                                matched += 1
                            if last.get("status") == "error":
                                errors += 1

                    except Exception:
                        errors += 1
                        rows.append(
                            dict(
                                catalog=catalog,
                                event_id="",
                                skymap_path=key,
                                status="error",
                                error="download_or_process_failed",
                            )
                        )

                    if pbar is not None:
                        pbar.update(1)
                        try:
                            pbar.set_postfix(matched=matched, errors=errors)
                        except Exception:
                            pass
                        try:
                            pbar.set_description(f"S3 skymaps ({catalog})")
                        except Exception:
                            pass

        finally:
            if pbar is not None:
                pbar.close()

        out_tsv_path = Path(out_events_tsv)
        _write_tsv(out_tsv_path, rows)

        if out_report_html:
            _write_search_skymaps_report(
             out_path=Path(out_report_html),
             rows=rows,
             tsv_path=Path(out_events_tsv),
        )

        return

    # -------------------------
    # Local / Galaxy
    # -------------------------

    if skymaps:
        # Galaxy: explicit list of files (collection elements)
        skymap_paths: list[Path] = []
        for s in skymaps:
            p = Path(s)
            if p.is_dir():
                skymap_paths.extend(list(_iter_skymaps(p)))
            else:
                skymap_paths.append(p)

        catalog_label = catalogs[0] if catalogs else ""

        for skymap_path in skymap_paths:
            if wf_filter_on and wf not in skymap_path.name:
                continue
            _process_one_skymap(
                rows=rows,
                catalog=catalog_label,
                skymap_path=skymap_path,
                selected_event_ids=selected_event_ids,
                ra_deg=ra_deg,
                dec_deg=dec_deg,
                prob=prob,
                requested_percent=requested_percent,
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
                if wf_filter_on and wf not in skymap_path.name:
                    continue
                _process_one_skymap(
                    rows=rows,
                    catalog=catalog,
                    skymap_path=skymap_path,
                    selected_event_ids=selected_event_ids,
                    ra_deg=ra_deg,
                    dec_deg=dec_deg,
                    prob=prob,
                    requested_percent=requested_percent,
                    plots_path=plots_path,
                )

    out_tsv_path = Path(out_events_tsv)
    _write_tsv(out_tsv_path, rows)

    if out_report_html:
        Path(out_report_html).write_text(
            "<html><body><p>search_skymaps: TSV written; plots in ./plots (if enabled).</p></body></html>\n",
            encoding="utf-8",
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
    plots_path: Path,
) -> None:
    event_id = _event_id_from_filename(skymap_path)

    if selected_event_ids and event_id not in selected_event_ids:
        return

    plot_png = ""
    try:
        cls_at_point = credible_level_at_radec_percent(skymap_path, ra_deg, dec_deg)
        inside = cls_at_point <= requested_percent

        # Plot ONLY for hits (prevents lots of irrelevant plots for small prob)
        if inside:
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

def _process_one_skymap_zenodo_member(
    *,
    rows: list[dict[str, Any]],
    catalog: str,
    member_name: str,
    member_bytes: bytes,
    selected_event_ids: set[str],
    ra_deg: float,
    dec_deg: float,
    prob: float,
    requested_percent: float,
    plots_path: Path,
) -> None:
    """Adapter for Zenodo skymaps.

    - writes member bytes to a temporary FITS/FITS.GZ file
    - suppresses noisy FITS/plot warnings during batch runs
    - reuses _process_one_skymap() unchanged
    - patches output row so TSV contains Zenodo member path (not deleted temp path)
    """
    tmp_dir = plots_path / ".tmp_skymaps"
    tmp_dir.mkdir(parents=True, exist_ok=True)

    tmp_name = Path(member_name).name
    tmp_path = tmp_dir / tmp_name

    try:
        tmp_path.write_bytes(member_bytes)

        with warnings.catch_warnings():
            # ligo.skymap / matplotlib noise
            warnings.filterwarnings(
                "ignore",
                message=r".*Setting the 'color' property.*",
                category=UserWarning,
            )
            warnings.filterwarnings(
                "ignore",
                message=r".*set_ticklabels\(\) should only be used.*FixedLocator.*",
                category=UserWarning,
            )

            # astropy WCS FITSFixedWarning: RADECSYS deprecated
            try:
                from astropy.wcs.wcs import FITSFixedWarning  # type: ignore

                warnings.filterwarnings(
                    "ignore",
                    category=FITSFixedWarning,
                    module=r"astropy\.wcs\.wcs",
                )
            except Exception:
                pass

            before_len = len(rows)

            _process_one_skymap(
                rows=rows,
                catalog=catalog,
                skymap_path=tmp_path,
                selected_event_ids=selected_event_ids,
                ra_deg=ra_deg,
                dec_deg=dec_deg,
                prob=prob,
                requested_percent=requested_percent,
                plots_path=plots_path,
            )

            # Patch the row added by _process_one_skymap() so it reflects Zenodo origin
            if len(rows) > before_len:
                rows[-1]["skymap_path"] = f"zenodo:{catalog}:{member_name}"
                rows[-1]["skymap_local_tmp"] = ""

    finally:
        try:
            tmp_path.unlink(missing_ok=True)
        except Exception:
            pass

def _process_one_skymap_s3_object(
    *,
    rows: list[dict[str, Any]],
    catalog: str,
    s3_key: str,
    s3_bytes: bytes,
    selected_event_ids: set[str],
    ra_deg: float,
    dec_deg: float,
    prob: float,
    requested_percent: float,
    plots_path: Path,
) -> None:
    """Adapter for S3 skymaps.

    - writes bytes to a temporary FITS/FITS.GZ file
    - suppresses noisy FITS/plot warnings during batch runs
    - reuses _process_one_skymap() unchanged
    - patches the output row so TSV contains the S3 key (not the deleted temp path)
    """
    tmp_dir = plots_path / ".tmp_skymaps"
    tmp_dir.mkdir(parents=True, exist_ok=True)

    tmp_name = Path(s3_key).name
    tmp_path = tmp_dir / tmp_name

    try:
        tmp_path.write_bytes(s3_bytes)

        with warnings.catch_warnings():
            # ligo.skymap / matplotlib noise
            warnings.filterwarnings(
                "ignore",
                message=r".*Setting the 'color' property.*",
                category=UserWarning,
            )
            warnings.filterwarnings(
                "ignore",
                message=r".*set_ticklabels\(\) should only be used.*FixedLocator.*",
                category=UserWarning,
            )

            # astropy WCS FITSFixedWarning: RADECSYS deprecated
            try:
                from astropy.wcs.wcs import FITSFixedWarning  # type: ignore

                warnings.filterwarnings(
                    "ignore",
                    category=FITSFixedWarning,
                    module=r"astropy\.wcs\.wcs",
                )
            except Exception:
                pass

            before_len = len(rows)

            _process_one_skymap(
                rows=rows,
                catalog=catalog,
                skymap_path=tmp_path,
                selected_event_ids=selected_event_ids,
                ra_deg=ra_deg,
                dec_deg=dec_deg,
                prob=prob,
                requested_percent=requested_percent,
                plots_path=plots_path,
            )

            # Patch the row added by _process_one_skymap() so it reflects S3 origin
            if len(rows) > before_len:
                rows[-1]["skymap_path"] = s3_key
                # Optional: keep a debug column (empty by default to avoid confusion)
                rows[-1]["skymap_local_tmp"] = ""

    finally:
        try:
            tmp_path.unlink(missing_ok=True)
        except Exception:
            pass


def _write_search_skymaps_report(*, out_path: Path, rows: list[dict[str, Any]], tsv_path: Path) -> None:
    import html

    n_total = len(rows)
    n_ok = sum(1 for r in rows if r.get("status") == "ok")
    n_err = sum(1 for r in rows if r.get("status") == "error")
    n_hit = sum(1 for r in rows if r.get("status") == "ok" and r.get("inside_requested_credible") is True)
    n_plots = sum(1 for r in rows if r.get("plot_png"))

    # Show only hits with plots
    hit_rows = [
        r for r in rows
        if r.get("status") == "ok"
        and r.get("inside_requested_credible") is True
        and r.get("plot_png")
    ]

    lines: list[str] = []
    lines.append("<html><head><meta charset='utf-8'><title>search_skymaps report</title></head><body>")
    lines.append("<h1>search_skymaps report</h1>")
    lines.append("<ul>")
    lines.append(f"<li>Total processed: {n_total}</li>")
    lines.append(f"<li>OK: {n_ok} | Errors: {n_err}</li>")
    lines.append(f"<li>Hits (inside requested credible region): {n_hit}</li>")
    lines.append(f"<li>Plots produced: {n_plots}</li>")
    lines.append(f"<li>TSV: {html.escape(str(tsv_path))}</li>")
    lines.append("</ul>")

    if not hit_rows:
        lines.append("<p><b>No hits</b> (no skymaps contained the requested position at the requested probability).</p>")
    else:
        lines.append("<h2>Hits</h2>")
        lines.append("<table border='1' cellspacing='0' cellpadding='6'>")
        lines.append("<tr><th>Catalog</th><th>Event</th><th>Skymap</th><th>Credible @ point (%)</th><th>Plot</th></tr>")

        for r in hit_rows:
            catalog = html.escape(str(r.get("catalog", "")))
            event_id = html.escape(str(r.get("event_id", "")))
            skymap = html.escape(str(r.get("skymap_path", "")))
            cap = r.get("credible_at_point_percent", "")
            cap_s = html.escape(f"{cap:.3g}" if isinstance(cap, (int, float)) else str(cap))
            plot_png = str(r.get("plot_png", ""))

            # Use relative paths if possible (nicer HTML portability)
            try:
                rel_plot = str(Path(plot_png).relative_to(out_path.parent))
            except Exception:
                rel_plot = plot_png
            rel_plot_esc = html.escape(rel_plot)

            lines.append(
                "<tr>"
                f"<td>{catalog}</td>"
                f"<td>{event_id}</td>"
                f"<td><code>{skymap}</code></td>"
                f"<td>{cap_s}</td>"
                f"<td><a href='{rel_plot_esc}'><img src='{rel_plot_esc}' style='max-width:480px; height:auto;'/></a></td>"
                "</tr>"
            )

        lines.append("</table>")

    lines.append("</body></html>\n")
    out_path.write_text("\n".join(lines), encoding="utf-8")


def credible_level_at_radec_percent(skymap_path: Path, ra_deg: float, dec_deg: float) -> float:
    """Exact greedy credible level at the HEALPix pixel containing (ra,dec).

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
    candidates = [produced_path, produced_path.with_suffix(".png")]
    found: Path | None = None
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
            data = json.loads(p.read_text(encoding="utf-8"))
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
        path.write_text("", encoding="utf-8")
        return

    keys: list[str] = []
    seen: set[str] = set()
    for r in rows:
        for k in r.keys():
            if k not in seen:
                keys.append(k)
                seen.add(k)

    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=keys, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow(r)
