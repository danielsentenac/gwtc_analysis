from __future__ import annotations
import base64
from pathlib import Path
from typing import List, Optional

# The options of the command-line run, set by cli.main and listed at the end of every report
_RUN_OPTIONS: Optional[dict] = None


def set_run_options(mode: str, rows: list[tuple[str, object, object, bool]]) -> None:
    """Record the options of this run: (option, value, default, set by the user) for each option of the mode,
    defaults included."""
    global _RUN_OPTIONS
    _RUN_OPTIONS = dict(mode=mode, rows=rows)


_PE_FILES: dict[str, str] = {}       # event -> PE file read for it
_PE_SOURCES: dict[str, str] = {}     # PE file -> where this run took it from


def record_pe_file(event: str, filename: str, source: Optional[str] = None) -> None:
    """Record the PE file read for an event (and, if known, where it came from) for the report."""
    _PE_FILES[str(event)] = Path(str(filename)).name
    if source:
        note_pe_source(filename, source)


def note_pe_source(filename: str, source: str) -> None:
    """Where this run took a PE file from: Zenodo, the S3 mirror, Galaxy, a PE cache, ..."""
    _PE_SOURCES[Path(str(filename)).name] = source


_DATA_FILES: dict[str, str] = {}     # public data file (URL) -> where this run took it from


def record_data_file(url: str, source: str) -> None:
    """Record a public data file the run used (injections, data releases, ...) and where it came from."""
    _DATA_FILES[str(url)] = source


def data_files_html() -> str:
    """The public data files used by the run, as an HTML table (empty when it used none)."""
    if not _DATA_FILES:
        return ""
    from html import escape

    body = "\n".join(f"<tr><td><code>{escape(u.rsplit('/', 1)[-1])}</code></td><td><a href='{escape(u)}'>"
                      f"{escape(u)}</a></td><td>{escape(src)}</td></tr>" for u, src in sorted(_DATA_FILES.items()))
    return ("<h2>Public data files of this run</h2>\n<p>Other public files read by the run, with the URL that "
            "identifies them; with <code>--data-repo s3</code> a file is downloaded from the S3 mirror when it holds "
            "the same file.</p>\n<table border='1' cellspacing='0' cellpadding='4'><tr><th>File</th><th>Source</th>"
            f"<th>Taken from</th></tr>\n{body}\n</table>\n")


def pe_files_html() -> str:
    """The PE files read by the run, with their type, as an HTML table (empty when the run read none)."""
    if not _PE_FILES:
        return ""
    from collections import Counter
    from html import escape

    from .pe_files import describe

    kinds = Counter(describe(f)[0] for f in _PE_FILES.values())
    body = "\n".join(
        f"<tr><td>{escape(ev)}</td><td><code>{escape(f)}</code></td><td>{escape(describe(f)[1])}</td>"
        f"<td>{escape(_PE_SOURCES.get(f, 'samples already extracted from it (PE cache)'))}</td></tr>"
        for ev, f in sorted(_PE_FILES.items()))
    count = ", ".join(f"{n} {k}" for k, n in kinds.most_common())
    return (f"<h2>PE files of this run</h2>\n<p>{len(_PE_FILES)} event(s): {escape(count)}. The file of each event is "
            "chosen in the Zenodo listing of its release by the rule of the mode, whatever --data-repo (see "
            "<a href='https://danielsentenac.github.io/gwtc_analysis/data-sources/#which-pe-file-is-read'>Which PE "
            "file is read</a>).</p>\n<table border='1' cellspacing='0' cellpadding='4'><tr><th>Event</th>"
            f"<th>File</th><th>Type</th><th>Taken from</th></tr>\n{body}\n</table>\n")


def options_html() -> str:
    """The PE files and the options of the run as HTML tables (the options are left out when the run did not
    come from the command line)."""
    if not _RUN_OPTIONS:
        return pe_files_html() + data_files_html()
    from html import escape

    def fmt(v):
        if v is None:
            return "—"
        if isinstance(v, (list, tuple)):
            return " ".join(map(str, v)) if v else "—"
        return str(v)

    body = "\n".join(
        f"<tr><td><code>{escape(o)}</code></td><td>{'<b>' if given else ''}{escape(fmt(v))}"
        f"{'</b>' if given else ''}</td><td>{escape(fmt(d))}</td></tr>"
        for o, v, d, given in _RUN_OPTIONS["rows"])
    return pe_files_html() + data_files_html() + (
        f"<h2>Options of this run</h2>\n<p><code>gwtc_analysis {escape(_RUN_OPTIONS['mode'])}</code>: every "
        "option with the value used; in bold, those set on the command line or in a settings file.</p>\n"
        "<table border='1' cellspacing='0' cellpadding='4'><tr><th>Option</th><th>Value</th><th>Default</th>"
        f"</tr>\n{body}\n</table>\n")


def _img_to_data_uri(path: str | Path) -> str:
    p = Path(path)
    data = p.read_bytes()
    suffix = p.suffix.lower().lstrip(".")
    mime = (
        "image/png"
        if suffix == "png"
        else "image/jpeg"
        if suffix in ("jpg", "jpeg")
        else "application/octet-stream"
    )
    b64 = base64.b64encode(data).decode("ascii")
    return f"data:{mime};base64,{b64}"

def write_simple_html_report(
    out_html: str | Path,
    title: str,
    paragraphs: List[str],
    images: Optional[List[str | Path]] = None,
    tables: Optional[List[tuple[str, str]]] = None,
) -> None:
    """
    Create a single-file HTML report suitable for Galaxy (no external assets).

    - paragraphs: list of HTML-escaped or simple text paragraphs
    - images: list of image paths; embedded as data URIs
    - tables: list of (caption, html_table_string)

    Plots are rendered in a responsive 2-column grid.
    """
    images = images or []
    tables = tables or []

    # --- Tables (full width) ---
    tbl_html = "\n".join(
        f"<h3>{cap}</h3>\n{tbl}"
        for cap, tbl in tables
    )

    # --- Images (2-column grid) ---
    if images:
        img_html = (
            '<h2>Plots</h2>\n'
            '<div class="plot-grid">\n'
            + "\n".join(
                f"""
<figure>
  <img src="{_img_to_data_uri(p)}" />
  <figcaption>{Path(p).name}</figcaption>
</figure>
"""
                for p in images
            )
            + "\n</div>"
        )
    else:
        img_html = ""

    body = f"""
<!doctype html>
<html>
<head>
<meta charset="utf-8"/>
<title>{title}</title>
<style>
body {{
  font-family: -apple-system, Segoe UI, Roboto, Helvetica, Arial, sans-serif;
  margin: 24px;
}}
h1, h2, h3 {{
  margin-top: 1.2em;
}}

code, pre {{
  background: #f6f8fa;
  padding: 2px 4px;
  border-radius: 4px;
}}

table {{
  border-collapse: collapse;
  width: 100%;
}}
th, td {{
  border: 1px solid #ddd;
  padding: 6px 8px;
}}
th {{
  background: #f3f3f3;
  text-align: left;
}}

/* ---- Plot grid ---- */
.plot-grid {{
  display: grid;
  grid-template-columns: 1fr 1fr;
  gap: 16px;
  margin-top: 1em;
}}
.plot-grid figure {{
  margin: 0;
}}
.plot-grid img {{
  width: 100%;
  height: auto;
  border: 1px solid #ddd;
  border-radius: 8px;
  display: block;
}}
.plot-grid figcaption {{
  color: #555;
  font-size: 0.9rem;
  margin-top: 4px;
  text-align: center;
}}

/* Responsive: 1 column on narrow screens */
@media (max-width: 900px) {{
  .plot-grid {{
    grid-template-columns: 1fr;
  }}
}}
</style>
</head>
<body>

<h1>{title}</h1>

{''.join(f'<p>{p}</p>' for p in paragraphs)}

{tbl_html}

{img_html}

{options_html()}

</body>
</html>
"""
    Path(out_html).write_text(body, encoding="utf-8")

