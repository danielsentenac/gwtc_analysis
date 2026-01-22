from __future__ import annotations
import base64
from pathlib import Path
from typing import List, Optional

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

</body>
</html>
"""
    Path(out_html).write_text(body, encoding="utf-8")

