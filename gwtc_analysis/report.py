from __future__ import annotations
import base64
from pathlib import Path
from typing import List, Optional

def _img_to_data_uri(path: str | Path) -> str:
    p = Path(path)
    data = p.read_bytes()
    # naive mime handling
    suffix = p.suffix.lower().lstrip(".")
    mime = "image/png" if suffix == "png" else "image/jpeg" if suffix in ("jpg","jpeg") else "application/octet-stream"
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
    - paragraphs: list of HTML-escaped or simple text paragraphs (we don't escape aggressively here).
    - images: list of image paths; embedded as data URIs.
    - tables: list of (caption, html_table_string)
    """
    images = images or []
    tables = tables or []
    img_html = "\n".join(
        f'<figure style="margin: 1.2rem 0;"><img style="max-width: 100%; height: auto; border: 1px solid #ddd; border-radius: 8px;" src="{_img_to_data_uri(p)}"/><figcaption style="color:#555;font-size:0.9rem;">{Path(p).name}</figcaption></figure>'
        for p in images
    )
    tbl_html = "\n".join(
        f"<h3>{cap}</h3>\n{tbl}"
        for cap, tbl in tables
    )
    body = f"""
<!doctype html>
<html>
<head>
<meta charset="utf-8"/>
<title>{title}</title>
<style>
body {{ font-family: -apple-system, Segoe UI, Roboto, Helvetica, Arial, sans-serif; margin: 24px; }}
h1, h2, h3 {{ margin-top: 1.2em; }}
code, pre {{ background: #f6f8fa; padding: 2px 4px; border-radius: 4px; }}
table {{ border-collapse: collapse; width: 100%; }}
th, td {{ border: 1px solid #ddd; padding: 6px 8px; }}
th {{ background: #f3f3f3; text-align: left; }}
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
