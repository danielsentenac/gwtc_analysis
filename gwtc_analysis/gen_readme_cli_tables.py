from __future__ import annotations

import argparse
import sys
from pathlib import Path

# --- Ensure package root is on sys.path ---
HERE = Path(__file__).resolve()
PKG_ROOT = HERE.parent.parent   # repo_root/
if str(PKG_ROOT) not in sys.path:
    sys.path.insert(0, str(PKG_ROOT))

import gwtc_analysis.cli as cli


CATALOGS_START = "<!-- CATALOG_COVERAGE_BEGIN -->"
CATALOGS_END = "<!-- CATALOG_COVERAGE_END -->"

START = "<!-- CLI_TABLES_BEGIN -->"
END = "<!-- CLI_TABLES_END -->"

def _parser_to_md_tables(parser: argparse.ArgumentParser) -> str:
    lines: list[str] = []
    sub_action = next(a for a in parser._actions if isinstance(a, argparse._SubParsersAction))
    for mode, sub in sub_action.choices.items():
        lines.append(f"### `{mode}`")
        lines.append("")
        lines.append("| Option | Default | Description |")
        lines.append("|---|---:|---|")
        for a in sub._actions:
            if not a.option_strings:
                continue
            if a.help is argparse.SUPPRESS:
                continue
            opt = ", ".join(a.option_strings)
            default = a.default
            if default is None or default is argparse.SUPPRESS:
                default_s = ""
            elif default is False and isinstance(a, argparse._StoreTrueAction):
                default_s = "False"
            else:
                default_s = str(default)
            help_ = (a.help or "").strip().replace("\n", " ")
            help_ = help_.replace("|", "\\|")
            lines.append(f"| `{opt}` | `{default_s}` | {help_} |")
        lines.append("")
    return "\n".join(lines).strip() + "\n"

def _fill_coverage(txt: str) -> str:
    """The block naming the catalogs of this version, generated from the catalog registry."""
    if CATALOGS_START not in txt:
        return txt
    from gwtc_analysis import __version__
    from gwtc_analysis import catalog_registry as reg
    before, rest = txt.split(CATALOGS_START, 1)
    after = rest.split(CATALOGS_END, 1)[1]
    return (before + CATALOGS_START + "\n> " + reg.coverage_text(__version__, markdown=True) + "\n"
            + CATALOGS_END + after)


def main() -> None:
    parser = cli.build_parser()
    tables = _parser_to_md_tables(parser)

    readme_path = Path("README.md")
    txt = readme_path.read_text(encoding="utf-8")
    if START not in txt or END not in txt:
        raise RuntimeError("README.md missing CLI table markers")

    before = txt.split(START)[0]
    after = txt.split(END)[1]
    new_txt = before + START + "\n\n" + tables + "\n" + END + after
    new_txt = _fill_coverage(new_txt)
    readme_path.write_text(new_txt, encoding="utf-8")
    print("✔ README.md updated from cli.py")

    docs_path = Path("docs") / "cli-reference.md"
    if docs_path.parent.is_dir():
        docs_path.write_text(
            "# CLI reference\n\n"
            "Every option of every mode, generated from `gwtc_analysis/cli.py` by\n"
            "`python gwtc_analysis/gen_readme_cli_tables.py`. Each mode also has its own help:\n"
            "`gwtc_analysis <MODE> -h`.\n\n"
            + tables.replace("### `", "## `"),
            encoding="utf-8",
        )
        print("✔ docs/cli-reference.md updated from cli.py")

    index_path = Path("docs") / "index.md"
    if index_path.is_file():
        index_path.write_text(_fill_coverage(index_path.read_text(encoding="utf-8")), encoding="utf-8")
        print("✔ docs/index.md catalog coverage updated from the registry")

if __name__ == "__main__":
    main()
