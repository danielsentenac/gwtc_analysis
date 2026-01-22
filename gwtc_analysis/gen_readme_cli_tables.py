from __future__ import annotations

import argparse
from pathlib import Path

import gwtc_analysis.cli as cli

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
    readme_path.write_text(new_txt, encoding="utf-8")
    print("✔ README.md updated from cli.py")

if __name__ == "__main__":
    main()
