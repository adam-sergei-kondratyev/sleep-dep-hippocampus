"""Fill README.template.md with numbers from results/summary.json and write README.md.

Usage: python -m src.render_readme
Fails if the template names a key that summary.json lacks, so README numbers can never be hand-typed.
"""
from __future__ import annotations

import json
import re
import sys

from src.common import ROOT, SUMMARY_JSON

TEMPLATE = ROOT / "README.template.md"
OUTPUT = ROOT / "README.md"
PLACEHOLDER = re.compile(r"\{\{\s*([a-z0-9_]+)\s*\}\}")


def fmt(value: object) -> str:
    if isinstance(value, bool):
        return str(value)
    if isinstance(value, int):
        return f"{value:,}"
    return str(value)


def render(template: str, summary: dict[str, object]) -> str:
    unknown = sorted({k for k in PLACEHOLDER.findall(template) if k not in summary})
    if unknown:
        raise KeyError(f"template keys missing from summary.json: {unknown}")
    return PLACEHOLDER.sub(lambda m: fmt(summary[m.group(1)]), template)


def main() -> int:
    if not SUMMARY_JSON.is_file():
        print("results/summary.json missing; run `make compare` first", file=sys.stderr)
        return 1
    summary = json.loads(SUMMARY_JSON.read_text())
    OUTPUT.write_text(render(TEMPLATE.read_text(), summary))
    print(f"wrote README.md from {TEMPLATE.name}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
