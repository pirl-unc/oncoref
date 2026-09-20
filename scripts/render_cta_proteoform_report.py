#!/usr/bin/env python3
"""Render verified current selections and every indexed figure."""

import argparse
from pathlib import Path

from cta_report_common import deterministic_zip
from cta_report_render import render_report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--out",
        type=Path,
        default=Path(__file__).resolve().parents[1] / "outputs/cta_proteoform_report_20260917",
    )
    out = parser.parse_args().out.resolve()
    index = render_report(out)
    deterministic_zip(out, out / "cta-analysis-bundle.zip")
    print(f"Rendered {len(index)} pages")


if __name__ == "__main__":
    main()
