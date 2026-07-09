#!/usr/bin/env python
"""Build a registry-backed Treehouse compendium source matrix.

This script is intentionally thin: Treehouse clinical routing, log2(TPM+1)
inverse transform, symbol canonicalization, parse diagnostics, sample QC, and
summary-row production live in :mod:`oncoref.expression_builders`.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from oncoref.expression_builders import (
    build_treehouse_source_matrices,
    treehouse_source_from_registry,
)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source_id_arg", nargs="?", help="Source id from expression_sources.yaml")
    parser.add_argument("--source-id", help="Source id from expression_sources.yaml")
    parser.add_argument(
        "--registry",
        type=Path,
        default=None,
        help="Optional expression_sources.yaml path; defaults to the packaged registry.",
    )
    parser.add_argument(
        "--cache-dir",
        type=Path,
        default=None,
        help="Treehouse cache directory. Defaults to ~/.cache/oncoref/expression/<source-id>/.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Directory for per-code parquet matrices and sidecars. Defaults to cache_dir/derived.",
    )
    parser.add_argument(
        "--tpm-path",
        type=Path,
        default=None,
        help="Use an already-downloaded Treehouse TPM matrix instead of registry/cache lookup.",
    )
    parser.add_argument(
        "--clinical-path",
        type=Path,
        default=None,
        help="Use an already-downloaded Treehouse clinical table instead of registry/cache lookup.",
    )
    parser.add_argument("--force-download", action="store_true")
    parser.add_argument("--high-expression-threshold", type=float, default=1.0)
    args = parser.parse_args(argv)

    source_id = args.source_id or args.source_id_arg
    if source_id is None:
        raise SystemExit("provide --source-id <id>")

    source = treehouse_source_from_registry(source_id, registry_path=args.registry)
    cache_dir = args.cache_dir or Path.home() / ".cache" / "oncoref" / "expression" / source_id
    result = build_treehouse_source_matrices(
        source,
        cache_dir=cache_dir,
        output_dir=args.output_dir,
        tpm_path=args.tpm_path,
        clinical_path=args.clinical_path,
        force_download=args.force_download,
        high_expression_threshold=args.high_expression_threshold,
    )
    summary = {
        "source_id": source_id,
        "source_cohort": source.source_cohort,
        "matrix_paths": {code: str(path) for code, path in result.matrix_paths.items()},
        "sidecar_paths": {name: str(path) for name, path in result.sidecar_paths.items()},
        "sample_counts": {
            code: len([c for c in matrix.columns if c not in {"Ensembl_Gene_ID", "Symbol"}])
            for code, matrix in result.matrices.items()
        },
    }
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
