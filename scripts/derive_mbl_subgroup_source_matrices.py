# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Derive four MBL molecular-subgroup source matrices from the parent MBL matrix.

The historical pirlygenes subgroup approximation assigns each sample to the
largest-TPM marker among WIF1 (WNT), GLI2 (SHH), MYC (Group 3), and KCNA1
(Group 4). The classifier implementation lives in ``oncoref.expression_builders``
so the release assets and regression tests share exactly one rule.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from oncoref import source_matrices
from oncoref.expression_builders import (
    _write_parquet_atomic,
    medulloblastoma_subgroup_matrices,
)


def derive(
    parent_path: Path,
    *,
    output_dir: Path,
    release_dir: Path | None = None,
) -> dict[str, Path]:
    """Write subgroup matrices and return their cache paths by cancer code."""
    parent = pd.read_parquet(parent_path)
    matrices = medulloblastoma_subgroup_matrices(parent)
    output_dir.mkdir(parents=True, exist_ok=True)
    if release_dir is not None:
        release_dir.mkdir(parents=True, exist_ok=True)

    paths: dict[str, Path] = {}
    for code, matrix in matrices.items():
        path = output_dir / f"{code}.parquet"
        _write_parquet_atomic(matrix, path)
        if release_dir is not None:
            _write_parquet_atomic(matrix, release_dir / f"{code}_per_sample_tpm.parquet")
        paths[code] = path
        print(f"{code}: {len(matrix.columns) - 2} samples -> {path}", flush=True)
    return paths


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--parent",
        type=Path,
        default=None,
        help="Parent MBL parquet (default: source_matrices.ensure('MBL'))",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Cache output directory (default: active source-matrix cache)",
    )
    parser.add_argument(
        "--release-dir",
        type=Path,
        default=None,
        help="Also write GitHub release asset filenames to this directory",
    )
    args = parser.parse_args(argv)

    parent_path = args.parent or source_matrices.ensure("MBL")
    output_dir = args.output_dir or source_matrices.cache_dir()
    derive(
        parent_path.expanduser(),
        output_dir=output_dir.expanduser(),
        release_dir=args.release_dir.expanduser() if args.release_dir else None,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
