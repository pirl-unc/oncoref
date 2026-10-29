"""Reproduce a small ORFdata regression dataset; no network or full vendor input."""

import gzip
import hashlib
import json
import sys
import tempfile
from pathlib import Path
from types import SimpleNamespace

import pytest


def main():
    root = Path(sys.argv[1]).resolve()
    if not (root / "tests/test_pfizer_native.py").is_file():
        raise SystemExit("Provide the pinned ORFdata checkout")
    sys.path.insert(0, str(root / "tests"))
    sys.path.insert(0, str(root))
    from orfdata.sources.pfizer_native import import_native
    from test_pfizer_evidence import source
    from test_pfizer_native import inputs

    out = Path(__file__).resolve().parent
    out.mkdir(parents=True, exist_ok=True)
    with (
        tempfile.TemporaryDirectory(prefix="peptide-fixture-") as td,
        pytest.MonkeyPatch.context() as patch,
    ):
        tmp = Path(td)
        slice_input = source.__wrapped__(tmp, patch, SimpleNamespace(param="synthetic-memberships"))
        paths, manifest, _ = inputs.__wrapped__(slice_input, tmp)
        dataset = import_native(paths, tmp / "dataset.sqlite", manifest)
        raw = dataset.path.read_bytes()
        # Paths inside this fixture name its temporary source inputs; there is no acquisition.
        with (
            (out / "dataset.sqlite.gz").open("wb") as handle,
            gzip.GzipFile(fileobj=handle, mode="wb", filename="", mtime=0) as stream,
        ):
            stream.write(raw)
        (out / "protein.tsv").write_bytes(paths["protein.tsv"].read_bytes())
        result = {
            "producer": "openvax/orfdata",
            "commit": "985a93e59c34b1c59dee76ac52b044af4a1c32c8",
            "construction": "test_pfizer_evidence.source(synthetic-memberships), test_pfizer_native.inputs, public import_native",
            "dataset_sha256": hashlib.sha256(raw).hexdigest(),
            "dataset_size": len(raw),
            "compressed_sha256": hashlib.sha256(
                (out / "dataset.sqlite.gz").read_bytes()
            ).hexdigest(),
            "counts": dataset.counts,
            "scope": "Audited factual native slices plus explicitly synthetic summary quantities, acquisition ordering and shared/dangling membership links; not the released full source",
        }
        (out / "receipt.json").write_text(json.dumps(result, indent=2) + "\n")
        (out / "PRODUCER-NOTICE.md").write_bytes(
            (root / "tests/data/pfizer/vendor/NOTICE.md").read_bytes()
        )
        print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
