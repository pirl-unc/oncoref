# Licensed under the Apache License, Version 2.0.
"""Replaying a source snapshot must not require a prepopulated destination."""

import hashlib
import json
import subprocess
import sys
from importlib.resources import files
from pathlib import Path

import build_canine_source_manifests as builder
import yaml


def test_replay_to_fresh_directory_uses_packaged_reference_seed(tmp_path):
    root = tmp_path / "snapshot"
    docs = root / "docs"
    docs.mkdir(parents=True)

    def source(relative, content):
        path = root / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content)
        return {
            "asset_id": path.name,
            "path": relative,
            "url": f"https://example.org/{path.name}",
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "size_bytes": path.stat().st_size,
            "source_version": "fixture-v1",
            "license": "Synthetic fixture",
        }

    def receipt(name, data):
        (docs / name).write_text(json.dumps(data))

    receipt("bgee-quantitative-assets.json", {"assets": []})
    bark = source(
        "data/processed/barkbase-tmm-gene.tsv",
        "gene\tname\tlength\tAdult_1_heart\nG1\tGENE1\t100\t2.0\n",
    )
    receipt("source-registry.json", {"assets": [bark]})
    source("data/report-reconciliation/barkbase-paper.xml", "<article>Fixture</article>")
    source(
        "data/inventory-2026-10-08/Canis_lupus_familiaris_RNA-Seq_libraries.tsv", "Run IDs\nR1\n"
    )
    source(
        "data/inventory-2026-10-08/barkbase-ena",
        "run_accession\ttax_id\tlibrary_strategy\nR1\t9615\tRNA-Seq\n",
    )
    protein = source(
        "data/report-reconciliation/dog-protein-abundance.tsv",
        "Run\tChannel\tBioReplicate\tCondition\tTechRepMixture\tMixture\tsample.type\n"
        "run1\t1\tpool1\theart\t1\t1\ttissue\n",
    )
    fasta = source("data/report-reconciliation/dog-protein-search.fasta", ">P1\nMPEPTIDE\n")
    receipt("report-protein-data-receipts.json", [protein, fasta])
    receipt(
        "real-receipts.json",
        [
            source(f"data/real/{name}", "Synthetic metadata\n")
            for name in ("doga-data.zip", "GSE325523.txt", "GSE304066.txt")
        ],
    )
    seed = files("oncoref").joinpath("data/species/dog/assets.json")
    seed_bytes = seed.read_bytes()
    seed_manifest = json.loads(seed_bytes)
    inputs_before = {path: path.read_bytes() for path in root.rglob("*") if path.is_file()}
    output = tmp_path / "new" / "manifests"
    assert not output.exists()

    # Run the public CLI away from the repository so neither cwd nor output
    # can accidentally supply the seed.
    result = subprocess.run(
        [
            sys.executable,
            str(Path(builder.__file__).resolve()),
            "--source-root",
            str(root),
            "--output",
            str(output),
        ],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        check=True,
    )
    assert len(json.loads(result.stdout)["audit"]) == 2
    generated = json.loads((output / "assets.json").read_text())
    assert generated["references"] == seed_manifest["references"]
    expected_references = [
        asset
        for asset in seed_manifest["assets"]
        if asset["group"] != "expression" and asset["id"] != "pfizer-dog-search-fasta"
    ]
    assert generated["assets"][: len(expected_references)] == expected_references
    assert (
        next(a for a in generated["assets"] if a["id"] == "barkbase-matrix")["sha256"]
        == bark["sha256"]
    )
    sources = yaml.safe_load((output / "expression_sources.yaml").read_text())["sources"]
    assert len(sources) == 5
    assert json.loads((output / "lineage_overlaps.json").read_text())["shared_run_ids"] == ["R1"]
    assert (output / "barkbase-samples.json").is_file()
    assert (output / "pfizer-dog-protein-samples.json").is_file()

    # Existing generated output is not an input: rebuilding replaces even a
    # broken destination manifest without mutating the package or snapshot.
    generated_bytes = {path.name: path.read_bytes() for path in output.iterdir()}
    (output / "assets.json").write_text("not a seed")
    builder.build(root, output)
    assert {path.name: path.read_bytes() for path in output.iterdir()} == generated_bytes
    assert seed.read_bytes() == seed_bytes
    assert {path: path.read_bytes() for path in root.rglob("*") if path.is_file()} == inputs_before
