#!/usr/bin/env python3
"""Export complete source mappings and species applicability, without imputing coverage."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from oncoref.normal_tissue_gene_filters import normal_tissue_gene_filter_profiles  # noqa: E402
from oncoref.normal_tissues import (  # noqa: E402
    normal_tissue_inventory,
    normal_tissue_mapping_inventory,
    normal_tissue_taxonomy,
)


def write_csv(path, rows):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    k: json.dumps(v) if isinstance(v, (tuple, list, dict)) else v
                    for k, v in row.items()
                }
            )


def build(out):
    if out.exists() and any(out.iterdir()):
        raise ValueError("Use a new empty output directory")
    out.mkdir(parents=True, exist_ok=True)
    concepts, mappings, counts = [], [], []
    for species in ("human", "dog", "mouse"):
        anatomy = normal_tissue_inventory(species=species)
        native = normal_tissue_mapping_inventory(species=species)
        concepts.extend({"species": species, **r} for r in anatomy)
        mappings.extend({"species": species, **r} for r in native)
        categories = Counter(
            r["tissue_category"] for r in anatomy if r["mapping_status"] == "mapped"
        )
        counts.append(
            {
                "species": species,
                **{
                    name: categories[name]
                    for name in (
                        "reproductive",
                        "disposable",
                        "manageable",
                        "protected",
                        "vital",
                        "unresolved",
                    )
                },
                "inapplicable": sum(
                    r["mapping_status"] == "not_applicable_to_species" for r in anatomy
                ),
                "native_mapping_rows": len(native),
                "native_sources": len({(r["source_id"], r["source_version"]) for r in native}),
                "mapping_availability": "bundled" if native else "not_bundled",
            }
        )
    write_csv(out / "canonical-anatomy.csv", concepts)
    write_csv(out / "source-mappings.csv", mappings)
    write_csv(out / "species-summary.csv", counts)
    taxonomy = normal_tissue_taxonomy()
    for name, obj in (
        ("taxonomy.json", taxonomy),
        ("gene-filter-profiles.json", normal_tissue_gene_filter_profiles()),
    ):
        (out / name).write_text(json.dumps(obj, indent=2) + "\n")
    table = "\n".join(
        "| " + " | ".join(map(str, row)) + " |"
        for row in [
            [{"manageable": "Manageable", "protected": "Protected"}.get(k, k) for k in counts[0]],
            ["---"] * len(counts[0]),
            *[row.values() for row in counts],
        ]
    )
    (out / "index.md").write_text(
        "# Complete normal-tissue registry — 20 October 2026\n\n"
        f"Policy `{taxonomy['policy_version']}`. Every canonical concept is listed for "
        "human, dog and mouse, including explicitly inapplicable concepts. Every "
        "bundled native source/version/tissue/cell mapping is exported separately.\n\n"
        + table
        + "\n\n"
        "Counts are concepts and mapping rows, **not samples or measured negative coverage**. "
        "The same source tissue can have several cell-type rows. Mouse source mappings are "
        "not bundled; no mouse expression data or absence calls are invented.\n\n"
        "## Exact inventories\n\n"
        "- [All species and canonical concepts](canonical-anatomy.csv)\n"
        "- [Every native dataset mapping, category and citation](source-mappings.csv)\n"
        "- [Taxonomy and source-artifact receipts](taxonomy.json)\n"
        "- [Patient-context gene-filter definitions](gene-filter-profiles.json)\n"
        "- [Checksums and implementation provenance](run-manifest.json)\n\n"
        "## Interpretation\n\n"
        "Reproductive, Disposable, Manageable, Protected and Vital are the five resolved categories. "
        "Unresolved is a QC state. Disposable denotes the narrowed conditional project "
        "allowance, whose compatibility scope ID is `permitted_normal`.\n\n"
        "Uterine horns are applicable to dog/mouse; coagulating gland is a mouse prostate "
        "descendant. Dog seminal vesicles/bulbourethral glands, mouse tonsil, and dog/mouse "
        "vermiform appendix are explicitly inapplicable. Applicability is checked through "
        "anatomical parents. Broad head, foot and soft-tissue labels stay unresolved.\n\n"
        "Native mapping rows do not supply sex or developmental metadata that the source "
        "does not report. At sample annotation time, explicit fetal/embryonic context "
        "overrides adult organ categories, while pregnancy of the mother does not. "
        "Sex and physiological state are preserved separately. Gene filters require a "
        "declared negative panel and cannot pool species, sources, versions or assays.\n\n"
        "See [definitions](../../normal-tissue-scopes.md), "
        "[manageability evidence](../../normal-tissue-manageability.md), and "
        "[patient-context gene filters](../../normal-tissue-gene-filters.md).\n\n"
        "Reproduce with `python scripts/normal_tissue_registry_report.py --out <empty-directory>` "
        "at this report's implementation revision.\n"
    )
    implementation = (
        "scripts/normal_tissue_registry_report.py",
        "oncoref/normal_tissues.py",
        "oncoref/normal_tissue_gene_filters.py",
        "oncoref/data/normal-tissue-taxonomy.yaml",
        "oncoref/data/normal-tissue-mappings.csv",
    )

    def digest(path):
        return hashlib.sha256(path.read_bytes()).hexdigest()

    (out / "run-manifest.json").write_text(
        json.dumps(
            {
                "schema": "oncoref.normal-tissue-registry-report.v1",
                "taxonomy_id": taxonomy["taxonomy_id"],
                "implementation_sha256": {p: digest(ROOT / p) for p in implementation},
                "outputs": {p.name: digest(p) for p in sorted(out.iterdir()) if p.is_file()},
            },
            indent=2,
        )
        + "\n"
    )
    print(table)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    build(parser.parse_args().out)
