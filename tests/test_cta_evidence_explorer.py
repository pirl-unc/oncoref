"""Presentation data must preserve released identities, values and evidence limits."""

import gzip
import hashlib
import importlib.util
import json
from pathlib import Path

import pandas as pd
import pytest

from oncoref import cta_selection

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "docs/audits/cta-evidence-guide-20261106"


@pytest.fixture(scope="module")
def report():
    return json.loads(gzip.decompress((OUT / "data.json.gz").read_bytes()))


def test_every_released_protein_and_exact_122_complement_are_preserved(report):
    current = cta_selection.profiles().set_index("entity_id")
    paper_entities = set(cta_selection.audit().memberships.entity_id)
    rows = {r["entity_id"]: r for r in report["rows"]}
    assert set(rows) == set(current.index)
    for entity, row in rows.items():
        assert row["classification"] == current.loc[entity, "classification"]
        assert tuple(row["gene_ids"]) == current.loc[entity, "gene_ids"]
        assert bool(row["papers"]) == (entity in paper_entities)
    subset = [r for r in rows.values() if r["category"] == "CTA" and not r["papers"]]
    assert len(subset) == 122
    assert sum(len(r["gene_ids"]) for r in subset) == 123
    assert {r["entity_id"] for r in subset} == {
        r["entity_id"] for r in json.loads((OUT / "unlisted-ctas.json").read_text())
    }
    assert report["comparison"]["unlisted"]["routes"] == {
        "Developmental only": 60,
        "Adult reproductive only": 49,
        "Adult reproductive + developmental": 13,
    }


def test_source_values_do_not_turn_missing_observations_into_zero(report):
    sources = {s["source_id"] for s in report["sources"]}
    assert set(report["columns"]) <= sources
    saw_zero = saw_missing = False
    for row in report["rows"]:
        assert set(row["values"]) == set(report["columns"])
        for source in report["numeric_columns"]:
            value = row["values"][source]
            assert value is None or isinstance(value, (float, int))
            saw_zero |= value == 0
            saw_missing |= value is None
        assert row["values"]["oncoref_tumor_union"] == row["tumor_max_tpm"]
        assert row["values"]["hpa_cancer"] == row["tumor_max_ptpm"]
    assert saw_zero and saw_missing


def test_top_contexts_are_unique_within_source_unit_and_scope():
    d = pd.read_csv(OUT / "normal-top-contexts.csv.gz", low_memory=False)
    keys = ["gene_id", "source", "unit", "scope"]
    assert not d.duplicated([*keys, "tissue", "cell_type"]).any()
    assert d.groupby(keys).size().max() == 5
    # Verify actual ranked contexts against a complete frozen native source,
    # including a shared protein from the user's concrete example.
    raw = pd.read_csv(
        ROOT / "docs/audits/cta-selection-20261028/alternative-evidence/normal-bulk.csv.gz",
        low_memory=False,
    )
    facts = pd.read_csv(
        ROOT / "docs/audits/cta-selection-20261028/protein-group-classifications.csv.gz",
        low_memory=False,
    ).set_index("protein_group_label")
    for label in ["CT45A5/6/7", "ENPP4", "MRGPRX1"]:
        gene = facts.loc[label, "representative_gene_id"]
        expected = raw[raw.gene_id.eq(gene) & raw.source.eq("hpa_brain")]
        expected = expected.sort_values(["value", "tissue"], ascending=[False, True]).head(5)
        actual = d[d.gene_id.eq(gene) & d.source.eq("hpa_brain")]
        assert list(zip(actual.tissue, actual.value, strict=True)) == list(
            zip(expected.tissue, expected.value, strict=True)
        )


def test_literature_annotations_do_not_pretend_all_122_were_reviewed(report):
    subset = {r["label"]: r for r in report["rows"] if r["category"] == "CTA" and not r["papers"]}
    reviewed = set()
    for finding in report["literature"]["findings"]:
        assert set(finding["proteins"]) <= set(subset)
        assert finding["url"].startswith("https://")
        assert finding["model"] and finding["limitation"] and finding["next_action"]
        reviewed.update(finding["proteins"])
    assert 0 < len(reviewed) < 122
    for label, row in subset.items():
        assert bool(row["literature_review"]) == (label in reviewed)
        assert row["values"]["gene-knowledge"] is None
    assert subset["MRGPRX1"]["classification"] == "CTA"
    assert "conflict" in subset["MRGPRX1"]["literature_review"][0]["category"]


def test_explorer_manifest_binds_outputs_and_source_artifacts():
    manifest = json.loads((OUT / "manifest.json").read_text())
    for category, base in [("inputs", ROOT), ("outputs", OUT)]:
        for name, digest in manifest[category].items():
            assert hashlib.sha256((base / name).read_bytes()).hexdigest() == digest, name


def test_top_five_reduction_is_valid_across_chunks_and_retains_raw_unsupported_values():
    spec = importlib.util.spec_from_file_location(
        "explorer_builder", ROOT / "scripts/build_cta_evidence_explorer.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    rows = pd.DataFrame(
        [
            {
                "gene_id": "g",
                "source": source,
                "unit": unit,
                "scope": "normal",
                "tissue": str(i),
                "cell_type": "bulk",
                "value": value,
                "raw_value": raw,
            }
            for source, unit in [("a", "TPM"), ("b", "nTPM")]
            for i, (value, raw) in enumerate(
                [(2, 2), (0, 0), (None, 900), (4, 4), (6, 6), (8, 8), (10, 10)]
            )
        ]
    )
    expected = module.top_contexts(rows)
    chunked = module.top_contexts(
        pd.concat([module.top_contexts(rows.iloc[:4]), module.top_contexts(rows.iloc[4:])])
    )
    pd.testing.assert_frame_equal(expected, chunked)
    assert set(expected.value) == {2, 4, 6, 8, 10}
    # Unsupported raw RNA cannot outrank supported expression.
    assert 900 not in set(expected.raw_value)
