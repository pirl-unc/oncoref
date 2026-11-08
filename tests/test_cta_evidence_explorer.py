"""Presentation data must preserve released identities, values and evidence limits."""

import gzip
import hashlib
import importlib.util
import json
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "docs/audits/cta-evidence-guide-20261106"


@pytest.fixture(scope="module")
def report():
    return json.loads(gzip.decompress((OUT / "data.json.gz").read_bytes()))


def test_frozen_131_proteins_and_exact_122_complement_are_preserved(report):
    assets = ROOT / "scripts/cta_evidence_assets"
    current = pd.read_csv(assets / "cta-selection-1.31.0.csv.gz").set_index("entity_id")
    current["gene_ids"] = current.gene_ids.map(json.loads)
    paper_entities = set(pd.read_csv(assets / "historical-paper-memberships.csv.gz").entity_id)
    rows = {r["entity_id"]: r for r in report["rows"]}
    assert set(rows) == set(current.index)
    for entity, row in rows.items():
        assert row["classification"] == current.loc[entity, "classification"]
        assert row["gene_ids"] == current.loc[entity, "gene_ids"]
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


def test_additional_catalog_matches_preserve_locus_and_source_scope(report):
    unlisted = {
        r["entity_id"]: r for r in report["rows"] if r["category"] == "CTA" and not r["papers"]
    }
    catalogs = {c["id"]: c for c in report["catalog_review"]["catalogs"]}
    expected = {
        "chew2019": {"DUXA", "DUXB", "LEUTX", "PRAMEF25/6"},
        "smith2023": {
            "C19orf85",
            "CPHXL",
            "DUXA",
            "LEUTX",
            "MBD3L2B",
            "OR10G3",
            "TAF11L12",
            "ZNF280A",
            "ZNF878",
        },
    }
    for name, catalog in catalogs.items():
        assert {m["protein"] for m in catalog["matches"]} == expected[name]
        assert catalog["matched_proteins"] == len(expected[name])
        for match in catalog["matches"]:
            protein = unlisted[match["entity_id"]]
            assert match["gene_id"] in protein["gene_ids"]
            assert match["row"] >= 4
            assert match["sheet"]
            if match["match_method"] == "Ensembl ID":
                assert match["source_gene_id"] == match["gene_id"]
            else:
                assert match["source_symbol"] in {g["symbol"] for g in protein["gene_identities"]}
    shared = next(m for m in catalogs["chew2019"]["matches"] if m["protein"] == "PRAMEF25/6")
    assert shared["gene_id"] == "ENSG00000229571"  # PRAMEF25 only, not every shared locus
    assert catalogs["chew2019"]["screen_table"]["overlapping_gene_rows"] == 74
    assert len(set.union(*expected.values())) == 11


def test_biology_categories_retain_evidence_and_overlap(report):
    ctas = [r for r in report["rows"] if r["category"] == "CTA"]
    assert len(ctas) == 602
    categories = {c["id"] for c in report["biology"]["categories"]}
    sources = {s["source_id"] for s in report["sources"]}
    for row in report["rows"]:
        assert bool(row["biology_categories"]) == (row["category"] == "CTA")
        assert set(row["biology_categories"]) <= categories
        evidence_categories = {e["category_id"] for e in row["biology_evidence"]}
        assert evidence_categories == set(row["biology_categories"]) - {"uncategorized"}
        for evidence in row["biology_evidence"]:
            assert evidence["basis"] and evidence["model"]
            assert (
                evidence.get("source_url", "").startswith("https://")
                or evidence.get("source_id") in sources
            )
            assert set(evidence["gene_ids"]) & set(row["gene_ids"])
    for category in report["biology_counts"]:
        assert category["all"] == sum(category["id"] in r["biology_categories"] for r in ctas)
        assert category["unlisted"] == sum(
            category["id"] in r["biology_categories"] for r in ctas if not r["papers"]
        )
        assert category["registered"] + category["unlisted"] == category["all"]
    by_label = {r["label"]: r for r in ctas}
    assert by_label["WEE2"]["route"] == "Developmental only"
    assert "oocyte_maternal" in by_label["WEE2"]["biology_categories"]
    for gene in ["PATE1", "WFDC8", "RNASE9", "SEMG2"]:
        assert "reproductive_secretion" in by_label[gene]["biology_categories"]
        assert any(
            e["kind"] == "Publication" and e["category_id"] == "reproductive_secretion"
            for e in by_label[gene]["biology_evidence"]
        )
    assert any("Mouse" in e["model"] for e in by_label["PABPN1L"]["biology_evidence"])
    exported = pd.read_csv(OUT / "cta-biology-categories.csv")
    assert set(exported.entity_id) == {r["entity_id"] for r in ctas}


def test_biology_rna_uses_all_contexts_before_the_top_five_cut():
    spec = importlib.util.spec_from_file_location(
        "explorer_biology_builder", ROOT / "scripts/build_cta_evidence_explorer.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    collector = module.BiologyCollector(pd.DataFrame({"representative_gene_id": ["g"]}))
    frame = pd.DataFrame(
        [
            {
                "gene_id": "g",
                "source": "hpa_single_cell_tissue",
                "unit": "nTPM",
                "tissue": "testis",
                "cell_type": f"cell {i}",
                "value": 100 - i,
            }
            for i in range(6)
        ]
        + [
            {
                "gene_id": "g",
                "source": "hpa_single_cell_tissue",
                "unit": "nTPM",
                "tissue": "ovary",
                "cell_type": "oocytes",
                "value": 6,
            },
            {
                "gene_id": "g",
                "source": "hpa_single_cell_tissue",
                "unit": "nTPM",
                "tissue": "ovary",
                "cell_type": "oocytes",
                "value": None,
                "raw_maximum": 900,
            },
        ]
    )
    collector.collect(frame, "HPA tissue clusters")
    selected = pd.concat(collector.parts)
    assert len(selected) == 7
    assert selected.loc[selected.cell_type.eq("oocytes"), "value"].tolist() == [6]


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


def test_explorer_manifest_binds_outputs_and_source_artifacts(historical_cta_input):
    manifest = json.loads((OUT / "manifest.json").read_text())
    for category, base in [("inputs", ROOT), ("outputs", OUT)]:
        for name, digest in manifest[category].items():
            path = historical_cta_input(name, digest) if category == "inputs" else base / name
            assert hashlib.sha256(path.read_bytes()).hexdigest() == digest, name


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
