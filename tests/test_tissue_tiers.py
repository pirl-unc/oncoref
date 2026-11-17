"""Registry updates add tissues and can never change a base assessment."""

import hashlib
from pathlib import Path

import pytest

from oncoref import tissue_tiers

DATA = Path(__file__).resolve().parents[1] / "oncoref/data"


def test_updates_add_seven_cited_tissues_over_the_frozen_registry():
    added = tissue_tiers.updates()
    assert dict(zip(added.id, added.tier)) == {
        "olfactory_epithelium": "manageable",
        "taste_bud": "protected",
        "tooth": "manageable",
        "sensory_ganglion": "vital",
        "nasal_mucosa": "protected",
        "lacrimal_gland": "protected",
        "bone": "vital",
    }
    assert added.confidence.isin(["high", "medium", "low"]).all()
    ledger = tissue_tiers.update_references().set_index("reference_id")
    assert added.references.map(
        lambda keys: set(keys) <= set(ledger.index) and len(keys) >= 2
    ).all()
    assert ledger.url.str.startswith("https://").all()
    # Each tier set on the oncology precedent records that precedent and its sources.
    precedent = added[added.oncology_precedent.notna()]
    assert set(precedent.id) == {"olfactory_epithelium", "taste_bud"}
    assert "ctcae5_anosmia" in added.set_index("id").references["olfactory_epithelium"]
    # The base file and its identity are untouched; released receipts bind them.
    base = hashlib.sha256((DATA / "normal-tissue-tiers.yaml").read_bytes()).hexdigest()
    assert tissue_tiers.registry()["registry_id"].endswith(base)
    assert tissue_tiers.policy_id().startswith(tissue_tiers.registry()["registry_id"] + "+")
    assert not set(added.id) & set(tissue_tiers.registry_assessments().id)


@pytest.mark.parametrize(
    ("tissue", "tier"),
    [
        ("olfactory epithelium", "manageable"),
        ("Olfactory mucosa", "manageable"),
        ("circumvallate papillae", "protected"),
        ("dorsal root ganglion", "vital"),
        ("dental pulp", "manageable"),
        ("bone", "vital"),
    ],
)
def test_added_tissues_resolve_only_through_the_updates(tissue, tier):
    assert tissue_tiers.registry_assess(tissue)["tier"] == "unresolved"
    found = tissue_tiers.assess(tissue)
    assert found["tier"] == tier and found["reason"] == "registry_update"
    assert found["references"] and found["updates_id"]
    cell = tissue_tiers.assess(tissue, "some cell", resolution="cell")
    assert cell["tier"] == tier and cell["reason"] == "registry_update_organ_fallback"


def test_base_assessments_pass_through_unchanged():
    for tissue, cell in [
        ("lung", None),
        ("skin", "melanocytes"),
        ("bone marrow", None),
        ("nose", None),
    ]:
        options = {"resolution": "cell"} if cell else {}
        assert tissue_tiers.assess(tissue, cell, **options) == tissue_tiers.registry_assess(
            tissue, cell, **options
        )
    # Outside the adult human policy an added tissue stays unresolved.
    assert tissue_tiers.assess("bone", stage="fetal")["tier"] == "unresolved"


def test_public_tables_cover_the_base_registry_and_its_updates():
    rules = tissue_tiers.assessments()
    assert rules.source.value_counts().to_dict() == {"registry": 107, "update": 7}
    assert rules.id.is_unique
    ledger = tissue_tiers.references()
    assert ledger.reference_id.is_unique
    assert len(ledger) == len(tissue_tiers.registry_references()) + len(
        tissue_tiers.update_references()
    )
    assert tissue_tiers.TIERS == ("vital", "protected", "manageable", "disposable")


def test_loss_consequences_preserve_scope_and_link_reviewed_evidence():
    testis = tissue_tiers.consequences("testis")
    assert testis["tier"] == "protected"
    assert "androgen" in testis["persistent_loss"]
    assert "cannot restore fertility" in testis["management"]
    assert all(r["url"].startswith("https://") for r in testis["evidence"])
    unknown_cell = tissue_tiers.consequences("testis", "unreviewed cell", resolution="cell")
    assert unknown_cell["reason"] == "organ_fallback"
    assert unknown_cell["matched_level"] == "bulk"
    assert tissue_tiers.consequences("breast")["tier"] == "manageable"
    assert tissue_tiers.consequences("dental pulp")["loss_assessed"]
    embryo = tissue_tiers.consequences("embryo", stage="embryonic")
    assert not embryo["loss_assessed"] and "persistent_loss" not in embryo


def test_frozen_release_inputs_match_their_manifest():
    import hashlib
    import json

    folder = Path(__file__).parent / "fixtures/release-1.37-inputs"
    recorded = json.loads((folder / "manifest.json").read_text())["files"]
    for name, entry in recorded.items():
        assert hashlib.sha256((folder / name).read_bytes()).hexdigest() == entry["sha256"], name
