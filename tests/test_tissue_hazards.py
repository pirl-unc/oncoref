from copy import deepcopy

import pandas as pd
import pytest

from oncoref import tissue_hazards as h
from oncoref.normal_tissues import normal_tissue_screening_tier, resolve_normal_tissue


@pytest.mark.parametrize(
    "tissue,cell,tier",
    [
        ("lung", "club cell", "vital"),
        ("lung", "respiratory ciliated cell", "protected"),
        ("lung", "goblet cells", "protected"),
        ("lung", "basal respiratory cells", "vital"),
        ("lung", "type II pneumocyte", "vital"),
        ("lung", "alveolar macrophage", "protected"),
        ("colon", "undifferentiated cells", "vital"),
        ("colon", "paneth cells", "protected"),
        ("kidney", "podocyte", "vital"),
        ("pancreas", "beta cells", "protected"),
        ("thyroid gland", "follicular cells", "manageable"),
        ("bone marrow", "mesenchymal stromal cells", "vital"),
        ("heart", "fibroblasts", "protected"),
        ("blood", "B cell", "manageable"),
        ("skin", "melanocytes", "disposable"),
        ("eye", "melanocytes", "vital"),
        ("blood", "naive B-cell", "manageable"),
        ("blood", "b_cell_naive", "manageable"),
        ("blood", "cd8_stim", "protected"),
        ("blood", "total PBMC", "protected"),
        ("blood", "Memory CD4 T-cell Th1/Th17", "protected"),
        ("blood", "Non-Vd2 gdTCR", "protected"),
        ("blood", "Vd2 gdTCR", "protected"),
        ("blood", "Plasmablast", "protected"),
        ("eye", "limbal superficial epithelium", "protected"),
        ("eye", "limbal progenitor cells", "vital"),
    ],
)
def test_compartment_hazards(tissue, cell, tier):
    r = h.assess(tissue, cell, resolution="cell")
    assert r["tier"] == tier
    assert r["references"]


def test_bulk_never_receives_cell_downgrade():
    assert h.assess("lung", "ciliated cells")["tier"] == "vital"
    assert h.assess("skin", "melanocytes")["tier"] == "protected"
    assert h.assess("blood", "B cell")["tier"] == "vital"
    assert h.assess("lung", "ciliated and basal cells", resolution="cell")["tier"] == "vital"
    assert h.assess("lung", "unrecognized cells", resolution="cell")["reason"] == "organ_fallback"


def test_no_species_age_or_reproductive_safety_inference():
    assert h.assess("lung", taxon_id=10090)["tier"] == "unresolved"
    assert h.assess("lung", stage="fetal")["reason"] == "outside_adult_human_policy"
    assert h.assess("ovary")["tier"] == "protected"
    assert "reproductive" in h.assess("ovary")["anatomical_membership"]
    assert h.assess("unknown")["tier"] == "unresolved"


def test_inventory_separates_reproductive_scope_from_missing_hazard_calls():
    inventory = h.anatomy_inventory().set_index("concept_id")
    for concept in ["placenta", "trophoblast", "amnion", "yolk_sac", "embryo", "fetus"]:
        row = inventory.loc[concept]
        assert row.selection_category == "Reproductive"
        assert "reproductive" in row.anatomical_membership
        assert not row.requires_hazard_review
    assert inventory.loc["embryo", "stage"] == "embryonic"
    assert inventory.loc["fetus", "stage"] == "fetal"
    assert inventory.loc["ovary", "tier"] == "protected"
    assert inventory.loc["ovary", "selection_category"] == "Reproductive"
    for concept in ["soft_tissue", "head", "foot"]:
        row = inventory.loc[concept]
        assert row.requires_localization
        assert row.context_kind == "broad_location"
        assert pd.isna(row.selection_category)
        assert row.tier == "unresolved"
    assert inventory.loc["skin", "selection_category"] == "Protected"
    assert inventory.loc["hair_follicle", "selection_category"] == "Disposable"


def test_registry_is_detached_and_rejects_uncited_or_duplicate_calls():
    doc = h.registry()
    for change in ("references", "id", "tier"):
        bad = deepcopy(doc)
        if change == "references":
            bad["assessments"][0][change] = ["invented_source"]
        elif change == "id":
            bad["assessments"][0][change] = bad["assessments"][1]["id"]
        else:
            bad["assessments"][0][change] = "safe"
        with pytest.raises(ValueError):
            h.validate_registry(bad)
    doc["assessments"][0]["tier"] = "changed"
    assert h.registry()["assessments"][0]["tier"] == "vital"


def test_every_existing_somatic_root_has_a_cited_bulk_call():
    legacy = h.registry()["legacy_nonhuman_screening_roots"]
    for root in {r for roots in legacy.values() for r in roots}:
        result = h.assess(root)
        assert result["tier"] in h.TIERS, (root, result)
        assert result["references"]
        a = resolve_normal_tissue(root, species="human")
        assert normal_tissue_screening_tier(a) == result["tier"]


def test_cell_override_is_consistent_for_rna_and_ihc_but_not_bulk():
    df = pd.DataFrame(
        {
            "tissue": ["lung"] * 3,
            "cell_type": ["ciliated cells"] * 3,
            "category": ["vital"] * 3,
            "mapping_status": ["mapped"] * 3,
            "source": ["hpa_single_cell_tissue", "hpa_normal_tissue", "hpa_rna_consensus"],
        }
    )
    assert h.context_tiers(df).tolist() == ["protected", "protected", "vital"]
    df["mapping_status"] = "unmapped"
    assert h.context_tiers(df).tolist() == ["vital"] * 3


def test_conflicting_equal_priority_rules_fail_closed(monkeypatch):
    doc = h.registry()
    rule = deepcopy(next(r for r in doc["assessments"] if r["id"] == "airway_ciliated"))
    rule.update(id="conflict", tier="manageable")
    doc["assessments"].append(rule)
    monkeypatch.setattr(h, "_load", lambda: doc)
    h._resolve.cache_clear()
    try:
        assert (
            h.assess("lung", "ciliated cells", resolution="cell")["reason"] == "ambiguous_cell_rule"
        )
    finally:
        h._resolve.cache_clear()


def test_purified_immune_assay_can_resolve_cells_but_bulk_blood_cannot():
    frame = pd.DataFrame(
        {
            "tissue": ["blood"] * 3,
            "cell_type": ["naive B-cell", "b_cell_naive", "B cell"],
            "source": ["hpa_immune_cell", "dice_b_cell_naive", "hpa_rna_consensus"],
            "category": ["protected"] * 3,
            "mapping_status": ["mapped"] * 3,
        }
    )
    assert h.context_tiers(frame).tolist() == ["manageable", "manageable", "vital"]


def test_removable_organ_does_not_make_shared_cells_disposable():
    assert h.assess("tonsil")["tier"] == "disposable"
    assert h.assess("tonsil", "B cells", resolution="cell")["tier"] == "manageable"
    assert h.assess("tonsil", "B cells and T cells", resolution="cell")["tier"] == "protected"
    frame = pd.DataFrame(
        {
            "tissue": ["tonsil"],
            "cell_type": ["B cells"],
            "source": ["hpa_normal_tissue"],
            "category": ["disposable"],
            "mapping_status": ["mapped"],
        }
    )
    assert h.context_tiers(frame).iloc[0] == "manageable"
