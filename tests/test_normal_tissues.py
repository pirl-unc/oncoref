# Licensed under the Apache License, Version 2.0.
import json
from dataclasses import replace
from pathlib import Path

import pandas as pd
import pytest

from oncoref import canine, cta_regen, human, mouse, reference_data
from oncoref.normal_tissues import (
    classify_normal_tissue_pattern,
    normal_tissue_category,
    normal_tissue_inventory,
    normal_tissue_mapping_inventory,
    normal_tissue_mappings,
    normal_tissue_scope,
    normal_tissue_scopes,
    normal_tissue_screening_details,
    normal_tissue_screening_policy,
    normal_tissue_screening_tier,
    normal_tissue_taxonomy,
    resolve_normal_tissue,
)


@pytest.mark.parametrize("species", ["human", "mouse", "dog"])
def test_restriction_is_subset_of_tissue_union_not_union_of_restricted_gene_lists(species):
    def classify(*labels):
        return classify_normal_tissue_pattern(
            [resolve_normal_tissue(label, species=species) for label in labels]
        )

    mixed = classify("testis", "oocyte")
    assert "reproductive" in mixed["matching_scopes"]
    assert not {"testis", "oocyte", "male_reproductive", "female_reproductive"} & set(
        mixed["matching_scopes"]
    )
    assert classify("testis", "mammary gland")["narrowest_scopes"] == ("reproductive",)
    assert classify("testis", "prostate")["narrowest_scopes"] == ("male_reproductive",)
    assert classify("testis", "cutaneous melanocytes")["narrowest_scopes"] == (
        "reproductive_or_disposable",
    )
    assert classify("hair follicle", "cutaneous melanocytes")["narrowest_scopes"] == (
        "disposable_normal",
    )
    assert classify("testis", "heart")["matching_scopes"] == ()
    assert classify()["status"] == "no_detections"
    assert classify("testis", "unrecognized tissue")["status"] == "unresolved"
    assert mixed["absence_evidence"] == "not_assessed"


@pytest.mark.parametrize("species", ["human", "mouse", "dog"])
def test_scopes_are_exact_unions_and_citations_are_present(species):
    taxonomy = normal_tissue_taxonomy()
    scopes = normal_tissue_scopes(species=species)
    assert scopes["reproductive"] == set.union(
        *[
            set(scopes[s])
            for s in (
                "male_reproductive",
                "female_reproductive",
                "extraembryonic",
                "embryonic",
                "fetal",
            )
        ]
    )
    assert (
        scopes["reproductive_or_disposable"] == scopes["reproductive"] | scopes["disposable_normal"]
    )
    assert not scopes["disposable_normal"] & scopes["reproductive"]
    assert normal_tissue_scope("male_strict", species=species) == scopes["testis"]
    assert normal_tissue_scope("female_strict", species=species) == scopes["oocyte"]
    for scope in scopes.values():
        for concept in scope:
            assert taxonomy["concepts"][concept]["references"]
    assert "thymus" not in scopes["reproductive_or_disposable"]
    assert scopes["reproductive_or_disposable"] < scopes["reproductive_or_permitted_normal"]
    assert (
        scopes["reproductive_or_permitted_normal"]
        == scopes["reproductive"] | scopes["permitted_normal"]
    )
    taxonomy["concepts"]["testis"]["parents"].append("bad")
    assert "bad" not in normal_tissue_taxonomy()["concepts"]["testis"]["parents"]


@pytest.mark.parametrize("species", ["human", "dog", "mouse"])
def test_expanded_permission_and_screening_tiers_are_species_aware(species):
    def tier(tissue, **kwargs):
        return normal_tissue_screening_tier(
            resolve_normal_tissue(tissue, species=species, **kwargs)
        )

    for tissue in (
        "testis",
        "placenta",
        "cutaneous melanocytes",
        "hair follicle",
        "gallbladder",
    ):
        assert tier(tissue) == "permitted"
    for tissue in (
        "heart",
        "left ventricle",
        "lung",
        "bronchus",
        "amygdala",
        "caudate",
        "cerebellum",
        "hypothalamus",
        "spinal cord",
        "liver",
        "renal cortex",
        "pancreas",
        "adrenal gland",
        "pituitary gland",
        "bone marrow",
        "duodenum",
        "skeletal muscle",
        "retina",
        "retinal pigment epithelium",
        "eye",
        "inner ear",
    ):
        assert tier(tissue) == "vital"
    for tissue in (
        "lymph node",
        "stomach",
        "cartilage",
        "esophagus",
        "salivary gland",
        "parotid gland",
    ):
        assert tier(tissue) == "protected"
    for tissue in ("skin", "nail bed", "thymus", "spleen", "thyroid gland", "adipose tissue"):
        assert tier(tissue) == "manageable"
    assert tier("tonsil") == ("unresolved" if species == "mouse" else "permitted")
    assert tier("appendix") == ("permitted" if species == "human" else "unresolved")
    for tissue in ("head", "foot", "soft tissue", "melanocytes", "unknown"):
        assert tier(tissue) == "unresolved"
    assert tier("liver", stage="fetal") == "permitted"
    assert tier("liver", context="tumor") == "unresolved"
    mixed = classify_normal_tissue_pattern(
        [resolve_normal_tissue(tissue, species=species) for tissue in ("testis", "thymus", "skin")]
    )
    assert mixed["matching_scopes"] == ()


@pytest.mark.parametrize("species", ["human", "dog", "mouse"])
def test_every_anatomical_concept_has_one_disjoint_screening_assignment(species):
    taxonomy = normal_tissue_taxonomy()
    policy = normal_tissue_screening_policy()
    taxon = {"human": 9606, "dog": 9615, "mouse": 10090}[species]

    def ancestors(cid):
        return {cid}.union(*(ancestors(p) for p in taxonomy["concepts"][cid]["parents"]))

    permitted = normal_tissue_scope(policy["permitted_scope"], species=species)
    for cid, concept in taxonomy["concepts"].items():
        if taxon not in concept.get("taxon_ids", taxonomy["taxon_ids"]):
            continue
        # Explicit narrow permitted compartments override the broader skin tier.
        matches = (
            ["permitted"]
            if cid in permitted
            else [
                name
                for name in ("vital", "protected", "manageable")
                if ancestors(cid) & set(policy[name]["roots"])
            ]
        )
        matches += ["unresolved"] if cid in policy["unresolved_roots"] else []
        assert len(matches) == 1, (cid, matches)
        assert (
            normal_tissue_screening_tier(resolve_normal_tissue(cid, species=species)) == matches[0]
        )


def test_acceptable_aliases_and_lymphoid_review_do_not_allow_all_immune_cells():
    assert normal_tissue_scope("acceptable_normal", species="human") == normal_tissue_scope(
        "permitted_normal", species="human"
    )
    spleen = normal_tissue_screening_details(human.resolve_normal_tissue("spleen"))
    assert spleen["screening_tier"] == "manageable"
    assert "asplenia_infection_risk_review" in spleen["review_flags"]
    assert "immune_cell_distribution_review" in spleen["review_flags"]
    assert any("cdc.gov" in url for url in spleen["references"])
    spleen["context_reviews"].clear()
    assert normal_tissue_screening_details(human.resolve_normal_tissue("spleen"))["context_reviews"]
    fetal = normal_tissue_screening_details(human.resolve_normal_tissue("spleen", stage="fetal"))
    assert fetal["screening_tier"] == "permitted" and not fetal["review_flags"]
    animal = normal_tissue_screening_details(canine.resolve_normal_tissue("spleen"))
    assert "not_validated" in animal["tolerability_evidence"]
    for tissue in ("blood", "bone marrow"):
        assert normal_tissue_screening_tier(human.resolve_normal_tissue(tissue)) == "vital"
    for cell in ("B-cells", "T-cells"):
        annotation = human.resolve_normal_tissue(
            cell_type=cell, source_id="hpa_single_cell", source_version="v23"
        )
        assert normal_tissue_screening_tier(annotation) == "unresolved"
    assert not normal_tissue_screening_details(mouse.resolve_normal_tissue("tonsil"))[
        "review_flags"
    ]


def test_policy_is_detached_and_stale_annotations_are_rejected():
    policy = normal_tissue_screening_policy()
    assert policy["rna_unit"] == "nTPM" and policy["rna_species"] == "human"
    policy["vital"]["maximum_rna"] = 999
    assert normal_tissue_screening_policy()["vital"]["maximum_rna"] == 1
    with pytest.raises(ValueError, match="Unknown normal-tissue screening policy"):
        normal_tissue_screening_policy("unknown")
    with pytest.raises(ValueError, match="current taxonomy"):
        normal_tissue_screening_tier(
            replace(human.resolve_normal_tissue("heart"), taxonomy_id="old")
        )


@pytest.mark.parametrize("species", ["human", "dog", "mouse"])
def test_five_categories_and_native_mapping_inventory_remain_distinct_from_coverage(species):
    inventory = normal_tissue_inventory(species=species)
    by_concept = {r["concept_id"]: r for r in inventory}
    assert by_concept["testis"]["tissue_category"] == "reproductive"
    assert by_concept["gallbladder"]["tissue_category"] == "disposable"
    assert by_concept["lymph_node"]["tissue_category"] == "protected"
    assert by_concept["skin"]["tissue_category"] == "manageable"
    assert by_concept["skin"]["tissue_category_label"] == "Manageable risk (conditional)"
    assert by_concept["lymph_node"]["tissue_category_label"] == "High consequence"
    assert by_concept["cutaneous_melanocyte"]["tissue_category"] == "disposable"
    assert by_concept["skin"]["category_policy_references"]
    assert by_concept["retina"]["tissue_category"] == "vital"
    for cid in ("head", "foot", "soft_tissue"):
        assert by_concept[cid]["tissue_category"] == "unresolved"
    for cid in (
        "uterine_horn",
        "coagulating_gland",
        "tonsil",
        "appendix",
        "seminal_vesicle",
        "bulbourethral_gland",
    ):
        row = by_concept[cid]
        assert (row["mapping_status"] == "mapped") == (
            row["taxon_id"] in row["applicable_taxon_ids"]
        )
    mappings = normal_tissue_mapping_inventory(species=species)
    for row in mappings:
        assert row["expression_coverage"] == "not_assessed"
        assert row["source_version"] and row["source_reference"]
    if species == "mouse":
        assert mappings == ()
    assert (
        normal_tissue_category(resolve_normal_tissue("liver", species=species, stage="fetal"))
        == "reproductive"
    )
    assert (
        normal_tissue_category(resolve_normal_tissue("liver", species=species, stage="adult"))
        == "vital"
    )


def test_unknown_stage_and_pregnant_maternal_context_are_not_reclassified():
    unknown = human.resolve_normal_tissue("testis", stage="unreviewed age label")
    assert unknown.mapping_status == "unresolved_stage" and not unknown.scopes
    maternal = human.resolve_normal_tissue(
        "liver", stage="adult", sex="female", physiological_state="pregnant"
    )
    assert maternal.sex == "female" and maternal.physiological_state == "pregnant"
    assert normal_tissue_category(maternal) == "vital"


def test_placenta_does_not_require_maternal_fetal_compartment_resolution():
    for species in (human, canine, mouse):
        placenta = species.resolve_normal_tissue("placenta")
        assert placenta.mapping_status == "mapped"
        assert "placental" in placenta.scopes
        assert "reproductive" in placenta.scopes


def test_exact_hpa_compartments_are_distinct_from_bulk_labels():
    def ihc(tissue, cell=None):
        return human.resolve_normal_tissue(
            tissue, cell_type=cell, source_id="hpa_normal_tissue", source_version="v23"
        )

    assert "mammary_glandular" in ihc("breast", "glandular cells").scopes
    assert "mammary_glandular" in ihc("breast", "myoepithelial cells").scopes
    assert "mammary_glandular" not in ihc("breast").scopes
    assert "mammary_glandular" not in ihc("breast", "adipocytes").scopes
    assert "female_reproductive" in ihc("breast").scopes
    assert ihc("lactating breast", "lactating glandular cells").physiological_state == "lactating"
    assert "oocyte" not in ihc("ovary", "follicle cells").scopes
    assert "oocyte" not in ihc("ovary").scopes
    assert "disposable_normal" in ihc("skin 1", "melanocytes").scopes
    assert "disposable_normal" in ihc("skin 2", "hair follicles").scopes
    assert "disposable_normal" in ihc("hair", "cells in external root sheath").scopes
    for tissue, cell in (
        ("skin 1", "keratinocytes"),
        ("skin 2", None),
        ("retina", "pigment epithelial cells"),
    ):
        assert "disposable_normal" not in ihc(tissue, cell).scopes
    assert "placental" in ihc("placenta", "decidual cells").scopes
    assert ihc("endometrium 1").concept_id == ihc("endometrium 2").concept_id == "endometrium"


def test_single_cell_source_resolves_only_supported_anatomical_sites():
    def sc(cell):
        return human.resolve_normal_tissue(
            cell_type=cell, source_id="hpa_single_cell", source_version="v23"
        )

    assert sc("Oocytes").concept_id == "oocyte"
    assert "disposable_normal" in sc("Melanocytes").scopes
    assert sc("Glandular and luminal cells").mapping_status == "unmapped"
    assert sc("Fibroblasts").scopes == ()
    assert human.resolve_normal_tissue("follicle cells").scopes == ()
    assert "testis" not in human.resolve_normal_tissue("spermatozoa").scopes
    assert "male_reproductive" in human.resolve_normal_tissue("spermatozoa").scopes
    # A generic, unsited melanocyte label cannot inherit the HPA skin mapping.
    assert human.resolve_normal_tissue("melanocytes").scopes == ()
    assert (
        human.resolve_normal_tissue("eye", cell_type="cutaneous melanocytes").mapping_status
        == "conflicting_anatomy"
    )


def test_stage_specific_expression_and_species_applicability():
    fetal = mouse.resolve_normal_tissue("liver", stage="fetal")
    adult = mouse.resolve_normal_tissue("liver", stage="adult")
    assert fetal.concept_id == adult.concept_id == "liver"
    assert "fetal" in fetal.scopes and "reproductive" in fetal.scopes
    assert adult.scopes == ()
    assert classify_normal_tissue_pattern([fetal, adult])["matching_scopes"] == ()
    assert (
        mouse.resolve_normal_tissue("whole embryo", stage="adult").mapping_status
        == "conflicting_stage"
    )
    assert "male_reproductive" in mouse.resolve_normal_tissue("coagulating gland").scopes
    for label in ("seminal vesicle", "bulbourethral gland", "coagulating gland"):
        assert canine.resolve_normal_tissue(label).mapping_status == "not_applicable_to_species"
    assert (
        canine.resolve_normal_tissue(
            "testis", source_id="hpa_rna_consensus", source_version="v23"
        ).mapping_status
        == "source_unavailable"
    )
    assert (
        human.resolve_normal_tissue(
            "testis", source_id="hpa_rna_consensus", source_version="latest"
        ).mapping_status
        == "source_unavailable"
    )
    assert human.resolve_normal_tissue("testis", context="tumor").scopes == ()
    with pytest.raises(ValueError, match="different species"):
        classify_normal_tissue_pattern(
            [human.resolve_normal_tissue("testis"), canine.resolve_normal_tissue("testis")]
        )
    with pytest.raises(ValueError, match="current tissue taxonomy"):
        classify_normal_tissue_pattern([replace(fetal, taxonomy_id="older-taxonomy")])
    assert "not_validated" in canine.resolve_normal_tissue("hair follicle").tolerability_evidence


def test_every_bundled_canine_normal_sample_tissue_maps_without_human_fallback():
    from oncoref.expression_registry import expression_source_registry_entries

    root = Path(__file__).parents[1] / "oncoref/data/species/dog"
    for source in expression_source_registry_entries(species="dog"):
        if "sample_metadata" not in source:
            continue
        samples = json.loads((root / source["sample_metadata"]).read_text())["samples"]
        for sample in samples:
            annotation = canine.resolve_normal_tissue(
                sample["tissue"],
                stage=sample.get("stage"),
                source_id=source["id"],
                source_version=source["source_version"],
            )
            if sample["tissue"] is not None:
                assert annotation.mapping_status == "mapped", sample
            if sample.get("stage") == "embryonic":
                assert "embryonic" in annotation.scopes
                assert "reproductive" in annotation.scopes
    a = canine.resolve_normal_tissue(
        "Breast/Mammary",
        source_id="pfizer-dog-protein",
        source_version="MSV000099744-update-2026-03-23",
    )
    assert "female_reproductive" in a.scopes
    assert "mammary_glandular" not in a.scopes


def test_hpa_mappings_match_pinned_source_vocabulary():
    for source in ("hpa_rna_consensus", "hpa_normal_tissue", "hpa_single_cell"):
        if not reference_data.is_cached(source, "v23"):
            pytest.skip("HPA v23 not cached")
        df = pd.read_csv(
            reference_data.local_path(source, "v23"),
            sep="\t",
            usecols=lambda c: c in {"Tissue", "Cell type"},
        ).fillna("")
        pairs = set(
            zip(
                df.get("Tissue", pd.Series("", index=df.index)),
                df.get("Cell type", pd.Series("", index=df.index)),
            )
        )
        mappings = normal_tissue_mappings(species="human", source_id=source, source_version="v23")
        mapped_pairs = {(r["source_tissue"], r["source_cell_type"]) for r in mappings}
        assert pairs - {("", "")} <= mapped_pairs
        for row in mappings:
            if row["source_cell_type"] or source == "hpa_rna_consensus":
                assert (row["source_tissue"], row["source_cell_type"]) in pairs
            assert row["source_reference"].startswith("https://v23.proteinatlas.org/")


def test_explicit_snapshot_annotation_preserves_frozen_archive():
    snapshot = canine.cta_snapshot()
    original = snapshot.as_dict()
    annotations = snapshot.tissue_annotations(candidate_id="CTCFL")
    assert annotations
    assert any("placental" in r["annotation"]["scopes"] for r in annotations)
    assert snapshot.as_dict() == original
    assert mouse.cta_tissue_annotations() == ()
    with pytest.raises(ValueError, match="Unknown CTA candidate"):
        snapshot.tissue_annotations(candidate_id="missing")


def test_human_snapshot_annotations_require_declared_hpa_version():
    from oncoref.species_antigens import CTAEvidenceSnapshot

    snapshot = human.cta_snapshot()
    candidate = next(c for c in snapshot.candidates() if c["label"] == "CGB1")["id"]
    annotations = snapshot.tissue_annotations(candidate_id=candidate)
    assert any("placental" in r["annotation"]["scopes"] for r in annotations)
    payload = snapshot.as_dict()
    payload.pop("snapshot_id")
    payload["sources"][0]["lineage"].remove("HPA:v23")
    undeclared = CTAEvidenceSnapshot.create(payload)
    assert all(
        not r["annotation"]["scopes"] for r in undeclared.tissue_annotations(candidate_id=candidate)
    )


@pytest.mark.parametrize("scope", ["core", "extended"])
def test_regenerated_restriction_accounts_for_accessory_and_thymus_detections(scope):
    seed = pd.DataFrame({"Ensembl_Gene_ID": ["T", "P", "B", "S"], "Symbol": ["T", "P", "B", "S"]})
    values = {
        "T": {"testis": 10, "prostate": 10},
        "P": {"placenta": 10, "breast": 10},
        "B": {"breast": 10},
        "S": {"testis": 10, "thymus": 10},
    }
    cta_regen._recompute_rna_columns(seed, values, tissue_scope=scope)
    assert list(seed.rna_restriction) == ["REPRODUCTIVE", "REPRODUCTIVE", "REPRODUCTIVE", "SOMATIC"]
    for tissues, restriction in (
        ("testis;prostate", "REPRODUCTIVE"),
        ("placenta;breast", "REPRODUCTIVE"),
        ("testis;thymus", "SOMATIC"),
        ("testis", "TESTIS"),
    ):
        assert (
            cta_regen._assign_protein_restriction(pd.Series({"protein_strict_expression": tissues}))
            == restriction
        )
