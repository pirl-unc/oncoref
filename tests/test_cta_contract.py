"""Consumer roles must not convert unresolved evidence into normal clearance."""

import pandas as pd
import pytest

from oncoref import cta_contract, cta_selection


def test_assessment_identity_and_exclusion_precedence_across_other_subtypes():
    rows = pd.DataFrame(
        [
            ("candidate", "CTA", "", "screened", "resolved", "allow", True),
            ("assessed", "Other", "Nonreproductive", "screened", "resolved", "allow", False),
            ("held", "Other", "Held", "normal_evidence_hold", "resolved", "review", False),
            (
                "missing",
                "Other",
                "Nonreproductive",
                "incomplete_evidence",
                "resolved",
                "review",
                False,
            ),
            ("unknown_id", "CTA", "", "screened", "ambiguous", "review", False),
            ("unknown_status", "CTA", "", "future_status", "resolved", "review", False),
            ("no_status", "Other", "Diffuse", None, "resolved", "review", False),
            ("no_identity", "CTA", "", "screened", None, "review", False),
            ("native", None, None, "unassessed", "resolved", "review", False),
            ("excluded", "Unsafe", "", "incomplete_evidence", "resolved", "exclude", False),
        ],
        columns=[
            "entity_id",
            "category",
            "subtype",
            "assessment_status",
            "identity_status",
            "expected_role",
            "expected_candidate",
        ],
    ).convert_dtypes()
    result = cta_contract.sequence_roles(rows)
    assert result.shared_sequence.tolist() == rows.expected_role.tolist()
    assert result.direct_target_candidate.tolist() == rows.expected_candidate.tolist()
    assert "assessment_requires_review" in result.iloc[-1].sequence_role_reasons
    assert "applicable_policy_exclusion" in result.iloc[-1].sequence_role_reasons
    reversed_rows = cta_contract.sequence_roles(rows.iloc[::-1]).sort_values("entity_id")
    pd.testing.assert_frame_equal(reversed_rows, result.sort_values("entity_id"))


def test_every_current_incomplete_other_record_requires_review():
    result = cta_selection.sequence_dispositions()
    other = result[result.category.eq("Other")]
    unresolved = other[~other.assessment_status.eq("screened")]
    assert not unresolved.empty
    assert not unresolved.subtype.eq("Held").all(), "Subtype alone misses unresolved cases"
    assert unresolved.shared_sequence.eq("review").all()
    assert not unresolved.direct_target_candidate.any()
    assert result.loc[result.category.eq("Unsafe"), "shared_sequence"].eq("exclude").all()
    assert result.role_policy_id.eq(cta_contract.ROLE_POLICY_ID).all()
    assert result.role_scope.eq("selected_representative_protein").all()


def test_roles_do_not_fill_missing_fields_or_accept_duplicate_identities():
    frame = cta_selection.results().iloc[:2]
    with pytest.raises(ValueError, match="Missing decision fields"):
        cta_contract.sequence_roles(frame.drop(columns="assessment_status"))
    with pytest.raises(ValueError, match="Unique"):
        cta_contract.sequence_roles(pd.concat([frame, frame]))
    with pytest.raises(ValueError, match="Unknown decision"):
        cta_contract.sequence_roles(frame.assign(category="Safe"))
    assert cta_contract.sequence_roles(frame.iloc[:0]).empty


def test_contract_is_detached_and_nonhuman_selection_remains_unavailable():
    first = cta_selection.contract()
    first["categories"]["Other"] = "mutated"
    assert cta_selection.contract()["categories"]["Other"] != "mutated"
    assert cta_selection.contract()["selection_policy_id"] == cta_selection.coverage()["policy_id"]
    mouse = cta_selection.contract(species="mouse")
    assert mouse["taxon_id"] == 10090
    assert mouse["selection_policy_id"] is None
    assert mouse["selection_status"] == "unavailable"
    assert cta_selection.sequence_dispositions(species="mouse").shared_sequence.eq("review").all()
    gates = cta_selection.gates()
    assert gates.gate_id.is_unique
    assert {"required_coverage", "normal_knowledge", "translation_qualification"} <= set(
        gates.gate_id
    )


def test_entire_inventory_is_classified_with_translation_support_separate():
    summary = cta_selection.coverage()
    counts = summary["universe_counts"]
    assert counts["classified_groups"] == len(cta_selection.results())
    assert counts["input_groups"] == counts["classified_groups"] + counts["out_of_universe_groups"]
    assert counts["classified_gene_loci"] <= counts["input_gene_loci"]
    assert counts["out_of_universe_groups"] == 0
    assert counts["out_of_universe_gene_loci"] == 0
    assert counts["input_groups"] == counts["classified_groups"]
    assert (
        counts["translation_support_insufficient_groups"]
        == (~cta_selection.profiles().translation_support_eligible).sum()
    )
    assert summary["universe"].startswith(str(counts["classified_groups"]))


def test_insufficient_translation_support_does_not_hide_normal_exclusion():
    import json
    from pathlib import Path

    import pandas as pd

    root = Path(__file__).resolve().parents[1]
    prior = pd.read_csv(root / "docs/audits/cta-inventory-20261105/formerly-out-of-universe.csv.gz")
    profiles = cta_selection.profiles().set_index("entity_id")
    frozen = pd.read_csv(
        root / "docs/audits/cta-release-1.21.0/protein-group-classifications.csv.gz"
    )
    members = frozen.set_index("protein_equivalence_group_id").member_gene_ids.map(json.loads)
    genes = {g for values in members.loc[prior.protein_equivalence_group_id] for g in values}
    restored = profiles[profiles.gene_ids.map(lambda ids: bool(set(ids) & genes))]
    assert genes <= {g for values in restored.gene_ids for g in values}
    assert not restored.translation_qualified.any()
    assert not restored.translation_support_eligible.any()
    assert set(restored.classification) <= {"Unsafe", "Other:Translation"}
    assert restored.gene_ids.map(len).gt(0).all()
    roles = cta_selection.sequence_dispositions().set_index("entity_id").loc[restored.index]
    assert roles.loc[restored.category.eq("Unsafe"), "shared_sequence"].eq("exclude").all()
    assert roles.loc[restored.category.eq("Other"), "shared_sequence"].eq("review").all()
    assert not roles.direct_target_candidate.any()


def test_expression_profile_never_labels_known_normal_production_as_restricted():
    profiles = cta_selection.profiles().set_index("label")
    examples = profiles.loc[["IFNW1", "IFNA4", "MRGPRX4", "CNGA2", "SCN10A", "GFRAL"]]
    assert examples.preliminary_tumor_restriction_profile.all()
    assert examples.category.eq("Unsafe").all()
    assert examples.tumor_restriction_status.eq("excluded_by_normal_evidence").all()
    assert not examples.evidence_flags.map(lambda flags: "tumor_restricted_lead" in flags).any()
