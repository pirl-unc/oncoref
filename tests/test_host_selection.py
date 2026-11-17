"""Protection follows evidence and declared anatomy, without named-gene exemptions."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parents[1] / "scripts"))
from types import SimpleNamespace

import pandas as pd
import pytest

from oncoref import antigen_sets
from oncoref import host_context as hc
from oncoref import target_filters as tf


def context(protect=("ovary",), **kwargs):
    return hc.HostContext(absent=tuple(c for c in hc.COMPARTMENTS if c not in protect), **kwargs)


def observation(value=0, **kwargs):
    return {
        "coverage": "measured",
        "max_tpm": value,
        "scope_complete": True,
        "max_ntpm": None,
        "protein_status": "not_detected",
        **kwargs,
    }


@pytest.fixture
def rows(monkeypatch):
    # These are already measured rows; evidence joins are tested on the real catalog elsewhere.
    monkeypatch.setattr(tf, "measures", lambda table=None: table.copy())
    return pd.DataFrame(
        [
            {
                "entity_id": "a",
                "protein_group_label": "PRAME",
                "member_gene_ids": '["ga"]',
                "category": "CTA",
                "classification": "CTA",
                "cta_tier": "core",
                "is_cta": True,
                "is_cancer_antigen": True,
                "is_reproductive_antigen": True,
                "is_target_candidate": True,
                "normal_reference_tpm": 1.0,
                "disposable_max_tpm": 0.0,
                "off_target_tpm": 1.0,
                "tumor_fold": 50.0,
                "tumor_max_tpm": 50.0,
                "reproductive_fold": 30.0,
                "curated_tier": "",
                "papers_counted": 10,
                "independent_patient_tumor_folds": "[50, 20]",
            }
        ]
    )


@pytest.mark.parametrize(
    "item,status",
    [
        (None, "unresolved"),
        ({}, "unresolved"),
        (observation(0), "eligible"),
        (observation(None), "unresolved"),
        (observation(float("nan")), "unresolved"),
        (observation(float("inf")), "unresolved"),
        (observation(-1), "unresolved"),
        (observation(True), "unresolved"),
        (observation(4), "excluded"),
        (observation(0, scope_complete=False), "unresolved"),
        (observation(0, protein_status="unresolved"), "unresolved"),
        (observation(0, max_ntpm=5), "excluded"),
    ],
)
def test_missing_is_not_zero_and_native_assays_remain_separate(rows, item, status):
    result = hc.assess(rows, context(), evidence={"a": {"ovary": item}}).iloc[0]
    assert result.host_status == status
    assert result.pass_host == (status == "eligible")


def test_anatomy_presets_never_accept_gonadal_loss():
    for profile, gonad in (("male-patient", "testis"), ("female-patient", "ovary")):
        host = hc.from_legacy(profile)
        assert gonad in host.protected and not host.accept_loss
        assert "resting_breast" in host.protected
        assert host.assumptions
    pregnant = hc.from_legacy("female-patient", pregnant=True)
    assert {"embryo", "placenta", "lactating_breast"} <= set(pregnant.protected)
    with pytest.raises(ValueError, match="biological origin"):
        tf.TargetFilter(host_profile="core")
    with pytest.raises(ValueError, match="both absent"):
        hc.HostContext(absent=("ovary",), accept_loss=("ovary",))


@pytest.mark.parametrize(
    "kwargs",
    [{"max_tpm": float("nan")}, {"max_ntpm": -1}, {"max_tpm": True}, {"absent": ["ovary"]}],
)
def test_invalid_contexts_fail(kwargs):
    with pytest.raises(ValueError):
        hc.HostContext(**kwargs)


def test_host_protection_can_exclude_positive_control_without_relabeling(rows):
    d = tf.apply(
        tf.TargetFilter(host=context()), rows, host_evidence={"a": {"ovary": observation(12.65)}}
    )
    assert not d.keep.item() and d.host_status.item() == "excluded"
    assert d.category.item() == d.filtered_category.item() == "CTA"
    assert d.is_cta.item() and d.filtered_cancer_antigen.item()


def test_off_target_and_patient_thresholds_include_host_expression(rows):
    ev = {"a": {"ovary": observation(2)}}
    d = tf.apply(
        tf.TargetFilter(host=context(max_tpm=20), max_off_target_tpm=2), rows, host_evidence=ev
    )
    assert not d.pass_off_target.item()
    # Two patients at 20x the baseline do not mean two at 20x the enlarged host reference.
    f = tf.TargetFilter(host=context(), min_tumor_fold=20, min_tumor_patients=2)
    d = tf.apply(f, rows, host_evidence=ev)
    assert d.filtered_cancer_antigen.item() and not d.pass_host.item()
    assert d.effective_tumor_fold.item() == 25
    assert not tf.apply(tf.TargetFilter(host=context()), rows, host_evidence={}).keep.item()


def test_protecting_more_tissues_or_lowering_limits_cannot_add_candidates(rows):
    ev = {"a": {"ovary": observation(3), "testis": observation(10)}}
    permissive = tf.apply(tf.TargetFilter(host=context()), rows, host_evidence=ev)
    assert permissive.keep.item()
    assert not tf.apply(
        tf.TargetFilter(host=context(max_tpm=2)), rows, host_evidence=ev
    ).keep.item()
    assert not tf.apply(
        tf.TargetFilter(host=context(("ovary", "testis"))), rows, host_evidence=ev
    ).keep.item()
    # Absence and accepted loss have distinct provenance even with the same local exposure mask.
    absent = context(())
    accepted = hc.HostContext(
        absent=tuple(c for c in hc.COMPARTMENTS if c != "ovary"), accept_loss=("ovary",)
    )
    assert absent.protected == accepted.protected == ()
    assert absent.states() != accepted.states()


def test_selection_views_and_host_background_share_one_decision(rows):
    second = rows.iloc[0].to_dict() | {
        "entity_id": "b",
        "member_gene_ids": '["gb"]',
        "protein_group_label": "second",
    }
    table = pd.concat([rows, pd.DataFrame([second])], ignore_index=True)
    selected = antigen_sets.select(
        host=context(),
        table=table,
        host_evidence={"a": {"ovary": observation(0)}, "b": {"ovary": observation(5)}},
    )
    proteins = selected.proteins()
    assert proteins.entity_id.tolist() == ["a"]
    assert proteins.is_cta.all() and proteins.classification.eq(proteins.category).all()
    assert selected.background_proteins().entity_id.tolist() == ["b"]
    assert selected.gene_ids() == frozenset({"ga"})

    def gene(g):
        return SimpleNamespace(
            transcripts=[
                SimpleNamespace(
                    id=g + "t",
                    biotype="protein_coding",
                    protein_sequence="LLLLLLLL" if g == "gb" else "IIIIIIII",
                )
            ]
        )

    genome = SimpleNamespace(
        release=115, species=SimpleNamespace(latin_name="homo_sapiens"), gene_by_id=gene
    )
    sequences = selected.sequences(genome=genome)
    assert set(sequences.gene_id) == selected.gene_ids()
    assert sequences.attrs["host_context"] == selected.metadata["host_context"]
    from oncoref.peptide_subtraction import MatchPolicy

    index = selected.background_index(
        genome=genome, policy=MatchPolicy(equivalence_groups=("ILV",)), lengths=(8,)
    )
    assert len(index.matches("IIIIIIII")) == 1
    # Returned frames/metadata must not mutate the captured decision.
    detached = selected.decisions()
    detached.loc[:, "keep"] = False
    meta = selected.metadata
    meta["host_context"]["max_tpm"] = 999
    assert len(selected.proteins()) == 1 and selected.metadata["host_context"]["max_tpm"] == 4


def test_empty_selection_keeps_schema(rows):
    selection = antigen_sets.select(table=rows.iloc[:0], host=context(), host_evidence={})
    assert selection.proteins().empty and selection.background_proteins().empty
    assert "host_status" in selection.decisions()


def test_exposure_builder_ignores_positive_ra_support_for_negative_reference():
    from build_host_exposure import summarize

    data = pd.DataFrame(
        {
            "identity_resolved": [True, True, False],
            "value": [0, 3, 99],
            "qc_pass": [False, False, False],
            "library_qc_review": ["", "", ""],
            "source_id": ["x", "x", "x"],
        }
    )
    adult = summarize(data, compartment="ovary")
    assert adult["max_tpm"] == 3 and adult["measured_observations"] == 2
    assert summarize(data, compartment="embryo")["scope_complete"] is False
    assert summarize(data.iloc[:0], compartment="ovary")["max_tpm"] is None


def test_tightened_biological_membership_updates_all_views_consistently(rows):
    selection = antigen_sets.select("RA", min_tumor_fold=60, table=rows)
    protein = selection.proteins().iloc[0]
    assert protein.category == protein.classification == "Reproductive Antigen only"
    assert not protein.is_cta and not protein.is_cancer_antigen
    assert protein.is_reproductive_antigen and protein.is_target_candidate
    assert protein.catalog_category == "CTA" and protein.catalog_is_cta
    assert protein.cta_tier == "" and protein.catalog_cta_tier == "core"
    assert selection.decisions().category.item() == "CTA"


def test_host_background_keeps_global_exclusions_but_not_paper_only_failures(rows):
    second = rows.iloc[0].to_dict() | {
        "entity_id": "other",
        "category": "Other",
        "classification": "Other:Review",
        "is_cta": False,
        "is_cancer_antigen": False,
        "is_reproductive_antigen": False,
        "is_target_candidate": False,
    }
    table = pd.concat([rows, pd.DataFrame([second])], ignore_index=True)
    selection = antigen_sets.select(table=table, min_papers=100)
    assert selection.proteins().empty
    assert selection.background_proteins().entity_id.tolist() == ["other"]
    assert selection.background_proteins("Unsafe").empty


def test_flagged_normal_library_cannot_lower_host_exposure_and_clear_coverage():
    from build_host_exposure import summarize

    data = pd.DataFrame(
        {
            "identity_resolved": [True, True],
            "value": [0, 20],
            "library_qc_review": ["", "unrelated tissue signatures"],
            "source_id": ["x", "x"],
        }
    )
    result = summarize(data, compartment="ovary")
    assert result["max_tpm"] == 20
    assert result["qc_review_observations"] == 1
    assert not result["scope_complete"]
