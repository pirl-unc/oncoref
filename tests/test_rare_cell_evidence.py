"""Donor-level rare-cell evidence, its mapping-risk flags and their provenance."""

import hashlib
import json
import sys
from pathlib import Path

import pandas as pd
import pytest

from oncoref import target_selection as policy

ROOT = Path(__file__).parents[1]
DATA = ROOT / "oncoref/data"
AUDIT = ROOT / "docs/audits/normal-cell-atlas-20261110/rare-cell"
sys.path.insert(0, str(ROOT / "scripts"))
from build_gene_mapping_risk import K, kmers  # noqa: E402
from build_rare_cell_evidence import panel  # noqa: E402


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


@pytest.fixture(scope="module")
def evidence():
    return pd.read_csv(DATA / "rare-cell-evidence.csv.gz").set_index("protein_group_label")


def test_evidence_and_flags_match_their_provenance(evidence):
    provenance = json.loads((DATA / "rare-cell-evidence.json").read_text())
    assert provenance["sha256"] == sha(DATA / "rare-cell-evidence.csv.gz")
    assert provenance["policy_id"] == policy.POLICY_ID
    assert provenance["rules"] == policy.definition()["rare_cell_rna"]
    assert provenance["mapping_risk"]["sha256"] == sha(DATA / "gene-mapping-risk.csv.gz")
    assert provenance["fractional_count_cells"] == 0
    assert provenance["status_counts"] == evidence.rare_cell_status.value_counts().to_dict()
    assert sha(AUDIT / "gene-mapping-risk.csv.gz") == sha(DATA / "gene-mapping-risk.csv.gz")


def test_actions_follow_the_grade_and_never_rest_on_single_molecules(evidence):
    unsafe = evidence[evidence.rare_cell_action.eq("unsafe")]
    assert (
        unsafe.rare_cell_status.eq("established").all() and unsafe.rare_cell_tier.eq("vital").all()
    )
    assert not (unsafe.paralog_ambiguous | unsafe.host_overlapped).any()
    assert unsafe.finding_replicated.all() and unsafe.finding_orthogonal_confirmation.all()
    assert unsafe.finding_bulk_tpm.ge(policy.DETECTION_FLOOR).all()
    assert unsafe.finding_studies_multi_molecule.ge(policy.RARE_CELL_REPLICATE_STUDIES).all()
    graded = evidence[evidence.rare_cell_status.isin(["established", "unresolved"])]
    assert graded.finding_top_multi_molecule_pct.ge(policy.RARE_CELL_DETECTION_PCT).all()
    assert graded.finding_top_n_multi.ge(policy.RARE_CELL_MIN_CELLS).all()
    single_study = graded[~graded.finding_replicated.astype(bool)]
    assert single_study.finding_top_multi_molecule_pct.ge(policy.RARE_CELL_SINGLE_STUDY_PCT).all()
    assert single_study.rare_cell_status.eq("unresolved").all()
    review = evidence[evidence.rare_cell_action.eq("review")]
    assert review.rare_cell_unresolved_because.str.len().gt(0).all()
    # Frequent detection without multi-molecule cells is reported, not graded.
    quiet = evidence[evidence.rare_cell_status.eq("below_screen")]
    assert quiet.max_detection_pct.max() > 5 and quiet.rare_cell_action.eq("none").all()


def test_known_mapping_confounders_are_flagged(evidence):
    assert evidence.loc["TAS2R1", "host_overlapped"] and evidence.loc["TMEM269", "host_overlapped"]
    assert evidence.loc["KRTAP5-11", "paralog_ambiguous"]
    assert "KRTAP5" in evidence.loc["KRTAP5-11", "sequence_partners"]
    risk = pd.read_csv(DATA / "gene-mapping-risk.csv.gz").set_index("symbol")
    assert (
        risk.loc["MAGEA12", "paralog_ambiguous"]
        and "MAGEA" in risk.loc["MAGEA12", "sequence_partners"]
    )
    # Alternate-scaffold copies of a gene are not paralogs of it.
    for name in ["PRAME", "NUTM1", "PLAC1"]:
        assert not risk.loc[name, "paralog_ambiguous"], name
        partners = str(risk.loc[name, "sequence_partners"]).split("; ")
        assert not any(p.startswith(name + " (") for p in partners), name
    assert kmers("A" * (K + 1)) == {"A" * K}
    assert kmers("A" * (K - 1)) == set()


def test_panel_keeps_adult_donor_resolved_named_populations():
    base = {
        "age_scope": "adult",
        "donor_identity_known": True,
        "count_unit": "UMI_protocol",
        "tissue_general": "brain",
        "cell_type": "cerebellar granule cell",
        "collection_id": "study",
        "donor_id": "d1",
        "n_cells": 5000,
    }
    rows = [
        base,
        dict(base, age_scope="prenatal"),
        dict(base, age_scope="unresolved"),
        dict(base, donor_identity_known=False),
        dict(base, count_unit="assay_review"),
        dict(base, cell_type="neuron"),
        dict(base, cell_type="microglial cell"),
        dict(base, tissue_general="liver", cell_type="hepatocyte"),
        dict(base, cell_type="oligodendrocyte", n_cells=10**6),
    ]
    populations = pd.DataFrame(rows).assign(population_id=range(len(rows)))
    kept = panel(populations)
    assert kept.population_id.tolist() == [0]
    assert kept.tier.tolist() == ["vital"]
