"""Release evidence, package decisions and presentation counts must reconcile."""

import hashlib
import json
from pathlib import Path

import pandas as pd

from oncoref import cta_selection

ROOT = Path(__file__).resolve().parents[1]
AUDIT = ROOT / "docs/audits/cta-release-1.21.0"


def test_release_receipts_verify_all_owned_evidence():
    for filename in ("manifest.json", "cluster-replay.json"):
        receipt = json.loads((AUDIT / filename).read_text())
        for role in ("inputs", "outputs"):
            for name, expected in receipt[role].items():
                path = (ROOT if role == "inputs" else AUDIT) / name
                assert hashlib.sha256(path.read_bytes()).hexdigest() == expected, str(path)


def test_frozen_release_counts_match_funnel_and_presentation():
    selected = pd.read_csv(AUDIT / "protein-group-classifications.csv.gz", low_memory=False)
    summary = json.loads((AUDIT / "summary.json").read_text())
    deck = json.loads((AUDIT / "presentation-data.json").read_text())
    assert selected.category.value_counts().to_dict() == summary["categories"]
    assert deck["summary"] == summary
    assert (
        deck["report_sha256"] == hashlib.sha256((AUDIT / "manifest.json").read_bytes()).hexdigest()
    )
    funnel = pd.read_csv(AUDIT / "filtering-funnel.csv")
    assert funnel[["CTA_eligible_so_far", "Unsafe", "Neither"]].sum(axis=1).eq(len(selected)).all()
    assert funnel.new_candidate_removals.sum() + summary["categories"]["CTA"] == len(selected)
    final = funnel.iloc[-1]
    assert [final.CTA_eligible_so_far, final.Unsafe, final.Neither] == [
        summary["categories"][category] for category in cta_selection.CATEGORIES
    ]
    known = pd.read_csv(AUDIT / "known-cta-decisions.csv.gz", low_memory=False)
    assert (
        known.loc[known.historical_CTA, "protein_group_category"].value_counts().to_dict()
        == deck["historical_loci"]
    )
    assert "Expression category" not in known
    assert "historical_expression_category" in known


def test_current_named_normal_evidence_keeps_units_and_margin():
    evidence = pd.read_csv(AUDIT / "named-normal-evidence.csv")
    assert not evidence.duplicated(["protein_group_label", "assay"]).any()
    assert evidence.unit.isin(["nTPM", "TPM", "CP10K"]).all()
    assert (evidence.strong_boundary == evidence.candidate_ceiling * 2).all()
    bronchial = evidence[
        (evidence.protein_group_label == "CTAG2") & (evidence.assay == "HPA tissue clusters")
    ].iloc[0]
    assert "basal" in bronchial.cell_type.lower()
    assert bronchial.value == 24.9 and bronchial.candidate_ceiling == 4
    assert bronchial.category_decision == "Unsafe"
    protein = pd.read_csv(AUDIT / "imported-protein-review-decisions.csv")
    adrenal = protein.loc[protein["Review ID"].eq("review-20261025-PRAME-adrenal")].iloc[0]
    assert adrenal.protein_level == "Low" and adrenal.protein_scale == "author_description"
    assert adrenal.low_normal_protein_flag and not adrenal.retained_protein_exclusion
