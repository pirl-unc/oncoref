"""Published candidate choices must reconcile to every locus and evidence layer."""

import hashlib
import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
AUDIT = ROOT / "docs/audits/cta-compartment-policy-20261027"


def test_compartment_candidates_reconcile_and_retain_prior_selections():
    d = pd.read_csv(AUDIT / "protein-group-classifications.csv.gz", low_memory=False)
    summary = json.loads((AUDIT / "summary.json").read_text())
    assert d.protein_equivalence_group_id.is_unique and len(d) == 19863
    assert d.category.value_counts().to_dict() == summary["categories"]
    assert d.loc[d.baseline_category.eq("CTA"), "category"].eq("CTA").all()
    selected = d[d.category.eq("CTA")]
    assert not selected[[c for c in d if c.startswith("hold:")]].any(axis=None)
    assert not selected[[c for c in d if c.startswith("unsafe:")]].any(axis=None)
    crosswalk = pd.read_csv(AUDIT / "gene-to-protein-group-decisions.csv.gz", low_memory=False)
    assert crosswalk.gene_id.is_unique and len(crosswalk) == 20097
    assert crosswalk.protein_group_category.eq("CTA").sum() == summary["selected_member_loci"]
    assert "gene_category" not in crosswalk  # No stale locus result under the new policy.
    assert not any(c.startswith(("hold:", "unsafe:")) for c in crosswalk)
    for field in [
        "historical_ensembl_gene_ids",
        "alternate_ensembl_gene_ids",
        "gene_name_synonyms",
    ]:
        assert crosswalk[field].map(json.loads).map(type).eq(list).all()


def test_named_calibration_target_keeps_its_evidence_flags():
    d = pd.read_csv(AUDIT / "named-antigen-decisions.csv").set_index("protein_group_label")
    assert d.loc["PRAME", "category"] == "CTA"
    assert d.loc["PRAME", "limited_protein_review_flag"]
    assert d.loc["PRAME", "low_normal_protein_flag"]
    assert d.loc["PRAME", "stem_progenitor_expression_flag"]
    assert d.loc["CTAG2", "category"] == "Unsafe"
    assert d.loc["TBXT", "category"] == "CTA"
    # These are output regression checks; the report identifies PRAME as a
    # calibration example, not held-out clinical validation.


def test_compartment_report_receipts_match_implementation_and_outputs():
    for manifest in [AUDIT / "manifest.json", AUDIT / "cluster-replay.json"]:
        data = json.loads(manifest.read_text())
        for role in ["inputs", "outputs"]:
            for name, expected in data[role].items():
                source = ROOT / name if role == "inputs" else manifest.parent / name
                assert hashlib.sha256(source.read_bytes()).hexdigest() == expected, str(source)
