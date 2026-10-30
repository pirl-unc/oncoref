"""Tumor restriction must use independent, resolved native-TPM observations."""

import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from build_tumor_restricted_support import sample_key  # noqa: E402


def test_duplicate_views_and_tcga_aliquots_do_not_supply_two_tumors():
    assert sample_key("view-a", "TCGA-AB-1234-01A") == sample_key("view-b", "TCGA-AB-1234-02B")
    assert sample_key("view-a", "TCGA-AB-1234-01A") != sample_key("view-a", "TCGA-AB-5678-01A")
    assert sample_key("study", "sample1") != sample_key("study", "sample2")


def test_packaged_leads_have_two_distinct_resolved_native_tpm_observations():
    folder = ROOT / "docs/audits/cta-selection-20261028"
    facts = pd.read_csv(folder / "protein-group-classifications.csv.gz", low_memory=False)
    obs = pd.read_csv(folder / "tumor-restricted-observations.csv.gz")
    leads = facts[facts.tumor_restricted_lead]
    assert len(leads) == 66
    assert leads.category.eq("Other").all()
    assert not leads.normal_restriction_reviewed.any()
    obs = obs[obs.protein_equivalence_group_id.isin(leads.protein_equivalence_group_id)]
    assert obs.unit.eq("TPM").all()
    assert obs.value.ge(10).all()
    assert obs.identity_resolved.all()
    assert obs.groupby("protein_equivalence_group_id").sample_key.nunique().eq(2).all()
    assert obs.protein_equivalence_group_id.nunique() == len(leads)
    receipt = json.loads((folder / "tumor-restricted-manifest.json").read_text())
    assert receipt["threshold_tpm"] == 10
    assert receipt["distinct_samples_required"] == 2
