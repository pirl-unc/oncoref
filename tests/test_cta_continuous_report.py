"""Pinned audit regressions protect clinical retention, gates and honest missingness."""

import hashlib
import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "docs/audits/cta-continuous-20261022"
BASE = ROOT / "docs/audits/cta-curation-20261021-r2"


def test_pilot_preserves_existing_candidates_and_never_invents_uncertainty():
    old = pd.read_csv(BASE / "candidate-decisions.csv.gz", keep_default_na=False)
    results = pd.read_csv(REPORT / "gene-results.csv", keep_default_na=False)
    passed = results.loc[results["Continuous point screen"].eq("Pass")]
    assert len(results) == 2540
    assert len(passed) == 758
    assert set(old.loc[old.retained, "Ensembl_Gene_ID"]) <= set(passed["Gene ID"])
    assert results.Change.eq("Gained").sum() == 12
    assert results.Change.eq("Lost").sum() == 0
    assert passed["Uncertainty result"].eq("Review").all()
    assert results["Joint mean-rule probability"].eq("").all()
    assert results["Joint probability including donor tail"].eq("").all()
    assert results["Posterior draws"].eq(0).all()
    assert not results["Uncertainty result"].eq("Included").any()
    clinical = pd.read_csv(BASE / "clinical-target-audit.csv")
    assert clinical.retained.sum() == 16
    assert set(clinical.loc[clinical.retained, "Ensembl_Gene_ID"]) <= set(passed["Gene ID"])
    held = set(old.loc[old.external_safety_hold, "Ensembl_Gene_ID"])
    assert results.loc[results["Gene ID"].isin(held), "Uncertainty result"].eq("Excluded").all()


def test_passing_candidates_obey_every_tissue_gate():
    results = pd.read_csv(REPORT / "gene-results.csv")
    evidence = pd.read_csv(REPORT / "tissue-evidence.csv.gz")
    ids = set(results.loc[results["Continuous point screen"].eq("Pass"), "Gene ID"])
    rows = evidence.loc[
        evidence.Ensembl_Gene_ID.isin(ids)
        & evidence.screening_tier.isin(["vital", "protected", "manageable"])
    ]
    assert rows.point_result.eq("Pass").all()
    assert rows.reproductive_margin_ntpm.gt(0).all()
    assert rows.absolute_margin_ntpm.ge(0).all()
    assert rows.continuous_selectivity_margin_ntpm.ge(0).all()
    assert not rows.disallowed_positive_ihc.any()
    assert not rows.unresolved_positive_ihc.any()
    assert rows.groupby("Ensembl_Gene_ID").size().eq(36).all()


def test_boundary_correction_and_full_report_provenance():
    sensitivity = pd.read_csv(REPORT / "sensitivity.csv")
    assert sensitivity["Pass point screen"].tolist() == [747, 752, 758, 769, 782, 790]
    manifest = json.loads((REPORT / "run-manifest.json").read_text())
    assert manifest["donor_model_fitted"] is False
    assert manifest["allowance_ntpm"] == 0.1
    for section, root in [("inputs", ROOT), ("implementation_sha256", ROOT), ("outputs", REPORT)]:
        for name, expected in manifest[section].items():
            assert hashlib.sha256((root / name).read_bytes()).hexdigest() == expected, name
