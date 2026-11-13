import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from audit_target_thresholds import (  # noqa: E402
    OUT,
    cohort_quantile,
    passes,
    patient_values,
    second_largest,
)


def test_two_independent_patients_are_not_two_aliquots_or_two_unknown_ids():
    values = pd.DataFrame({"a": [10, 12, 1, 100, 200]})
    samples = pd.DataFrame({"participant_id": ["p1", "p1", "p2", "", ""]})
    donors = patient_values(values, samples)
    assert donors.a.to_dict() == {"p1": 12, "p2": 1}
    assert second_largest(donors.to_numpy()).tolist() == [1]


def test_percentile_does_not_interpolate_one_outlier_into_a_supported_signal():
    values = pd.DataFrame({"a": [0] * 19 + [100], "sparse": [np.nan] + [100] * 19})
    samples = pd.DataFrame({"participant_id": [f"p{i}" for i in range(20)], "cancer_code": "A"})
    q = cohort_quantile(values, samples, 0.95)
    assert q.a == 0 and pd.isna(q["sparse"])
    values.loc[18, "a"] = 100
    assert cohort_quantile(values, samples, 0.95).a == 100


def test_enrichment_floor_and_equality_boundaries():
    assert passes(np.array([1, 2, 3]), np.array([0, 2, 2]), 1, 1).tolist() == [True, False, True]
    assert passes(np.array([9.9, 10, 19.9, 20]), np.array([0, 0, 2, 2]), 10, 1).tolist() == [
        False,
        True,
        False,
        True,
    ]


def test_all_scenarios_keep_set_algebra_and_unsafe_membership():
    grid = pd.read_csv(OUT / "threshold-grid.csv")
    assert len(grid) == 516 and grid.id.nunique() == 86
    assert grid.Unsafe.eq(18259).all()
    assert grid["Cancer Antigens"].eq(grid.CTA + grid["Cancer Antigen only"]).all()
    assert (
        grid["Target candidates"]
        .eq(grid.CTA + grid["Cancer Antigen only"] + grid["Reproductive Antigen only"])
        .all()
    )
    assert (grid["Target candidates"] + grid.Other + grid.Unsafe).eq(19893).all()
    selected = grid[grid.id.eq("fold_L1_R1_S10") & grid.disposable_cap_tpm.eq(32)].iloc[0]
    assert selected["Cancer Antigens"] == 880 and selected.CTA == 630
    assert selected["Target candidates"] == 910 and selected.Other == 724
    changes = pd.read_csv(OUT / "recommended-category-changes.csv")
    assert changes.groupby(["category", "proposed_category"]).size().to_dict() == {
        ("CTA", "Reproductive Antigen only"): 8,
        ("Cancer Antigen only", "Other"): 10,
    }
    moved_to_other = changes[changes.proposed_category.eq("Other")]
    assert moved_to_other.proposed_classification.eq("Other:Review").all()
    assert moved_to_other.proposed_cancer_evidence.eq("unassessed").all()


def test_tightening_cap_or_enrichment_never_adds_qualifiers():
    grid = pd.read_csv(OUT / "threshold-grid.csv")
    for _, rows in grid[grid.disposable_cap_tpm.notna()].groupby("id"):
        assert rows.sort_values("disposable_cap_tpm")["Target candidates"].is_monotonic_increasing
    rules = grid[
        grid.disposable_cap_tpm.isna()
        & grid.kind.eq("fold")
        & grid.single_fold.gt(0)
        & grid.repeated_fold.gt(0)
    ]
    for _, rows in rules.groupby(["floor", "repeated_fold"]):
        assert rows.sort_values("single_fold")["Cancer Antigens"].is_monotonic_decreasing


def test_disposable_caps_keep_native_tpm_and_ntpm_separate():
    evidence = pd.read_csv(OUT / "protein-threshold-evidence.csv.gz").set_index(
        "protein_group_label"
    )
    assert evidence.loc["LCE3C", "disposable_max_tpm"] == 21.4
    assert evidence.loc["LCE3C", "disposable_max_ntpm"] == 23.4
    native = pd.read_csv(OUT / "disposable-bulk-ntpm.csv")
    assert set(native.unit) == {"nTPM"}
    assert set(native.Tissue) == {"appendix", "gallbladder", "tonsil"}
    assert native.entity_id.nunique() == 957


def test_threshold_provenance_matches_committed_implementation_and_outputs():
    manifest = json.loads((OUT / "manifest.json").read_text())
    for section, base in [("implementation_sha256", ROOT), ("outputs_sha256", OUT)]:
        for path, expected in manifest[section].items():
            assert hashlib.sha256((base / path).read_bytes()).hexdigest() == expected, path
