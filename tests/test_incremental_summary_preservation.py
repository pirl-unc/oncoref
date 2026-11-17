"""A targeted rebuild must preserve other cohorts and its hardlinked base."""

import hashlib
import json

import pandas as pd
from scripts import merge_expression_artifact_update as update


def test_incremental_update_preserves_shared_source_and_base(tmp_path, monkeypatch):
    base, rebuilt, out = [tmp_path / name for name in ("base", "rebuilt", "out")]
    for root, codes in [(base, ["KEEP", "TGCT"]), (rebuilt, ["TGCT"])]:
        (root / "cancer-reference-expression").mkdir(parents=True)
        (root / update.PROVENANCE_PATH).parent.mkdir(parents=True)
        pd.DataFrame({"cancer_code": codes, "value": range(len(codes))}).to_csv(
            root / update.PROVENANCE_PATH, index=False
        )
        pd.DataFrame({"cancer_code": codes}).to_csv(root / update.SAMPLE_QC_PATH, index=False)
        pd.DataFrame(
            {
                "cancer_code": codes,
                "sample_qc_policy_version": ["v4"] * len(codes),
                "n_source_samples": [1] * len(codes),
                "n_cohort_samples": [1] * len(codes),
                "n_negative_values_clipped": [0] * len(codes),
                "sample_qc_fallback_reason": [""] * len(codes),
            }
        ).to_csv(root / update.BUILD_METADATA_CSV, index=False)
        (root / update.BUILD_METADATA_JSON).write_text(json.dumps({}))
    path = "cancer-reference-expression/SHARED.csv.gz"
    pd.DataFrame({"cancer_code": ["KEEP", "TGCT"], "value": [7.0, 9.0]}).to_csv(
        base / path, index=False
    )
    pd.DataFrame({"cancer_code": ["TGCT"], "value": [3.0]}).to_csv(rebuilt / path, index=False)
    hashes = {
        p.relative_to(base): hashlib.sha256(p.read_bytes()).hexdigest()
        for p in base.rglob("*")
        if p.is_file()
    }
    monkeypatch.setattr(update, "assign_representative_partitions", lambda frame: frame)
    monkeypatch.setattr(
        update, "representative_partition_cohort_metadata", lambda frame: frame[["cancer_code"]]
    )
    monkeypatch.setattr(update, "representative_partition_build_metadata", lambda frame: {})
    assert update.merge_expression_artifact_update(base, rebuilt, out) == {"TGCT"}
    result = pd.read_csv(out / path).set_index("cancer_code").value.to_dict()
    assert result == {"KEEP": 7.0, "TGCT": 3.0}
    assert hashes == {
        p.relative_to(base): hashlib.sha256(p.read_bytes()).hexdigest()
        for p in base.rglob("*")
        if p.is_file()
    }
