import hashlib
import json
import sys
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from audit_target_expression_maxima import OUT, evaluate, maxima  # noqa: E402


def test_strict_target_maximum_and_independent_inclusive_cap():
    assert evaluate(100, 2, 16, cap=16) == {"maximum_status": "passes", "cap_status": "passes"}
    assert evaluate(100, 2, 17, cap=16) == {"maximum_status": "passes", "cap_status": "fails"}
    assert evaluate(16, 2, 16)["maximum_status"] == "fails"
    assert evaluate(16, 16, 2)["maximum_status"] == "fails"
    assert evaluate(0.9, 0, 0)["maximum_status"] == "fails"
    assert evaluate(1, 0, 0)["maximum_status"] == "passes"


def test_missing_values_are_not_zeros_and_invalid_measurements_are_rejected():
    for values in [(float("nan"), 2, 3), (4, float("nan"), 3), (4, 2, float("nan"))]:
        assert evaluate(*values, cap=16)["maximum_status"] == "unassessed"
    assert evaluate(4, 2, float("nan"), cap=16)["cap_status"] == "unassessed"
    for values in [(-1, 0, 0), (1, float("inf"), 0)]:
        with pytest.raises(ValueError):
            evaluate(*values)


def test_tied_witnesses_are_retained_and_maxima_not_added():
    rows = pd.DataFrame(
        [
            {
                "entity_id": "a",
                "scope": "disposable",
                "value": value,
                "source_id": "source",
                "source_context": tissue,
                "statistic": "mean TPM",
            }
            for value, tissue in [(2, "appendix"), (2, "tonsil"), (1, "gallbladder")]
        ]
    )
    result = maxima(rows, "disposable").loc["a"]
    assert result.disposable_max_tpm == 2
    assert len(json.loads(result.disposable_witnesses)) == 2


def test_counterfactual_report_counts_and_fingerprints():
    summary = json.loads((OUT / "summary.json").read_text())
    assert summary["maximum_status"] == {"passes": 920}
    table = pd.read_csv(OUT / "cap-scenarios.csv").set_index("disposable_cap_tpm")
    assert table.loc[4, "moved_to_other"] == 2
    assert table.loc[16, "moved_to_other"] == 1
    assert table.loc[32, "moved_to_other"] == 0
    records = pd.read_csv(OUT / "protein-maxima.csv")
    assert len(records) == records.entity_id.nunique() == 920
    assert records.disposable_max_tpm.notna().all()
    assert records[records.disposable_max_tpm.gt(16)].protein.tolist() == ["LCE3C"]
    assert summary["specific_registered_contexts_needing_hazard_review"] == []
    manifest = json.loads((OUT / "manifest.json").read_text())
    for section, base in [
        ("implementation_sha256", ROOT),
        ("inputs_sha256", ROOT),
        ("outputs_sha256", OUT),
    ]:
        for path, expected in manifest[section].items():
            if path == "oncoref/data/target-selection.csv.gz":
                # The audit read the 2026-11-09.1 draft catalog, preserved beside it.
                path = "docs/audits/target-selection-20261109/rebuild/packaged-catalog-2026-11-09.1.csv.gz"
            assert hashlib.sha256((base / path).read_bytes()).hexdigest() == expected, path
