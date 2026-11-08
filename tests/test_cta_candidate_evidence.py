"""Normal findings and prevalence preserve assay, copy and sample identity."""

import copy
import gzip
import hashlib
import json
from pathlib import Path

import pandas as pd
import pytest
from build_cta_candidate_evidence import count_samples, replay_updates, tumor_tables

from oncoref import gene_knowledge

ROOT = Path(__file__).resolve().parents[1]
AUDIT = ROOT / "docs/audits/cta-candidate-evidence-20261108"


def test_pending_curations_are_explicit_and_ambiguous_assays_cannot_exclude():
    assert gene_knowledge.findings("PRAMEF18").empty
    reviewed = gene_knowledge.findings("PRAMEF18", include_updates=True)
    assert len(reviewed) == 2 and reviewed.effect.eq("review").all()
    assert reviewed.interpretation_status.eq("unresolved").all()
    for record in reviewed.to_dict("records"):
        wrong = copy.deepcopy(record)
        wrong["effect"] = "exclude"
        with pytest.raises(ValueError, match="Unresolved or rejected"):
            gene_knowledge.validate([wrong])
    rejected = gene_knowledge.findings("UBTFL1", include_updates=True).iloc[0]
    assert rejected.effect == "context" and rejected.interpretation_status == "rejected"
    assert gene_knowledge.findings("TRIM75", taxonomy_id=10090, include_updates=True).empty


def test_complete_funnel_replay_holds_seven_proteins_without_new_unsafe_calls():
    baseline = pd.read_csv(
        ROOT / "docs/audits/cta-selection-20261028/protein-group-classifications.csv.gz",
        low_memory=False,
    )
    result = replay_updates(baseline).classifications
    changed = result[result.category.ne(baseline.category)]
    assert set(changed.protein_group_label) == {
        "CALHM3",
        "PRAMEF18",
        "TRIM75",
        "TAF11L11",
        "TAF11L12",
        "TAF11L13",
        "TAF11L14",
    }
    assert changed.classification.eq("Other:Review").all()
    assert changed.assessment_status.eq("normal_evidence_hold").all()
    pd.testing.assert_series_equal(
        result.category.eq("Unsafe"), baseline.category.eq("Unsafe"), check_dtype=False
    )
    family = result[result.protein_group_label.str.startswith("TAF11L")]
    assert family.knowledge_normal_hold.all()
    assert family[family.category.ne("Unsafe")].classification.eq("Other:Review").all()
    # High developmental/tumor evidence remains present and cannot erase holds.
    assert changed.positive_context.all() and changed.cancer_association.all()


def test_sample_namespaces_aliases_quality_and_missingness():
    raw = pd.DataFrame(
        [
            ["g", "OLD", "x", 10.0, "a", "source", "shared", "1"],
            ["g", "NEW", "x", 10.0, "b", "source", "shared", "1"],
            ["g", "NEW", "bad", 100.0, "a", "source", "shared", "1"],
            ["g", "NEW", "zero", 0.0, "a", "source", "shared", "1"],
        ],
        columns=[
            "Ensembl_Gene_ID",
            "Symbol",
            "sample",
            "tpm",
            "cohort",
            "source_cohort",
            "namespace",
            "matrix_version",
        ],
    )
    qc = pd.DataFrame(
        [
            ["a", "x", True, True],
            ["b", "x", True, True],
            ["a", "bad", False, True],
            ["a", "zero", True, True],
        ],
        columns=[
            "cancer_code",
            "sample_id",
            "recommended_for_absolute_tpm_floor",
            "linear_tpm_comparable",
        ],
    )
    coverage = pd.DataFrame(
        [
            ["g", "a", "measured", "a"],
            ["g", "b", "measured", "b"],
            ["g", "c", "gene_absent", "c"],
        ],
        columns=["gene_id", "cohort", "status", "matrix_asset"],
    )
    candidates = [{"gene_id": "g", "protein": "NEW"}]
    counts, cohorts, _ = tumor_tables(raw, qc, candidates, coverage)
    r = counts.iloc[0]
    assert r.measured_samples == 2 and r.ge_10 == 1 and r.maximum_tpm == 10
    assert r.observed_sample_records == 3 and r.excluded_sample_records == 1
    absent = cohorts[cohorts.cohort.eq("c")].iloc[0]
    assert absent.measured_samples == 0 and pd.isna(absent.ge_1)
    conflict = raw.copy()
    conflict.loc[1, "tpm"] = 11
    with pytest.raises(ValueError, match="disagree"):
        tumor_tables(conflict, qc, candidates, coverage)
    with pytest.raises(ValueError, match="Duplicate"):
        tumor_tables(pd.concat([raw, raw.iloc[:1]]), qc, candidates, coverage)
    with pytest.raises(ValueError, match="finite"):
        count_samples(pd.Series([1.0, float("nan")]))


def test_ssu72l2_counts_and_no_false_protein_negatives():
    counts = pd.read_csv(AUDIT / "candidates.csv").set_index("protein")
    row = counts.loc["SSU72L2"]
    assert (row.measured_samples, row.ge_1, row.ge_5, row.ge_10) == (14143, 50, 6, 3)
    assert row.maximum_tpm == pytest.approx(990.91)
    assert counts.hpa_ihc.eq("Not assayed").sum() == 9
    assert counts.hpa_ms.eq("Not available").all()
    positives = pd.read_csv(AUDIT / "ssu72l2-positive-samples.csv")
    assert positives["sample"].nunique() == 50 and positives.eligible.all()
    peak = positives.sort_values("tpm").iloc[-1]
    assert peak["sample"] == "THR24_2154_S01"


def test_primer_conflict_and_audit_receipts():
    inputs = AUDIT / "inputs"
    receipt = json.loads((inputs / "UBTFL1-primer-check.json").read_text())
    seq = json.loads((inputs / "UBTFL1.cdna.json").read_text())["seq"]
    complement = str.maketrans("ACGT", "TGCA")
    for primer in receipt["primers"]:
        assert primer not in seq and primer.translate(complement)[::-1] not in seq
    mouse = json.loads((inputs / "Sod1-mouse-primer-targets.json").read_text())
    for seq in mouse.values():
        for primer in receipt["primers"]:
            assert primer in seq or primer.translate(complement)[::-1] in seq
    assert b"Sod1" in gzip.decompress((inputs / "XM_021185540.2.gb.gz").read_bytes())
    manifest = json.loads((AUDIT / "manifest.json").read_text())
    for section, root in [("inputs", ROOT), ("outputs", AUDIT)]:
        for name, digest in manifest[section].items():
            assert hashlib.sha256((root / name).read_bytes()).hexdigest() == digest, name
    # Saved compressed responses must be actual HPA XML, not retrieval error pages.
    for path in inputs.glob("*.xml.gz"):
        assert b"<proteinAtlas" in gzip.decompress(path.read_bytes())
