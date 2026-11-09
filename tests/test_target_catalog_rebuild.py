"""Check the rebuilt inventory against its measurements and retained provenance."""

import gzip
import hashlib
import json
import sys
from pathlib import Path

import pandas as pd
import pytest

from oncoref import target_selection as policy

ROOT = Path(__file__).parents[1]
OUT = ROOT / "docs/audits/target-selection-20261109/rebuild"
sys.path.insert(0, str(ROOT / "scripts"))
from rebuild_target_catalog import group_normal_rows, participant_id  # noqa: E402


@pytest.fixture(scope="module")
def catalog():
    return policy.results()


def test_complete_inventory_set_algebra_and_release_migrations(catalog):
    assert len(catalog) == 19893 and catalog.entity_id.is_unique
    genes = [g for ids in catalog.member_gene_ids for g in json.loads(ids)]
    assert len(genes) == len(set(genes)) == 20099
    assert catalog.category.value_counts().to_dict() == {
        "Unsafe": 18259,
        "Other": 724,
        "CTA": 632,
        "Cancer Antigen only": 250,
        "Reproductive Antigen only": 28,
    }
    assert catalog.is_cta.equals(catalog.is_cancer_antigen & catalog.is_reproductive_antigen)
    assert catalog.is_target_candidate.equals(
        catalog.is_cancer_antigen | catalog.is_reproductive_antigen
    )
    assert catalog.loc[catalog.is_target_candidate, list(policy.SHARED_CHECKS)].all().all()
    assert catalog.category.eq("Unsafe").equals(catalog.released_category.eq("Unsafe"))
    assert len(policy.results(changed_only=True)) == 343
    assert set(policy.results(category="Cancer Antigen only").category) == {"Cancer Antigen only"}
    with pytest.raises(ValueError, match="Unknown category"):
        policy.results(category="Cancer-only")


def test_named_assignments_and_accepted_review_holds(catalog):
    by_name = catalog.set_index("protein_group_label")
    for name in ["H2BW2", "HUS1B", "OR2B2", "TMEM269", "SSU72L2"]:
        assert by_name.loc[name, "category"] == "CTA"
    assert by_name.loc["EGFL6", "category"] == "Unsafe"
    assert set(
        catalog.loc[
            catalog.released_category.eq("CTA") & catalog.category.eq("Other"),
            "protein_group_label",
        ]
    ) == {"CALHM3", "PRAMEF18", "TRIM75", "TAF11L11", "TAF11L12", "TAF11L13", "TAF11L14"}
    assert by_name.loc["CDY1B", "cancer_evidence"] == "does_not_qualify"
    assert by_name.loc["CDY1B", "tumor_max_tpm"] < 1
    assert by_name.loc["CDY1B", "category"] == "Reproductive Antigen only"
    assert by_name.loc["FAM237A", "cancer_evidence"] == "unassessed"
    assert by_name.loc["FAM237A", "unresolved_patient_samples_above_normal"] >= 2
    assert by_name.loc["FAM237A", "category"] == "Reproductive Antigen only"
    gene = by_name.loc["H2BW2", "representative_gene_id"]
    assert policy.results(gene).protein_group_label.tolist() == ["H2BW2"]


def test_comparison_has_required_contexts_and_witnesses(catalog):
    summary = policy.metadata()["summary"]
    normal = pd.read_csv(OUT / "normal-comparison.csv.gz")
    assert set(normal.unit) == {"TPM"}
    normal = normal[normal.scope.eq("normal") & normal.qc_pass & normal.identity_resolved]
    required = set(summary["comparison_required_contexts"])
    maxima = normal.groupby("entity_id").value.max()
    assert len(maxima) == summary["shared_eligible"] == 957
    for _, group in normal.groupby("entity_id"):
        assert required <= set(group.context_id)
    assessed = catalog[catalog.normal_reference_tpm.notna()].set_index("entity_id")
    pd.testing.assert_series_equal(
        assessed.normal_reference_tpm.sort_index(), maxima.sort_index(), check_names=False
    )
    witnesses = pd.read_csv(OUT / "tumor-witnesses.csv.gz").fillna({"participant_id": ""})
    assert set(witnesses.unit) == {"TPM"}
    by_entity = dict(iter(witnesses.groupby("entity_id")))
    for row in catalog[catalog.is_cancer_antigen].itertuples():
        observations = by_entity[row.entity_id]
        above = observations[
            observations.value.ge(1) & observations.value.gt(row.normal_reference_tpm)
        ]
        patients = above.loc[above.participant_id.ne(""), "participant_id"].nunique()
        single = observations.value.ge(5 * max(row.normal_reference_tpm, 1)).any()
        assert patients >= 2 or single, row.protein_group_label
    samples = pd.read_csv(OUT / "tumor-samples.csv.gz")
    assert samples.sample_key.is_unique
    assert not samples.duplicated(["namespace", "sample_id"]).any()


def test_missing_comparison_is_not_reported_as_zero(catalog):
    failed = ~catalog[list(policy.SHARED_CHECKS)].all(axis=1) | catalog.category.eq("Unsafe")
    assert catalog.loc[failed, "normal_reference_tpm"].isna().all()
    assert catalog.loc[failed, "tumor_max_tpm"].isna().all()
    coverage = pd.read_csv(OUT / "tumor-source-coverage.csv")
    assert coverage.status.eq("assessed").sum() == 130
    nutm = coverage[coverage.cancer_code.eq("NUTM")]
    assert len(nutm) == 1 and not nutm.status.eq("assessed").any()


def test_paper_nominations_preserve_exact_members_and_full_claim_ledger(catalog):
    with gzip.open(OUT / "paper-attestations.json.gz", "rt") as handle:
        links = json.load(handle)
    for row in catalog.itertuples():
        ids = json.loads(row.paper_attestation_indices)
        genes = set(json.loads(row.member_gene_ids))
        claims = [links[i] for i in ids]
        assert all(genes & set(c["gene_ids"]) for c in claims)
        nominees = [c for c in claims if c["published_nomination"]]
        assert bool(nominees) == row.published_nomination
        assert {c["source_id"] for c in nominees} == set(json.loads(row.nomination_source_ids))
        covered = set().union(*(set(c["gene_ids"]) for c in nominees)) & genes
        assert covered == set(json.loads(row.nomination_member_ids))
        assert row.nomination_scope == (
            "all_members" if covered == genes else "some_members" if covered else "none"
        )
    assert catalog[catalog.is_cta].published_nomination.sum() == 466


def test_packaged_snapshot_matches_report_metadata(catalog):
    meta = policy.metadata()
    path = ROOT / "oncoref/data/target-selection.csv.gz"
    assert hashlib.sha256(path.read_bytes()).hexdigest() == meta["data_sha256"]
    assert meta["summary"] == json.loads((OUT / "summary.json").read_text())
    assert meta["comparison"]["normal_sources"] == json.loads(
        (OUT / "downloaded-hpa-sources.json").read_text()
    )
    assert catalog.policy.eq(policy.POLICY_ID).all()
    for category in policy.CATEGORIES:
        subset = pd.read_csv(OUT / (category.lower().replace(" ", "-") + ".csv.gz"))
        assert set(subset.entity_id) == set(catalog.loc[catalog.category.eq(category), "entity_id"])


def test_report_provenance_manifest():
    manifest = json.loads((OUT / "manifest.json").read_text())
    for paths, base in [
        (manifest["implementation_sha256"], ROOT),
        (manifest["outputs_sha256"], OUT),
    ]:
        for path, digest in paths.items():
            assert hashlib.sha256((base / path).read_bytes()).hexdigest() == digest, path


def test_same_context_sum_and_missing_member_coverage():
    facts = pd.DataFrame([{"entity_id": "AB", "member_gene_ids": '["A", "B"]', "member_count": 2}])
    raw = pd.DataFrame(
        [
            ["A", "lung", "lung", 9],
            ["B", "lung", "lung", 1],
            ["A", "brain", "brain", 1],
            ["B", "brain", "brain", 9],
            ["A", "heart", "heart", 9],
        ],
        columns=["gene_id", "tissue", "source_context", "value"],
    )
    result = group_normal_rows(raw, facts, source="test", statistic="mean").set_index("tissue")
    assert result.loc["lung", "value"] == result.loc["brain", "value"] == 10
    assert pd.isna(result.loc["heart", "value"])
    with pytest.raises(ValueError, match="Duplicate normal"):
        group_normal_rows(pd.concat([raw, raw.iloc[:1]]), facts, source="test", statistic="mean")


@pytest.mark.parametrize(
    "sample,expected",
    [
        ("TCGA-AB-1234-01A", "TCGA:TCGA-AB-1234"),
        ("TCGA-AB-1234-02B", "TCGA:TCGA-AB-1234"),
        ("TARGET-10-ABC123-01A", "TARGET:TARGET-10-ABC123"),
        ("MMRF_1234_1_BM", "MMRF:1234"),
        ("SRR123456", ""),
        ("Patient1_sample2", ""),
    ],
)
def test_documented_patient_identifiers_only(sample, expected):
    assert participant_id("source", sample) == expected
