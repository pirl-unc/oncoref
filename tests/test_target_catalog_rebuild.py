"""Check the rebuilt inventory against its measurements and retained provenance."""

import gzip
import hashlib
import json
import sys
from pathlib import Path

import pandas as pd
import pytest

from oncoref import samples as sample_api
from oncoref import target_selection as policy

ROOT = Path(__file__).parents[1]
OUT = ROOT / "docs/audits/target-selection-20261110/rebuild"
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
        "Unsafe": 18253,
        "Other": 1116,
        "CTA": 327,
        "Cancer Antigen only": 70,
        "Reproductive Antigen only": 127,
    }
    assert catalog.is_cta.equals(catalog.is_cancer_antigen & catalog.is_reproductive_antigen)
    assert catalog.is_target_candidate.equals(
        catalog.is_cancer_antigen | catalog.is_reproductive_antigen
    )
    assert catalog.loc[catalog.is_target_candidate, list(policy.SHARED_CHECKS)].all().all()
    # Only the rare-cell regrade and the literature review change Unsafe membership.
    unsafe, released = catalog.category.eq("Unsafe"), catalog.released_category.eq("Unsafe")
    assert (released & ~unsafe).sum() == 86
    assert catalog.loc[released & ~unsafe, "legacy_exclusion_reasons"].str.contains("Rare").all()
    assert catalog.loc[unsafe & ~released, "Unsafe_reasons"].value_counts().to_dict() == {
        "Established rare normal-cell RNA": 67,
        "Published normal-tissue expression": 13,
    }
    # Strong unresolved signals are held; accepted exclusions still take precedence.
    assert set(catalog.rare_cell_status) <= {*policy.RARE_CELL_STATUSES, "not_extracted"}
    assert set(catalog.rare_cell_action) == {"none", "unsafe", "review"}
    assert not catalog.loc[catalog.rare_cell_action.eq("review"), "is_target_candidate"].any()
    assert catalog.loc[catalog.rare_cell_status.eq("established"), "category"].eq("Unsafe").all()
    assert len(policy.results(changed_only=True)) == 461
    assert set(policy.results(category="Cancer Antigen only").category) == {"Cancer Antigen only"}
    with pytest.raises(ValueError, match="Unknown category"):
        policy.results(category="Cancer-only")


def test_named_assignments_and_accepted_review_holds(catalog):
    by_name = catalog.set_index("protein_group_label")
    # Reproductive coexpression is advisory, not evidence against tumor origin.
    for name in [
        "ZNRF4",
        "PRM2",
        "CCDC70",
        "BPIFA3",
        "CYLC2",
        "SPATA19",
        "PGK2",
        "CETN1",
        "TEX38",
        "PRR30",
        "TSSK2",
        "IQCF5",
        "BOD1L2",
    ]:
        assert by_name.loc[name, "category"] == "CTA"
        assert by_name.loc[name, "tumor_top_sample"] == "TREEHOUSE_POLYA_25_01:TCGA-2G-AAFL-01"
    for name in [
        "MAGEC2",
        "PRDM14",
        "TMEM202",
        "SSU72L2",
        "SSX2/B",
        "PRAMEF1",
        "PRAMEF2",
        "CSAG1",
    ]:
        assert by_name.loc[name, "category"] == "CTA"
    assert 681 < by_name.loc["MAGEC2", "tumor_fold"] < 682
    assert by_name.loc["GOLGA6D", "classification"] == "CTA"
    assert by_name.loc["DAZ1", "category"] == by_name.loc["DAZ4", "category"] == "Unsafe"
    assert (
        by_name.loc["DAZ2", "classification"]
        == by_name.loc["DAZ3", "classification"]
        == "Other:Sequence"
    )
    assert by_name.loc["MAGEA4", "category"] == "CTA"
    for name in ["HUS1B", "ZNF705D"]:
        assert by_name.loc[name, "category"] == "Cancer Antigen only"
        assert by_name.loc[name, "reproductive_fold"] < policy.REPRODUCTIVE_FOLD
    assert by_name.loc["OR2B2", "classification"] == "Other:Unscreened tissue"
    for name in ["DAZL", "CTAG2", "DCST1", "NKX2-6"]:
        assert by_name.loc[name, "classification"] == "Other:Review"
        assert by_name.loc[name, "rare_cell_status"] == "unresolved"
    # No tumor reaches 10x the normal maximum, so the cancer claim is not made.
    for name in ["H2BW2", "EXD1"]:
        assert by_name.loc[name, "category"] == (
            "Other" if name == "H2BW2" else "Reproductive Antigen only"
        )
        assert by_name.loc[name, "tumor_fold"] < policy.TUMOR_FOLD
    # PRAME is in the eye of two donors in two studies, 1.2-1.5% of cells in different
    # populations. Two lone donors that disagree on the cell type are a recorded warning.
    for name in ["PRAME", "EXD1", "LRRIQ4"]:
        assert by_name.loc[name, "is_reproductive_antigen"]
        assert by_name.loc[name, "rare_cell_status"] == "dismissed"
        assert "different cell types" in by_name.loc[name, "rare_cell_basis"]
        assert by_name.loc[name, "finding_donors_multi_molecule"] == 2
    assert by_name.loc["PRAME", "finding_tissue"] == "eye"
    # Replicated multi-molecule RNA at a locus without mapping risk.
    for name in ["PLAC1", "ADAM7", "CBY3"]:
        assert by_name.loc[name, "category"] == "Unsafe"
        assert by_name.loc[name, "rare_cell_status"] == "established"
    # A mapping-risk flag alone cannot dismiss strong normal-cell RNA.
    for name, note in [("TAS2R1", "same-strand gene"), ("TMEM269", "same-strand gene")]:
        assert by_name.loc[name, "rare_cell_status"] == "unresolved"
        assert note in by_name.loc[name, "rare_cell_mapping_notes"]
        assert by_name.loc[name, "literature_effect"] == "none"
        assert by_name.loc[name, "classification"] == "Other:Review"
    # Published brain RNA turns MAGEA12's single-molecule detection into a finding.
    assert by_name.loc["MAGEA12", "max_multi_molecule_pct"] < 1
    assert "published" in by_name.loc["MAGEA12", "rare_cell_basis"]
    assert by_name.loc["MAGEA12", "category"] == "Unsafe"
    assert by_name.loc["MYH13", "literature_effect"] == "exclude"
    assert by_name.loc["MYH13", "category"] == "Unsafe"
    # The candidate sweep: a quoted report in a Vital or Protected tissue excludes.
    for name, tissue in [("MRGPRX1", "dorsal root ganglion"), ("SCGB1C1", "nasal mucosa")]:
        assert by_name.loc[name, "literature_effect"] == "exclude"
        assert by_name.loc[name, "literature_decision_tissue"] == tissue
        assert by_name.loc[name, "Unsafe_reasons"] == "Published normal-tissue expression"
    # Olfactory epithelium is Manageable, so its finding is recorded, and the protein has
    # no reproductive evidence for a home the panels do not resolve.
    assert by_name.loc["BPIFB3", "literature_effect"] == "none"
    assert by_name.loc["BPIFB3", "literature_recorded_tier"] == "manageable"
    assert by_name.loc["BPIFB3", "classification"] == "Other:Unscreened tissue"
    # A confirmed finding below the rubric, and a family rule, are recorded only.
    assert by_name.loc["LRRC52", "literature_recorded_tier"] == "vital"
    assert by_name.loc["LRRC52", "classification"] == "CTA"
    assert by_name.loc["OR2B2", "literature_family"] == "olfactory_receptor"
    assert by_name.loc["OR2B2", "literature_family_tier"] == "manageable"
    dismissed = catalog[catalog.is_target_candidate & catalog.rare_cell_status.eq("dismissed")]
    released = catalog[catalog.is_target_candidate & catalog.released_category.eq("Unsafe")]
    assert pd.concat([dismissed, released]).literature_effect.eq("none").all()
    held = set(
        catalog.loc[
            catalog.released_category.eq("CTA") & catalog.category.eq("Other"),
            "protein_group_label",
        ]
    )
    # New enrichment failures and normal-evidence holds coexist with read-through exclusions.
    readthrough = set(
        catalog.loc[
            catalog.released_category.eq("CTA") & catalog.readthrough_locus, "protein_group_label"
        ]
    )
    assert len(held) == 98 and len(held - readthrough) == 92
    assert {"ATF7-NPFF", "ENSG00000273155"} <= readthrough <= held
    assert held >= {"CALHM3", "PRAMEF18", "TRIM75", "TAF11L11", "TAF11L12", "TAF11L13", "TAF11L14"}
    assert by_name.loc["CDY1B", "cancer_evidence"] == "does_not_qualify"
    assert by_name.loc["CDY1B", "tumor_max_tpm"] < 1
    assert by_name.loc["CDY1B", "classification"] == "Other:Expression"
    # Unresolved patient identity no longer holds a protein: one sample decides.
    assert by_name.loc["FAM237A", "cancer_evidence"] == "does_not_qualify"
    assert by_name.loc["FAM237A", "unresolved_patient_samples_above_normal"] == 2
    assert by_name.loc["FAM237A", "tumor_fold"] < policy.TUMOR_FOLD
    assert by_name.loc["FAM237A", "classification"] == "Other:Expression"
    gene = by_name.loc["H2BW2", "representative_gene_id"]
    assert policy.results(gene).protein_group_label.tolist() == ["H2BW2"]


def test_comparison_has_required_contexts_and_witnesses(catalog):
    summary = policy.metadata()["summary"]
    normal = pd.read_csv(OUT / "normal-comparison.csv.gz")
    assert set(normal.unit) == {"TPM"}
    normal = normal[normal.scope.eq("normal") & normal.qc_pass & normal.identity_resolved]
    required = set(summary["comparison_required_contexts"])
    maxima = normal.groupby("entity_id").value.max()
    assert len(maxima) == summary["assessed_for_cancer_evidence"] == 918
    assert summary["shared_eligible"] == 841
    assert summary["expression_bounds"]["fail"] == 22
    failed = catalog[catalog.expression_bounds_status.eq("fail")]
    assert (
        not failed.is_target_candidate.any()
        and failed.classification.isin(["Other:Expression", "Other:Review"]).all()
    )
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
        floor = max(row.normal_reference_tpm, row.disposable_max_tpm, 1)
        assert observations.value.ge(policy.TUMOR_FOLD * floor).any(), row.protein_group_label
    samples = pd.read_csv(OUT / "tumor-samples.csv.gz")
    assert samples.sample_key.is_unique
    assert not samples.duplicated(["namespace", "sample_id"]).any()


def test_missing_comparison_is_not_reported_as_zero(catalog):
    prior = [c for c in policy.SHARED_CHECKS if c != "expression_bounds_pass"]
    failed = ~catalog[prior].all(axis=1) | catalog.category.eq("Unsafe")
    # Newly unresolved rare-cell records retain their measured counterfactuals.
    failed &= ~catalog.rare_cell_action.eq("review")
    assert catalog.loc[failed, "normal_reference_tpm"].isna().all()
    assert catalog.loc[failed, "tumor_max_tpm"].isna().all()
    coverage = pd.read_csv(OUT / "tumor-source-coverage.csv")
    assert coverage.status.eq("assessed").sum() == 130
    germ = coverage[coverage.cancer_code.eq("TGCT")]
    assert len(germ) == 1 and germ.status.eq("assessed").all()
    assert germ.samples.iloc[0] == 154
    nutm = coverage[coverage.cancer_code.eq("NUTM")]
    assert len(nutm) == 1 and not nutm.status.eq("assessed").any()


def test_tenfold_witness_tables_are_complete_and_reproductive_evidence_is_symmetric(catalog):
    witness = pd.read_csv(OUT / "tumor-witnesses.csv.gz", keep_default_na=False)
    passing = witness[witness.single_sample_pass]
    assert not passing.duplicated(["entity_id", "source_id", "sample_id"]).any()
    sample_counts = passing.groupby("entity_id").size()
    known = passing[passing.participant_id.ne("")].groupby("entity_id").participant_id.nunique()
    anonymous = passing[passing.participant_id.eq("")].groupby("entity_id").size()
    measured = catalog[catalog.normal_reference_tpm.notna()].set_index("entity_id")
    for column, counts in [
        ("samples_at_tenfold", sample_counts),
        ("independent_patients_at_tenfold", known),
        ("anonymous_samples_at_tenfold", anonymous),
    ]:
        assert measured[column].eq(counts.reindex(measured.index, fill_value=0)).all()
    ra = catalog[catalog.is_reproductive_antigen]
    assert ra.reproductive_fold.ge(policy.REPRODUCTIVE_FOLD).all()
    denominator = ra[["normal_reference_tpm", "disposable_max_tpm"]].max(axis=1).clip(lower=1)
    assert (ra.reproductive_supported_tpm / denominator - ra.reproductive_fold).abs().lt(1e-8).all()
    reproductive = pd.read_csv(OUT / "reproductive-witnesses.csv.gz")
    assert set(reproductive.unit) == {"TPM"}
    assert reproductive.single_sample_pass.all()
    assert set(ra.entity_id) <= set(reproductive.entity_id)
    development = reproductive[reproductive.scope.eq("developmental")]
    small = development[development.required_embryos.eq(2)]
    well_sampled = development[development.required_embryos.eq(1)]
    assert len(small) + len(well_sampled) == len(development)
    assert small.cells.ge(5).all()
    assert (small.positive_cells / small.cells).ge(0.5).all()
    keys = ["entity_id", "source_id", "developmental_stage", "source_context"]
    donors = development.groupby(keys).participant_id.nunique().rename("supporting_donors")
    assert small.join(donors, on=keys).supporting_donors.ge(2).all()
    assert well_sampled.cells.ge(20).all() and well_sampled.positive_cells.ge(5).all()
    assert (well_sampled.positive_cells / well_sampled.cells).ge(0.1).all()
    assert development.participant_id.notna().all()
    normal = pd.read_csv(OUT / "normal-comparison.csv.gz")
    individual = normal[normal.observation_kind.eq("individual_tissue_sample")]
    assert individual.sample_id.nunique() == 54
    assert set(individual.scope) == {"reproductive"}
    percentiles = normal[normal.observation_kind.eq("biopsy_percentile")]
    assert set(percentiles.scope) == {"normal", "disposable"}
    assert percentiles.sample_id.isna().all()
    assert individual.participant_id.isna().all()  # HPA does not supply donor identity.


def test_release_flag_means_scientific_review_completion():
    metadata = policy.metadata()
    assert metadata["release_ready"] is False
    assert "not software publication" in metadata["release_ready_meaning"]


def test_historical_threshold_generator_rejects_live_catalog_before_writing(monkeypatch, tmp_path):
    import build_cta_threshold_review as historical

    output = tmp_path / "historical"
    monkeypatch.setattr(historical, "OUT", output)
    with pytest.raises(RuntimeError, match=r"dated 1\.38\.0"):
        historical.main()
    assert not output.exists()


def test_paper_nominations_preserve_exact_members_and_full_claim_ledger(catalog):
    with gzip.open(OUT / "paper-attestations.json.gz", "rt") as handle:
        links = json.load(handle)
    for row in catalog.itertuples():
        ids = json.loads(row.paper_attestation_indices)
        genes = set(json.loads(row.member_gene_ids))
        claims = [links[i] for i in ids]
        assert all(genes & set(c["gene_ids"]) for c in claims)
        nominees = [c for c in claims if genes & set(c["nomination_gene_ids"])]
        assert bool(nominees) == row.published_nomination
        assert {c["source_id"] for c in nominees} == set(json.loads(row.nomination_source_ids))
        covered = set().union(*(set(c["nomination_gene_ids"]) for c in nominees)) & genes
        assert covered == set(json.loads(row.nomination_member_ids))
        assert row.nomination_scope == (
            "all_members" if covered == genes else "some_members" if covered else "none"
        )
    assert catalog[catalog.is_cta].published_nomination.sum() == 252
    assert catalog[catalog.category.eq("Cancer Antigen only")].published_nomination.sum() == 28


def test_packaged_snapshot_matches_report_metadata(catalog):
    meta = policy.metadata()
    path = ROOT / "oncoref/data/target-selection.csv.gz"
    assert hashlib.sha256(path.read_bytes()).hexdigest() == meta["data_sha256"]
    assert meta["summary"] == json.loads((OUT / "summary.json").read_text())
    from oncoref import reference_data

    sources = meta["comparison"]["normal_sources"]
    assert {row["source"] for row in sources} == {
        "rna_tissue_hpa",
        "rna_tissue_gtex",
        "rna_brain_hpa",
        "rna_brain_gtex",
    }
    for row in sources:
        pin = reference_data.REFERENCE_SOURCES[row["source"]]["pins"][row["version"]]
        assert all(row[key] == pin[key] for key in ("url", "sha256", "archive_sha256", "bytes"))
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


def test_expression_bounds_use_current_reproductive_observations(monkeypatch):
    import rebuild_target_catalog as rebuild

    selected = pd.DataFrame(
        {"entity_id": ["embryo-lineage", "legacy-only"], "developmental_mean_tpm": [0, 100]}
    )
    evidence = pd.DataFrame(
        {
            "entity_id": selected.entity_id,
            "normal_reference_tpm": [1, 1],
            "tumor_max_tpm": [0.5, 0.5],
            "reproductive_max_sample_tpm": [20, 0.5],
        }
    )
    normal = pd.DataFrame(
        {
            "entity_id": list(selected.entity_id) * 2,
            "scope": ["disposable", "disposable", "reproductive", "reproductive"],
            "value": [0, 0, 0.5, 0.5],
            "qc_pass": True,
            "identity_resolved": True,
        }
    )
    monkeypatch.setattr(
        rebuild, "disposable_ntpm", lambda _: pd.Series(0.0, index=selected.entity_id)
    )
    result = rebuild.expression_bounds(selected, evidence, normal).set_index("entity_id")
    assert result.loc["embryo-lineage", "expression_bounds_status"] == "pass"
    assert result.loc["embryo-lineage", "target_peak_tpm"] == 20
    assert result.loc["legacy-only", "expression_bounds_status"] == "fail"
    assert result.loc["legacy-only", "target_peak_tpm"] == 0.5


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


def test_rebuild_patient_metadata_matches_shared_api():
    metadata = pd.read_csv(OUT / "tumor-samples.csv.gz", keep_default_na=False)
    expected = sample_api.sample_identities(metadata[["source_cohort", "sample_id"]])
    pd.testing.assert_frame_equal(metadata[expected.columns], expected)
    assert metadata.participant_id.ne("").sum() == 13164
    assert not metadata.cancer_code.isin(policy.GERM_CELL_TUMORS).any()
    assert metadata.patient_identity_status.eq("conflict").sum() == 0


def test_lactating_mammary_separation_and_cta_tiers(catalog):
    meta = policy.metadata()
    lactation = meta["comparison"]["mammary_lactation"]
    assert "breast.410" in lactation["lactating_samples"]
    assert len(lactation["lactating_samples"]) == 6
    assert len(lactation["resting_samples"]) == 512
    assert {"breast.373", "breast.390", "breast.405"} <= set(lactation["resting_samples"])

    reproductive = pd.read_csv(
        OUT / "reproductive-comparison.csv.gz",
        usecols=["tissue", "compartment", "sample_id"],
        low_memory=False,
    )
    assert set(reproductive.compartment) == set(policy.REPRODUCTIVE_COMPARTMENTS)
    breast = reproductive[reproductive.tissue.eq("breast")]
    assert set(breast.loc[breast.compartment.eq("lactating_breast"), "sample_id"]) == set(
        lactation["lactating_samples"]
    )
    assert set(breast.loc[breast.compartment.eq("resting_breast"), "sample_id"]) == set(
        lactation["resting_samples"]
    )

    ctas = catalog[catalog.is_cta]
    non_ctas = catalog[~catalog.is_cta]
    assert non_ctas.cta_tier.fillna("").eq("").all()
    assert ctas.cta_tier.value_counts().to_dict() == {
        "core": 206,
        "extended_only": 21,
        "core_with_extended_exposure": 100,
    }
    by_name = ctas.set_index("protein_group_label")
    assert by_name.loc["LALBA", "cta_tier"] == "extended_only"
    assert set(by_name.loc["LALBA", "reproductive_qualifying_compartments"].split(" | ")) == set(
        policy.MAMMARY_COMPARTMENTS
    )
    assert by_name.loc["BTN1A1", "reproductive_qualifying_compartments"] == "lactating_breast"
    assert set(by_name.loc["TGM4", "reproductive_qualifying_compartments"].split(" | ")) == {
        "prostate",
        "resting_breast",
    }
    assert by_name.loc["MAGEA4", "cta_tier"] == "core"
    assert by_name.loc["MAGEA11", "cta_tier"] == "core_with_extended_exposure"
