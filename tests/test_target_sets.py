"""Separate weak/variable evidence from the forbidden protein background."""

import hashlib
import json
from pathlib import Path

import pandas as pd
import pytest

from oncoref import target_sets
from oncoref.reference_identity import ReferenceIdentity

ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "docs/audits/target-sets-20261023"


def inputs(r=4.9, n=0.2, status="Retained", ihc=None):
    categories = pd.DataFrame(
        [
            {
                "Ensembl_Gene_ID": "g",
                "Symbol": "G",
                "category_policy_id": "oncoref.normal-expression-categories.2026-10-21.2",
                "reproductive_rna_max_ntpm": r,
                "rna_exclusions": "",
                "ihc_exclusions": "",
                "selectivity_exclusions": "",
                "other_subgroup": "",
                "quality_flags": "",
                "rna_routine_missing": 0,
            }
        ]
    )
    screen = pd.DataFrame(
        [
            {
                "Gene ID": "g",
                "Category": status,
                "Policy": "oncoref.cta-screen.2026-10-23.1",
                "Donor review": "Not assessed",
                "Donor findings": "",
                "Clinical evidence": "",
                "Reason": "Original screen",
            }
        ]
    )
    rows = []
    for tissue, tier, value in [("testis", "permitted", r), ("heart muscle", "vital", n)]:
        rows.append(
            {
                "Ensembl_Gene_ID": "g",
                "source_name": "hpa_rna_consensus",
                "source_version": "v23",
                "taxonomy_id": 9606,
                "modality": "bulk_rna",
                "tissue": tissue,
                "screening_tier": tier,
                "tissue_category": "reproductive" if tissue == "testis" else "vital",
                "rna_ntpm": value,
                "valid_measurement": True,
                "strict_positive": value > 0,
                "mapping_status": "mapped",
                "ihc_reliability": "",
                "ihc_level": "",
            }
        )
    for tissue, tier, level in ihc or []:
        rows.append(
            {
                "Ensembl_Gene_ID": "g",
                "source_name": "hpa_normal_tissue",
                "source_version": "v23",
                "taxonomy_id": 9606,
                "modality": "ihc",
                "tissue": tissue,
                "screening_tier": tier,
                "tissue_category": "reproductive" if tissue == "testis" else tier,
                "rna_ntpm": float("nan"),
                "valid_measurement": True,
                "strict_positive": level != "Not detected",
                "mapping_status": "mapped",
                "ihc_reliability": "Enhanced",
                "ihc_level": level,
            }
        )
    return categories, screen, pd.DataFrame(rows)


def decide(**kwargs):
    return target_sets.classify_genes(*inputs(**kwargs)).iloc[0]


def test_admission_floor_is_inclusive_and_does_not_make_low_genes_forbidden():
    low = decide(r=4.9)
    assert (low.Category, low.Detail) == ("Neither", "Low")
    assert not low.Target and not low.Forbidden
    assert decide(r=5).Category == "CTA"
    assert decide(r=0, n=0).Category == "Neither"


def test_negative_ihc_is_not_positive_support_and_not_missing():
    negative = decide(ihc=[("testis", "permitted", "Not detected")])
    assert negative.Category == "Neither" and negative.Detail == "Low"
    assert negative["IHC coverage"] == "Negative where measured"
    assert negative["IHC measurements"] == 1
    assert decide().get("IHC coverage") == "Missing"
    assert decide(r=1, ihc=[("testis", "permitted", "High")]).Category == "CTA"


def test_five_ntpm_is_not_a_blanket_normal_tissue_allowance():
    row = decide(r=4, n=3)
    assert row.Category == "Unsafe" and row.Detail == "RNA" and row.Forbidden


def test_positive_normal_ihc_overrides_low_rna():
    row = decide(r=1, n=0, ihc=[("heart muscle", "vital", "Low")])
    assert row.Category == "Unsafe" and row.Detail == "Protein"


def test_selectivity_only_failure_is_not_forbidden():
    categories, screen, observations = inputs(r=10, status="Excluded")
    categories["selectivity_exclusions"] = "Fails ratio"
    row = target_sets.classify_genes(categories, screen, observations).iloc[0]
    assert row.Category == "Neither" and row.Detail == "Selectivity"
    assert not row.Forbidden


def test_variable_is_neither_but_protein_and_safety_vetoes_win():
    values = inputs(
        r=10.1,
        n=5.1,
        status="Variable",
        ihc=[("testis", "permitted", "High"), ("heart muscle", "vital", "Not detected")],
    )
    _, screen, _ = values
    screen["Donor review"] = "Reviewed"
    screen["Donor findings"] = "0 and failing donor"
    row = target_sets.classify_genes(*values).iloc[0]
    assert row.Category == "Neither" and row.Detail == "Variable"
    assert not row.Target and not row.Forbidden
    assert row["Donor findings"] == "0 and failing donor"
    held = target_sets.classify_genes(*values, safety_holds={"g"}).iloc[0]
    assert held.Category == "Unsafe" and held.Detail == "Safety"
    positive = decide(r=10.1, n=5.1, status="Variable", ihc=[("heart muscle", "vital", "Low")])
    assert positive.Category == "Unsafe"


def test_missing_coverage_cannot_become_low_or_a_target():
    values = inputs(r=50)
    values[0]["rna_routine_missing"] = 1
    row = target_sets.classify_genes(*values).iloc[0]
    assert row.Category == "Neither" and row.Detail == "Unassessed"
    assert not row.Forbidden


def test_above_ceiling_donor_is_unsafe_without_qualifying_variable_pattern():
    values = inputs(r=10, n=0.2, status="Excluded")
    donors = pd.DataFrame(
        {"Gene ID": ["g"], "Highest RNA (nTPM)": [13.5], "Absolute ceiling (nTPM)": [12.0]}
    )
    row = target_sets.classify_genes(*values, donor_results=donors).iloc[0]
    assert row.Category == "Unsafe" and row["Normal donor RNA above ceiling"]


def test_wrong_species_or_policy_is_rejected():
    values = inputs()
    values[2]["taxonomy_id"] = 10090
    with pytest.raises(ValueError, match="another species"):
        target_sets.classify_genes(*values)
    values = inputs()
    values[0]["category_policy_id"] = "other"
    with pytest.raises(ValueError, match="pinned"):
        target_sets.classify_genes(*values)


@pytest.fixture
def reference():
    return ReferenceIdentity(
        taxon_id=9606,
        assembly_accession="fixture.1",
        annotation_source="fixture",
        annotation_release="1",
        source_version="1",
        identifier_namespaces={"gene": "fixture", "transcript": "fixture", "protein": "fixture"},
        content_provenance={"proteins": "a" * 64},
    )


def occurrence(reference, gene):
    return {
        "reference_id": reference.reference_id,
        "occurrence_id": gene,
        "gene_id": gene,
        "transcript_id": gene + "t",
        "protein_id": gene + "p",
        "source": "fixture",
        "sequence": "MACDEF",
        "complete": True,
    }


@pytest.mark.parametrize(
    "other,expected,forbidden",
    [
        ("unsafe", "Unsafe", True),
        ("neither", "Neither", False),
        ("unknown", "Neither", False),
        ("cta2", "CTA", False),
    ],
)
def test_shared_full_protein_assignment(reference, other, expected, forbidden):
    genes = pd.DataFrame(
        {
            "Gene ID": ["cta", "cta2", "unsafe", "neither"],
            "Category": ["CTA", "CTA", "Unsafe", "Neither"],
        }
    )
    proteins = target_sets.partition_proteins(
        [occurrence(reference, "cta"), occurrence(reference, other)],
        genes,
        reference=reference,
        expression_reference_id=reference.reference_id,
    )
    assert len(proteins) == 1
    assert proteins[0]["group"] == expected
    assert proteins[0]["forbidden"] is forbidden
    assert len(proteins[0]["occurrences"]) == 2
    inventory = pd.DataFrame(
        {
            "Protein sequence": [proteins[0]["protein_sequence_id"]],
            "Genes": [f"cta;{other}"],
            "Group": ["Reproductive"],
            "Reference": [reference.reference_id],
            "Exclude by default": [False],
        }
    )
    classified = target_sets.repartition_protein_table(
        inventory, genes, reference_id=reference.reference_id
    )
    assert classified.iloc[0].Group == expected
    assert bool(classified.iloc[0].Forbidden) is forbidden


def test_published_gene_groups_and_important_targets():
    genes = target_sets.results()
    assert len(genes) == 20097 and not genes["Gene ID"].duplicated().any()
    groups = [target_sets.gene_ids(group) for group in target_sets.GROUPS]
    assert len(set.union(*groups)) == len(genes)
    assert not any(groups[i] & groups[j] for i in range(3) for j in range(i + 1, 3))
    assert len(groups[0]) == 399
    indexed = genes.set_index("Gene")
    for name, detail in [("CTAG2", "Variable"), ("ZNF730", "Cellular"), ("FOXR2", "Antibody")]:
        assert indexed.loc[name, "Category"] == "Neither"
        assert indexed.loc[name, "Detail"] == detail
        assert not indexed.loc[name, "Forbidden"]
    assert indexed.loc["MAGEA12", "Category"] == "Unsafe"
    for name in ["IL5", "IL13", "ABCB5", "PAPPA"]:
        assert indexed.loc[name, "Category"] == "Unsafe"
    assert indexed.loc["CT55", "Category"] == "CTA"
    for name in ["HORMAD1", "TDRD1", "BRDT", "DAZL"]:
        assert indexed.loc[name, "Category"] == "Neither"
    for name in ["ACRBP", "TEX101", "GAGE12F", "PAGE2B", "TRIP13", "MAGEA3", "ROPN1"]:
        assert indexed.loc[name, "Category"] == "Unsafe"
        assert indexed.loc[name, "Detail"] == "High RNA"
    benchmark = {
        "AKAP4",
        "CT83",
        "CTAG1A",
        "CTAG1B",
        "CTCFL",
        "GAGE1",
        "MAGEA1",
        "MAGEA10",
        "MAGEA4",
        "MAGEC1",
        "MAGEC2",
        "PRAME",
        "SAGE1",
        "SPANXB1",
        "SSX2",
        "SSX4",
        "SSX5",
        "XAGE1A",
        "XAGE1B",
    }
    assert benchmark <= target_sets.gene_names("CTA")
    removed = pd.read_csv(REPORT / "removed-cta-candidates.csv")
    assert len(removed) == 151 and removed.Category.eq("Neither").all()
    assert removed["Positive IHC measurements"].eq(0).all()
    assert removed["Maximum RNA (nTPM)"].lt(5).all()


def test_current_report_and_dataset_checksums():
    manifest = json.loads((REPORT / "run-manifest.json").read_text())
    assert manifest["forbidden_categories"] == ["Unsafe"]
    for section, root in [("inputs", ROOT), ("implementation_sha256", ROOT), ("outputs", REPORT)]:
        for name, expected in manifest[section].items():
            assert hashlib.sha256((root / name).read_bytes()).hexdigest() == expected, name
    historical = REPORT / "gene-results.csv.gz"
    assert (
        hashlib.sha256(historical.read_bytes()).hexdigest() == manifest["published_dataset_sha256"]
    )
    recipe = json.loads((REPORT / "exclusion-basis.json").read_text())
    assert recipe["forbidden"] == ["unsafe.fasta.gz", "additional-forbidden.fasta.gz"]
    assert recipe["neither_forbidden"] is False


def test_historical_normal_cell_report_checksums():
    report = ROOT / "docs/audits/cta-normal-cells-20261023"
    manifest = json.loads((report / "run-manifest.json").read_text())
    assert manifest["forbidden_categories"] == ["Unsafe"]
    for section, root in [("implementation_and_inputs", ROOT), ("outputs", report)]:
        for name, expected in manifest[section].items():
            assert hashlib.sha256((root / name).read_bytes()).hexdigest() == expected, name


def test_current_neither_review_and_dataset_checksums():
    report = ROOT / "docs/audits/cta-neither-review-20261024"
    manifest = json.loads((report / "run-manifest.json").read_text())
    assert manifest["forbidden_categories"] == ["Unsafe"]
    for section, root in [("implementation_and_inputs", ROOT), ("outputs", report)]:
        for name, expected in manifest[section].items():
            assert hashlib.sha256((root / name).read_bytes()).hexdigest() == expected, name
    assert (ROOT / "oncoref/data/target-gene-sets.csv.gz").read_bytes() == (
        report / "gene-results.csv.gz"
    ).read_bytes()
    assert (ROOT / "oncoref/data/cta-followup-reviews.csv").read_bytes() == (
        report / "follow-up.csv"
    ).read_bytes()
    assert (ROOT / "oncoref/data/sequence-gene-rna-reviews.csv.gz").read_bytes() == (
        report / "sequence-gene-reviews.csv.gz"
    ).read_bytes()


def test_no_remaining_neither_or_cta_has_high_nonpermitted_measured_rna():
    report = ROOT / "docs/audits/cta-neither-review-20261024"
    genes = target_sets.results()
    maxima = pd.read_csv(report / "rna-context-maxima.csv.gz", low_memory=False)
    high = maxima.loc[~maxima.permitted & maxima.maximum.gt(30)]
    from oncoref.normal_rna_review import forbidden_gene_ids, sequence_gene_reviews

    extra = sequence_gene_reviews()
    all_genes = pd.concat([genes, extra], ignore_index=True)
    remaining = all_genes.loc[all_genes.Category.isin(["CTA", "Neither"])]
    assert remaining["Normal RNA assessed"].all()
    assert not set(high.gene_id) & set(remaining["Gene ID"])
    observations = pd.read_csv(report / "high-normal-rna.csv.gz", low_memory=False)
    assert set(high.gene_id) == set(observations.gene_id)
    assert not observations.finding_id.duplicated().any()
    assert observations.unit.isin(["TPM", "nTPM"]).all()
    assert observations.value.gt(30).all()
    assert not observations.permitted.any()
    assert maxima.loc[maxima.category.eq("reproductive"), "maximum"].max() > 1000
    assert len(extra) == 1721 and not extra.Target.any()
    assert set(extra.loc[extra.Category.eq("Unsafe"), "Gene ID"]) <= forbidden_gene_ids()
    assert not set(extra.loc[extra.Category.eq("Neither"), "Gene ID"]) & forbidden_gene_ids()


def test_followup_records_preserve_questions_without_gene_exceptions():
    from oncoref.normal_rna_review import followup_reviews

    records = followup_reviews()
    assert len(records) == 15 and not records.Gene.duplicated().any()
    for col in ["Review question", "Reconsideration evidence", "Limitations", "Evidence report"]:
        assert records[col].notna().all() and records[col].ne("").all()
    indexed = records.set_index("Gene")
    assert indexed.loc["CTAG2", "Detail"] == "Variable"
    assert indexed.loc["FOXR2", "Detail"] == "Antibody"
    assert indexed.loc["ACRBP", "Forbidden"]
    assert indexed.loc["MAGEA3", "Forbidden"]
    assert indexed.loc["ROPN1", "Forbidden"]
