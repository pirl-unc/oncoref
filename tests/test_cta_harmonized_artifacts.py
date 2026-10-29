"""Published research counts, counting units and provenance must remain coherent."""

import hashlib
import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
AUDIT = ROOT / "docs/audits/cta-harmonized-20261027"


def test_every_reference_locus_has_one_auditable_group_decision():
    genes = pd.read_csv(AUDIT / "gene-to-protein-group-decisions.csv.gz", low_memory=False)
    groups = pd.read_csv(AUDIT / "protein-group-classifications.csv.gz", low_memory=False)
    assert genes.gene_id.is_unique and len(genes) == 20097
    assert groups.protein_equivalence_group_id.is_unique
    assert genes.protein_equivalence_group_id.isin(groups.protein_equivalence_group_id).all()
    assert genes.identity_catalog_sha256.notna().all()
    for field in [
        "historical_ensembl_gene_ids",
        "alternate_ensembl_gene_ids",
        "gene_name_synonyms",
    ]:
        assert genes[field].map(json.loads).map(type).eq(list).all()
    # Group exports cannot misrepresent a representative gene's old flags as
    # group-level measurements or an unchanged gene identity.
    assert not {"gene_id", "selected", "tumor_maximum_max_ptpm"} & set(groups)
    summary = json.loads((AUDIT / "summary.json").read_text())
    assert groups.category.value_counts().to_dict() == summary["protein_group_categories"]
    assert genes.protein_group_category.eq("CTA").sum() == summary["selected_member_loci"]


def test_selected_groups_pass_all_layers_and_all_historical_ctas_are_retained():
    groups = pd.read_csv(AUDIT / "protein-group-classifications.csv.gz", low_memory=False)
    holds = [c for c in groups if c.startswith("hold:")]
    strong = [c for c in groups if c.startswith("unsafe:")]
    assert not groups.loc[groups.category.eq("CTA"), holds].any(axis=None)
    assert groups.loc[groups.category.eq("Unsafe"), strong].any(axis=1).all()
    historical = pd.read_csv(AUDIT / "historical-cta-decisions.csv")
    genes = pd.read_csv(AUDIT / "gene-classifications.csv.gz", low_memory=False)
    assert set(historical.gene_id) == set(
        genes.loc[genes["Expression category"].eq("CTA"), "gene_id"]
    )
    named = pd.read_csv(AUDIT / "named-antigen-decisions.csv")
    assert {"NY-ESO-1", "XAGE1A/B", "PRAME", "CTAG2", "CT83", "MAGEA4", "SSX2/B", "TBXT"} <= set(
        named.protein_group_label
    )


def test_active_report_receipts_pin_implementation_and_artifacts(historical_cta_input):
    for path in [AUDIT / "manifest.json", AUDIT / "group-normal/manifest.json"]:
        manifest = json.loads(path.read_text())
        for role in ["inputs", "outputs"]:
            for name, expected in manifest[role].items():
                assert not Path(name).is_absolute(), "Published receipts must be portable"
                source = (
                    historical_cta_input(name, expected) if role == "inputs" else path.parent / name
                )
                assert hashlib.sha256(source.read_bytes()).hexdigest() == expected, str(source)
