#!/usr/bin/env python3
"""Explain current coverage gaps against pinned source rows and identity audits.

Matching uses the shared human registry, including registered historical IDs and
unique synonyms. This audit does not impute measurements or merge conflicting
loci. Presence in a count matrix does not supply a missing TPM measurement.
"""

from __future__ import annotations

import json
import sys
import zipfile
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from cta_report_common import input_hashes, sha256, write_csv  # noqa: E402

from oncoref import cell_rna_policy, hpa_cancer, normal_cell_screen, reference_data  # noqa: E402
from oncoref.gene_matching import gene_identity_catalog, match_reference_genes  # noqa: E402

OUT = ROOT / "docs/audits/cta-source-walkthrough-20261028"
CURRENT = ROOT / "docs/audits/cta-selection-20261028"
EXPANDED = ROOT / "docs/audits/cta-expanded-funnel-20261027"
PANELS = [
    ("hpa_rna_consensus", "HPA v23 consensus bulk", ["Tissue"], "nTPM", 50),
    (
        "hpa_single_cell_tissue",
        "HPA v23 tissue-cell clusters",
        ["Tissue", "Cell type", "Cluster"],
        "nTPM",
        426,
    ),
    ("hpa_immune_cell", "HPA v23 immune cells", ["Immune cell"], "nTPM", 19),
    ("hpa_immune_monaco", "Monaco immune cells via HPA v23", ["Immune cell"], "TPM", 30),
]


def encoded(values) -> str:
    """Serialize distinct identifiers or contexts with deterministic ordering."""
    return json.dumps(sorted(set(values)), separators=(",", ":"))


def build() -> None:
    """Export source-specific group gaps, member identities and cohort coverage."""
    facts = pd.read_csv(CURRENT / "protein-group-classifications.csv.gz", low_memory=False)
    target = facts[facts.category.isna() | ~facts.required_coverage].copy()
    target["member_gene_ids"] = target.member_gene_ids.map(json.loads)
    target_ids = {g for ids in target.member_gene_ids for g in ids}
    catalog = gene_identity_catalog().set_index("canonical_gene_id")
    catalog_hash = catalog.identity_catalog_sha256.iloc[0]
    catalog = catalog.reindex(catalog.index.union(sorted(target_ids)))
    for column in [
        "historical_ensembl_gene_ids",
        "alternate_ensembl_gene_ids",
        "gene_name_synonyms",
    ]:
        catalog[column] = catalog[column].fillna("[]")
    inputs = [
        Path(__file__),
        CURRENT / "protein-group-classifications.csv.gz",
        CURRENT / "hpa-protein-evidence.csv.gz",
        EXPANDED / "tumor-union/manifest.json",
    ]
    source_receipts, matches, normal_rows, members = [], [], [], []
    for source, title, context_cols, unit, expected in PANELS:
        receipt = reference_data.provenance(source, "v23", verify_content=True)
        source_receipts.append({k: receipt[k] for k in ["name", "version", "url", "sha256"]})
        d = pd.read_csv(
            receipt["path"], sep="\t", usecols=["Gene", "Gene name", *context_cols, unit]
        )
        if source == "hpa_single_cell_tissue":
            contexts = (
                d[["Tissue", "Cell type"]]
                .drop_duplicates()
                .rename(columns={"Tissue": "tissue", "Cell type": "cell_type"})
            )
            contexts = normal_cell_screen.contexts(
                contexts.assign(source=source, version="v23", taxonomy_id=9606, stage="adult")
            )
            accepted = contexts[~contexts.category.isin(cell_rna_policy.PERMITTED)][
                ["tissue", "cell_type"]
            ]
            d = d.merge(
                accepted.rename(columns={"tissue": "Tissue", "cell_type": "Cell type"}),
                on=["Tissue", "Cell type"],
                validate="many_to_one",
            )
        pairs = d[["Gene", "Gene name"]].drop_duplicates()
        identity = match_reference_genes(
            pairs, id_column="Gene", symbol_column="Gene name", source=source, source_version="v23"
        )
        relevant = identity.source_gene_id.isin(target_ids) | identity.candidate_gene_ids.map(
            lambda x: bool(set(json.loads(x)) & target_ids)
        )
        matches.append(identity[relevant])
        identity_candidates = {}
        for row in identity[relevant].itertuples():
            for gid in (
                set(json.loads(row.candidate_gene_ids)) | {row.source_gene_id}
            ) & target_ids:
                identity_candidates.setdefault(gid, []).append(row)
        accepted_path = CURRENT / f"{source}-accepted-identities.csv.gz"
        inputs.append(accepted_path)
        good = pd.read_csv(accepted_path)
        lookup = good.set_index("source_gene_id").canonical_gene_id
        if not lookup.index.is_unique:
            raise ValueError(f"Conflicting source names for one ID: {source}")
        d["canonical_gene_id"] = d.Gene.map(lookup)
        expected_contexts = set(
            d[context_cols].drop_duplicates().astype(str).agg(" / ".join, axis=1)
        )
        if len(expected_contexts) != expected:
            raise ValueError(f"Unexpected context universe in {source}")
        relevant = d[d.canonical_gene_id.isin(target_ids)].copy()
        relevant["context"] = relevant[context_cols[0]].astype(str)
        for column in context_cols[1:]:
            relevant["context"] = relevant.context + " / " + relevant[column].astype(str)
        seen = relevant[relevant[unit].notna()].groupby("canonical_gene_id").context.agg(set)
        source_ids = relevant.groupby("canonical_gene_id").Gene.agg(set)
        source_names = relevant.groupby("canonical_gene_id")["Gene name"].agg(
            lambda x: set(x.dropna())
        )
        for gid in sorted(target_ids):
            contexts = seen.get(gid, set())
            missing = expected_contexts - contexts
            members.append(
                {
                    "source_id": source,
                    "source": title,
                    "gene_id": gid,
                    "symbol": catalog.loc[gid, "canonical_symbol"],
                    "historical_gene_ids": catalog.loc[gid, "historical_ensembl_gene_ids"],
                    "alternate_gene_ids": catalog.loc[gid, "alternate_ensembl_gene_ids"],
                    "synonyms": catalog.loc[gid, "gene_name_synonyms"],
                    "source_gene_ids": encoded(source_ids.get(gid, set())),
                    "source_gene_names": encoded(source_names.get(gid, set())),
                    "candidate_source_ids": encoded(
                        r.source_gene_id for r in identity_candidates.get(gid, [])
                    ),
                    "identity_statuses": encoded(
                        r.identity_status for r in identity_candidates.get(gid, [])
                    ),
                    "observed_contexts": len(contexts),
                    "expected_contexts": expected,
                    "missing_contexts": encoded(missing),
                    "coverage_status": "complete"
                    if not missing
                    else "partial"
                    if contexts
                    else "absent",
                    "source_url": receipt["url"],
                }
            )
        for r in target.itertuples():
            absent = [g for g in r.member_gene_ids if not seen.get(g, set())]
            incomplete = [g for g in r.member_gene_ids if seen.get(g, set()) != expected_contexts]
            if not incomplete:
                continue
            missing = set().union(
                *(expected_contexts - seen.get(g, set()) for g in r.member_gene_ids)
            )
            normal_rows.append(
                {
                    "entity_id": r.protein_equivalence_group_id,
                    "label": r.protein_group_label,
                    "category": r.category,
                    "assessment_status": r.assessment_status,
                    "source_id": source,
                    "source": title,
                    "member_gene_ids": encoded(r.member_gene_ids),
                    "incomplete_member_gene_ids": encoded(incomplete),
                    "absent_member_gene_ids": encoded(absent),
                    "missing_contexts": encoded(missing),
                    "expected_contexts_per_member": expected,
                    "minimum_observed_contexts": min(
                        len(seen.get(g, set())) for g in r.member_gene_ids
                    ),
                    "source_url": receipt["url"],
                }
            )
        print(f"Audited {title}", flush=True)
        del d, relevant

    atlas = ROOT / "tmp/cta-quantitative/atlas_tpm"
    with zipfile.ZipFile(atlas) as archive:
        name = next(n for n in archive.namelist() if n.endswith("_rows"))
        rows = pd.read_csv(archive.open(name), sep="\t", header=None, names=["gene_id", "label"])
    identity = match_reference_genes(
        rows,
        id_column="gene_id",
        source="tyser_E-MTAB-9388",
        source_version="Ensembl106_filtered_TPM",
    )
    matches.append(
        identity[identity.candidate_gene_ids.map(lambda x: bool(set(json.loads(x)) & target_ids))]
    )
    lineage_rows = pd.read_csv(
        CURRENT / "developmental-lineages.csv.gz", usecols=["Gene ID", "source"]
    )
    measured = set(lineage_rows["Gene ID"])
    developmental_sources = {
        source: set(part["Gene ID"]) for source, part in lineage_rows.groupby("source")
    }
    developmental_sources["combined_developmental_TPM"] = measured
    developmental = []
    for r in target.itertuples():
        for source, available in developmental_sources.items():
            missing = set(r.member_gene_ids) - available
            developmental.append(
                {
                    "entity_id": r.protein_equivalence_group_id,
                    "label": r.protein_group_label,
                    "assessment_status": r.assessment_status,
                    "missing_member_gene_ids": encoded(missing),
                    "source_id": source,
                    "complete": not missing,
                    "source_version": "Ensembl106_filtered_TPM",
                }
            )
    source_receipts.append(
        {
            "name": "tyser_E-MTAB-9388",
            "sha256": sha256(atlas),
            "url": "https://www.ebi.ac.uk/gxa/sc/experiments/E-MTAB-9388/Results",
        }
    )

    unassigned = target[target.category.isna()]
    tumor_missing = unassigned[
        unassigned[["oncoref_maximum_tpm", "hpa_tumor_group_lower_bound_ptpm"]].isna().all(axis=1)
    ]
    tumor_missing_ids = {g for ids in tumor_missing.member_gene_ids for g in ids}
    cohorts = []
    tumor_identities = []
    for path in sorted((EXPANDED / "tumor-union/cohorts").glob("*.csv.gz")):
        repaired_path = CURRENT / "tumor-identity-repair/cohorts" / path.name
        if repaired_path.exists():
            path = repaired_path
        d = pd.read_csv(path, usecols=["gene_id", "maximum_tpm", "source_cohort"])
        observed = set(d.loc[d.maximum_tpm.notna(), "gene_id"])
        cohort = str(d.source_cohort.iloc[0]) if len(d) else ""
        identity_path = path.parent.parent / "identity" / path.name
        inputs.extend([path, identity_path])
        match = pd.read_csv(identity_path, low_memory=False).fillna("")
        relevant = match.source_gene_id.isin(tumor_missing_ids) | match.candidate_gene_ids.map(
            lambda x: bool(set(json.loads(x)) & tumor_missing_ids)
        )
        tumor_identities.append(match[relevant].assign(cohort=path.stem.removesuffix(".csv")))
        for r in unassigned.itertuples():
            missing = set(r.member_gene_ids) - observed
            cohorts.append(
                {
                    "entity_id": r.protein_equivalence_group_id,
                    "label": r.protein_group_label,
                    "provider": "OncoRef",
                    "cohort": path.stem.removesuffix(".csv"),
                    "source_cohort": cohort,
                    "unit": "TPM",
                    "member_count": r.member_count,
                    "measured_members": r.member_count - len(missing),
                    "missing_member_gene_ids": encoded(missing),
                }
            )
    cancer_receipt = reference_data.provenance("hpa_cancer_rna", verify_content=True)
    source_receipts.append({k: cancer_receipt[k] for k in ["name", "version", "url", "sha256"]})
    for (code, cohort), d in hpa_cancer.hpa_cancer_rna_prevalence().groupby(
        ["cancer_code", "cohort"]
    ):
        observed = set(d.loc[d.samples.gt(0) & d.max_ptpm.notna(), "gene_id"])
        for r in unassigned.itertuples():
            missing = set(r.member_gene_ids) - observed
            cohorts.append(
                {
                    "entity_id": r.protein_equivalence_group_id,
                    "label": r.protein_group_label,
                    "provider": "HPA 25.1",
                    "cohort": code,
                    "source_cohort": cohort,
                    "unit": "pTPM",
                    "member_count": r.member_count,
                    "measured_members": r.member_count - len(missing),
                    "missing_member_gene_ids": encoded(missing),
                }
            )
    outputs = {
        "normal-panel-coverage-gaps.csv.gz": pd.DataFrame(normal_rows),
        "normal-panel-member-coverage.csv.gz": pd.DataFrame(members),
        "coverage-identity-matches.csv.gz": pd.concat(matches, ignore_index=True),
        "embryonic-coverage-gaps.csv.gz": pd.DataFrame(developmental),
        "combined-rna-coverage.csv.gz": pd.read_csv(CURRENT / "combined-rna-coverage.csv.gz"),
        "tumor-cohort-coverage.csv.gz": pd.DataFrame(cohorts),
        "tumor-missing-identity-candidates.csv.gz": pd.concat(tumor_identities, ignore_index=True),
    }
    proteins = pd.read_csv(CURRENT / "hpa-protein-evidence.csv.gz", low_memory=False)
    outputs["unresolved-protein-contexts.csv.gz"] = proteins[
        proteins.unresolved_evidence
        & proteins.representative_gene_id.isin(unassigned.representative_gene_id)
    ]
    for name, frame in outputs.items():
        write_csv(frame, OUT / name)
    summary = []
    normal = outputs["normal-panel-coverage-gaps.csv.gz"]
    incomplete = normal[normal.assessment_status.eq("incomplete_evidence")]
    for source, title, _, _, _ in PANELS:
        s = incomplete[incomplete.source_id.eq(source)]
        summary.append(
            {
                "source_id": source,
                "source": title,
                "groups": s.entity_id.nunique(),
                "example_genes": "; ".join(s.label.dropna().head(5)),
            }
        )
    combined = outputs["combined-rna-coverage.csv.gz"]
    combined = combined[
        combined.entity_id.isin(
            target.loc[
                target.assessment_status.eq("incomplete_evidence"), "protein_equivalence_group_id"
            ]
        )
    ]
    missing = combined[~combined.bulk_complete]
    summary.insert(
        1,
        {
            "source_id": "combined_normal_bulk",
            "source": "Combined normal bulk: HPA consensus + HPA brain + direct GTEx",
            "groups": len(missing),
            "example_genes": "; ".join(missing.label.fillna(missing.entity_id).head(5)),
        },
    )
    gtex_missing = combined[~combined.gtex_core_complete]
    summary.append(
        {
            "source_id": "gtex_core_panel",
            "source": "Required GTEx v11 bulk panel (52 contexts)",
            "groups": len(gtex_missing),
            "example_genes": "; ".join(gtex_missing.label.head(5)),
        }
    )
    pd.DataFrame(summary).to_csv(OUT / "normal-panel-gap-counts.csv", index=False)
    (OUT / "coverage-audit-receipt.json").write_text(
        json.dumps(
            {
                "inputs": input_hashes(inputs),
                "sources": source_receipts,
                "identity_catalog_sha256": catalog_hash,
                "outputs": {
                    name: sha256(OUT / name) for name in [*outputs, "normal-panel-gap-counts.csv"]
                },
            },
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    build()
