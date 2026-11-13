"""Cited CTA literature reviews and distinct source roles.

This registry annotates the frozen 1.31.0 selection. A paper's placement in the
core or supplementary reading set never determines whether its findings apply.
This frozen registry does not itself edit classifications. Current transcript
decisions and reviewed coding holds are exposed by ``transcript_selection``.
"""

from __future__ import annotations

import gzip
import json
from collections import defaultdict
from copy import deepcopy
from importlib.resources import files

import pandas as pd

NOMINATION_SOURCE_TAGS = frozenset(
    {
        "Wang2016_CT",
        "Bruggeman2018_GC",
        "daSilva2017_CT",
        "Jamin2021_CT",
        "Jamin2021_core",
        "Chang2019_TGCT",
        "Carter2023_CT",
        "Seager2024_CTA",
        "Loriot2025_S1",
        "Loriot2025_S2",
        "Bai2016_EGFL6",
        "Boel1995_BAGE",
        "Bradley2020_CPA",
        "Garg2007_SPAG9",
        "Sasaroli2011_vascular",
    }
)


def registry() -> dict:
    """Return an independent copy of the packaged, versioned evidence registry."""
    payload = files("oncoref").joinpath("data/cta-literature.json.gz").read_bytes()
    return json.loads(gzip.decompress(payload))


def sources(*, collection: str | None = None, role: str | None = None) -> pd.DataFrame:
    """Source citations, coverage, overlap and retention reasons.

    `collection` is core or supplementary. `role` is one of registry()['roles'].
    All roles remain available regardless of whether a source adds unique names.
    """
    data = registry()
    if collection not in {None, "core", "supplementary"}:
        raise ValueError("collection must be core, supplementary or None")
    if role is not None and role not in data["roles"]:
        raise ValueError(f"Unknown literature role: {role}")
    result = pd.DataFrame(data["sources"])
    if collection is not None:
        result = result[result.collection.eq(collection)]
    if role is not None:
        result = result[result.roles.map(lambda roles: role in roles)]
    return result.reset_index(drop=True)


def reviews(gene: str | None = None) -> pd.DataFrame:
    """Individual reviews for the 41 unmatched and 55 broad-table-only proteins.

    An optional exact protein label or Ensembl gene ID selects its review. Source
    scope and unresolved assay identity are retained alongside every citation.
    """
    result = pd.DataFrame(registry()["reviews"])
    if gene is not None:
        result = result[result.protein.eq(gene) | result.gene_ids.map(lambda ids: gene in ids)]
    return result.reset_index(drop=True)


def memberships() -> pd.DataFrame:
    """All paper-to-protein links, including unresolved and noncoverage links."""
    return pd.DataFrame(registry()["memberships"])


def attestations(gene: str | None = None) -> pd.DataFrame:
    """Protein-paper occurrences with their citations and original claim scope.

    Includes current nominations, the frozen literature collection and subsequent
    corrections, for every category. Nothing is dropped because a paper is
    supplementary, redundant, adverse, or cannot resolve a related copy.
    ``published_nomination`` is a curated catalog membership, not arbitrary paper
    coverage. Source-level roles are explicitly distinguished from the claim on
    the protein-paper link. Original records retain assay/identity limitations.
    """
    from . import cta_selection
    from .cta_provenance import paper_membership

    literature = registry()
    followup = followup_findings()
    audit = cta_selection.audit()
    catalog_sources = {r["source_id"]: r for r in audit.sources.to_dict("records")}
    paper_sources = {r["source_id"]: r for r in literature["sources"]}
    followup_sources = {r["source_id"]: r for r in followup["sources"]}
    # A DOI can contain a reproductive background, a tumor-positive catalog,
    # and proteomics. Preserve and evaluate the exact table, not the DOI's role.
    claims = defaultdict(list)
    publication_rows = paper_membership().fillna("")
    mapped_papers = set(publication_rows.doi)
    claim_columns = [
        "source_tag",
        "Ensembl_Gene_ID",
        "source_file",
        "source_sheet",
        "source_rows",
        "source_subset",
        "validation_scope",
        "assay_resolution",
    ]
    for record in publication_rows.to_dict("records"):
        if record["Ensembl_Gene_ID"]:
            claims[(record["doi"], record["Ensembl_Gene_ID"])].append(
                {key: record[key] for key in claim_columns}
            )
    collections = [
        ("current_nomination", audit.memberships.to_dict("records"), catalog_sources),
        ("literature_registry", literature["memberships"], paper_sources),
        ("literature_followup", followup["memberships"], followup_sources),
    ]
    rows = []
    for origin, links, source_records in collections:
        for link in links:
            genes = list(link["gene_ids"])
            protein = link.get("protein", link.get("label", ""))
            if gene is not None and gene != protein and gene not in genes:
                continue
            sid = link["source_id"]
            source = source_records[sid]
            paper = paper_sources.get(sid, {})
            nomination = (
                (origin == "current_nomination" or link.get("kind") == "historical_catalog")
                and link.get("eligible_for_coverage", True)
                and "catalog_nomination" in paper.get("roles", [])
            )
            exact_claims = [claim for g in genes for claim in claims[(sid, g)]]
            nomination_genes = (
                (
                    sorted(
                        {
                            claim["Ensembl_Gene_ID"]
                            for claim in exact_claims
                            if claim["source_tag"] in NOMINATION_SOURCE_TAGS
                        }
                    )
                    if sid in mapped_papers
                    else genes
                )
                if nomination
                else []
            )
            rows.append(
                {
                    "origin": origin,
                    "source_id": sid,
                    "source_url": source.get("url", source.get("source_url", "")),
                    "citation": source.get("citation", ""),
                    "title": source.get("title", paper.get("title", "")),
                    "source_roles": list(source.get("roles", paper.get("roles", []))),
                    "collection": paper.get("collection", ""),
                    "protein": protein,
                    "gene_ids": genes,
                    "entity_id": link.get("entity_id", ""),
                    "published_nomination": bool(nomination_genes),
                    "nomination_gene_ids": nomination_genes,
                    "publication_claims": deepcopy(exact_claims),
                    "eligible_for_coverage": bool(
                        link.get("eligible_for_coverage", origin == "current_nomination")
                    ),
                    "claim_scope": link.get("scope", source.get("evidence_scope", "")),
                    "source_scope": source.get("scope", source.get("evidence_scope", "")),
                    "membership_record": deepcopy(link),
                    "source_record": deepcopy(source),
                }
            )
    columns = [
        "origin",
        "source_id",
        "source_url",
        "citation",
        "title",
        "source_roles",
        "collection",
        "protein",
        "gene_ids",
        "entity_id",
        "published_nomination",
        "nomination_gene_ids",
        "publication_claims",
        "eligible_for_coverage",
        "claim_scope",
        "source_scope",
        "membership_record",
        "source_record",
    ]
    return pd.DataFrame(rows, columns=columns)


def cta_coverage() -> pd.DataFrame:
    """Paper coverage for all 602 released CTAs; empty lists are coverage gaps."""
    return pd.DataFrame(registry()["cta_coverage"])


def followup_findings() -> dict:
    """Dated literature corrections; recommendations do not change released labels."""
    return json.loads(files("oncoref").joinpath("data/cta-literature-followup.json").read_text())


def _reconcile_coverage(profiles, nominations, links, additions):
    """Join exact loci, never transfer an old protein label to its other copies."""
    nominated, covered, related = (defaultdict(set) for _ in range(3))
    for row in nominations:
        for gene in row["gene_ids"]:
            nominated[gene].add(row["source_id"])
            covered[gene].add(row["source_id"])
    for row in [*links, *additions]:
        target = covered if row["eligible_for_coverage"] else related
        for gene in row["gene_ids"]:
            target[gene].add(row["source_id"])
    result = []
    for row in profiles:
        if row["category"] != "CTA":
            continue
        genes = sorted(row["gene_ids"])
        matched = [g for g in genes if covered[g]]
        result.append(
            {
                "entity_id": row["entity_id"],
                "protein": row["label"],
                "gene_ids": genes,
                "nomination_sources": sorted(set().union(*(nominated[g] for g in genes))),
                "exact_paper_sources": sorted(set().union(*(covered[g] for g in genes))),
                "related_sources": sorted(set().union(*(related[g] for g in genes))),
                "covered_gene_ids": matched,
                "uncovered_gene_ids": [g for g in genes if g not in matched],
                "paper_coverage": (
                    "all_members"
                    if len(matched) == len(genes)
                    else "some_members"
                    if matched
                    else "no_paper"
                ),
                "per_gene_sources": {g: sorted(covered[g]) for g in genes},
            }
        )
    return pd.DataFrame(result)


def current_coverage() -> pd.DataFrame:
    """Paper coverage of current CTAs, preserving exact member-locus scope.

    Includes current nomination memberships, the frozen literature registry and
    dated follow-up findings. A covered paper may concern normal function or
    counterevidence; coverage does not establish CTA nomination or safety.
    ``some_members`` means the shared protein has paper coverage for only some
    source loci. The original 602-protein ``cta_coverage()`` remains reproducible.
    """
    from . import cta_selection

    return _reconcile_coverage(
        cta_selection.profiles().to_dict("records"),
        cta_selection.audit().memberships.to_dict("records"),
        registry()["memberships"],
        followup_findings()["memberships"],
    )


def _protein_rows(key: str, gene: str | None) -> pd.DataFrame:
    result = pd.DataFrame(registry()[key])
    if gene is not None:
        result = result[result.protein.eq(gene) | result.gene_ids.map(lambda ids: gene in ids)]
    return result.reset_index(drop=True)


def gap_assessments(gene: str | None = None) -> pd.DataFrame:
    """All 33 paper gaps: biology, released measurements and review recommendations."""
    return _protein_rows("gap_assessments", gene)


def concerns(gene: str | None = None) -> pd.DataFrame:
    """The 29 flagged reviews, separated by evidence type and next action."""
    return _protein_rows("concerns", gene)


def localization(gene: str | None = None) -> pd.DataFrame:
    """Gene-mapped locations for all 602 CTAs, with shared and partial scope.

    Facets are annotation summaries, not validated tumor-surface epitopes or
    selected-sequence localization. Native evidence, isoform notes and primary
    paper additions are retained in each row.
    """
    return _protein_rows("localization", gene)


def surface_domains(gene: str | None = None) -> pd.DataFrame:
    """37 membrane CTAs: inclusive extracellular coordinates and sequence attribution.

    UniProt canonical dimensions and selected-sequence projections are separate.
    Missing boundaries are unknown, not zero. Predicted topology is not measured
    tumor-surface accessibility.
    """
    return _protein_rows("surface_domains", gene)


def surface_cancer_rna(gene: str | None = None, *, cohort: str | None = None) -> pd.DataFrame:
    """HPA 25.1 native-pTPM sample counts for the 37 membrane CTAs.

    TCGA and validation cohorts remain separate. Numerators use >= thresholds;
    denominators count measured samples, not the population prevalence of cancer.
    """
    if cohort not in {None, "TCGA", "validation"}:
        raise ValueError("cohort must be TCGA, validation or None")
    result = _protein_rows("surface_cancer_rna", gene)
    if cohort is not None:
        result = result[result.cohort.eq(cohort)]
    return result.reset_index(drop=True)


def surface_cancer_ihc(gene: str | None = None) -> pd.DataFrame:
    """HPA scored-patient staining counts; retain unmeasured records explicitly.

    These counts do not establish antibody specificity, subcellular localization,
    the fraction of malignant cells positive, or selected-proteoform expression.
    """
    return _protein_rows("surface_cancer_ihc", gene)
