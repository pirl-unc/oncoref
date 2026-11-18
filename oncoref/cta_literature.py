"""Cited CTA literature reviews and distinct source roles.

This registry annotates the frozen 1.31.0 selection. A paper's placement in the
core or supplementary reading set never determines whether its findings apply.
This frozen registry does not itself edit classifications. Current transcript
decisions and reviewed coding holds are exposed by ``transcript_selection``.
"""

from __future__ import annotations

import gzip
import json
import re
from collections import defaultdict
from copy import deepcopy
from importlib.resources import files

import pandas as pd

# Uses describe the claim's role in curation, never automatic set admission.
PAPER_USES = {
    "cta": "cta_nomination",
    "oncofetal": "oncofetal_support",
    "tumor": "tumor_support",
    "reproductive": "reproductive_support",
}
ALLOWED_USES = frozenset(PAPER_USES.values()) | {
    "translation_support",
    "translation_context",
    "functional_context",
    "normal_tissue_review",
    "identity_review",
    "immune_context",
}
EVIDENCE_TYPES = frozenset(
    {
        "rna_expression",
        "induced_rna_expression",
        "immunostaining",
        "ribosome_profiling",
        "regulatory_assay",
        "combined_gene_perturbation",
        "mass_spectrometry",
        "assay_identity_review",
        "cta_nomination",
        "cancer_embryo_program",
    }
)

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


def current_attestations(
    gene: str | None = None, *, catalog: pd.DataFrame | None = None
) -> pd.DataFrame:
    """Literature occurrences mapped onto the current catalog by exact member loci.

    Historical group IDs and follow-ups without a group ID cannot serve as join
    keys. Original records remain intact; ``gene_ids`` and nomination scope are
    restricted to the current group's intersecting loci. Related-copy and rejected
    assay records remain visible with ``eligible_for_coverage=False``.
    """
    from . import target_selection

    if catalog is None:
        catalog = target_selection.results(gene)
    elif gene is not None:
        catalog = catalog[
            catalog.protein_group_label.eq(gene)
            | catalog.member_gene_ids.map(lambda ids: gene in json.loads(ids))
        ]
    by_gene = {
        member: (row.entity_id, row.protein_group_label)
        for row in catalog.itertuples()
        for member in json.loads(row.member_gene_ids)
    }
    links = attestations()
    reviewed = {
        (source["source_id"], member): claim
        for source in publication_registry()["sources"]
        for claim in source["gene_claims"]
        for member in claim["gene_ids"]
    }
    rows = []
    for record in links.to_dict("records"):
        groups = defaultdict(list)
        for member in record["gene_ids"]:
            if member in by_gene:
                groups[by_gene[member]].append(member)
        for (entity, label), members in groups.items():
            gene_claims = []
            for member in sorted(set(members)):
                claim = deepcopy(reviewed.get((record["source_id"], member)))
                if claim is None:
                    claim = _inherited_claim(record, member)
                else:
                    claim["curation_status"] = "reviewed_gene_claim"
                gene_claims.append(claim | {"gene_id": member})
            nominations = [
                c["gene_id"] for c in gene_claims if "cta_nomination" in c["allowed_uses"]
            ]
            rows.append(
                record
                | {
                    "recorded_entity_id": record["entity_id"],
                    "recorded_gene_ids": record["gene_ids"],
                    "entity_id": entity,
                    "protein": label,
                    "gene_ids": sorted(set(members)),
                    "nomination_gene_ids": nominations,
                    "published_nomination": bool(nominations),
                    "gene_claims": gene_claims,
                    "allowed_uses": sorted({u for c in gene_claims for u in c["allowed_uses"]}),
                    "eligible_for_coverage": any(c["eligible_for_coverage"] for c in gene_claims),
                    "recorded_claim_scope": record["claim_scope"],
                    "claim_scope": " | ".join(dict.fromkeys(c["scope"] for c in gene_claims)),
                }
            )
    return pd.DataFrame(
        rows,
        columns=[
            *links.columns,
            "recorded_entity_id",
            "recorded_gene_ids",
            "recorded_claim_scope",
            "gene_claims",
            "allowed_uses",
        ],
    )


def _inherited_claim(record, member):
    """Preserve legacy scope; never pretend an imported role is a new review."""
    roles = record["source_roles"]
    uses = []
    if record["eligible_for_coverage"]:
        for role, use in (
            ("tumor_expression", "tumor_support"),
            ("cancer_embryo_program", "oncofetal_support"),
            ("reproductive_developmental", "reproductive_support"),
        ):
            if role in roles:
                uses.append(use)
        if member in record["nomination_gene_ids"]:
            uses.append("cta_nomination")
    for role, use in (
        ("normal_tissue", "normal_tissue_review"),
        ("function", "functional_context"),
        ("translation_identity", "translation_context"),
        ("immune_recognition", "immune_context"),
    ):
        if role in roles:
            uses.append(use)
    native = [c for c in record["publication_claims"] if c["Ensembl_Gene_ID"] == member]
    scope = (
        " | ".join(
            dict.fromkeys(
                c["validation_scope"]
                + ("; source subset: " + c["source_subset"] if c["source_subset"] else "")
                for c in native
            )
        )
        if native
        else record["claim_scope"]
    )
    return {
        "gene_ids": [member],
        "protein": record["protein"],
        "scope": scope,
        "allowed_uses": sorted(uses),
        "eligible_for_coverage": record["eligible_for_coverage"],
        "curation_status": "imported_scope",
        "reported_as": scope,
        "context": "; ".join(dict.fromkeys(c["source_subset"] for c in native)) if native else None,
        "source_location": deepcopy(native),
        "limitations": record["membership_record"].get(
            "identity_limit",
            record["source_record"].get("limitations", "See original source record"),
        ),
    }


def validate_publications(data):
    """Validate the current paper-centric ledger and each exact gene claim."""
    if data.get("schema") != "oncoref.cta-literature-followup" or data.get("schema_version") != 2:
        raise ValueError("Unsupported publication curation schema")
    seen = set()
    for source in data["sources"]:
        sid = source["source_id"]
        if sid in seen or not all(
            source.get(k) for k in ("citation", "title", "url", "gene_claims")
        ):
            raise ValueError("Publications require unique source IDs, citations and gene claims")
        seen.add(sid)
        genes = set()
        for claim in source["gene_claims"]:
            ids = claim.get("gene_ids", [])
            if not ids or any(not re.fullmatch(r"ENSG\d{11}", gene) for gene in ids):
                raise ValueError("Publication claims require exact Ensembl gene IDs")
            if genes.intersection(ids) or len(set(ids)) != len(ids):
                raise ValueError("Consolidate each paper's claims for a gene in one record")
            genes.update(ids)
            required = (
                "protein",
                "reported_as",
                "scope",
                "context",
                "assay",
                "source_location",
                "limitations",
                "evidence_types",
                "reviewed_on",
            )
            if not all(claim.get(k) for k in required) or claim.get("taxonomy_id") != 9606:
                raise ValueError("Gene claims require human identity, assay, context and scope")
            uses, types = claim.get("allowed_uses"), claim["evidence_types"]
            if not isinstance(uses, list) or not set(uses) <= ALLOWED_USES:
                raise ValueError("Unknown publication use")
            if not isinstance(types, list) or not set(types) <= EVIDENCE_TYPES:
                raise ValueError("Unknown publication evidence type")
            if not isinstance(claim.get("eligible_for_coverage"), bool):
                raise ValueError("Coverage eligibility must be explicit")
            if not claim["eligible_for_coverage"] and set(uses) & set(PAPER_USES.values()):
                raise ValueError("Ineligible claims cannot supply positive antigen support")
            if "cta_nomination" in uses and "cta_nomination" not in types:
                raise ValueError("CTA nomination requires an explicit nomination claim")
            if "translation_support" in uses and not set(types) & {
                "ribosome_profiling",
                "mass_spectrometry",
            }:
                raise ValueError("Translation support requires translation/protein evidence")


def publication_registry() -> dict:
    """Central current curation: papers containing exact gene lists and allowed uses.

    Frozen nomination tables remain source evidence; new reviews and corrections
    belong here. A gene-specific review supersedes inherited scope for that
    paper/gene in current views, without editing the historical record.
    """
    data = json.loads(files("oncoref").joinpath("data/cta-literature-followup.json").read_text())
    validate_publications(data)
    return data


def followup_findings() -> dict:
    """Compatibility tables derived from the single current publication ledger."""
    data = publication_registry()
    sources, memberships = [], []
    for source in data["sources"]:
        sources.append({k: v for k, v in source.items() if k != "gene_claims"})
        memberships.extend(c | {"source_id": source["source_id"]} for c in source["gene_claims"])
    return data | {"schema_version": 1, "sources": sources, "memberships": memberships}


def paper_gene_claims(gene: str | None = None, *, use: str | None = None) -> pd.DataFrame:
    """Current exact-gene claims, citations, native contexts and permitted uses.

    ``reviewed_gene_claim`` is centrally adjudicated; ``imported_scope`` preserves
    prior source/table curation. Neither status asserts clinical targetability.
    """
    if use is not None and use not in ALLOWED_USES:
        raise ValueError(f"Unknown publication use: {use}")
    rows = {}
    for record in current_attestations(gene).to_dict("records"):
        for claim in record["gene_claims"]:
            # The reviewed claim also appears on historical occurrences of this DOI.
            if (
                claim["curation_status"] == "reviewed_gene_claim"
                and record["origin"] != "literature_followup"
            ):
                continue
            if use is not None and use not in claim["allowed_uses"]:
                continue
            row = {
                k: record[k] for k in ("entity_id", "source_id", "source_url", "citation", "title")
            }
            row.update(claim)
            rows[json.dumps(row, sort_keys=True)] = row
    return pd.DataFrame(
        list(rows.values()),
        columns=[
            "entity_id",
            "source_id",
            "source_url",
            "citation",
            "title",
            "gene_id",
            "protein",
            "reported_as",
            "reported_name",
            "context",
            "assay",
            "evidence_types",
            "allowed_uses",
            "eligible_for_coverage",
            "curation_status",
            "source_location",
            "limitations",
            "scope",
            "taxonomy_id",
            "reviewed_on",
            "provenance",
        ],
    )


def papers(gene: str | None = None) -> pd.DataFrame:
    """One row per current paper, with its exact gene list and per-gene claims."""
    claims = paper_gene_claims(gene)
    rows = []
    for source, group in claims.groupby("source_id", sort=True):
        row = group.iloc[0]
        rows.append(
            {
                "source_id": source,
                "citation": row.citation,
                "title": row.title,
                "url": row.source_url,
                "gene_ids": sorted(group.gene_id.unique()),
                "gene_claims": group.to_dict("records"),
            }
        )
    return pd.DataFrame(
        rows, columns=["source_id", "citation", "title", "url", "gene_ids", "gene_claims"]
    )


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
    from . import target_selection

    profiles = target_selection.results(category="CTA").rename(
        columns={"protein_group_label": "label"}
    )
    profiles["gene_ids"] = profiles.member_gene_ids.map(json.loads)
    links = current_attestations().to_dict("records")
    nominations = [
        row | {"gene_ids": row["nomination_gene_ids"]}
        for row in links
        if row["published_nomination"]
    ]
    gene_links = [
        {
            "source_id": row["source_id"],
            "gene_ids": [claim["gene_id"]],
            "eligible_for_coverage": claim["eligible_for_coverage"],
        }
        for row in links
        for claim in row["gene_claims"]
    ]
    return _reconcile_coverage(profiles.to_dict("records"), nominations, gene_links, [])


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
