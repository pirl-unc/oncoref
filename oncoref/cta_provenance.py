"""Paper-backed CTA intake and an exact, explicitly bounded source cover.

The optimization covers the screened panel, then takes *all* entries from the
selected papers upstream of mapping and filtering. A citation is nomination
provenance, never proof of locus-specific protein, HLA presentation or safety.
"""

from itertools import combinations

import pandas as pd

from .cta_sources import publication_membership, publication_sources
from .load_dataset import get_data


def minimum_covers(target, source_sets):
    """Enumerate every minimum-cardinality cover, deterministically.

    Appropriate for the small, curated publication registry. Unlike a greedy
    approximation this establishes minimality; uncovered targets fail closed.
    Keys must identify papers (DOIs), not individual supplementary tables.
    """
    target = set(target)
    sets = {key: set(value) & target for key, value in sorted(source_sets.items())}
    uncovered = target - set().union(*sets.values())
    if uncovered:
        raise ValueError(f"Retained genes lack paper provenance: {sorted(uncovered)}")
    if not target:
        return [()]
    keys = [key for key, value in sets.items() if value]
    for size in range(1, len(keys) + 1):
        covers = [
            combo
            for combo in combinations(keys, size)
            if target <= set().union(*(sets[key] for key in combo))
        ]
        if covers:
            return covers
    raise AssertionError("A covered finite universe must have a cover")


def paper_membership():
    """Complete source entries, with paper identifiers and exact source locations."""
    sources = publication_sources()
    if sources.source_tag.duplicated().any() or sources.doi.isna().any():
        raise ValueError("Publication registry requires unique source tags and DOIs")
    refs = publication_membership().merge(
        sources[["source_tag", "doi", "citation", "source_url", "source_table", "input_sha256"]],
        on="source_tag",
        how="left",
        validate="many_to_one",
    )
    if refs.doi.isna().any():
        raise ValueError("Publication membership has no registered paper")
    return refs


def source_cover():
    """Minimum cover of the current default among fully imported candidate lists.

    Targeted citations whose complete candidate lists were not retrieved are
    available as evidence but cannot masquerade as complete source universes.
    Ties choose lexicographically by DOI and all alternatives remain recorded.
    """
    from .cta import cta_gene_ids

    refs = paper_membership()
    coding = refs[refs.protein_candidate_eligible]
    sets = {doi: set(group.Ensembl_Gene_ID) for doi, group in coding.groupby("doi")}
    target = set(cta_gene_ids())
    covers = minimum_covers(target, sets)
    chosen = covers[0]
    union = set().union(*(sets[doi] for doi in chosen))
    witnesses = {
        doi: sorted(
            (sets[doi] & target)
            - set().union(*(members for other, members in sets.items() if other != doi))
        )
        for doi in chosen
    }
    return {
        "scope": "Exact minimum among the registered complete candidate-source papers; not all literature",
        "target_definition": "oncoref.cta.cta_gene_ids() default after screening full imported intake",
        "papers_considered": sorted(sets),
        "minimum_papers": len(chosen),
        "selected_dois": list(chosen),
        "all_minimum_covers": [list(combo) for combo in covers],
        "indispensable_witness_gene_ids": witnesses,
        "retained_genes": len(target),
        "candidate_genes": len(union),
        "candidate_gene_ids": sorted(union),
    }


def selected_membership():
    """Full lists of the selected papers, including noncoding/unmapped entries."""
    refs = paper_membership()
    return refs[refs.doi.isin(source_cover()["selected_dois"])].copy()


def candidate_evidence():
    """Active, paper-backed coding candidate pool, before any HPA filtering."""
    ids = set(source_cover()["candidate_gene_ids"])
    raw = get_data("cancer-testis-antigens")
    if not ids <= set(raw.Ensembl_Gene_ID):
        raise ValueError("Unassessed genes in the complete paper union")
    return raw[raw.Ensembl_Gene_ID.isin(ids)].copy()


def legacy_only_candidates():
    """Historical rows outside the selected paper union, retained for audit only."""
    raw = get_data("cancer-testis-antigens")
    result = raw[~raw.Ensembl_Gene_ID.isin(source_cover()["candidate_gene_ids"])].copy()
    result["provenance_status"] = "outside_selected_complete_paper_union"
    return result


def candidate_provenance():
    """One row per active coding candidate with DOI and exact row-level links.

    The accompanying selected_membership() table preserves assay resolution,
    original identity, mapping method, source rows and validation limitations.
    """
    from .cta import cta_gene_ids

    refs = selected_membership()
    refs = refs[refs.protein_candidate_eligible].fillna("")
    default = cta_gene_ids()
    rows = []
    normal_sources = {
        "daSilva2017_testis_biased",
        "Gong2021_reproductive_PC",
        "Gong2021_reproductive_ncRNA",
        "Gong2021_placenta_PC",
        "Gong2021_placenta_ncRNA",
    }
    for gid, group in refs.groupby("Ensembl_Gene_ID", sort=True):
        normal_only = set(group.source_tag) <= normal_sources
        rows.append(
            {
                "Ensembl_Gene_ID": gid,
                "Symbol": group.Symbol.iloc[0],
                "paper_dois": ";".join(sorted(set(group.doi))),
                "source_tags": ";".join(sorted(set(group.source_tag))),
                "source_locations": " | ".join(
                    sorted(
                        {
                            f"{r.doi}: {r.source_file or r.source_table} / "
                            f"{r.source_sheet or r.source_table} / rows {r.source_rows or 'see named figure'}"
                            for r in group.itertuples(index=False)
                        }
                    )
                ),
                "default_panel": gid in default,
                "nomination_scope": (
                    "normal_reproductive_expression_only"
                    if normal_only
                    else "includes_cancer_associated_source"
                ),
                "legacy_retention_evidence_caveat": gid in default and normal_only,
                "evidence_scope": "Paper nomination/expression; antigen validation must be assessed separately",
            }
        )
    return pd.DataFrame(rows)


def paper_intake_counts():
    """Unique-identity funnel per selected paper, never summing nested tables."""
    from . import cta

    refs = selected_membership().fillna("")
    raw = candidate_evidence()
    family = cta.cta_unfiltered_gene_ids()
    hpa = set(raw.loc[cta.passes_filters_mask(raw), "Ensembl_Gene_ID"]) & family
    default = cta.cta_gene_ids()
    rows = []
    for doi, group in refs.groupby("doi", sort=True):
        identities = {
            r.Ensembl_Gene_ID
            or f"unmapped:{r.source_id_type}:{r.source_gene_id or r.source_symbol}"
            for r in group.itertuples(index=False)
        }
        mapped = set(group.Ensembl_Gene_ID) - {""}
        coding = set(group.loc[group.protein_candidate_eligible, "Ensembl_Gene_ID"])
        for stage, members in (
            ("published", identities),
            ("mapped", mapped),
            ("protein_coding", coding),
            ("family_eligible", coding & family),
            ("hpa_restriction", coding & hpa),
            ("default_panel", coding & default),
        ):
            rows.append({"source_tag": doi, "stage": stage, "remaining": len(members)})
    return pd.DataFrame(rows)


def historical_tag_audit():
    """Reconcile legacy labels without silently rewriting historical evidence."""
    raw = get_data("cancer-testis-antigens")
    refs = paper_membership()
    comparisons = {
        "CTexploreR_CT": {"Loriot2025_S1"},
        "CTexploreR_CTP": {"Loriot2025_S2"},
        "daSilva2017_protein": {"daSilva2017_tumor_proteomics"},
        "CTpedia": set(),
    }
    rows = []
    for tag, source_tags in comparisons.items():
        historical = raw.source_databases.fillna("").str.split(";").map(lambda x, tag=tag: tag in x)
        verified = set(refs.loc[refs.source_tag.isin(source_tags), "Ensembl_Gene_ID"])
        for r in raw[historical].itertuples(index=False):
            rows.append(
                {
                    "Ensembl_Gene_ID": r.Ensembl_Gene_ID,
                    "Symbol": r.Symbol,
                    "historical_tag": tag,
                    "comparison_sources": ";".join(sorted(source_tags)),
                    "status": (
                        "historical_catalog_tag_not_primary_paper_provenance"
                        if not source_tags
                        else "reproduced"
                        if r.Ensembl_Gene_ID in verified
                        else "not_reproduced_by_complete_import"
                    ),
                    "effect": "Audit only; active provenance and plots use exact imported paper memberships",
                }
            )
    return pd.DataFrame(rows)


def gene_citation_evidence_report():
    """Separate nomination and assay modalities for every active candidate.

    Empty evidence columns mean no imported link in that category, not a
    biological negative. Full row-level assay limitations accompany this index.
    Additional targeted/unselected papers can supply evidence without changing
    the minimum-cover starting pool or the current retention policy.
    """
    from .cta_sources import gene_publication_evidence

    report = candidate_provenance()
    refs = paper_membership().fillna("")
    targeted = gene_publication_evidence().fillna("")
    normal = {
        "daSilva2017_testis_biased",
        "Gong2021_reproductive_PC",
        "Gong2021_reproductive_ncRNA",
        "Gong2021_placenta_PC",
        "Gong2021_placenta_ncRNA",
    }
    cancer = {
        "Wang2016_CT",
        "Bruggeman2018_GC",
        "daSilva2017_CT",
        "Jamin2021_CT",
        "Jamin2021_core",
        "Chang2019_TGCT",
        "Carter2023_CT",
        "Loriot2025_S1",
        "Loriot2025_S2",
        "Bai2016_EGFL6",
        "Bradley2020_CPA",
    }
    categories = {
        "published_list_nomination_sources": set(refs.source_tag),
        "normal_reproductive_expression_sources": normal
        | (cancer - {"Bai2016_EGFL6", "Bradley2020_CPA"}),
        "cancer_expression_nomination_sources": cancer,
        "protein_evidence_sources": {"daSilva2017_tumor_proteomics"},
    }
    for column, tags in categories.items():
        grouped = refs[refs.source_tag.isin(tags)].groupby("Ensembl_Gene_ID").doi.agg(set)
        report[column] = report.Ensembl_Gene_ID.map(grouped).map(
            lambda value: ";".join(sorted(value)) if isinstance(value, set) else ""
        )
    # This is the sole newly validated target in Bradley Figure 3. Its peptide
    # and T-cell findings must not propagate to the other nine CPA candidates.
    validated = refs[refs.source_tag.eq("Bradley2020_CPA") & refs.source_symbol.eq("VGLL1")]
    links = validated.groupby("Ensembl_Gene_ID").doi.agg(lambda x: ";".join(sorted(set(x))))
    for column in ("hla_peptide_evidence_sources", "antigen_specific_t_cell_evidence_sources"):
        report[column] = report.Ensembl_Gene_ID.map(links).fillna("")
    extra = {
        "normal_reproductive_expression_sources": {
            "Amoushahi2020_targeted",
            "Rull2005_CGB_placenta",
        },
        "cancer_expression_nomination_sources": {
            "Park2010_targeted",
            "Wollenzien2023_targeted",
            "Traynor2022_targeted",
            "Gure2002_targeted",
            "Song2022_targeted",
            "Shukla2018_targeted",
            "Caballero2014_targeted",
            "Kristensen2008_targeted",
            "Bialas2020_CGB_cancer",
            "Kubiczak2013_CGB_ovarian",
            "McKellar2025_CGB7_cancer",
        },
        "protein_evidence_sources": {"Song2022_targeted", "Amoushahi2020_targeted"},
    }
    for column, tags in extra.items():
        links = targeted[targeted.source_tag.isin(tags)].groupby("Ensembl_Gene_ID").doi.agg(set)
        report[column] = [
            ";".join(sorted((set(value.split(";")) - {""}) | links.get(gid, set())))
            for gid, value in zip(report.Ensembl_Gene_ID, report[column])
        ]
    extra_links = targeted.groupby("Ensembl_Gene_ID").doi.agg(lambda x: ";".join(sorted(set(x))))
    report["additional_targeted_citations"] = report.Ensembl_Gene_ID.map(extra_links).fillna("")
    report["interpretation"] = (
        "Imported links only; empty is not assessed. RNA/protein labels may be shared across paralogs; "
        "see source-row assay_resolution and targeted evidence limitations. No locus-specific validation inferred."
    )
    return report
