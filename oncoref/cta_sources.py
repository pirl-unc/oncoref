"""Published CTA/placental nominations, separate from antigen validation.

Gong's complete placenta-enriched S5/S6 lists are retained, including historical
noncoding and unmapped rows. Only currently protein-coding mapped genes enter
the HPA candidate table; source membership never overrides specificity filters.
"""

from __future__ import annotations

import pandas as pd

from .gene_ids import canonical_gene_id, canonical_gene_space
from .load_dataset import get_data

GONG_PC = "Gong2021_placenta_PC"
GONG_NC = "Gong2021_placenta_ncRNA"
BRADLEY = "Bradley2020_CPA"
SOURCE_LABELS = {
    GONG_PC: "Gong 2021 · S5 coding",
    GONG_NC: "Gong 2021 · S6 noncoding",
    BRADLEY: "Bradley 2020 · CPA",
}
BRADLEY_GENES = (
    "VGLL1",
    "PLAC1",
    "CGB3",
    "CGB5",
    "IGF2BP3",
    "DEPDC1B",
    "ADAM12",
    "SLC38A9",
    "CAPN6",
    "MMP11",
)


def publication_membership() -> pd.DataFrame:
    """One row per source/gene, including entries ineligible for protein CTAs."""
    return get_data("cta-publication-membership")


def publication_sources() -> pd.DataFrame:
    """Publication, exact table/figure, input checksum and evidence scope."""
    return get_data("cta-publication-sources")


def gene_publication_evidence() -> pd.DataFrame:
    """Targeted primary evidence, with assay resolution and validation limits.

    These are curated evidence rows for specific genes, not complete nomination
    lists from the papers. In particular, a combined CGB1/CGB2 result is never
    counted as a separately measured positive for each gene.
    """
    return get_data("cta-gene-publication-evidence").copy()


def add_gene_evidence_tags(table: pd.DataFrame, evidence=None) -> pd.DataFrame:
    """Annotate already nominated genes without changing evidence/filter values."""
    evidence = gene_publication_evidence() if evidence is None else evidence
    result = table.copy()
    tags = evidence.groupby("Ensembl_Gene_ID").source_tag.agg(set)
    for idx, row in result.iterrows():
        if row.Ensembl_Gene_ID in tags.index:
            old = set(str(row.source_databases).split(";")) - {"", "nan"}
            result.at[idx, "source_databases"] = ";".join(sorted(old | tags[row.Ensembl_Gene_ID]))
    return result


def placental_source_coverage() -> pd.DataFrame:
    """Coverage of the prior placental nominations, independent of CTA admission."""
    from .cta import cta_gene_ids

    raw = get_data("cancer-testis-antigens")
    prior = raw[
        raw.source_databases.fillna("").str.split(";").map(lambda tags: "placental_antigen" in tags)
    ][["Symbol", "Ensembl_Gene_ID"]].copy()
    lists = publication_membership()
    focused = gene_publication_evidence()
    sources = pd.concat(
        [lists[["Ensembl_Gene_ID", "source_tag"]], focused[["Ensembl_Gene_ID", "source_tag"]]]
    )
    tags = sources.groupby("Ensembl_Gene_ID").source_tag.agg(lambda x: ";".join(sorted(set(x))))
    prior["publication_sources"] = prior.Ensembl_Gene_ID.map(tags).fillna("")
    prior["covered_by_publication"] = prior.publication_sources.ne("")
    prior["default_panel"] = prior.Ensembl_Gene_ID.isin(cta_gene_ids())
    prior["coverage_meaning"] = "Nomination/expression provenance; not antigen validation"
    return prior


def extract_publications(gong_workbook) -> pd.DataFrame:
    """Extract exactly 71 S5 and 74 S6 placenta rows; transcribe Bradley Fig. 3.

    Source IDs/labels/biotypes remain unchanged. Resolve Gong by its source
    Ensembl ID, never by a potentially ambiguous historical symbol fallback.
    """
    import openpyxl

    workbook = openpyxl.load_workbook(gong_workbook, read_only=True, data_only=True)
    rows = []
    for sheet, tag, biotype, expected in (
        ("Data 5 - tissue-enriched PC", GONG_PC, "protein_coding", 71),
        ("Data 6 tissue-enriched lncR", GONG_NC, "noncoding", 74),
    ):
        selected = [r for r in workbook[sheet].values if str(r[0]).strip().lower() == "placenta"]
        if len(selected) != expected or len({r[1] for r in selected}) != expected:
            raise ValueError(f"Unexpected {sheet} membership: {len(selected)} rows")
        for row in selected:
            rows.append(
                {
                    "source_tag": tag,
                    "source_gene_id": row[1],
                    "source_symbol": row[2],
                    "source_biotype": biotype,
                    "evidence_type": "placental_RNA_enrichment",
                    "validation_scope": "No tumor-antigen validation in this source",
                }
            )
    workbook.close()
    for symbol in BRADLEY_GENES:
        rows.append(
            {
                "source_tag": BRADLEY,
                "source_gene_id": "",
                "source_symbol": symbol,
                "source_biotype": "not specified",
                "evidence_type": "cancer_placenta_candidate",
                "validation_scope": (
                    "HLA peptide and antigen-specific T-cell validation"
                    if symbol == "VGLL1"
                    else "Expression-nominated CPA; not newly validated by this study"
                ),
            }
        )
    return resolve_membership(pd.DataFrame(rows))


def resolve_membership(rows: pd.DataFrame) -> pd.DataFrame:
    """Canonicalize identities without turning unmapped/noncoding rows into CTAs."""
    result = rows.copy()
    space = canonical_gene_space().set_index("ensembl_gene_id")
    gids, symbols, biotypes = [], [], []
    for row in result.itertuples(index=False):
        identifier = (
            row.source_gene_id
            if pd.notna(row.source_gene_id) and row.source_gene_id
            else row.source_symbol
        )
        gid = canonical_gene_id(identifier)
        gids.append(gid or "")
        symbols.append(space.loc[gid, "symbol"] if gid else "")
        biotypes.append(space.loc[gid, "biotype"] if gid else "")
    result["Ensembl_Gene_ID"] = gids
    result["Symbol"] = symbols
    result["biotype"] = biotypes
    result["mapping_status"] = ["mapped" if gid else "unmapped" for gid in gids]
    result["protein_candidate_eligible"] = result.biotype.eq("protein_coding")
    if result.duplicated(["source_tag", "source_symbol"]).any():
        raise ValueError("Duplicate source/gene nomination")
    return result


def add_publication_candidates(table: pd.DataFrame, membership: pd.DataFrame) -> pd.DataFrame:
    """Union coding nominations, annotate existing rows, regenerate only new rows.

    Existing annotations/filter decisions are preserved. Absent HPA RNA remains
    missing and fails the restriction gate; it is never treated as zero signal.
    New canonical transcript IDs are left unassigned, not guessed.
    """
    from .cta_regen import regenerate_cta_columns

    result = table.copy()
    if result.Ensembl_Gene_ID.duplicated().any():
        raise ValueError("Duplicate candidate Ensembl IDs")
    existing = set(result.Ensembl_Gene_ID)
    eligible = membership[membership.biotype.eq("protein_coding")]
    new = []
    for gid, group in eligible.groupby("Ensembl_Gene_ID", sort=True):
        if not gid or gid in existing:
            continue
        record = dict.fromkeys(result.columns, None)
        record.update(
            Symbol=group.iloc[0].Symbol,
            Ensembl_Gene_ID=gid,
            source_databases=";".join(sorted(set(group.source_tag))),
            biotype="protein_coding",
            Aliases="",
            Full_Name="",
            Function="",
            Canonical_Transcript_ID="",
        )
        new.append(record)
    if new:
        seed = pd.DataFrame(new, columns=result.columns)
        # New rows start without RNA values; the existing regenerator retains
        # these missing values when HPA has no observation, so they fail closed.
        for col in seed.columns:
            if col.startswith("rna_") and (col.endswith("ntpm") or "frac" in col):
                seed[col] = pd.to_numeric(seed[col], errors="coerce")
        regenerated = regenerate_cta_columns(seed)
        result = pd.concat([result, regenerated], ignore_index=True)
    tags_by_id = (
        membership[membership.Ensembl_Gene_ID.fillna("").ne("")]
        .groupby("Ensembl_Gene_ID")["source_tag"]
        .agg(set)
    )
    for idx, row in result.iterrows():
        added = tags_by_id.get(row.Ensembl_Gene_ID, set())
        if added:
            old = {t for t in str(row.source_databases or "").split(";") if t and t != "nan"}
            result.at[idx, "source_databases"] = ";".join(sorted(old | added))
    return result


def intake_membership(membership=None, table=None) -> pd.DataFrame:
    """Source-row audit through mapping, coding, family, HPA and public defaults."""
    from . import cta

    refs = publication_membership() if membership is None else membership.copy()
    raw = get_data("cancer-testis-antigens") if table is None else table
    # Public sets retain reviewed specificity decisions and expression rescues.
    unfiltered = cta.cta_unfiltered_gene_ids()
    default = cta.cta_gene_ids()
    hpa_pass = set(raw.loc[cta.passes_filters_mask(raw), "Ensembl_Gene_ID"])
    refs["mapped"] = refs.Ensembl_Gene_ID.fillna("").ne("")
    refs["protein_coding"] = refs.mapped & refs.biotype.eq("protein_coding")
    refs["family_eligible"] = refs.protein_coding & refs.Ensembl_Gene_ID.isin(unfiltered)
    refs["hpa_restriction"] = refs.family_eligible & refs.Ensembl_Gene_ID.isin(hpa_pass)
    refs["default_panel"] = refs.family_eligible & refs.Ensembl_Gene_ID.isin(default)
    if (refs.default_panel & ~refs.hpa_restriction).any():
        raise ValueError("Specificity overrides require a branching source funnel")
    return refs


def intake_counts(membership=None) -> pd.DataFrame:
    refs = intake_membership(membership)
    rows = []
    for tag, group in refs.groupby("source_tag", sort=False):
        previous = len(group)
        for stage in (
            "published",
            "mapped",
            "protein_coding",
            "family_eligible",
            "hpa_restriction",
            "default_panel",
        ):
            remaining = len(group) if stage == "published" else int(group[stage].sum())
            rows.append(
                {
                    "source_tag": tag,
                    "stage": stage,
                    "remaining": remaining,
                    "dropped": previous - remaining,
                }
            )
            previous = remaining
    return pd.DataFrame(rows)
