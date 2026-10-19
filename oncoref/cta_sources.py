"""Published CTA/placental nominations, separate from antigen validation.

Gong's complete placenta-enriched S5/S6 lists are retained, including historical
noncoding and unmapped rows. Only currently protein-coding mapped genes enter
the HPA candidate table; source membership never overrides specificity filters.
"""

from __future__ import annotations

import hashlib

import pandas as pd

from .cta_landscape import PAPERS
from .cta_tissues import cta_dataset_name
from .gene_ids import canonical_gene_id, canonical_gene_space
from .load_dataset import get_data

GONG_PC = "Gong2021_placenta_PC"
GONG_NC = "Gong2021_placenta_ncRNA"
BRADLEY = "Bradley2020_CPA"
SOURCE_LABELS = {
    GONG_PC: "Gong 2021 · S5 (historically coding)",
    GONG_NC: "Gong 2021 · S6 (historically noncoding)",
    BRADLEY: "Bradley 2020 · CPA",
}
SOURCE_LABELS.update({tag: paper[0] for tag, paper in PAPERS.items()})

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
    lists from the papers; use publication_membership() for complete lists and
    cta_provenance.candidate_provenance() for panel-wide citation coverage.
    A combined CGB1/CGB2 result is never counted as separately measured positives.
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


def placental_source_coverage(*, tissue_scope="core") -> pd.DataFrame:
    """Coverage of the prior placental nominations, independent of CTA admission."""
    from . import cta

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
    cta_dataset_name(tissue_scope)
    default = cta.cta_gene_ids() if tissue_scope == "core" else cta.cta_extended_gene_ids()
    prior["default_panel"] = prior.Ensembl_Gene_ID.isin(default)
    prior["coverage_meaning"] = "Nomination/expression provenance; not antigen validation"
    return prior


def extract_publications(gong_workbook) -> pd.DataFrame:
    """Extract exactly 71 S5 and 74 S6 placenta rows; transcribe Bradley Fig. 3.

    Source IDs/labels/biotypes remain unchanged. Resolve Gong by its source
    Ensembl ID, never by a potentially ambiguous historical symbol fallback.
    """
    # The row-count checks alone cannot detect a different or edited workbook.
    # Pin the original published input before interpreting any of its cells.
    from pathlib import Path

    sources = publication_sources()
    expected = set(sources.loc[sources.source_tag.isin([GONG_PC, GONG_NC]), "input_sha256"])
    actual = hashlib.sha256(Path(gong_workbook).read_bytes()).hexdigest()
    if len(expected) != 1 or actual not in expected:
        raise ValueError("Gong workbook checksum differs from the curated publication source")

    import openpyxl

    workbook = openpyxl.load_workbook(gong_workbook, read_only=True, data_only=True)
    rows = []
    try:
        for sheet, tag, biotype, expected in (
            ("Data 5 - tissue-enriched PC", GONG_PC, "protein_coding", 71),
            ("Data 6 tissue-enriched lncR", GONG_NC, "noncoding", 74),
        ):
            selected = [
                r for r in workbook[sheet].values if str(r[0]).strip().lower() == "placenta"
            ]
            distinct = len({r[1] for r in selected})
            if len(selected) != expected:
                raise ValueError(
                    f"Unexpected {sheet} membership: {len(selected)} placenta rows, expected {expected}"
                )
            if distinct != expected:
                raise ValueError(
                    f"Unexpected {sheet} membership: {len(selected)} placenta rows carry only "
                    f"{distinct} distinct source Ensembl IDs, expected {expected}"
                )
            for row in selected:
                if not row[1] or not row[2]:
                    raise ValueError(f"Missing source Ensembl ID or symbol in {sheet}: {row[:3]}")
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
    finally:
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
    # Distinct historical symbols within one source can still canonicalize onto a
    # single gene; that would double-count it in the intake funnel and make the
    # Symbol chosen for a new candidate row arbitrary, so reject it explicitly.
    mapped = result.loc[result.Ensembl_Gene_ID.astype(str).ne("")]
    collisions = mapped.loc[mapped.duplicated(["source_tag", "Ensembl_Gene_ID"], keep=False)]
    if not collisions.empty:
        detail = "; ".join(
            f"{tag}: {sorted(grp.source_symbol)} -> {gid}"
            for (tag, gid), grp in collisions.groupby(["source_tag", "Ensembl_Gene_ID"])
        )
        raise ValueError(f"Distinct source symbols canonicalize to one gene ({detail})")
    return result


def _with_tissue_scope(table: pd.DataFrame, tissue_scope: str) -> pd.DataFrame:
    """Identify legacy core evidence and reject incomplete or conflicting labels."""
    result = table.copy()
    column = "rna_reproductive_tissue_scope"
    if column not in result:
        # Scope-free snapshots predate the extended panel and contain core RNA
        # fractions. Normalize before concatenation, even when no rows are new.
        result[column] = "core"
    scopes = result[column].astype("string")
    if not scopes.eq(tissue_scope).fillna(False).all():
        observed = sorted(scopes.fillna("<missing>").unique())
        raise ValueError(
            f"Expected RNA tissue_scope={tissue_scope!r}; found {observed}. "
            "Legacy tables without a scope column contain core evidence. "
            "Use regenerate_cta_columns(table, tissue_scope=...) to compute a different "
            "scope or repair incomplete labels; changing labels alone does not recompute RNA."
        )
    return result


def add_publication_candidates(table: pd.DataFrame, membership: pd.DataFrame) -> pd.DataFrame:
    """Union coding nominations, annotate existing rows, regenerate only new rows.

    Existing annotations/filter decisions are preserved. Absent HPA RNA remains
    missing and fails the restriction gate; it is never treated as zero signal.
    New canonical transcript IDs are left unassigned, not guessed.
    Legacy input without ``rna_reproductive_tissue_scope`` is core evidence;
    the returned table explicitly labels every row ``core``, even on a no-op
    import. Explicit missing, mixed or non-core scope labels are rejected.
    """
    from .cta_regen import regenerate_cta_columns

    result = _with_tissue_scope(table, "core")
    result["Symbol"] = result["Symbol"].fillna(result["Ensembl_Gene_ID"])
    ids = result.Ensembl_Gene_ID.astype(str).str.split(".").str[0]
    if ids.duplicated().any():
        raise ValueError("Duplicate candidate Ensembl IDs")
    existing = set(ids)
    membership = membership.copy()
    membership["Ensembl_Gene_ID"] = (
        membership.Ensembl_Gene_ID.fillna("").astype(str).str.split(".").str[0]
    )
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
        # New rows start without RNA values. Regeneration leaves measurements
        # missing and fails the restriction gate when HPA has no observation.
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
        added = tags_by_id.get(str(row.Ensembl_Gene_ID).split(".")[0], set())
        if added:
            old = (
                {t.strip() for t in row.source_databases.split(";") if t.strip()}
                if pd.notna(row.source_databases)
                else set()
            )
            result.at[idx, "source_databases"] = ";".join(sorted(old | added))
    return result


def intake_membership(membership=None, table=None, *, tissue_scope="core") -> pd.DataFrame:
    """Source-row audit through mapping, coding, family, HPA and public defaults.

    A supplied table must contain evidence for ``tissue_scope``. Legacy tables
    without a scope column are core only; the argument selects the matching
    review policy, and never recalculates or relabels RNA fractions.
    """
    from . import cta
    from .gene_families import gene_family_ids

    refs = publication_membership() if membership is None else membership.copy()
    dataset = cta_dataset_name(tissue_scope)
    raw = _with_tissue_scope(get_data(dataset) if table is None else table, tissue_scope)
    raw["Ensembl_Gene_ID"] = raw.Ensembl_Gene_ID.astype(str).str.split(".").str[0]
    refs["Ensembl_Gene_ID"] = refs.Ensembl_Gene_ID.fillna("").astype(str).str.split(".").str[0]
    # Apply the same family/expression/specificity rules to this table, including
    # dry-run candidates. Mixing it with global shipped sets hides new rows.
    family = ~raw.Ensembl_Gene_ID.isin(gene_family_ids("histone")) & ~raw.Symbol.map(
        cta._alpha_tubulin_symbol
    )
    curated = cta._with_specificity_columns(
        raw.loc[family].reset_index(drop=True), tissue_scope=tissue_scope
    )
    unfiltered = set(curated.Ensembl_Gene_ID)
    default = set(curated.loc[cta._canonical_default_mask(curated), "Ensembl_Gene_ID"])
    hpa_pass = set(raw.loc[cta.passes_filters_mask(raw), "Ensembl_Gene_ID"])
    refs["mapped"] = refs.Ensembl_Gene_ID.ne("")
    refs["protein_coding"] = refs.mapped & refs.biotype.eq("protein_coding")
    refs["family_eligible"] = refs.protein_coding & refs.Ensembl_Gene_ID.isin(unfiltered)
    refs["hpa_restriction"] = refs.family_eligible & refs.Ensembl_Gene_ID.isin(hpa_pass)
    refs["default_panel"] = refs.family_eligible & refs.Ensembl_Gene_ID.isin(default)
    if (refs.default_panel & ~refs.hpa_restriction).any():
        raise ValueError("Specificity overrides require a branching source funnel")
    return refs


def intake_counts(membership=None, table=None, *, tissue_scope="core") -> pd.DataFrame:
    refs = intake_membership(membership, table, tissue_scope=tissue_scope)
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
            remaining = (
                len(group)
                if stage == "published"
                else group.loc[group[stage], "Ensembl_Gene_ID"].nunique()
            )
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
