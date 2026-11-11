#!/usr/bin/env python3
"""Gene-by-gene sanity review of every Target candidate against independent sources.

For each candidate this collects the highest non-reproductive value in GTEx v11,
HPA single-cell clusters, HPA immunohistochemistry and every Tabula Sapiens cell
population (immune and stromal cells included, which the donor-level rule leaves
out), then attaches the hand-written review note where one exists. It reports;
it does not classify.
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
CATALOG = ROOT / "docs/audits/target-selection-20261110/rebuild/protein-decisions.csv.gz"
OUT = ROOT / "docs/audits/target-selection-20261110/sanity-review"
SOURCES = Path.home() / ".cache/oncoref/sources"
GTEX = ROOT / "tmp/cta-source-consolidation/gtex-v11-median-tpm.gct.gz"
TABULA = ROOT / "tmp/normal-atlas/tabula-sapiens-2025-11-08/donor-populations.parquet"
TABULA_COLLECTION = "e5f58829-1a66-40b5-a624-9046778e74f5"
MIN_CELLS = 200
GTEX_REPRODUCTIVE = ("Testis", "Ovary", "Uterus", "Vagina", "Cervix", "Fallopian", "Prostate", "Breast", "Cells_")  # fmt: skip
TABULA_REPRODUCTIVE = {"testis", "ovary", "uterus", "prostate gland"}
HPA_REPRODUCTIVE_CELLS = {
    "Spermatogonia", "Spermatocytes", "Early spermatids", "Late spermatids", "Oocytes",
    "Sertoli cells", "Leydig cells", "Peritubular cells", "Granulosa cells", "Theca cells",
    "Cytotrophoblasts", "Syncytiotrophoblasts", "Extravillous trophoblasts", "Hofbauer cells",
    "Endometrial stromal cells", "Ovarian stromal cells", "Prostatic glandular cells",
    "Glandular and luminal cells", "Breast glandular cells", "Breast myoepithelial cells",
    "Basal prostatic cells",
}  # fmt: skip
HPA_REPRODUCTIVE_TISSUES = {
    "testis", "ovary", "placenta", "endometrium 1", "endometrium 2", "fallopian tube", "cervix",
    "vagina", "epididymis", "seminal vesicle", "prostate", "breast", "lactating breast",
}  # fmt: skip
#: Gene families whose normal tissue is absent from every bulk and single-cell panel used.
FAMILIES = [
    ("Olfactory receptors", r"^OR\d", "Olfactory epithelium is not sampled."),
    ("Taste, vomeronasal and trace-amine receptors", r"^(TAS[12]R|VN1R|TAAR)", "Taste buds and nasal epithelium are not sampled."),
    ("Hair keratins and keratin-associated proteins", r"^(KRTAP|KRT3\d|TCHHL)", "Hair follicle; treated as Disposable, replay still open (#687)."),
    ("Cornified-envelope proteins", r"^LCE", "Upper epidermis."),
    ("Enamel and dentin proteins", r"^(AMEL|AMBN|MMP20|DSPP|ENAM)", "Tooth; adult ameloblasts are gone, odontoblasts are not."),
    ("Male reproductive-tract secretions", r"^(DEFB|SPAG11|EDDM|WFDC|PATE|CST8|CST9L|CST11|LCN9|RNASE9|RNASE1[0-3]|SEMG|LYZL)", "Epididymis, seminal vesicle and prostate."),
    ("Unnamed gene records", r"^ENSG", "No symbol; several lie inside another gene."),
]  # fmt: skip
CONCERNS = {
    "normal_cells": "Normal non-reproductive cells express it in an independent atlas",
    "immune_cells": "Blood or immune cells express it in Tabula Sapiens (outside the donor panel by design)",
    "parent_gene": "Counts are likely those of a parent or host gene, in tumors as much as in normal tissue",
    "unsampled_tissue": "Known normal tissue is absent or diluted in every panel",
    "antibody_only": "Normal IHC without any RNA; probably antibody cross-reaction",
    "accepted_context": "Explained by a permitted tissue or a known technical effect",
}


def family(symbol):
    return next((name for name, pattern, _ in FAMILIES if re.search(pattern, symbol)), "")


def dossier(proteins):
    proteins = proteins.copy()
    proteins["genes"] = proteins.member_gene_ids.map(json.loads)
    genes = {g for members in proteins.genes for g in members}
    gtex = pd.read_csv(GTEX, sep="\t", skiprows=2)
    gtex["gene"] = gtex.Name.str.split(".").str[0]
    gtex = gtex[gtex.gene.isin(genes)].set_index("gene").drop(columns=["Name", "Description"])
    gtex = gtex[[c for c in gtex.columns if not c.startswith(GTEX_REPRODUCTIVE)]].astype(float)
    cells = pd.read_csv(SOURCES / "hpa_single_cell/v23/rna_single_cell_type.tsv", sep="\t")
    cells = cells[cells.Gene.isin(genes) & ~cells["Cell type"].isin(HPA_REPRODUCTIVE_CELLS)]
    ihc = pd.read_csv(SOURCES / "hpa_normal_tissue/v23/normal_tissue.tsv", sep="\t")
    ihc = ihc[
        ihc.Gene.isin(genes)
        & ihc.Level.isin(["Low", "Medium", "High"])
        & ~ihc.Tissue.isin(HPA_REPRODUCTIVE_TISSUES)
    ]
    tabula = pd.read_parquet(TABULA)
    tabula = tabula[
        tabula.gene_id.isin(genes)
        & tabula.n_cells.ge(MIN_CELLS)
        & ~tabula.tissue_general.isin(TABULA_REPRODUCTIVE)
    ]
    rows = []
    for p in proteins.itertuples():
        row = {
            "protein": p.protein_group_label,
            "entity_id": p.entity_id,
            "classification": p.classification,
            "published_nomination": bool(p.published_nomination),
            "family": family(p.protein_group_label),
        }
        bulk = gtex.reindex(p.genes).max().dropna()
        row["gtex_top_tissue"] = bulk.idxmax() if len(bulk) else ""
        row["gtex_top_tpm"] = round(float(bulk.max()), 2) if len(bulk) else np.nan
        row["gtex_tissues_ge1"] = int(bulk.ge(1).sum())
        top = cells[cells.Gene.isin(p.genes)].sort_values("nTPM", ascending=False).head(1)
        row["hpa_cell_type"] = top["Cell type"].iloc[0] if len(top) else ""
        row["hpa_cell_ntpm"] = float(top.nTPM.iloc[0]) if len(top) else np.nan
        stain = ihc[ihc.Gene.isin(p.genes)]
        row["ihc_positive_contexts"] = len(stain)
        row["ihc_tissues"] = "; ".join(sorted(set(stain.Tissue))[:6])
        seen = tabula[tabula.gene_id.isin(p.genes)].sort_values("pct_ge2", ascending=False)
        signal = seen[seen.pct_ge2.ge(1)]
        best = seen.head(1)
        row["tabula_population"] = (
            f"{best.tissue_general.iloc[0]}: {best.cell_type.iloc[0]}" if len(best) else ""
        )
        row["tabula_detection_pct"] = round(float(best.pct.iloc[0]), 1) if len(best) else 0.0
        row["tabula_multi_molecule_pct"] = (
            round(float(best.pct_ge2.iloc[0]), 1) if len(best) else 0.0
        )
        row["tabula_donors_with_signal"] = int(signal.donor_id.nunique())
        for column in [
            "tumor_max_tpm", "tumor_top_cancer_code", "independent_patients_above_normal",
            "normal_reference_tpm", "reproductive_ntpm", "paralog_ambiguous", "host_overlapped",
            "rare_cell_status", "rare_cell_basis", "finding_tissue", "finding_top_cell_type",
        ]:  # fmt: skip
            row[column] = getattr(p, column)
        rows.append(row)
    return pd.DataFrame(rows)


def table(frame, columns):
    head = "| " + " | ".join(columns.values()) + " |\n|" + "---|" * len(columns) + "\n"
    return head + "\n".join(
        "| " + " | ".join("" if pd.isna(v) else str(v) for v in row) + " |"
        for row in frame[list(columns)].itertuples(index=False)
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.parse_args()
    catalog = pd.read_csv(CATALOG, low_memory=False)
    notes = pd.read_csv(OUT / "manual-notes.csv")
    if notes.symbol.duplicated().any() or not set(notes.concern) <= set(CONCERNS):
        raise ValueError("Review notes need unique symbols and known concerns")
    missing = set(notes.symbol) - set(catalog.protein_group_label)
    if missing:
        raise ValueError(f"Review notes name unknown proteins: {sorted(missing)}")
    reviewed = catalog[
        catalog.is_target_candidate.astype(bool) | catalog.protein_group_label.isin(notes.symbol)
    ]
    data = dossier(reviewed).merge(notes, left_on="protein", right_on="symbol", how="left")
    data = data.drop(columns="symbol").sort_values("protein")
    data.to_csv(OUT / "candidates.csv", index=False)
    candidate = data.classification.isin(
        ["CTA", "Cancer Antigen only", "Reproductive Antigen only"]
    )
    current = data[candidate]
    quiet = (
        current.gtex_top_tpm.fillna(0).lt(1)
        & current.hpa_cell_ntpm.fillna(0).lt(3)
        & current.ihc_positive_contexts.eq(0)
        & current.tabula_multi_molecule_pct.lt(1)
    )
    families = (
        current[current.family.ne("")]
        .groupby("family")
        .agg(
            proteins=("protein", "size"),
            cta=("classification", lambda s: int(s.eq("CTA").sum())),
            cancer_only=("classification", lambda s: int(s.eq("Cancer Antigen only").sum())),
            median_tumor_tpm=("tumor_max_tpm", "median"),
            median_patients=("independent_patients_above_normal", "median"),
        )
        .round(1)
        .reset_index()
    )
    families["why"] = families.family.map({name: why for name, _, why in FAMILIES})
    flagged = data[data.concern.notna()].copy()
    flagged["now"] = flagged.classification
    sections = []
    for concern, title in CONCERNS.items():
        part = flagged[flagged.concern.eq(concern)]
        still = int(
            part.now.isin(["CTA", "Cancer Antigen only", "Reproductive Antigen only"]).sum()
        )
        sections.append(
            f"### {title}\n\n{len(part)} reviewed; {still} are candidates in the current catalog.\n\n"
            + table(part, {"protein": "Protein", "now": "Current class", "note": "What I saw"})
        )
    weak = current[
        current.classification.ne("Reproductive Antigen only")
        & (current.tumor_max_tpm.lt(5) | current.independent_patients_above_normal.le(2))
    ]
    (OUT / "index.md").write_text(
        "# Gene-by-gene sanity review of the candidate sets\n\n"
        f"Every one of the {len(current)} Target candidates was checked against four sources "
        "that the classification does not use in this form: GTEx v11 tissue medians, HPA v23 "
        "single-cell clusters, HPA v23 immunohistochemistry and all "
        "1,093,048 Tabula Sapiens UMI cells (23 donors, 30 tissue groups, immune and stromal "
        "cells included). Review notes are hand-written, one per protein, in "
        "[manual-notes.csv](manual-notes.csv). [Full table](candidates.csv). "
        "These are outputs of screening rules, not validated tumor antigens.\n\n"
        "## Result\n\n"
        f"- **{int(quiet.sum())} of {len(current)}** candidates show nothing in any of the four "
        "sources: under 1 TPM in every non-reproductive GTEx tissue, under 3 nTPM in every "
        "non-reproductive HPA cell cluster, no positive normal IHC and no Tabula Sapiens "
        "population with 1% of cells carrying two molecules.\n"
        f"- **{len(flagged)}** proteins have a review note; "
        f"{int(flagged.now.isin(['CTA', 'Cancer Antigen only', 'Reproductive Antigen only']).sum())} "
        "of them are candidates in the current catalog.\n"
        f"- **{len(weak)}** Cancer Antigens rest on a tumor maximum under 5 TPM or on at most two "
        "patients.\n\n"
        "## Families whose normal tissue no panel samples\n\n"
        "Absence from the sampled tissues is not absence from the body. These families pass "
        "because their tissue is missing, and most qualify as Cancer Antigens on sparse tumor "
        "counts.\n\n"
        + table(
            families,
            {
                "family": "Family",
                "proteins": "Candidates",
                "cta": "CTA",
                "cancer_only": "Cancer Antigen only",
                "median_tumor_tpm": "Median tumor max (TPM)",
                "median_patients": "Median patients above normal",
                "why": "Normal tissue",
            },
        )
        + "\n\n## Reviewed proteins\n\n"
        + "\n\n".join(sections)
        + "\n\n## Limits\n\n"
        "- Tabula Sapiens is one study; its donors share libraries and a reference, so a signal "
        "there is replicated across donors but not across studies.\n"
        "- A monocyte population of one donor (TSP21 lymph node, 237 cells) shows low counts for "
        "many testis genes; single-population signals there were not treated as findings.\n"
        "- The notes use general knowledge of each gene. They are a reviewer's reading, not a "
        "literature search, and they change no category.\n"
    )
    print(
        json.dumps({"candidates": len(current), "quiet": int(quiet.sum()), "noted": len(flagged)})
    )


if __name__ == "__main__":
    main()
