#!/usr/bin/env python3
"""Normal nose and skin expression of Target candidates, from two references.

The expression panels behind the catalog hold no nasal tissue, and their bulk skin
dilutes single cell types. This script measures every Target candidate in:

- bulk RNA-seq of human olfactory mucosa (four people) and one nasal respiratory
  sample, Olender et al. 2016, BMC Genomics, Table S2 (FPKM);
- single cells of adult normal nose, skin, scalp and tongue in the pinned CELLxGENE
  Census, counted by ``extract_nose_skin_cell_counts.py``.

It records a tissue for a protein when the reference shows it there in several
people or several datasets. The catalog rebuild carries that record on the row. A
recorded nose tissue is a home tissue for the catalog's Other:Tissue rule.

    python scripts/build_nose_skin_reference.py tmp/normal-atlas/nose-skin-2025-11-08
"""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from oncoref import target_curation  # noqa: E402

DATA = ROOT / "oncoref/data"
OLFACTORY = ["OE7", "OE12", "OE15", "OE17"]
RULES = {
    "olfactory_min_fpkm": 1.0,
    "olfactory_min_samples": 3,
    "min_grade": "medium",
    "cell_min_cells": 200,
    "cell_min_detected_pct": 5.0,
    "cell_min_multi_molecule_pct": 1.0,
    "cell_min_datasets": 2,
    "excluded_cell_types": "melanocyte|unknown",
}
SITES = {"nose": "nasal mucosa", "skin": "skin"}  # sites with at least two datasets


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def candidates():
    from pyensembl import EnsemblRelease

    genome = EnsemblRelease(115)
    catalog = pd.read_csv(DATA / "target-selection.csv.gz", low_memory=False)
    rows = []
    for r in catalog[catalog.is_target_candidate].itertuples():
        for gene_id in json.loads(r.member_gene_ids):
            try:
                symbol = genome.gene_by_id(gene_id).gene_name
            except ValueError:
                symbol = ""
            rows.append((r.entity_id, r.protein_group_label, gene_id, symbol))
    return pd.DataFrame(rows, columns=["entity_id", "protein", "gene_id", "symbol"])


def olfactory(genes, table_path):
    table = pd.read_excel(table_path, sheet_name="Table S2")
    table = table[["symbol", *OLFACTORY, "respiratory"]]
    for column in [*OLFACTORY, "respiratory"]:
        table[column] = pd.to_numeric(table[column], errors="coerce")
    table = table.groupby("symbol", as_index=False).max()
    joined = genes.merge(table, on="symbol", how="left")
    joined["median"] = joined[OLFACTORY].median(axis=1)
    joined["maximum"] = joined[OLFACTORY].max(axis=1)
    joined["n"] = (joined[OLFACTORY] >= RULES["olfactory_min_fpkm"]).sum(axis=1)
    joined["measured"] = joined[OLFACTORY[0]].notna()
    return joined.groupby("entity_id").agg(
        olfactory_measured=("measured", "any"),
        olfactory_fpkm_median=("median", "max"),
        olfactory_fpkm_max=("maximum", "max"),
        olfactory_samples_ge1=("n", "max"),
        nasal_respiratory_fpkm=("respiratory", "max"),
    )


def cells(genes, folder):
    counts = pd.read_parquet(folder / "counts.parquet")
    pops = pd.read_parquet(folder / "populations.parquet")
    pops["site"] = pops.tissue_general.replace({"skin of body": "skin"})
    pops = pops[~pops.cell_type.str.contains(RULES["excluded_cell_types"], case=False)]
    key = ["site", "dataset_id", "cell_type"]
    groups = (
        pops.groupby(key)
        .agg(cells=("cells", "sum"), donors=("donor_id", "nunique"), depth=("raw_sum", "sum"))
        .reset_index()
    )
    groups = groups[groups.cells >= RULES["cell_min_cells"]]
    measured = (
        counts.merge(pops[["population_id", *key]], on="population_id")
        .groupby(["gene_id", *key])
        .agg(positive=("n_positive", "sum"), multi=("n_ge2", "sum"), total=("count_sum", "sum"))
        .reset_index()
        .merge(groups, on=key)
        .merge(genes, on="gene_id")
    )
    measured["pct"] = 100 * measured.positive / measured.cells
    measured["multi_pct"] = 100 * measured.multi / measured.cells
    measured["cpm"] = 1e6 * measured.total / measured.depth
    measured["detected"] = measured.pct.ge(RULES["cell_min_detected_pct"]) & measured.multi_pct.ge(
        RULES["cell_min_multi_molecule_pct"]
    )
    out = {}
    for site in ["nose", "skin", "scalp", "tongue"]:
        rows = measured[measured.site.eq(site)]
        top = rows.sort_values(["detected", "pct"], ascending=False).groupby("entity_id").head(1)
        top = top.set_index("entity_id")
        out[f"{site}_datasets_detected"] = (
            rows[rows.detected].groupby("entity_id").dataset_id.nunique()
        )
        out[f"{site}_top_cell_type"] = top.cell_type
        out[f"{site}_top_pct"] = top.pct.round(2)
        out[f"{site}_top_cpm"] = top.cpm.round(2)
    summary = groups.groupby("site").agg(
        datasets=("dataset_id", "nunique"),
        cell_types=("cell_type", "nunique"),
        cells=("cells", "sum"),
    )
    donors = pops.groupby("site").donor_id.nunique()
    return pd.DataFrame(out), summary.assign(donors=donors).reset_index().to_dict("records")


def main():
    folder = Path(sys.argv[1])
    genes = candidates()
    table = folder / "olender-2016-table-s2.xlsx"
    frame = (
        genes.groupby("entity_id")
        .agg(protein=("protein", "first"), symbols=("symbol", lambda s: "; ".join(sorted(set(s)))))
        .join(olfactory(genes, table))
    )
    cell_columns, coverage = cells(genes, folder)
    frame = frame.join(cell_columns)
    for column in [c for c in frame if c.endswith("_datasets_detected")]:
        frame[column] = frame[column].fillna(0).astype(int)
    grade = lambda v: target_curation.rna_grade(v) if pd.notna(v) else ""  # noqa: E731
    frame["olfactory_grade"] = frame.olfactory_fpkm_median.map(grade)
    for site in ["nose", "skin", "scalp", "tongue"]:
        frame[f"{site}_top_grade"] = frame[f"{site}_top_cpm"].map(grade)
    strong = {"medium", "high"}
    tissues, basis = [], []
    for r in frame.itertuples():
        found, why = [], []
        if (
            r.olfactory_samples_ge1 >= RULES["olfactory_min_samples"]
            and r.olfactory_grade in strong
        ):
            found.append("olfactory epithelium")
            why.append(
                f"olfactory mucosa bulk RNA: {r.olfactory_samples_ge1} of 4 people at 1 FPKM "
                f"or more, median {r.olfactory_fpkm_median:.1f} ({r.olfactory_grade})"
            )
        for site, tissue in SITES.items():
            n = getattr(r, f"{site}_datasets_detected")
            if n >= RULES["cell_min_datasets"] and getattr(r, f"{site}_top_grade") in strong:
                found.append(tissue)
                why.append(
                    f"{site} single cells: {getattr(r, f'{site}_top_pct'):.0f}% of "
                    f"{getattr(r, f'{site}_top_cell_type')} in the top dataset, "
                    f"detected in {n} datasets"
                )
        tissues.append("; ".join(found))
        basis.append("; ".join(why))
    frame["recorded_tissue"], frame["recorded_basis"] = tissues, basis
    frame = frame.reset_index().sort_values("protein")
    frame.to_csv(DATA / "nose-skin-reference.csv", index=False, float_format="%.4g")
    provenance = {
        "schema": "oncoref.nose-skin-reference",
        "schema_version": 1,
        "scope": "Target candidates of the catalog this was built against; one row each.",
        "catalog_sha256": digest(DATA / "target-selection.csv.gz"),
        "sources": {
            "olfactory_bulk": {
                "citation": "Olender T, et al. The human olfactory transcriptome. "
                "BMC Genomics 2016;17:619.",
                "doi": "10.1186/s12864-016-2960-3",
                "url": "https://static-content.springer.com/esm/art%3A10.1186%2Fs12864-016-2960-3/"
                "MediaObjects/12864_2016_2960_MOESM2_ESM.xlsx",
                "table": "Table S2, columns OE7, OE12, OE15, OE17 (olfactory mucosa) and "
                "respiratory (one nasal respiratory sample); FPKM",
                "sha256": digest(table),
                "matched_proteins": int(frame.olfactory_measured.sum()),
            },
            "cells": {
                "source": "CELLxGENE Census 2025-11-08, primary normal tissue, adult, UMI protocols",
                "extraction": "scripts/extract_nose_skin_cell_counts.py",
                "counts_sha256": digest(folder / "counts.parquet"),
                "populations_sha256": digest(folder / "populations.parquet"),
                "coverage": coverage,
            },
        },
        "rules": RULES,
        "grades": {
            "calibration": "oncoref/data/rna-ihc-calibration.json",
            "floors": target_curation.rna_ihc_calibration()["grade_floors_ntpm"],
            "note": "A tissue is recorded only at grade medium or high. FPKM and single-cell "
            "counts per million are read on the nTPM scale as an approximation.",
        },
        "olfactory_grades": frame.olfactory_grade.value_counts().to_dict(),
        "recorded": {
            tissue: int(frame.recorded_tissue.str.contains(tissue).sum())
            for tissue in ["olfactory epithelium", "nasal mucosa", "skin"]
        },
        "limits": [
            "Olfactory biopsies mix olfactory and respiratory epithelium. The table is from 2016 "
            "gene symbols; proteins without a matching symbol are unmeasured, not absent.",
            "One nasal respiratory bulk sample exists. Its value is reported and never recorded "
            "on its own.",
            "Census nasal cells are brushings of respiratory epithelium. They hold no olfactory "
            "neurons.",
            "Scalp and tongue have one dataset each, so nothing is recorded from them. The scalp "
            "cells are follicular keratinocytes, not hair-shaft cortex, where hair keratins and "
            "keratin-associated proteins are made.",
            "One skin dataset (e7c5ba7e) shows many testis genes at similar rates in every cell "
            "type. Requiring two datasets keeps that from recording a tissue.",
            "Melanocytes are left out by request, as are cells without a type.",
        ],
    }
    (DATA / "nose-skin-reference.json").write_text(json.dumps(provenance, indent=2) + "\n")
    print(provenance["recorded"], len(frame))
    print(frame[frame.recorded_tissue.ne("")][["protein", "recorded_basis"]].to_string(index=False))


if __name__ == "__main__":
    main()
