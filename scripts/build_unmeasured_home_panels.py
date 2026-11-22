#!/usr/bin/env python3
"""Build quantitative normal-tissue panels and protein adjudications for unmeasured home tissues (#735).

Processes checksum-pinned adult human bulk RNA-seq (linear TPM), donor-resolved single-cell
RNA-seq (UMI), and auxiliary targeted receptor copy-number references across all five
``UNMEASURED_HOME_TISSUES`` (``olfactory epithelium``, ``hair follicle``, ``taste bud``,
``dental pulp``, ``bone``) plus ``nasal mucosa``.

Outputs:
- ``oncoref/data/unmeasured-home-tissue-panels.csv`` (compact dataset/panel provenance & rules)
- ``oncoref/data/unmeasured-home-protein-adjudications.csv`` (all 568 home_tissue protein groups)
- ``oncoref/data/unmeasured-home-gene-expression.csv.gz`` (genome-wide Ensembl 115 gene table,
  externalized on GitHub Release ``reference-unmeasured-homes-20261122``)
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
import re
import sqlite3
import sys
import tarfile
from pathlib import Path

import numpy as np
import pandas as pd
from pyensembl import EnsemblRelease
from scipy.sparse import csr_matrix

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from oncoref import cta_selection_policy, gene_ids, target_selection, tissue_tiers  # noqa: E402

DATA_DIR = ROOT / "oncoref" / "data"
MIN_CELL_UMI = 500

BULK_TIER_CEILINGS_TPM: dict[str, float] = {
    "vital": float(cta_selection_policy.TIER_BASES["vital"]),
    "protected": float(cta_selection_policy.TIER_BASES["protected"]),
    "manageable": float(cta_selection_policy.TIER_BASES["manageable"]),
    "disposable": float(target_selection.DISPOSABLE_CAP),
}

SC_GE2_CEILINGS_PCT: dict[str, tuple[float, float]] = {
    # (replicated_ge2_donors_pct, single_donor_pct)
    "vital": (1.0, 5.0),
    "protected": (2.0, 10.0),
    "manageable": (5.0, 10.0),
    "disposable": (100.0, 100.0),
}


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(4 * 1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def build_ensembl115_gene_index() -> tuple[
    pd.DataFrame, dict[int, str], dict[str, str], dict[str, int]
]:
    """Return Ensembl 115 protein-coding gene universe, Entrez->ENSG, Symbol->ENSG, and exonic bp."""
    genome = EnsemblRelease(115)
    ts = pd.read_csv(DATA_DIR / "target-selection.csv.gz", low_memory=False)
    catalog_genes: set[str] = set()
    for members_json in ts["member_gene_ids"]:
        catalog_genes.update(json.loads(members_json))

    conn = sqlite3.connect(genome.db.local_db_path)
    genes_sql = pd.read_sql_query("SELECT gene_id, gene_name FROM gene", conn)
    exons_sql = pd.read_sql_query(
        "SELECT DISTINCT gene_id, start, end FROM exon ORDER BY gene_id, start, end",
        conn,
    )
    conn.close()

    sym_by_gid: dict[str, str] = {}
    symbol_counts: dict[str, list[str]] = {}
    for gid, gname in zip(genes_sql["gene_id"], genes_sql["gene_name"]):
        sym = str(gname or "")
        sym_by_gid[str(gid)] = sym
        if sym:
            symbol_counts.setdefault(sym.upper(), []).append(str(gid))

    exonic_bp: dict[str, int] = {}
    cur_gid = ""
    cur_start = 0
    cur_end = -1
    cur_total = 0
    for gid, s, e in zip(
        exons_sql["gene_id"].to_numpy(),
        exons_sql["start"].to_numpy(dtype=np.int64),
        exons_sql["end"].to_numpy(dtype=np.int64),
    ):
        if gid != cur_gid:
            if cur_gid:
                exonic_bp[cur_gid] = int(cur_total + (cur_end - cur_start + 1))
            cur_gid = str(gid)
            cur_start = int(s)
            cur_end = int(e)
            cur_total = 0
        elif int(s) <= cur_end + 1:
            if int(e) > cur_end:
                cur_end = int(e)
        else:
            cur_total += cur_end - cur_start + 1
            cur_start = int(s)
            cur_end = int(e)
    if cur_gid:
        exonic_bp[cur_gid] = int(cur_total + (cur_end - cur_start + 1))

    gene_rows = []
    for gid in sorted(catalog_genes):
        sym = sym_by_gid.get(gid) or gene_ids.canonical_gene_symbol(gid) or ""
        gene_rows.append({"gene_id": gid, "symbol": sym, "exonic_bp": exonic_bp.get(gid, 0)})

    symbol_to_ensg: dict[str, str] = {}
    for sym_upper, gids in symbol_counts.items():
        cat_hits = [x for x in gids if x in catalog_genes]
        if len(cat_hits) == 1:
            symbol_to_ensg[sym_upper] = cat_hits[0]
        elif len(gids) == 1:
            symbol_to_ensg[sym_upper] = gids[0]

    entrez_df = gene_ids.entrez_gene_mappings()
    entrez_to_ensg: dict[int, str] = {}
    for row in entrez_df.itertuples(index=False):
        ensg = str(row.canonical_ensembl_gene_id).split(".")[0]
        for raw_eid in (row.entrez_id, row.live_entrez_id):
            try:
                eid = int(raw_eid)
            except (TypeError, ValueError):
                continue
            if eid not in entrez_to_ensg or ensg in catalog_genes:
                entrez_to_ensg[eid] = ensg

    return pd.DataFrame(gene_rows), entrez_to_ensg, symbol_to_ensg, exonic_bp


def map_feature_label(raw_label: str, symbol_to_ensg: dict[str, str]) -> str:
    label = str(raw_label).strip().strip('"')
    if not label:
        return ""
    if label.startswith("ENSG"):
        return label.split(".")[0]
    cleaned = re.sub(r"^(GRCh38|hg38|hg19)[_:-]+", "", label)
    if cleaned.startswith("ENSG"):
        return cleaned.split(".")[0]
    hit = symbol_to_ensg.get(cleaned.upper(), "")
    if hit:
        return hit
    resolved = gene_ids.resolve_symbol(cleaned)
    return resolved if resolved and resolved.startswith("ENSG") else ""


def summarize_bulk_matrix(
    tpm_by_gene: pd.DataFrame, sample_cols: list[str], prefix: str
) -> pd.DataFrame:
    vals = tpm_by_gene[sample_cols].to_numpy(dtype=float)
    return pd.DataFrame(
        {
            "gene_id": tpm_by_gene["gene_id"].values,
            f"{prefix}_measured": True,
            f"{prefix}_p95_tpm": np.percentile(vals, 95, axis=1).round(4),
            f"{prefix}_max_tpm": np.max(vals, axis=1).round(4),
            f"{prefix}_median_tpm": np.median(vals, axis=1).round(4),
            f"{prefix}_donors_ge1_tpm": (vals >= 1.0).sum(axis=1).astype(int),
        }
    )


def load_ncbi_tpm(
    path: Path, sample_cols: list[str], entrez_to_ensg: dict[int, str], prefix: str
) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    df["gene_id"] = pd.to_numeric(df["GeneID"], errors="coerce").map(
        lambda x: entrez_to_ensg.get(int(x), "") if pd.notna(x) else ""
    )
    df = df[df["gene_id"].ne("")].copy()
    for c in sample_cols:
        df[c] = pd.to_numeric(df[c], errors="coerce").fillna(0.0)
    df = df.groupby("gene_id", as_index=False)[sample_cols].sum()
    return summarize_bulk_matrix(df, sample_cols, prefix)


def load_gse305125_dental_pulp_bulk(path: Path, exonic_bp: dict[str, int]) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    df["gene_id"] = df["ID"].astype(str).str.split(".").str[0]
    ctrl_cols = [c for c in df.columns if c.endswith("_ctr")]
    df = df.groupby("gene_id", as_index=False)[ctrl_cols].sum()
    lengths_kb = df["gene_id"].map(lambda g: exonic_bp.get(g, 0) / 1000.0).to_numpy(dtype=float)
    valid = lengths_kb > 0.05
    df = df.loc[valid].reset_index(drop=True)
    lengths_kb = lengths_kb[valid]

    # Convert raw counts -> exon-union TPM within each untreated library
    tpm_cols = {}
    for c in ctrl_cols:
        rpk = pd.to_numeric(df[c], errors="coerce").fillna(0.0).to_numpy(dtype=float) / lengths_kb
        denom = rpk.sum()
        tpm_cols[c] = (1e6 * rpk / denom) if denom > 0 else np.zeros_like(rpk)
    tpm_df = pd.DataFrame({"gene_id": df["gene_id"].values, **tpm_cols})

    # Aggregate paired central (DPC) and dentin-adherent (DAC) untreated libraries per donor (D1..D6)
    donor_cols = []
    for d_idx in range(1, 7):
        pair = [c for c in ctrl_cols if c.startswith(f"D{d_idx}_")]
        col_name = f"Donor_{d_idx}"
        tpm_df[col_name] = tpm_df[pair].max(axis=1)
        donor_cols.append(col_name)
    return summarize_bulk_matrix(tpm_df, donor_cols, "dental_pulp_bulk")


def load_gse344036_bone_bulk(path: Path, symbol_to_ensg: dict[str, str]) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", low_memory=False)
    ctrl_cols = [f"Control-{i}" for i in range(1, 7)]
    for c in ctrl_cols:
        df[c] = pd.to_numeric(df[c], errors="coerce").fillna(0.0)
    # Within-sample FPKM -> linear TPM over all quantified transcripts/genes before mapping filter
    for c in ctrl_cols:
        tot_fpkm = float(df[c].sum())
        if tot_fpkm > 0:
            df[c] = 1e6 * df[c] / tot_fpkm
    df["gene_id"] = df["id"].map(lambda s: map_feature_label(s, symbol_to_ensg))
    df = df[df["gene_id"].ne("")].groupby("gene_id", as_index=False)[ctrl_cols].sum()
    return summarize_bulk_matrix(df, ctrl_cols, "bone_bulk")


def summarize_umi_donor_matrix(
    mat_csr: csr_matrix,
    gene_ids_per_row: list[str],
    donor_id: str,
    study_id: str,
) -> tuple[pd.DataFrame, int]:
    """Filter droplets >= MIN_CELL_UMI, collapse duplicate ENSG rows, return per-gene donor UMI stats."""
    cell_totals = np.asarray(mat_csr.sum(axis=0)).ravel()
    keep_cells = cell_totals >= MIN_CELL_UMI
    n_qc_cells = int(keep_cells.sum())
    if n_qc_cells == 0:
        return pd.DataFrame(), 0
    sub = mat_csr[:, keep_cells].tocsr()
    total_umi = float(cell_totals[keep_cells].sum())

    # Map valid ENSG rows and collapse duplicate ENSG indices via sparse indicator multiplication
    valid_idx = [i for i, g in enumerate(gene_ids_per_row) if g]
    sub = sub[valid_idx, :]
    valid_genes = [gene_ids_per_row[i] for i in valid_idx]
    unique_genes, inv = np.unique(valid_genes, return_inverse=True)
    collapse = csr_matrix(
        (np.ones(len(valid_genes), dtype=np.float32), (inv, np.arange(len(valid_genes)))),
        shape=(len(unique_genes), len(valid_genes)),
    )
    sub = (collapse @ sub).tocsr()

    # Compute n_pos (>=1 UMI), n_ge2 (>=2 UMI), umi_sum per gene in vectorized sparse ops
    umi_sum = np.asarray(sub.sum(axis=1)).ravel()
    n_pos = np.diff(sub.indptr)
    ge2_data = (sub.data >= 2).astype(np.int32)
    ge2_mat = csr_matrix((ge2_data, sub.indices, sub.indptr), shape=sub.shape)
    n_ge2 = np.asarray(ge2_mat.sum(axis=1)).ravel()

    out = pd.DataFrame(
        {
            "gene_id": unique_genes,
            "study_id": study_id,
            "donor_id": donor_id,
            "n_cells": n_qc_cells,
            "n_pos": n_pos,
            "n_ge2": n_ge2,
            "pct_detected": 100.0 * n_pos / n_qc_cells,
            "pct_ge2": 100.0 * n_ge2 / n_qc_cells,
            "cpm": (1e6 * umi_sum / total_umi) if total_umi > 0 else 0.0,
            "umi_sum": umi_sum,
            "total_umi": total_umi,
        }
    )
    return out, n_qc_cells


def read_10x_mtx_pair(
    mtx_bytes: bytes,
    features_bytes: bytes,
    symbol_to_ensg: dict[str, str],
    donor_id: str,
    study_id: str,
) -> tuple[pd.DataFrame, int]:
    feat_txt = gzip.decompress(features_bytes).decode("utf-8", errors="replace")
    row_genes = []
    for ln in feat_txt.splitlines():
        if not ln.strip():
            continue
        parts = ln.split("\t")
        gid = map_feature_label(parts[0], symbol_to_ensg)
        if not gid and len(parts) > 1:
            gid = map_feature_label(parts[1], symbol_to_ensg)
        row_genes.append(gid)
    with gzip.open(io.BytesIO(mtx_bytes), "rt", encoding="ascii") as gz_in:
        dim_line = ""
        for line in gz_in:
            if not line.startswith("%"):
                dim_line = line
                break
        n_rows, n_cols, _nnz = [int(x) for x in dim_line.split()]
        triplets = pd.read_csv(
            gz_in,
            sep=r"\s+",
            header=None,
            names=["r", "c", "v"],
            dtype={"r": np.int32, "c": np.int32, "v": np.float32},
            engine="c",
        )
    mat = csr_matrix(
        (
            triplets["v"].to_numpy(),
            (triplets["r"].to_numpy() - 1, triplets["c"].to_numpy() - 1),
        ),
        shape=(n_rows, n_cols),
    )
    return summarize_umi_donor_matrix(mat, row_genes, donor_id, study_id)


def read_dense_dge_gz(
    dge_bytes: bytes,
    symbol_to_ensg: dict[str, str],
    donor_id: str,
    study_id: str,
) -> tuple[pd.DataFrame, int]:
    with gzip.open(io.BytesIO(dge_bytes), "rt", encoding="utf-8", errors="replace") as fh:
        df = pd.read_csv(fh, sep="\t", index_col=0)
    row_genes = [map_feature_label(idx, symbol_to_ensg) for idx in df.index]
    mat = csr_matrix(df.to_numpy(dtype=np.float32))
    return summarize_umi_donor_matrix(mat, row_genes, donor_id, study_id)


def aggregate_sc_studies(donor_frames: list[pd.DataFrame], prefix: str) -> pd.DataFrame:
    cat = pd.concat(donor_frames, ignore_index=True)
    # Multi-molecule donor call: >=2 multi-molecule cells in that donor
    cat["is_donor_ge2"] = cat["n_ge2"] >= 2
    cat_ge2 = cat[cat["is_donor_ge2"]]
    studies_ge2 = cat_ge2.groupby("gene_id")["study_id"].nunique()
    donors_ge2 = cat_ge2.groupby("gene_id")["donor_id"].nunique()

    agg = (
        cat.groupby("gene_id", as_index=False)
        .agg(
            max_pct_detected=("pct_detected", "max"),
            max_pct_ge2=("pct_ge2", "max"),
            umi_sum=("umi_sum", "sum"),
        )
        .copy()
    )
    total_depth = float(
        cat[["study_id", "donor_id", "total_umi"]].drop_duplicates()["total_umi"].sum()
    )
    agg[f"{prefix}_measured"] = True
    agg[f"{prefix}_studies_ge2"] = agg["gene_id"].map(studies_ge2).fillna(0).astype(int)
    agg[f"{prefix}_donors_ge2"] = agg["gene_id"].map(donors_ge2).fillna(0).astype(int)
    agg[f"{prefix}_max_pct_detected"] = agg["max_pct_detected"].round(4)
    agg[f"{prefix}_max_pct_ge2"] = agg["max_pct_ge2"].round(4)
    agg[f"{prefix}_pooled_cpm"] = (
        (1e6 * agg["umi_sum"] / total_depth).round(4) if total_depth > 0 else 0.0
    )
    return agg[
        [
            "gene_id",
            f"{prefix}_measured",
            f"{prefix}_studies_ge2",
            f"{prefix}_donors_ge2",
            f"{prefix}_max_pct_detected",
            f"{prefix}_max_pct_ge2",
            f"{prefix}_pooled_cpm",
        ]
    ]


def load_cellxgene_nose_umi(nose_skin_dir: Path) -> tuple[pd.DataFrame, int, int]:
    counts = pd.read_parquet(nose_skin_dir / "counts.parquet")
    pops = pd.read_parquet(nose_skin_dir / "populations.parquet")
    nose_pops = pops[
        pops["tissue_general"].eq("nose")
        & ~pops["cell_type"].str.contains("melanocyte|unknown", case=False, na=False)
    ].copy()
    donor_pop = (
        counts.merge(
            nose_pops[["population_id", "dataset_id", "donor_id", "cells", "raw_sum"]],
            on="population_id",
        )
        .groupby(["gene_id", "dataset_id", "donor_id"], as_index=False)
        .agg(
            n_cells=("cells", "sum"),
            n_pos=("n_positive", "sum"),
            n_ge2=("n_ge2", "sum"),
            umi_sum=("count_sum", "sum"),
        )
    )
    pop_totals = nose_pops.groupby(["dataset_id", "donor_id"], as_index=False).agg(
        total_cells=("cells", "sum"), total_umi=("raw_sum", "sum")
    )
    pop_totals = pop_totals[pop_totals["total_cells"] >= 50]
    donor_pop = donor_pop.merge(pop_totals, on=["dataset_id", "donor_id"])
    donor_pop["study_id"] = donor_pop["dataset_id"]
    donor_pop["pct_detected"] = 100.0 * donor_pop["n_pos"] / donor_pop["total_cells"]
    donor_pop["pct_ge2"] = 100.0 * donor_pop["n_ge2"] / donor_pop["total_cells"]
    agg = aggregate_sc_studies([donor_pop], "nasal_sc")
    return (
        agg,
        int(pop_totals["donor_id"].nunique()),
        int(pop_totals["total_cells"].sum()),
    )


def load_verbeurgt_2014(path: Path, symbol_to_ensg: dict[str, str]) -> pd.DataFrame:
    raw = pd.read_excel(path, sheet_name="Copies number")
    # Row index 3 holds "Well | Detector Name", data starts at row index 4
    data = raw.iloc[4:].copy()
    detectors = data.iloc[:, 1].astype(str)
    sample_cols = list(data.columns[2:28])
    for c in sample_cols:
        data[c] = pd.to_numeric(data[c], errors="coerce").fillna(0.0)
    symbols = detectors.str.split("-").str[0].str.strip()
    data["gene_id"] = symbols.map(lambda s: map_feature_label(s, symbol_to_ensg))
    data = data[data["gene_id"].ne("")].groupby("gene_id", as_index=False)[sample_cols].max()
    vals = data[sample_cols].to_numpy(dtype=float)
    return pd.DataFrame(
        {
            "gene_id": data["gene_id"].values,
            "verbeurgt_measured": True,
            "verbeurgt_median_copies": np.median(vals, axis=1).round(2),
            "verbeurgt_max_copies": np.max(vals, axis=1).round(2),
            "verbeurgt_donors_ge10_copies": (vals >= 10.0).sum(axis=1).astype(int),
        }
    )


def load_olender_respiratory_fpkm(path: Path, symbol_to_ensg: dict[str, str]) -> pd.DataFrame:
    tbl = pd.read_excel(path, sheet_name="Table S2")[["symbol", "respiratory"]].copy()
    tbl["respiratory"] = pd.to_numeric(tbl["respiratory"], errors="coerce")
    tbl = tbl.dropna(subset=["respiratory"])
    tbl["gene_id"] = tbl["symbol"].map(lambda s: map_feature_label(s, symbol_to_ensg))
    tbl = tbl[tbl["gene_id"].ne("")].groupby("gene_id", as_index=False)["respiratory"].max()
    return tbl.rename(columns={"respiratory": "nasal_respiratory_fpkm"})


def adjudicate_row(row: pd.Series) -> str:
    if not bool(row["is_ca_cta_home_gated"]):
        return "already_excluded_upstream"

    tissue = str(row["home_tissue"])
    tier = str(row["home_tissue_tier"])
    bulk_ceil = BULK_TIER_CEILINGS_TPM[tier]
    repl_sc_ceil, single_sc_ceil = SC_GE2_CEILINGS_PCT[tier]

    bulk_meas = bool(row["home_bulk_measured"])
    bulk_p95 = float(row["home_bulk_p95_tpm"]) if pd.notna(row["home_bulk_p95_tpm"]) else 0.0
    bulk_max = float(row["home_bulk_max_tpm"]) if pd.notna(row["home_bulk_max_tpm"]) else 0.0

    sc_donors_ge2 = int(row["home_sc_donors_ge2_umi"])
    sc_max_pct_ge2 = (
        float(row["home_sc_max_donor_pct_ge2"])
        if pd.notna(row["home_sc_max_donor_pct_ge2"])
        else 0.0
    )

    if tissue == "hair follicle":
        # Disposable compartment: somatic skin bulk already passed in target-selection
        if sc_donors_ge2 >= 1:
            return "home_disposable_follicle_expressed_tenfold"
        return "home_undetected_in_sampled_panel"

    # Non-disposable home tissue or cross-panel non-disposable violation:
    if str(row.get("cross_panel_ceiling_violation", "")).strip():
        return "home_expressed_above_tier_ceiling"
    sc_exceeds = (sc_donors_ge2 >= 2 and sc_max_pct_ge2 >= repl_sc_ceil) or (
        sc_donors_ge2 == 1 and sc_max_pct_ge2 >= single_sc_ceil
    )
    if (bulk_meas and bulk_p95 > bulk_ceil) or sc_exceeds:
        return "home_expressed_above_tier_ceiling"

    norm_ref = float(row["normal_reference_tpm"]) if pd.notna(row["normal_reference_tpm"]) else 0.0
    tum_max = float(row["tumor_max_tpm"]) if pd.notna(row["tumor_max_tpm"]) else 0.0
    combined_baseline = max(norm_ref, bulk_p95 if bulk_meas else 0.0, 1.0)
    if tum_max < 10.0 * combined_baseline:
        return "home_expressed_insufficient_tumor_fold"

    aux_pos = "copies>=" in str(row["auxiliary_assay_summary"]) and not str(
        row["auxiliary_assay_summary"]
    ).endswith("0/26 donors >=10 copies")
    if (bulk_meas and bulk_max >= 0.5) or sc_donors_ge2 >= 1 or aux_pos:
        return "home_cleared_panel_passes_tenfold"

    return "home_undetected_in_sampled_panel"


def build_all(raw_dir: Path, nose_skin_dir: Path, durante_xlsx: Path) -> dict:
    print("1. Indexing Ensembl 115 genes & exonic union lengths...")
    genes_df, entrez_to_ensg, symbol_to_ensg, exonic_bp = build_ensembl115_gene_index()

    print(
        "2. Loading linear bulk RNA-seq TPM panels (GSE80249, GSE221229, GSE305125, GSE344036)..."
    )
    oe_bulk_path = raw_dir / "GSE80249_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz"
    oe_bulk = load_ncbi_tpm(
        oe_bulk_path,
        ["GSM2122728", "GSM2122729", "GSM2122730", "GSM2122731"],
        entrez_to_ensg,
        "olfactory_bulk",
    )

    taste_bulk_path = raw_dir / "GSE221229_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz"
    taste_cols = [
        "GSM6856456",
        "GSM6856457",
        "GSM6856458",
        "GSM6856459",
        "GSM6856460",
        "GSM6856461",
        "GSM6856463",
    ]
    taste_bulk = load_ncbi_tpm(taste_bulk_path, taste_cols, entrez_to_ensg, "taste_bulk")

    pulp_bulk_path = raw_dir / "GSE305125_counts_raw.txt.gz"
    pulp_bulk = load_gse305125_dental_pulp_bulk(pulp_bulk_path, exonic_bp)

    bone_fpkm_path = raw_dir / "GSE344036_FPKM_anno.xls.gz"
    bone_bulk = load_gse344036_bone_bulk(bone_fpkm_path, symbol_to_ensg)

    print("3. Loading donor-resolved single-cell UMI panels...")
    # 3a. Olfactory neuroepithelium (GSE139522 Patient 2 & Patient 3)
    oe_sc_frames = []
    oe_sc_cells = 0
    for gsm, p_lbl in [("GSM4142871", "Patient2"), ("GSM4142872", "Patient3")]:
        f_df, n_c = read_10x_mtx_pair(
            (raw_dir / f"{gsm}_{p_lbl}_matrix.mtx.gz").read_bytes(),
            (raw_dir / f"{gsm}_{p_lbl}_features.tsv.gz").read_bytes(),
            symbol_to_ensg,
            donor_id=f"{gsm}_{p_lbl}",
            study_id="GSE139522",
        )
        oe_sc_frames.append(f_df)
        oe_sc_cells += n_c
    oe_sc = aggregate_sc_studies(oe_sc_frames, "olfactory_sc")

    # 3b. Nasal respiratory mucosa (CELLxGENE Census 2025-11-08)
    nasal_sc, nasal_donors, nasal_cells = load_cellxgene_nose_umi(nose_skin_dir)

    # 3c. Hair follicle: GSE129611 (5 donors) + GSE193269 (3 normal scalp hair follicle donors)
    hf_frames = []
    hf_129611_cells = 0
    with tarfile.open(raw_dir / "GSE129611_RAW.tar") as tf:
        for member in tf.getmembers():
            raw_bytes = tf.extractfile(member).read()
            donor_tag = member.name.split("_")[0]
            f_df, n_c = read_dense_dge_gz(raw_bytes, symbol_to_ensg, donor_tag, "GSE129611")
            hf_frames.append(f_df)
            hf_129611_cells += n_c

    hf_193269_cells = 0
    for gsm, hf_lbl in [("GSM6138434", "HF1"), ("GSM6138435", "HF2"), ("GSM6138436", "HF3")]:
        f_df, n_c = read_10x_mtx_pair(
            (raw_dir / f"{gsm}_{hf_lbl}_matrix.mtx.gz").read_bytes(),
            (raw_dir / f"{gsm}_{hf_lbl}_features.tsv.gz").read_bytes(),
            symbol_to_ensg,
            donor_id=f"{gsm}_{hf_lbl}",
            study_id="GSE193269",
        )
        hf_frames.append(f_df)
        hf_193269_cells += n_c
    hf_sc = aggregate_sc_studies(hf_frames, "hair_sc")

    # 3d. Dental pulp: GSE164157 (Pulp1..Pulp5)
    pulp_sc_frames = []
    pulp_sc_cells = 0
    with tarfile.open(raw_dir / "GSE164157_RAW.tar") as tf:
        by_name = {m.name: tf.extractfile(m).read() for m in tf.getmembers()}
        for idx, gsm in enumerate(
            ["GSM4998457", "GSM4998458", "GSM4998459", "GSM4998460", "GSM4998461"], start=1
        ):
            f_df, n_c = read_10x_mtx_pair(
                by_name[f"{gsm}_Pulp{idx}_matrix.mtx.gz"],
                by_name[f"{gsm}_Pulp{idx}_genes.tsv.gz"],
                symbol_to_ensg,
                donor_id=f"{gsm}_Pulp{idx}",
                study_id="GSE164157",
            )
            pulp_sc_frames.append(f_df)
            pulp_sc_cells += n_c
    pulp_sc = aggregate_sc_studies(pulp_sc_frames, "dental_pulp_sc")

    # 3e. Bone primary osteoblasts: GSE147390 (GSM4429895)
    bone_sc_frames = []
    bone_sc_cells = 0
    with tarfile.open(raw_dir / "GSE147390_RAW.tar") as tf:
        by_name = {m.name: tf.extractfile(m).read() for m in tf.getmembers()}
        f_df, n_c = read_10x_mtx_pair(
            by_name["GSM4429895_matrix.mtx.gz"],
            by_name["GSM4429895_features.tsv.gz"],
            symbol_to_ensg,
            donor_id="GSM4429895",
            study_id="GSE147390",
        )
        bone_sc_frames.append(f_df)
        bone_sc_cells += n_c
    bone_sc = aggregate_sc_studies(bone_sc_frames, "bone_sc")

    print("4. Loading auxiliary receptor qPCR / respiratory references...")
    verb_path = raw_dir / "verbeurgt_2014_s004_copies.xlsx"
    verb_df = load_verbeurgt_2014(verb_path, symbol_to_ensg)
    olender_xlsx = nose_skin_dir / "olender-2016-table-s2.xlsx"
    olender_resp = load_olender_respiratory_fpkm(olender_xlsx, symbol_to_ensg)

    print("5. Assembling genome-wide Ensembl 115 home-tissue expression matrix...")
    gw = genes_df.copy()
    for part in (
        oe_bulk,
        oe_sc,
        verb_df,
        nasal_sc,
        olender_resp,
        hf_sc,
        taste_bulk,
        pulp_bulk,
        pulp_sc,
        bone_bulk,
        bone_sc,
    ):
        gw = gw.merge(part, on="gene_id", how="left")

    for col in [c for c in gw.columns if c.endswith("_measured")]:
        gw[col] = gw[col].fillna(False).astype(bool)
    for col in [c for c in gw.columns if "_donors_" in c or "_studies_" in c]:
        gw[col] = gw[col].fillna(0).astype(int)

    gw_path = DATA_DIR / "unmeasured-home-gene-expression.csv.gz"
    gw.sort_values("gene_id").to_csv(gw_path, index=False, compression="gzip")

    print("6. Building panel provenance table & 568-protein adjudication table...")
    panels_rows = [
        {
            "panel_id": "gse80249_olfactory_bulk_tpm",
            "home_tissue": "olfactory epithelium",
            "tissue_tier": tissue_tiers.assess("olfactory epithelium")["tier"],
            "modality": "bulk_rnaseq_tpm",
            "accession": "GSE80249",
            "citation": "Olender et al. 2016 (BMC Genomics)",
            "pmid": "27515280",
            "url": "https://www.ncbi.nlm.nih.gov/geo/download/?type=rnaseq_counts&acc=GSE80249&format=file&file=GSE80249_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz",
            "sha256": sha256_file(oe_bulk_path),
            "n_donors": 4,
            "n_libraries_or_populations": 4,
            "n_qc_cells": 0,
            "native_unit": "TPM",
            "genes_mapped_ensembl115": int(gw["olfactory_bulk_measured"].sum()),
            "tier_ceiling_tpm": BULK_TIER_CEILINGS_TPM["manageable"],
            "within_assay_rule": "Manageable bulk P95 <= 16 TPM and tumor_max_tpm >= 10 x max(normal_reference_tpm, home_bulk_p95_tpm, 1.0)",
            "limitations": "Adult olfactory mucosal biopsies (ages 58-65, 2F/2M) admix sensory neuroepithelium with respiratory epithelium; monogenic OR expression is sparse in bulk.",
        },
        {
            "panel_id": "gse139522_olfactory_scrnaseq_umi",
            "home_tissue": "olfactory epithelium",
            "tissue_tier": tissue_tiers.assess("olfactory epithelium")["tier"],
            "modality": "scrnaseq_umi",
            "accession": "GSE139522",
            "citation": "Durante et al. 2020 (Nature Neuroscience)",
            "pmid": "32066986",
            "url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE139522",
            "sha256": sha256_file(raw_dir / "GSM4142871_Patient2_matrix.mtx.gz"),
            "n_donors": 2,
            "n_libraries_or_populations": 2,
            "n_qc_cells": oe_sc_cells,
            "native_unit": "UMI",
            "genes_mapped_ensembl115": int(gw["olfactory_sc_measured"].sum()),
            "tier_ceiling_tpm": np.nan,
            "within_assay_rule": "Evaluated on native UMI multi-molecule donor detection (never pooled into TPM); >=5% replicated or >=10% single-donor multi-molecule cells exceeds Manageable ceiling",
            "limitations": "Each mature olfactory sensory neuron expresses one allele of ~1 of ~400 OR genes (~0.2-0.5% of OSNs per receptor), requiring large OSN counts for exhaustive OR capture.",
        },
        {
            "panel_id": "verbeurgt_2014_olfactory_qpcr_copies",
            "home_tissue": "olfactory epithelium",
            "tissue_tier": tissue_tiers.assess("olfactory epithelium")["tier"],
            "modality": "targeted_qpcr_copies",
            "accession": "doi:10.1371/journal.pone.0096333.s004",
            "citation": "Verbeurgt et al. 2014 (PLoS ONE)",
            "pmid": "24800820",
            "url": "https://journals.plos.org/plosone/article/file?id=10.1371/journal.pone.0096333.s004&type=supplementary",
            "sha256": sha256_file(verb_path),
            "n_donors": 26,
            "n_libraries_or_populations": 26,
            "n_qc_cells": 0,
            "native_unit": "TaqMan_copies",
            "genes_mapped_ensembl115": int(gw["verbeurgt_measured"].sum()),
            "tier_ceiling_tpm": np.nan,
            "within_assay_rule": "Targeted TaqMan OR receptor presence across 26 adult human whole-olfactory-mucosa donors; auxiliary receptor detection only, never converted to TPM",
            "limitations": "Targeted to olfactory receptor (OR) gene family; TaqMan copy numbers are not transcriptome-normalized TPM.",
        },
        {
            "panel_id": "cellxgene_adult_nose_umi",
            "home_tissue": "nasal mucosa",
            "tissue_tier": tissue_tiers.assess("nasal mucosa")["tier"],
            "modality": "scrnaseq_umi",
            "accession": "CELLxGENE-2025-11-08-nose",
            "citation": "CELLxGENE Census 2025-11-08 adult nasal mucosa + Olender 2016 respiratory",
            "pmid": "27515280",
            "url": "https://chanzuckerberg.github.io/cellxgene-census/",
            "sha256": sha256_file(nose_skin_dir / "counts.parquet"),
            "n_donors": nasal_donors,
            "n_libraries_or_populations": 405,
            "n_qc_cells": nasal_cells,
            "native_unit": "UMI",
            "genes_mapped_ensembl115": int(gw["nasal_sc_measured"].sum()),
            "tier_ceiling_tpm": np.nan,
            "within_assay_rule": "Protected respiratory mucosa single-cell UMI multi-molecule screen (>=2% replicated multi-molecule cells exceeds ceiling)",
            "limitations": "Nasal brushings sample respiratory epithelium rather than olfactory cleft neuroepithelium.",
        },
        {
            "panel_id": "gse129611_hair_follicle_scrnaseq_umi",
            "home_tissue": "hair follicle",
            "tissue_tier": tissue_tiers.assess("hair follicle")["tier"],
            "modality": "scrnaseq_umi",
            "accession": "GSE129611",
            "citation": "Takahashi et al. 2020 (J Invest Dermatol)",
            "pmid": "31676413",
            "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE129nnn/GSE129611/suppl/GSE129611_RAW.tar",
            "sha256": sha256_file(raw_dir / "GSE129611_RAW.tar"),
            "n_donors": 5,
            "n_libraries_or_populations": 5,
            "n_qc_cells": hf_129611_cells,
            "native_unit": "UMI",
            "genes_mapped_ensembl115": int(gw["hair_sc_measured"].sum()),
            "tier_ceiling_tpm": BULK_TIER_CEILINGS_TPM["disposable"],
            "within_assay_rule": "Disposable hair-follicle single-cell UMI confirmation across 5 anagen scalp donors (Drop-seq + 10x); somatic non-follicular skin bulk already gated in catalog",
            "limitations": "Enzymatic single-cell dissociation under-samples terminally cornified hair-shaft cortex cells where ultra-high-sulfur KRTAPs peak.",
        },
        {
            "panel_id": "gse193269_hair_follicle_scrnaseq_umi",
            "home_tissue": "hair follicle",
            "tissue_tier": tissue_tiers.assess("hair follicle")["tier"],
            "modality": "scrnaseq_umi",
            "accession": "GSE193269",
            "citation": "Shim et al. 2022 (J Invest Dermatol)",
            "pmid": "35853486",
            "url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE193269",
            "sha256": sha256_file(raw_dir / "GSM6138434_HF1_matrix.mtx.gz"),
            "n_donors": 3,
            "n_libraries_or_populations": 3,
            "n_qc_cells": hf_193269_cells,
            "native_unit": "UMI",
            "genes_mapped_ensembl115": int(gw["hair_sc_measured"].sum()),
            "tier_ceiling_tpm": BULK_TIER_CEILINGS_TPM["disposable"],
            "within_assay_rule": "Independent 3-donor normal adult scalp hair follicle 10x UMI replication (HF1-HF3) paired with GSE129611",
            "limitations": "Single-cell UMI captures matrix, cuticle, IRS and ORS keratinocytes; cornified shaft interior remains partially depleted.",
        },
        {
            "panel_id": "gse221229_fungiform_papillae_bulk_tpm",
            "home_tissue": "taste bud",
            "tissue_tier": tissue_tiers.assess("taste bud")["tier"],
            "modality": "bulk_rnaseq_tpm",
            "accession": "GSE221229",
            "citation": "aye/GEO GSE221229 Human Fungiform Papilla Transcriptome",
            "pmid": "",
            "url": "https://www.ncbi.nlm.nih.gov/geo/download/?type=rnaseq_counts&acc=GSE221229&format=file&file=GSE221229_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz",
            "sha256": sha256_file(taste_bulk_path),
            "n_donors": 7,
            "n_libraries_or_populations": 7,
            "n_qc_cells": 0,
            "native_unit": "TPM",
            "genes_mapped_ensembl115": int(gw["taste_bulk_measured"].sum()),
            "tier_ceiling_tpm": BULK_TIER_CEILINGS_TPM["protected"],
            "within_assay_rule": "Protected bulk P95 <= 8 TPM and tumor_max_tpm >= 10 x max(normal_reference_tpm, home_bulk_p95_tpm, 1.0)",
            "limitations": "Whole fungiform papilla biopsies dilute rare Type II bitter-receptor (TAS2R) taste receptor cells within surrounding lingual epithelium.",
        },
        {
            "panel_id": "gse305125_dental_pulp_bulk_tpm",
            "home_tissue": "dental pulp",
            "tissue_tier": tissue_tiers.assess("dental pulp")["tier"],
            "modality": "bulk_rnaseq_tpm",
            "accession": "GSE305125",
            "citation": "GEO GSE305125 Peripheral and Central Human Dental Pulp Bulk RNA-Seq",
            "pmid": "",
            "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE305nnn/GSE305125/suppl/GSE305125_counts_raw.txt.gz",
            "sha256": sha256_file(pulp_bulk_path),
            "n_donors": 6,
            "n_libraries_or_populations": 12,
            "n_qc_cells": 0,
            "native_unit": "TPM",
            "genes_mapped_ensembl115": int(gw["dental_pulp_bulk_measured"].sum()),
            "tier_ceiling_tpm": BULK_TIER_CEILINGS_TPM["manageable"],
            "within_assay_rule": "Manageable bulk donor P95 <= 16 TPM across untreated central & dentin-adherent odontoblast pulp (6 donors) and 10x tumor-over-home fold",
            "limitations": "Ex-vivo untreated control pulp explants (6h control arm paired with central and dentin-adherent layers); paired with in-vivo GSE164157 scRNA-seq.",
        },
        {
            "panel_id": "gse164157_dental_pulp_scrnaseq_umi",
            "home_tissue": "dental pulp",
            "tissue_tier": tissue_tiers.assess("dental pulp")["tier"],
            "modality": "scrnaseq_umi",
            "accession": "GSE164157",
            "citation": "Pagella et al. 2021 (iScience)",
            "pmid": "33997688",
            "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE164nnn/GSE164157/suppl/GSE164157_RAW.tar",
            "sha256": sha256_file(raw_dir / "GSE164157_RAW.tar"),
            "n_donors": 5,
            "n_libraries_or_populations": 5,
            "n_qc_cells": pulp_sc_cells,
            "native_unit": "UMI",
            "genes_mapped_ensembl115": int(gw["dental_pulp_sc_measured"].sum()),
            "tier_ceiling_tpm": np.nan,
            "within_assay_rule": "In-vivo healthy adult molar dental pulp 10x UMI across 5 donors (Pulp1-Pulp5); >=5% replicated multi-molecule cells exceeds Manageable ceiling",
            "limitations": "Mineralized dentin-embedded odontoblast processes can be partially sheared during pulp extraction.",
        },
        {
            "panel_id": "gse344036_mandible_bone_bulk_tpm",
            "home_tissue": "bone",
            "tissue_tier": tissue_tiers.assess("bone")["tier"],
            "modality": "bulk_rnaseq_tpm",
            "accession": "GSE344036",
            "citation": "GEO GSE344036 Healthy Adult Human Mandible Bone Controls",
            "pmid": "",
            "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE344nnn/GSE344036/suppl/GSE344036_FPKM_anno.xls.gz",
            "sha256": sha256_file(bone_fpkm_path),
            "n_donors": 6,
            "n_libraries_or_populations": 6,
            "n_qc_cells": 0,
            "native_unit": "TPM",
            "genes_mapped_ensembl115": int(gw["bone_bulk_measured"].sum()),
            "tier_ceiling_tpm": BULK_TIER_CEILINGS_TPM["vital"],
            "within_assay_rule": "Vital bulk donor P95 <= 4 TPM across 6 healthy adult bone controls (Control-1..6) and 10x tumor-over-home fold",
            "limitations": "Mandibular cortical/trabecular bone biopsies; paired with primary femur ALPL+/CD45- osteoblast scRNA-seq (GSE147390).",
        },
        {
            "panel_id": "gse147390_primary_osteoblast_scrnaseq_umi",
            "home_tissue": "bone",
            "tissue_tier": tissue_tiers.assess("bone")["tier"],
            "modality": "scrnaseq_umi",
            "accession": "GSE147390",
            "citation": "Gong et al. 2021 (Aging)",
            "pmid": "37105556",
            "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE147nnn/GSE147390/suppl/GSE147390_RAW.tar",
            "sha256": sha256_file(raw_dir / "GSE147390_RAW.tar"),
            "n_donors": 1,
            "n_libraries_or_populations": 1,
            "n_qc_cells": bone_sc_cells,
            "native_unit": "UMI",
            "genes_mapped_ensembl115": int(gw["bone_sc_measured"].sum()),
            "tier_ceiling_tpm": np.nan,
            "within_assay_rule": "Vital primary human femur ALPLhigh/CD45low osteoblast 10x UMI (>=5% single-study multi-molecule osteoblast detection exceeds Vital ceiling)",
            "limitations": "Single 31-year-old male femur donor with osteopenia; corroborates bulk healthy bone RNA-seq (GSE344036).",
        },
    ]
    panels_df = pd.DataFrame(panels_rows)
    panels_df.to_csv(DATA_DIR / "unmeasured-home-tissue-panels.csv", index=False)
    return build_protein_adjudications(gw)


def build_protein_adjudications(gw: pd.DataFrame) -> dict:
    """Adjudicate all 568 home_tissue protein groups from the genome-wide panel matrix."""
    gw_path = DATA_DIR / "unmeasured-home-gene-expression.csv.gz"
    ts = pd.read_csv(DATA_DIR / "target-selection.csv.gz", low_memory=False)
    ht = ts[ts["home_tissue"].fillna("").ne("")].copy()
    gw_by_gene = gw.set_index("gene_id")

    tissue_spec = {
        "olfactory epithelium": (
            "gse80249_olfactory_bulk_tpm",
            "olfactory_bulk",
            "gse139522_olfactory_scrnaseq_umi",
            "olfactory_sc",
        ),
        "nasal mucosa": (
            "",
            "",
            "cellxgene_adult_nose_umi",
            "nasal_sc",
        ),
        "hair follicle": (
            "",
            "",
            "gse129611_hair_follicle_scrnaseq_umi;gse193269_hair_follicle_scrnaseq_umi",
            "hair_sc",
        ),
        "taste bud": (
            "gse221229_fungiform_papillae_bulk_tpm",
            "taste_bulk",
            "",
            "",
        ),
        "dental pulp": (
            "gse305125_dental_pulp_bulk_tpm",
            "dental_pulp_bulk",
            "gse164157_dental_pulp_scrnaseq_umi",
            "dental_pulp_sc",
        ),
        "bone": (
            "gse344036_mandible_bone_bulk_tpm",
            "bone_bulk",
            "gse147390_primary_osteoblast_scrnaseq_umi",
            "bone_sc",
        ),
    }

    adj_rows = []
    for r in ht.itertuples(index=False):
        tissue = str(r.home_tissue)
        tier = tissue_tiers.assess(tissue)["tier"]
        members = [g for g in json.loads(r.member_gene_ids) if g in gw_by_gene.index]
        sub = gw_by_gene.loc[members] if members else pd.DataFrame()

        b_id, b_pfx, sc_id, sc_pfx = tissue_spec[tissue]
        b_meas = bool(sub[f"{b_pfx}_measured"].any()) if (b_pfx and not sub.empty) else False
        b_p95 = float(sub[f"{b_pfx}_p95_tpm"].max()) if b_meas else np.nan
        b_max = float(sub[f"{b_pfx}_max_tpm"].max()) if b_meas else np.nan
        b_med = float(sub[f"{b_pfx}_median_tpm"].max()) if b_meas else np.nan
        b_ge1 = int(sub[f"{b_pfx}_donors_ge1_tpm"].max()) if b_meas else 0

        sc_meas = bool(sub[f"{sc_pfx}_measured"].any()) if (sc_pfx and not sub.empty) else False
        sc_studies = int(sub[f"{sc_pfx}_studies_ge2"].max()) if sc_meas else 0
        sc_donors = int(sub[f"{sc_pfx}_donors_ge2"].max()) if sc_meas else 0
        sc_pct_det = float(sub[f"{sc_pfx}_max_pct_detected"].max()) if sc_meas else np.nan
        sc_pct_ge2 = float(sub[f"{sc_pfx}_max_pct_ge2"].max()) if sc_meas else np.nan
        sc_cpm = float(sub[f"{sc_pfx}_pooled_cpm"].max()) if sc_meas else np.nan

        aux = ""
        if (
            tissue == "olfactory epithelium"
            and not sub.empty
            and bool(sub["verbeurgt_measured"].any())
        ):
            v_med = float(sub["verbeurgt_median_copies"].max())
            v_max = float(sub["verbeurgt_max_copies"].max())
            v_ge10 = int(sub["verbeurgt_donors_ge10_copies"].max())
            aux = f"Verbeurgt2014 WHOM qPCR: median={v_med:.1f}, max={v_max:.1f} copies>=10 in {v_ge10}/26 donors >=10 copies"
        elif (
            tissue == "nasal mucosa"
            and not sub.empty
            and sub["nasal_respiratory_fpkm"].notna().any()
        ):
            rfpkm = float(sub["nasal_respiratory_fpkm"].max())
            aux = f"Olender2016 nasal respiratory bulk={rfpkm:.2f} FPKM"

        s10x = int(r.samples_at_tenfold) if pd.notna(r.samples_at_tenfold) else 0
        is_gated = bool(
            r.classification == "Other:Unscreened tissue"
            or (r.classification == "Reproductive Antigen only" and s10x >= 1)
        )
        tum_max = float(r.tumor_max_tpm) if pd.notna(r.tumor_max_tpm) else np.nan
        norm_ref = float(r.normal_reference_tpm) if pd.notna(r.normal_reference_tpm) else 0.0
        fold_vs_home = (
            round(tum_max / max(norm_ref, b_p95, 1.0), 3)
            if (b_meas and pd.notna(tum_max))
            else np.nan
        )

        cross_violations = []
        if not sub.empty:
            for other_tissue, pfx, ceil_val in (
                ("bone", "bone_bulk", BULK_TIER_CEILINGS_TPM["vital"]),
                ("taste bud", "taste_bulk", BULK_TIER_CEILINGS_TPM["protected"]),
                ("olfactory epithelium", "olfactory_bulk", BULK_TIER_CEILINGS_TPM["manageable"]),
                ("dental pulp", "dental_pulp_bulk", BULK_TIER_CEILINGS_TPM["manageable"]),
            ):
                if other_tissue == tissue:
                    continue
                if bool(sub[f"{pfx}_measured"].any()):
                    other_p95 = float(sub[f"{pfx}_p95_tpm"].max())
                    if other_p95 > ceil_val:
                        cross_violations.append(
                            f"{other_tissue}:{other_p95:.2f}>{ceil_val:.1f}_TPM"
                        )

        rec = {
            "entity_id": r.entity_id,
            "protein_group_label": r.protein_group_label,
            "representative_gene_id": r.representative_gene_id,
            "home_tissue": tissue,
            "home_tissue_tier": tier,
            "active_category": r.category,
            "active_classification": r.classification,
            "is_ca_cta_home_gated": is_gated,
            "is_reproductive_antigen": bool(r.is_reproductive_antigen),
            "normal_reference_tpm": round(norm_ref, 4),
            "disposable_max_tpm": (
                round(float(r.disposable_max_tpm), 4) if pd.notna(r.disposable_max_tpm) else np.nan
            ),
            "tumor_max_tpm": round(tum_max, 4) if pd.notna(tum_max) else np.nan,
            "tumor_top_cancer_code": r.tumor_top_cancer_code
            if pd.notna(r.tumor_top_cancer_code)
            else "",
            "samples_at_tenfold": s10x,
            "home_bulk_panel_id": b_id,
            "home_bulk_measured": b_meas,
            "home_bulk_p95_tpm": round(b_p95, 4) if pd.notna(b_p95) else np.nan,
            "home_bulk_max_tpm": round(b_max, 4) if pd.notna(b_max) else np.nan,
            "home_bulk_median_tpm": round(b_med, 4) if pd.notna(b_med) else np.nan,
            "home_bulk_donors_ge1_tpm": b_ge1,
            "home_bulk_tier_ceiling_tpm": BULK_TIER_CEILINGS_TPM[tier] if b_id else np.nan,
            "tumor_over_home_bulk_fold": fold_vs_home,
            "home_sc_panel_ids": sc_id,
            "home_sc_measured": sc_meas,
            "home_sc_studies_ge2_umi": sc_studies,
            "home_sc_donors_ge2_umi": sc_donors,
            "home_sc_max_donor_pct_detected": round(sc_pct_det, 4)
            if pd.notna(sc_pct_det)
            else np.nan,
            "home_sc_max_donor_pct_ge2": round(sc_pct_ge2, 4) if pd.notna(sc_pct_ge2) else np.nan,
            "home_sc_pooled_cpm": round(sc_cpm, 4) if pd.notna(sc_cpm) else np.nan,
            "cross_panel_ceiling_violation": "; ".join(cross_violations),
            "auxiliary_assay_summary": aux,
        }
        rec["home_panel_disposition"] = adjudicate_row(pd.Series(rec))
        adj_rows.append(rec)

    adj_df = pd.DataFrame(adj_rows).sort_values(
        ["is_ca_cta_home_gated", "home_tissue", "protein_group_label"],
        ascending=[False, True, True],
    )
    adj_path = DATA_DIR / "unmeasured-home-protein-adjudications.csv"
    adj_df.to_csv(adj_path, index=False)

    summary = {
        "total_home_proteins": len(adj_df),
        "ca_cta_gated_proteins": int(adj_df["is_ca_cta_home_gated"].sum()),
        "dispositions_all": adj_df["home_panel_disposition"].value_counts().to_dict(),
        "dispositions_gated": (
            adj_df.loc[adj_df["is_ca_cta_home_gated"], "home_panel_disposition"]
            .value_counts()
            .to_dict()
        ),
        "gene_expression_bytes": gw_path.stat().st_size,
        "gene_expression_sha256": sha256_file(gw_path),
    }
    print(json.dumps(summary, indent=2))
    return summary


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--raw-dir",
        type=Path,
        default=None,
        help="Directory containing checksum-verified GEO/supplementary files",
    )
    parser.add_argument(
        "--from-gene-table",
        action="store_true",
        help="Recompute protein adjudications directly from unmeasured-home-gene-expression.csv.gz",
    )
    parser.add_argument(
        "--nose-skin-dir",
        type=Path,
        default=ROOT / "tmp/normal-atlas/nose-skin-2025-11-08",
    )
    parser.add_argument(
        "--durante-xlsx",
        type=Path,
        default=ROOT / "tmp/cta-other-literature/sources/human-olfactory-supplement.xlsx",
    )
    args = parser.parse_args(argv)
    if args.from_gene_table:
        gw = pd.read_csv(DATA_DIR / "unmeasured-home-gene-expression.csv.gz")
        build_protein_adjudications(gw)
        return
    if args.raw_dir is None:
        parser.error("--raw-dir is required unless --from-gene-table is passed")
    build_all(args.raw_dir, args.nose_skin_dir, args.durante_xlsx)


if __name__ == "__main__":
    main()
