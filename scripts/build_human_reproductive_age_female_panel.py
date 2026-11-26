#!/usr/bin/env python3
# Licensed under the Apache License, Version 2.0.
"""Build GTEx v11 donor-age metadata and the GSE107746 adult human oocyte/follicle RNA panel.

Generates:
1. ``oncoref/data/gtex-donor-age.csv`` from open-access
   ``GTEx_Analysis_v11_Annotations_SubjectPhenotypesDS.txt`` joined with
   ``oncoref.gtex.donor_sex()`` (recording ``20-29``..``70-79`` age brackets and
   ``is_reproductive_age`` for ages ``20-49``).
2. ``oncoref/data/human-adult-oocyte-follicle-rna.csv`` and ``.json`` from
   ``GSE107746`` (Zhang et al. 2018 *Mol Cell*, PMID:30420494), profiling adult
   human reproductive-age oocytes and granulosa cells across primordial, primary,
   secondary, antral, and preovulatory stages:
   - Inverts ``log2(FPKM + 1)`` to linear ``FPKM`` and normalizes within each
     library to linear ``TPM``.
   - Enforces the uniform ``<= 2.5%`` ectopic olfactory-receptor (``>= 5 TPM``)
     genomic-DNA contamination filter across all libraries.
   - Restricts to unambiguous 1:1 Ensembl 115 protein-coding symbols with low
     transcriptome cross-mapping risk (``shared_kmer_fraction < 0.05`` and no
     same-strand coding-exon overlap) so Cufflinks symbol-level FPKM never
     mis-attributes multicopy germline/embryo paralog families.
   - Summarizes each stage by the finite-sample replicated upper bound
     ``min(quantile(0.95), second_highest_sample)`` across stages with ``>= 2``
     QC-passing libraries.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import re
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from oncoref import gtex, library_qc  # noqa: E402
from oncoref.load_dataset import get_data  # noqa: E402

DEFAULT_SCRATCH_DIR = Path(
    os.environ.get("ONCOREF_REPRO_INPUTS_DIR", str(ROOT / "tmp" / "repro_age_inputs"))
)
PANEL_ID = "oncoref.human-adult-oocyte-follicle-rna.2026-11-26.1"
REPRODUCTIVE_AGE_BRACKETS = frozenset({"20-29", "30-39", "40-49"})
AGE_MIDPOINTS = {
    "20-29": 25,
    "30-39": 35,
    "40-49": 45,
    "50-59": 55,
    "60-69": 65,
    "70-79": 75,
}
OR_SYMBOL_RE = re.compile(r"^OR\d+[A-Z]+\d+$")


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(4 * 1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def build_gtex_donor_age(pheno_path: Path, out_csv: Path) -> dict:
    raw = pd.read_csv(pheno_path, sep="\t", dtype=str).fillna("")
    sex_df = gtex.donor_sex()
    chrom_by_donor = dict(zip(sex_df["donor_id"], sex_df["chromosomal_sex"]))

    rows = []
    for r in raw.itertuples(index=False):
        did = str(r.SUBJID).strip()
        if not did:
            continue
        sex_code = str(r.SEX).strip()
        rep_sex = "male" if sex_code == "1" else ("female" if sex_code == "2" else "unknown")
        c_sex = chrom_by_donor.get(did, "")
        bracket = str(r.AGE).strip()
        midpoint = AGE_MIDPOINTS.get(bracket, "")
        is_repro = bracket in REPRODUCTIVE_AGE_BRACKETS
        hardy = str(r.DTHHRDY).strip()
        rows.append(
            {
                "donor_id": did,
                "reported_sex": rep_sex,
                "chromosomal_sex": c_sex,
                "age_bracket": bracket,
                "age_midpoint": midpoint,
                "is_reproductive_age": is_repro,
                "hardy_scale": hardy,
            }
        )
    out_df = pd.DataFrame(rows).sort_values("donor_id").reset_index(drop=True)
    out_df.to_csv(out_csv, index=False, lineterminator="\n")
    return {
        "source_filename": pheno_path.name,
        "source_sha256": _sha256(pheno_path),
        "donors_total": len(out_df),
        "rnaseq_donors_matched": int(out_df["chromosomal_sex"].ne("").sum()),
        "xx_reproductive_age_donors": int(
            (out_df["chromosomal_sex"].eq("XX") & out_df["is_reproductive_age"]).sum()
        ),
        "xx_postmenopausal_donors": int(
            (out_df["chromosomal_sex"].eq("XX") & ~out_df["is_reproductive_age"]).sum()
        ),
        "output_sha256": _sha256(out_csv),
    }


def _parse_gse107746_sample_stages(
    series_matrix_path: Path, matrix_columns: list[str]
) -> dict[str, tuple[str, str]]:
    """Map every GSE107746 sample column to (cell_population, follicle_stage)."""
    with gzip.open(series_matrix_path, "rt", encoding="utf-8", errors="replace") as fh:
        lines = [ln.rstrip("\n") for ln in fh]
    titles = [
        x.strip('"')
        for x in next(ln for ln in lines if ln.startswith("!Sample_title")).split("\t")[1:]
    ]
    sources = [
        x.strip('"')
        for x in next(ln for ln in lines if ln.startswith("!Sample_source_name_ch1")).split("\t")[
            1:
        ]
    ]
    stage_by_col: dict[str, tuple[str, str]] = {}
    for title, src in zip(titles, sources):
        s_low = src.lower()
        pop = "oocyte" if "oocyte" in s_low else "granulosa"
        if "primordial" in s_low:
            st = "primordial"
        elif "primary" in s_low:
            st = "primary"
        elif "secondary" in s_low:
            st = "secondary"
        elif "antral" in s_low:
            st = "antral"
        elif "preovulatory" in s_low:
            st = "preovulatory"
        else:
            raise ValueError(f"Unrecognized GSE107746 source_name_ch1: {src!r}")
        stage_by_col[title] = (pop, st)

    for col in matrix_columns:
        if col in stage_by_col:
            continue
        if col.startswith("Preovulatory_follicle_MII_") and col.endswith("O"):
            stage_by_col[col] = ("oocyte", "preovulatory")
        elif col.startswith("Primary_follicle_") and col.endswith("G"):
            stage_by_col[col] = ("granulosa", "primary")
        else:
            raise ValueError(f"Unmapped GSE107746 matrix column: {col!r}")
    return stage_by_col


def _stage_replicated_bound(values_2d: np.ndarray) -> np.ndarray:
    """Return min(quantile(0.95), second_highest) per row across >=2 biological replicates."""
    if values_2d.shape[1] < 2:
        return np.zeros(values_2d.shape[0], dtype=np.float64)
    p95 = np.nanquantile(values_2d, 0.95, axis=1)
    sorted_vals = np.sort(np.where(np.isnan(values_2d), -np.inf, values_2d), axis=1)
    second_hi = sorted_vals[:, -2]
    second_hi = np.where(np.isneginf(second_hi), 0.0, second_hi)
    return np.minimum(p95, second_hi)


def build_adult_oocyte_follicle_panel(
    fpkm_log2_path: Path,
    series_matrix_path: Path,
    out_csv: Path,
    out_json: Path,
) -> dict:
    raw = pd.read_csv(fpkm_log2_path, sep="\t")
    gene_col = raw.columns[0]
    sample_cols = list(raw.columns[1:])
    stage_by_col = _parse_gse107746_sample_stages(series_matrix_path, sample_cols)

    # Collapse duplicate symbol rows if any, convert log2(FPKM + 1) -> linear FPKM -> linear TPM
    raw[gene_col] = raw[gene_col].astype(str).str.strip()
    raw = raw[raw[gene_col].ne("")].groupby(gene_col, sort=False)[sample_cols].max()
    log2_mat = raw.to_numpy(dtype=np.float64)
    fpkm_mat = np.clip(np.exp2(log2_mat) - 1.0, 0.0, None)
    col_sums = fpkm_mat.sum(axis=0)
    col_sums = np.where(col_sums > 0.0, col_sums, np.nan)
    tpm_mat = (fpkm_mat / col_sums) * 1e6
    tpm_df = pd.DataFrame(tpm_mat, index=raw.index, columns=sample_cols)

    # Apply uniform <= 2.5% ectopic olfactory receptor (>= 5.0 TPM) genomic-DNA contamination QC
    or_mask = np.array([bool(OR_SYMBOL_RE.match(sym)) for sym in tpm_df.index], dtype=bool)
    n_or_assayed = int(or_mask.sum())
    if n_or_assayed < 200:
        raise ValueError(f"Expected >=200 olfactory receptors in GSE107746, got {n_or_assayed}")
    or_shares = (tpm_mat[or_mask, :] >= library_qc.BACKGROUND_TPM).mean(axis=0)
    qc_pass_cols = [col for col, sh in zip(sample_cols, or_shares) if float(sh) <= 0.025]
    qc_fail_cols = [col for col in sample_cols if col not in set(qc_pass_cols)]

    # Group QC-passing columns by (population, stage)
    cols_by_group: dict[tuple[str, str], list[str]] = defaultdict(list)
    for col in qc_pass_cols:
        cols_by_group[stage_by_col[col]].append(col)

    # Build strict 1:1 symbol -> Ensembl 115 PEG crosswalk with <5% 31-mer cDNA paralog guard
    base_ts = pd.read_csv(
        ROOT / "tests/fixtures/release-1.50-inputs/oncoref/data/target-selection.csv.gz",
        usecols=["entity_id", "protein_group_label", "member_gene_ids"],
    )
    risk_df = get_data("gene-mapping-risk")
    high_risk_genes = set(
        risk_df.loc[
            (pd.to_numeric(risk_df["shared_kmer_fraction"], errors="coerce").fillna(0.0) >= 0.05)
            | risk_df["host_overlapped"].fillna(False).astype(bool)
            | risk_df["mapping_risk"].fillna(False).astype(bool),
            "gene_id",
        ].astype(str)
    )

    sym_to_entities: dict[str, set[str]] = defaultdict(set)
    entity_meta: dict[str, tuple[str, str, bool]] = {}
    for r in base_ts.itertuples(index=False):
        eid = str(r.entity_id)
        lbl = str(r.protein_group_label)
        gids = json.loads(r.member_gene_ids)
        is_single_low_risk = len(gids) == 1 and str(gids[0]) not in high_risk_genes
        entity_meta[eid] = (lbl, str(gids[0]) if gids else "", is_single_low_risk)
        if "/" not in lbl:
            sym_to_entities[lbl].add(eid)

    # Add unambiguous 1:1 NCBI synonym fallback ONLY when official symbol is absent from matrix
    syn_df = get_data("ncbi-symbol-synonyms")
    off_to_eid = {sym: next(iter(eids)) for sym, eids in sym_to_entities.items() if len(eids) == 1}
    syn_candidates: dict[str, set[str]] = defaultdict(set)
    for r in syn_df.itertuples(index=False):
        syn = str(r.alias).strip()
        off = str(r.official_symbol).strip()
        if syn and off in off_to_eid and syn not in tpm_df.index:
            syn_candidates[syn].add(off_to_eid[off])

    unique_row_for_eid: dict[str, str] = {}
    for sym, eids in sym_to_entities.items():
        if len(eids) == 1 and sym in tpm_df.index:
            eid = next(iter(eids))
            if entity_meta[eid][2]:
                unique_row_for_eid[eid] = sym
    for syn, eids in syn_candidates.items():
        if len(eids) == 1 and syn in tpm_df.index:
            eid = next(iter(eids))
            if eid not in unique_row_for_eid and entity_meta[eid][2]:
                unique_row_for_eid[eid] = syn

    ordered_eids = sorted(unique_row_for_eid)
    ordered_syms = [unique_row_for_eid[eid] for eid in ordered_eids]
    sub_tpm = tpm_df.loc[ordered_syms, qc_pass_cols]

    oocyte_stages = [
        st for (pop, st), cs in sorted(cols_by_group.items()) if pop == "oocyte" and len(cs) >= 2
    ]
    gran_stages = [
        st for (pop, st), cs in sorted(cols_by_group.items()) if pop == "granulosa" and len(cs) >= 2
    ]

    oocyte_bounds = {
        st: _stage_replicated_bound(sub_tpm[cols_by_group[("oocyte", st)]].to_numpy(dtype=float))
        for st in oocyte_stages
    }
    gran_bounds = {
        st: _stage_replicated_bound(sub_tpm[cols_by_group[("granulosa", st)]].to_numpy(dtype=float))
        for st in gran_stages
    }
    all_oocyte_cols = [c for st in oocyte_stages for c in cols_by_group[("oocyte", st)]]
    all_gran_cols = [c for st in gran_stages for c in cols_by_group[("granulosa", st)]]
    oocyte_max_arr = sub_tpm[all_oocyte_cols].max(axis=1).to_numpy(dtype=float)
    gran_max_arr = sub_tpm[all_gran_cols].max(axis=1).to_numpy(dtype=float)

    rows = []
    for idx, eid in enumerate(ordered_eids):
        lbl, gid, _ = entity_meta[eid]
        src_sym = ordered_syms[idx]
        best_oo_st = max(oocyte_stages, key=lambda st: float(oocyte_bounds[st][idx]))
        oo_supp = float(oocyte_bounds[best_oo_st][idx])
        best_gr_st = max(gran_stages, key=lambda st: float(gran_bounds[st][idx]))
        gr_supp = float(gran_bounds[best_gr_st][idx])
        foll_supp = max(oo_supp, gr_supp)
        if foll_supp < 0.5:
            continue
        dominant_cell = "oocyte" if oo_supp >= gr_supp else "granulosa"
        dominant_stage = best_oo_st if oo_supp >= gr_supp else best_gr_st
        rows.append(
            {
                "entity_id": eid,
                "protein_group_label": lbl,
                "gene_id": gid,
                "source_symbol": src_sym,
                "adult_follicle_supported_tpm": round(foll_supp, 6),
                "oocyte_supported_tpm": round(oo_supp, 6),
                "top_oocyte_stage": best_oo_st,
                "oocyte_max_sample_tpm": round(float(oocyte_max_arr[idx]), 6),
                "granulosa_supported_tpm": round(gr_supp, 6),
                "top_granulosa_stage": best_gr_st,
                "granulosa_max_sample_tpm": round(float(gran_max_arr[idx]), 6),
                "dominant_compartment_cell": dominant_cell,
                "dominant_follicle_stage": dominant_stage,
            }
        )

    out_df = (
        pd.DataFrame(rows)
        .sort_values(
            ["adult_follicle_supported_tpm", "protein_group_label"], ascending=[False, True]
        )
        .reset_index(drop=True)
    )
    out_df.to_csv(
        out_csv,
        index=False,
        compression={"method": "gzip", "mtime": 0},
        lineterminator="\n",
    )
    csv_sha = _sha256(out_csv)

    stage_summary = {f"{pop}:{st}": len(cs) for (pop, st), cs in sorted(cols_by_group.items())}
    receipt = {
        "schema": "oncoref.human-adult-oocyte-follicle-rna.v1",
        "panel_id": PANEL_ID,
        "study_accession": "GSE107746",
        "citation": "Zhang et al. 2018 Mol Cell (PMID:30420494)",
        "anatomical_compartment": "ovary",
        "developmental_stage": "adult (reproductive-age human ovarian follicles)",
        "unit": "TPM (normalized within-library from linear Cufflinks FPKM)",
        "inputs": {
            fpkm_log2_path.name: {
                "bytes": fpkm_log2_path.stat().st_size,
                "sha256": _sha256(fpkm_log2_path),
            },
            series_matrix_path.name: {
                "bytes": series_matrix_path.stat().st_size,
                "sha256": _sha256(series_matrix_path),
            },
        },
        "library_qc": {
            "rule": "olfactory_share_ge5 <= 0.025",
            "assayed_olfactory_genes": n_or_assayed,
            "libraries_total": len(sample_cols),
            "libraries_qc_pass": len(qc_pass_cols),
            "libraries_qc_excluded": len(qc_fail_cols),
            "qc_pass_by_stage": stage_summary,
        },
        "identity_mapping": {
            "uniquely_mapped_low_risk_pegs": len(ordered_eids),
            "paralog_kmer_ceiling": 0.05,
            "reported_rows_ge_0_5_tpm": len(out_df),
        },
        "sha256": csv_sha,
    }
    out_json.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    return receipt


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scratch-dir", type=Path, default=DEFAULT_SCRATCH_DIR)
    parser.add_argument("--root", type=Path, default=ROOT)
    args = parser.parse_args(argv)

    data_dir = args.root / "oncoref" / "data"
    age_info = build_gtex_donor_age(
        args.scratch_dir / "GTEx_v11_SubjectPhenotypesDS.txt",
        data_dir / "gtex-donor-age.csv",
    )
    foll_info = build_adult_oocyte_follicle_panel(
        args.scratch_dir / "GSE107746_Folliculogenesis_FPKM.log2.txt.gz",
        args.scratch_dir / "GSE107746_series_matrix.txt.gz",
        data_dir / "human-adult-oocyte-follicle-rna.csv.gz",
        data_dir / "human-adult-oocyte-follicle-rna.json",
    )
    print(json.dumps({"gtex_donor_age": age_info, "adult_follicle_panel": foll_info}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
