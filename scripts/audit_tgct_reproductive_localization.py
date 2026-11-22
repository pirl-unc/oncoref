#!/usr/bin/env python
"""Audit somatic testicular admixture vs malignant germ-cell expression in TCGA-TGCT (#713).

Computes orthogonal somatic niche signatures (Leydig steroidogenic interstitial cells
and Sertoli intratubular nurse cells), post-meiotic spermatid stoichiometry against
GTEx v11 Testis, and malignant TGCT lineage markers across all 154 TCGA-TGCT libraries,
then evaluates independent non-AAFL and non-R-TGCT tumor witness support across all
132 curation tumor matrices for every protein peaking in TCGA-2G-AAFL-01.
"""

from __future__ import annotations

import argparse
import gzip
import json
from pathlib import Path

import numpy as np
import pandas as pd

from oncoref import (
    curation,
    expression,
    library_qc,
    reference_data,
    source_matrices,
    target_selection,
)

LEYDIG_GENES: tuple[str, ...] = (
    "INSL3",
    "CYP17A1",
    "STAR",
    "CYP11A1",
    "HSD3B2",
    "LHCGR",
)
SERTOLI_GENES: tuple[str, ...] = (
    "INHA",
    "CLU",
    "SOX9",
    "WT1",
    "AMH",
    "FSHR",
    "DMRT1",
)
SPERMATID_GENES: tuple[str, ...] = (
    "PRM1",
    "PRM2",
    "TNP1",
    "LELP1",
    "SMCP",
    "PHF7",
    "ACR",
    "SPATA19",
    "PGK2",
    "TSSK2",
)
MALIGNANT_LINEAGE_PANELS: dict[str, tuple[str, ...]] = {
    "seminoma_gcnis": ("POU5F1", "NANOG", "SOX17", "TFAP2C", "KIT"),
    "embryonal_carcinoma": ("POU5F1", "NANOG", "SOX2", "DNMT3B"),
    "yolk_sac_tumor": ("AFP", "GPC3", "FOXA2"),
    "trophoblastic": ("CGB3", "CGA", "GATA3"),
}


def _load_gtex_v11_testis_tpm() -> tuple[dict[str, float], dict[str, float]]:
    """Return (gene_id -> Testis median TPM, symbol -> Testis median TPM) from GTEx v11."""
    path = reference_data.pinned_path("gtex_tissue_rna", "v11")
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        _version = handle.readline()
        _dims = handle.readline()
        df = pd.read_csv(handle, sep="\t")
    testis_col = next(c for c in df.columns if c.lower() == "testis")
    df["gene_id"] = df["Name"].astype(str).str.split(".").str[0]
    df["Testis"] = pd.to_numeric(df[testis_col], errors="coerce").fillna(0.0)
    by_gene = df.groupby("gene_id")["Testis"].max().to_dict()
    by_sym = df.groupby("Description")["Testis"].max().to_dict()
    return by_gene, by_sym


def _panel_ratios_to_testis(
    by_symbol: pd.DataFrame,
    genes: tuple[str, ...],
    gtex_by_symbol: dict[str, float],
    *,
    min_testis_tpm: float = 1.0,
) -> tuple[pd.Series, pd.Series, pd.Series]:
    usable = [
        g for g in genes if g in by_symbol.index and gtex_by_symbol.get(g, 0.0) >= min_testis_tpm
    ]
    if not usable:
        nan = pd.Series(np.nan, index=by_symbol.columns)
        return nan, nan, pd.Series(0, index=by_symbol.columns, dtype=int)
    ref = pd.Series({g: gtex_by_symbol[g] for g in usable}, dtype=float)
    ratios = by_symbol.loc[usable].div(ref, axis=0)
    return (
        ratios.min(axis=0),
        ratios.median(axis=0),
        (ratios >= 0.01).sum(axis=0).astype(int),
    )


def build_tgct_library_audit(
    tgct_raw: pd.DataFrame,
    gtex_by_symbol: dict[str, float],
) -> pd.DataFrame:
    sample_cols = [c for c in tgct_raw.columns if c not in ("Ensembl_Gene_ID", "Symbol")]
    sigs = library_qc.signatures(tgct_raw, sample_cols)
    qc_df = expression.sample_expression_qc_from_matrix(tgct_raw, cancer_type="TGCT").set_index(
        "sample_id"
    )

    numeric = tgct_raw[sample_cols].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    by_symbol = (
        pd.concat([tgct_raw[["Symbol"]], numeric], axis=1)
        .groupby("Symbol", sort=False)[sample_cols]
        .sum()
    )

    leydig_min, leydig_med, leydig_n1pct = _panel_ratios_to_testis(
        by_symbol, LEYDIG_GENES, gtex_by_symbol, min_testis_tpm=1.0
    )
    sertoli_min, sertoli_med, sertoli_n1pct = _panel_ratios_to_testis(
        by_symbol, SERTOLI_GENES, gtex_by_symbol, min_testis_tpm=1.0
    )
    sperm_min, sperm_med, sperm_n1pct = _panel_ratios_to_testis(
        by_symbol, SPERMATID_GENES, gtex_by_symbol, min_testis_tpm=10.0
    )

    sperm_usable = [
        g for g in SPERMATID_GENES if g in by_symbol.index and gtex_by_symbol.get(g, 0.0) >= 10.0
    ]
    ref_log = np.log10(np.array([gtex_by_symbol[g] for g in sperm_usable], dtype=float) + 1.0)

    rows = []
    for sample_id in sample_cols:
        sig_row = sigs.loc[sample_id]
        notes = library_qc.notes(sig_row)
        has_r_note = "R" in notes
        sample_vec = np.log10(by_symbol.loc[sperm_usable, sample_id].to_numpy(dtype=float) + 1.0)
        if np.std(sample_vec) > 1e-9 and np.std(ref_log) > 1e-9:
            sperm_r = float(np.corrcoef(sample_vec, ref_log)[0, 1])
        else:
            sperm_r = 0.0

        pou5f1 = float(by_symbol.loc["POU5F1", sample_id]) if "POU5F1" in by_symbol.index else 0.0
        sox17 = float(by_symbol.loc["SOX17", sample_id]) if "SOX17" in by_symbol.index else 0.0
        sox2 = float(by_symbol.loc["SOX2", sample_id]) if "SOX2" in by_symbol.index else 0.0
        afp = float(by_symbol.loc["AFP", sample_id]) if "AFP" in by_symbol.index else 0.0
        cgb3 = float(by_symbol.loc["CGB3", sample_id]) if "CGB3" in by_symbol.index else 0.0

        if pou5f1 >= 200.0 and sox2 >= 15.0 and sox2 > sox17:
            inferred_program = "embryonal_carcinoma_dominant"
        elif pou5f1 >= 200.0 and sox17 >= 10.0:
            inferred_program = "seminoma_dominant"
        elif afp >= 100.0:
            inferred_program = "yolk_sac_dominant"
        elif cgb3 >= 100.0:
            inferred_program = "trophoblastic_dominant"
        else:
            inferred_program = "teratoma_or_mixed_differentiated"

        l_med = float(leydig_med[sample_id])
        s_med = float(sertoli_med[sample_id])
        sp_med = float(sperm_med[sample_id])
        t_content = float(sig_row["testis_tissue_content"])

        if has_r_note and l_med >= 0.05 and s_med >= 0.05 and sp_med >= 0.02:
            localization_call = "admixed_seminiferous_and_interstitial_parenchyma"
        elif has_r_note and (l_med >= 0.01 or s_med >= 0.02):
            localization_call = "trace_admixed_testicular_parenchyma"
        elif has_r_note:
            localization_call = "reproductive_note_somatic_niche_low"
        else:
            localization_call = "no_coordinated_spermatid_program"

        rows.append(
            {
                "sample_id": sample_id,
                "case_id": sample_id[:12],
                "sample_qc_status": str(qc_df.loc[sample_id, "sample_qc_status"]),
                "library_qc_note_R": bool(has_r_note),
                "testis_tissue_content": round(t_content, 6),
                "spermatid_median_ratio_to_gtex": round(sp_med, 6),
                "spermatid_min_ratio_to_gtex": round(float(sperm_min[sample_id]), 6),
                "spermatid_markers_ge1pct_gtex": int(sperm_n1pct[sample_id]),
                "spermatid_gtex_log10_r": round(sperm_r, 4),
                "leydig_median_ratio_to_gtex": round(l_med, 6),
                "leydig_min_ratio_to_gtex": round(float(leydig_min[sample_id]), 6),
                "leydig_markers_ge1pct_gtex": int(leydig_n1pct[sample_id]),
                "sertoli_median_ratio_to_gtex": round(s_med, 6),
                "sertoli_min_ratio_to_gtex": round(float(sertoli_min[sample_id]), 6),
                "sertoli_markers_ge1pct_gtex": int(sertoli_n1pct[sample_id]),
                "PRM1_tpm": round(float(by_symbol.loc["PRM1", sample_id]), 3),
                "PRM2_tpm": round(float(by_symbol.loc["PRM2", sample_id]), 3),
                "TNP1_tpm": round(float(by_symbol.loc["TNP1", sample_id]), 3),
                "INSL3_tpm": round(float(by_symbol.loc["INSL3", sample_id]), 3),
                "CYP17A1_tpm": round(float(by_symbol.loc["CYP17A1", sample_id]), 3),
                "STAR_tpm": round(float(by_symbol.loc["STAR", sample_id]), 3),
                "CLU_tpm": round(float(by_symbol.loc["CLU", sample_id]), 3),
                "INHA_tpm": round(float(by_symbol.loc["INHA", sample_id]), 3),
                "SOX9_tpm": round(float(by_symbol.loc["SOX9", sample_id]), 3),
                "POU5F1_tpm": round(pou5f1, 3),
                "SOX17_tpm": round(sox17, 3),
                "SOX2_tpm": round(sox2, 3),
                "AFP_tpm": round(afp, 3),
                "inferred_malignant_program": inferred_program,
                "localization_call": localization_call,
            }
        )
    out = pd.DataFrame(rows).sort_values(
        ["library_qc_note_R", "testis_tissue_content", "sample_id"],
        ascending=[False, False, True],
    )
    return out.reset_index(drop=True)


def _load_hpa_v23_testis_tpm() -> dict[str, float]:
    """Return gene_id -> Testis TPM from HPA v23 rna_tissue_hpa."""
    path = reference_data.pinned_path("rna_tissue_hpa", "v23")
    df = pd.read_csv(path, sep="\t")
    sub = df[df["Tissue"].astype(str).str.lower().eq("testis")].copy()
    sub["gene_id"] = sub["Gene"].astype(str).str.split(".").str[0]
    sub["TPM"] = pd.to_numeric(sub["TPM"], errors="coerce").fillna(0.0)
    return sub.groupby("gene_id")["TPM"].max().to_dict()


def build_aafl_protein_support(
    root: Path,
    library_audit: pd.DataFrame,
    gtex_by_gene: dict[str, float],
    gtex_by_sym: dict[str, float] | None = None,
    hpa_testis_by_gene: dict[str, float] | None = None,
) -> pd.DataFrame:
    ts = target_selection.results()
    aafl_rows = ts[
        ts["tumor_top_sample"].astype(str).str.contains("TCGA-2G-AAFL-01", na=False)
    ].copy()

    target_entities = []
    needed_genes: set[str] = set()
    for row in aafl_rows.itertuples(index=False):
        members = [str(g).split(".")[0] for g in json.loads(str(row.member_gene_ids))]
        needed_genes.update(members)
        target_entities.append((row, members))

    r_tgct_samples = set(
        library_audit.loc[library_audit["library_qc_note_R"], "sample_id"].astype(str)
    )

    # Scan all 132 curation tumor matrices across QC-passing samples.
    qc_manifest = expression.source_matrix_sample_qc_manifest()
    pass_by_code: dict[str, set[str]] = {}
    for code, sub in qc_manifest[qc_manifest["sample_qc_status"].eq("pass")].groupby("cancer_code"):
        pass_by_code[str(code)] = set(sub["sample_id"].astype(str))

    per_entity_records: dict[str, list[tuple[str, str, float]]] = {
        row.entity_id: [] for row, _ in target_entities
    }

    for code, version, _sha in curation.required_tumor_matrices(root):
        path = source_matrices.local_path(code, version=version)
        df = pd.read_parquet(path)
        sample_cols = [c for c in df.columns if c not in ("Ensembl_Gene_ID", "Symbol")]
        allowed_samples = pass_by_code.get(code)
        if allowed_samples is not None:
            sample_cols = [c for c in sample_cols if c in allowed_samples]
        if not sample_cols:
            continue
        gene_ids = df["Ensembl_Gene_ID"].astype(str).str.split(".").str[0]
        mask = gene_ids.isin(needed_genes)
        if not mask.any():
            continue
        sub = df.loc[mask, ["Ensembl_Gene_ID", *sample_cols]].copy()
        sub["gene_id"] = sub["Ensembl_Gene_ID"].astype(str).str.split(".").str[0]
        grouped = sub.groupby("gene_id", sort=False)[sample_cols].sum()

        for row, members in target_entities:
            present = [m for m in members if m in grouped.index]
            if not present:
                continue
            vals = grouped.loc[present].sum(axis=0)
            for sid, val in vals.items():
                fval = float(val)
                if fval > 0.0:
                    per_entity_records[row.entity_id].append((code, str(sid), fval))

    out_rows = []
    for row, members in target_entities:
        eid = str(row.entity_id)
        threshold = float(row.single_tumor_threshold_tpm)
        gtex_present = any(m in gtex_by_gene for m in members)
        gtex_testis = sum(gtex_by_gene.get(m, 0.0) for m in members) if gtex_present else np.nan
        if pd.notna(gtex_testis) and gtex_testis > 0:
            ref_testis = float(gtex_testis)
            ref_source = "GTEx_v11_median"
        else:
            hpa_val = sum((hpa_testis_by_gene or {}).get(m, 0.0) for m in members)
            ref_testis = float(hpa_val) if hpa_val > 0 else np.nan
            ref_source = "HPA_v23_testis_tpm" if hpa_val > 0 else "unmeasured"
        obs = per_entity_records[eid]

        aafl_tpm = 0.0
        other_r_max, other_r_sample = 0.0, ""
        non_r_tgct_max, non_r_tgct_sample = 0.0, ""
        non_tgct_max, non_tgct_code, non_tgct_sample = 0.0, "", ""
        n_non_aafl_ge_thresh = 0
        n_clean_witnesses_ge_thresh = 0
        clean_cohorts_ge_thresh: set[str] = set()

        for code, sid, val in obs:
            if code == "TGCT" and sid == "TCGA-2G-AAFL-01":
                aafl_tpm = max(aafl_tpm, val)
                continue
            if val >= threshold:
                n_non_aafl_ge_thresh += 1
            if code == "TGCT" and sid in r_tgct_samples:
                if val > other_r_max:
                    other_r_max, other_r_sample = val, sid
            elif code == "TGCT":
                if val > non_r_tgct_max:
                    non_r_tgct_max, non_r_tgct_sample = val, sid
                if val >= threshold:
                    n_clean_witnesses_ge_thresh += 1
                    clean_cohorts_ge_thresh.add(code)
            else:
                if val > non_tgct_max:
                    non_tgct_max, non_tgct_code, non_tgct_sample = val, code, sid
                if val >= threshold:
                    n_clean_witnesses_ge_thresh += 1
                    clean_cohorts_ge_thresh.add(code)

        if str(row.category) != "CTA":
            support_tier = "below_cancer_threshold"
        elif n_clean_witnesses_ge_thresh >= 1:
            support_tier = "confirmed_outside_r_tgct"
        elif n_non_aafl_ge_thresh >= 1:
            support_tier = "r_tgct_multi_library_only"
        else:
            support_tier = "aafl_sole_witness_admixed_parenchyma_candidate"

        ratio_to_testis = (
            aafl_tpm / ref_testis if pd.notna(ref_testis) and ref_testis > 0 else np.nan
        )
        out_rows.append(
            {
                "entity_id": eid,
                "protein_group_label": str(row.protein_group_label),
                "member_gene_ids": ";".join(members),
                "released_category": str(row.category),
                "single_tumor_threshold_tpm": round(threshold, 4),
                "gtex_v11_testis_median_tpm": (
                    round(float(gtex_testis), 3) if pd.notna(gtex_testis) else np.nan
                ),
                "reference_testis_tpm": (
                    round(float(ref_testis), 3) if pd.notna(ref_testis) else np.nan
                ),
                "reference_testis_source": ref_source,
                "aafl_tpm": round(aafl_tpm, 4),
                "aafl_to_reference_testis_ratio": (
                    round(float(ratio_to_testis), 6) if pd.notna(ratio_to_testis) else np.nan
                ),
                "max_other_r_tgct_tpm": round(other_r_max, 4),
                "top_other_r_tgct_sample": other_r_sample,
                "max_non_r_tgct_tpm": round(non_r_tgct_max, 4),
                "top_non_r_tgct_sample": non_r_tgct_sample,
                "max_non_tgct_tpm": round(non_tgct_max, 4),
                "top_non_tgct_cancer_code": non_tgct_code,
                "top_non_tgct_sample": non_tgct_sample,
                "released_samples_at_tenfold": int(row.samples_at_tenfold),
                "non_aafl_witnesses_at_threshold": int(n_non_aafl_ge_thresh),
                "non_r_tgct_or_non_tgct_witnesses_at_threshold": int(n_clean_witnesses_ge_thresh),
                "clean_cancer_codes_at_threshold": ";".join(sorted(clean_cohorts_ge_thresh)),
                "support_tier": support_tier,
            }
        )

    order = {
        "aafl_sole_witness_admixed_parenchyma_candidate": 0,
        "r_tgct_multi_library_only": 1,
        "confirmed_outside_r_tgct": 2,
        "below_cancer_threshold": 3,
    }
    df_out = pd.DataFrame(out_rows)
    df_out["_ord"] = df_out["support_tier"].map(order)
    df_out = df_out.sort_values(
        ["_ord", "aafl_tpm", "protein_group_label"],
        ascending=[True, False, True],
    ).drop(columns=["_ord"])
    return df_out.reset_index(drop=True)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--root",
        type=Path,
        default=Path(__file__).resolve().parents[1],
    )
    parser.add_argument(
        "--library-output",
        type=Path,
        default=None,
    )
    parser.add_argument(
        "--protein-output",
        type=Path,
        default=None,
    )
    args = parser.parse_args(argv)
    root = args.root.resolve()
    lib_out = args.library_output or (
        root / "oncoref/data/tgct-reproductive-localization-audit.csv"
    )
    prot_out = args.protein_output or (root / "oncoref/data/tgct-aafl-protein-support.csv")

    gtex_by_gene, gtex_by_sym = _load_gtex_v11_testis_tpm()
    hpa_testis_by_gene = _load_hpa_v23_testis_tpm()
    tgct_raw = pd.read_parquet(source_matrices.local_path("TGCT"))

    lib_df = build_tgct_library_audit(tgct_raw, gtex_by_sym)
    prot_df = build_aafl_protein_support(
        root, lib_df, gtex_by_gene, gtex_by_sym, hpa_testis_by_gene
    )

    lib_out.parent.mkdir(parents=True, exist_ok=True)
    prot_out.parent.mkdir(parents=True, exist_ok=True)
    lib_df.to_csv(lib_out, index=False)
    prot_df.to_csv(prot_out, index=False)

    summary = {
        "n_tgct_libraries": len(lib_df),
        "n_r_libraries": int(lib_df["library_qc_note_R"].sum()),
        "localization_calls": lib_df["localization_call"].value_counts().to_dict(),
        "n_aafl_peak_proteins": len(prot_df),
        "support_tiers": prot_df["support_tier"].value_counts().to_dict(),
    }
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
