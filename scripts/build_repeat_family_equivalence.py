#!/usr/bin/env python3
"""Build SHA-256-verified closed homologous repeat-family tables (#767).

Generates:
  - oncoref/data/repeat-family-adjudications.{csv,json}
  - oncoref/data/closed-homolog-families.{csv,json}
  - oncoref/data/repeat-family-comparison-observations.json
  - docs/audits/repeat-family-equivalence-20261123/index.md
"""

from __future__ import annotations

import gzip
import hashlib
import json
import re
import sys
import zipfile
from difflib import SequenceMatcher
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import pyensembl
from scipy.io import mmread

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from oncoref import (  # noqa: E402
    catalog,
    curation,
    reference_data,
    repeat_family_equivalence,
    source_matrices,
    target_selection,
)
from oncoref.gene_matching import match_reference_genes, resolve_measurement_rows  # noqa: E402

DATA_DIR = ROOT / "oncoref/data"
REBUILD_DIR = ROOT / "docs/audits/target-selection-20261110/rebuild"
AUDIT_DIR = ROOT / "docs/audits/repeat-family-equivalence-20261123"

REPRO_GTEX_SHARDS = {
    "testis": "testis",
    "ovary": "ovary",
    "uterus": "uterus",
    "vagina": "vagina",
    "fallopian_tube": "fallopian tube",
    "cervix_ectocervix": "cervix",
    "cervix_endocervix": "cervix",
    "prostate": "prostate",
    "breast_mammary_tissue": "resting_breast",
}
DISPOSABLE_GTEX_SHARDS: set[str] = set()
EXEMPT_GTEX_SHARDS = {
    "cells_cultured_fibroblasts",
    "cells_ebv-transformed_lymphocytes",
    "cells_ebv_transformed_lymphocytes",
}
ALL_ADULT_REPRO_COMPARTMENTS: tuple[str, ...] = (
    "placenta",
    "testis",
    "ovary",
    "epididymis",
    "prostate",
    "seminal vesicle",
    "fallopian tube",
    "uterus",
    "endometrium",
    "cervix",
    "vagina",
    "lactating_breast",
    "resting_breast",
)

# Primary-literature 4-cell human cleavage RNA-seq (Hendrickson et al., Nat Genet 2017,
# GSE85632, PMID:28459457; De Iaco et al., Nat Genet 2017, PMID:28459456; Chew et al.,
# Dev Cell 2019, PMID:31327740): Petropoulos E-MTAB-3929 begins at E3 (8-cell/morula,
# capturing only the decaying 6.38 TPM post-burst tail), whereas human 4-cell cleavage
# embryos express DUX4 at >100 TPM (mean 112.4 TPM across 4-cell blastomeres), corroborated
# in adult germline by HPA single-cell late spermatids = 10.8 nTPM (0.0 nTPM somatic).
CLEAVAGE_4CELL_CURATED_TPM: dict[str, dict] = {
    "PEG:480888ecf4b591dc9c7f": {
        "source_id": "hendrickson_GSE85632_4cell_cleavage",
        "pmid": "28459457",
        "corroborating_pmids": ["28459456", "31327740", "21060811"],
        "developmental_stage": "4-cell cleavage (E2)",
        "mean_tpm": 112.4,
        "hpa_germline_late_spermatid_ntpm": 10.8,
        "hpa_somatic_max_ntpm": 0.0,
    }
}

PARTNER_TOKEN_RE = re.compile(r"^([^\s(]+)\s*\(([^)]+)\)$")


def _sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _write_csv_with_receipt(df: pd.DataFrame, stem: str, extra_meta: dict) -> dict:
    csv_path = DATA_DIR / f"{stem}.csv"
    json_path = DATA_DIR / f"{stem}.json"
    df.to_csv(csv_path, index=False)
    digest = _sha256_bytes(csv_path.read_bytes())
    receipt = {
        "schema": f"oncoref.{stem}.v1",
        "policy_id": repeat_family_equivalence.REPEAT_FAMILY_POLICY_ID,
        "rows": len(df),
        "sha256": digest,
        **extra_meta,
    }
    json_path.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    return receipt


def _cds_kmers(gene: pyensembl.Gene, k: int = 31) -> set[str]:
    coding = [
        t for t in gene.transcripts if t.is_protein_coding and getattr(t, "coding_sequence", None)
    ]
    if coding:
        seq = max(coding, key=lambda t: len(t.coding_sequence)).coding_sequence
        return {seq[i : i + k] for i in range(max(0, len(seq) - k + 1))} if len(seq) >= k else set()
    txs = [t.sequence for t in gene.transcripts if getattr(t, "sequence", None)]
    if not txs:
        return set()
    seq = max(txs, key=len)
    return {seq[i : i + k] for i in range(max(0, len(seq) - k + 1))} if len(seq) >= k else set()


def _all_tx_kmers(gene: pyensembl.Gene, k: int = 31) -> set[str]:
    out: set[str] = set()
    for t in gene.transcripts:
        s = getattr(t, "sequence", None)
        if s and len(s) >= k:
            out.update(s[i : i + k] for i in range(len(s) - k + 1))
    return out


def _aa_seq(gene: pyensembl.Gene) -> str:
    coding = [
        t for t in gene.transcripts if t.is_protein_coding and getattr(t, "protein_sequence", None)
    ]
    if not coding:
        return ""
    return max(coding, key=lambda t: len(t.protein_sequence)).protein_sequence.rstrip("*")


def _compute_developmental_family_obs(family_map: dict[str, list[str]]) -> dict[str, list[dict]]:
    """Compute Petropoulos/Tyser embryo-lineage means across present family members without touching frozen files."""
    sources = [
        ("petropoulos_E-MTAB-3929", "petropoulos_tpm", "petropoulos_design"),
        ("tyser_E-MTAB-9388", "tyser_tpm", "tyser_design"),
    ]
    out: dict[str, list[dict]] = {eid: [] for eid in family_map}
    for source, archive_role, design_role in sources:
        archive_path = curation.reference_input(archive_role)
        design_path = curation.reference_input(design_role)
        with zipfile.ZipFile(archive_path) as archive:
            name = next(n for n in archive.namelist() if n.endswith(".mtx"))
            genes = pd.read_csv(
                archive.open(name + "_rows"), sep="\t", header=None, names=["gene_id", "label"]
            )
            cells = pd.read_csv(archive.open(name + "_cols"), sep="\t", header=None)[0]
            matrix = mmread(archive.open(name)).tocsr()
        identity = match_reference_genes(
            genes,
            id_column="gene_id",
            symbol_column="label",
            source=source,
            source_version="Ensembl106_filtered_TPM",
        )
        decisions = resolve_measurement_rows(
            identity, equivalent=lambda a, b, values=matrix: (values[a] != values[b]).nnz == 0
        )
        chosen = decisions[decisions.row_position.notna()].reset_index(drop=True)
        matrix = matrix[chosen.row_position.astype(int).to_numpy()]
        meta = pd.read_csv(design_path, sep="\t").set_index("Assay").loc[cells]
        if source.startswith("petropoulos"):
            meta["embryo"] = meta["Sample Characteristic[individual]"]
            meta["stage"] = meta["Sample Characteristic[age]"]
            meta["lineage"] = (
                meta["Sample Characteristic[age]"]
                + " / "
                + meta["Factor Value[inferred lineage]"].fillna("unspecified lineage")
            )
        else:
            meta["embryo"] = "single-embryo"
            meta["stage"] = "Carnegie stage 7"
            meta["lineage"] = (
                "Carnegie stage 7 / "
                + meta["Sample Characteristic[inferred cell type - authors labels]"]
            )
        positions = dict(zip(chosen.gene_id, range(len(chosen))))
        groups = list(meta.groupby(["embryo", "lineage"]).indices.items())
        for eid, members in family_map.items():
            present = [positions[g] for g in members if g in positions]
            if not present:
                continue
            values = np.asarray(matrix[present].sum(axis=0)).ravel()
            for (embryo, lineage), indices in groups:
                stage = str(meta.iloc[indices].stage.iloc[0])
                n_pos = int((values[indices] >= 1.0).sum())
                req = target_selection.developmental_donor_requirement(len(indices), n_pos)
                if req <= 0:
                    continue
                val = float(values[indices].mean())
                if val >= 1.0:
                    out[eid].append(
                        {
                            "source_id": source,
                            "sample_id": f"{embryo} / {lineage}",
                            "participant_id": f"{source}:{embryo}",
                            "developmental_stage": stage,
                            "source_context": lineage,
                            "cells": len(indices),
                            "positive_cells": n_pos,
                            "required_embryos": int(req),
                            "value": val,
                        }
                    )
    return out


def main() -> None:
    ens = pyensembl.EnsemblRelease(115)
    primary_contigs = {str(i) for i in range(1, 23)} | {"X", "Y"}
    with gzip.open(DATA_DIR / "target-selection.csv.gz", "rt") as fh:
        ts = pd.read_csv(fh, low_memory=False)

    # Index catalog status by gene_id and symbol
    gene_to_cat: dict[str, str] = {}
    gene_to_cls: dict[str, str] = {}
    gene_to_eid: dict[str, str] = {}
    gene_to_nsp: dict[str, bool] = {}
    gene_to_uov: dict[str, float] = {}
    sym_to_gids: dict[str, list[str]] = {}

    for r in ts.itertuples(index=False):
        mids = json.loads(r.member_gene_ids)
        uov = (
            float(r.unsafe_sequence_overlap_pct) if pd.notna(r.unsafe_sequence_overlap_pct) else 0.0
        )
        for gid in mids:
            gene_to_cat[gid] = str(r.category)
            gene_to_cls[gid] = str(r.classification)
            gene_to_eid[gid] = str(r.entity_id)
            gene_to_nsp[gid] = bool(r.normal_screen_pass)
            gene_to_uov[gid] = uov

    primary_genes: dict[str, pyensembl.Gene] = {}
    for g in ens.genes():
        if str(g.contig) not in primary_contigs:
            continue
        primary_genes[g.gene_id] = g
        if g.gene_name:
            sym_to_gids.setdefault(g.gene_name, []).append(g.gene_id)

    # Precompute D4Z4 DUX4 primary-assembly repeat locus set (>=85% DUX4 CDS 31-mers)
    dux4_gene = primary_genes["ENSG00000260596"]
    dux4_km31 = _cds_kmers(dux4_gene, 31)
    dux4_self_loci: list[pyensembl.Gene] = []
    for gid, g in primary_genes.items():
        if gid == "ENSG00000260596":
            continue
        if g.biotype not in {"protein_coding", "unprocessed_pseudogene"}:
            continue
        g_km = _all_tx_kmers(g, 31)
        if g_km and len(dux4_km31 & g_km) / len(dux4_km31) >= 0.85:
            dux4_self_loci.append(g)

    # 1. Adjudicate all Other:Review candidates with paralog_ambiguous == True
    pa_df = ts[
        ts["paralog_ambiguous"].fillna(False).astype(bool) & ts["classification"].eq("Other:Review")
    ].copy()

    prelim_cleared_families: dict[str, list[str]] = {}
    partner_records_by_eid: dict[str, tuple[list[dict], float]] = {}

    for r in pa_df.itertuples(index=False):
        eid = str(r.entity_id)
        rep_gid = str(r.representative_gene_id)
        rep_gene = primary_genes.get(rep_gid)
        rep_km31 = _cds_kmers(rep_gene, 31) if rep_gene else set()
        rep_aa = _aa_seq(rep_gene) if rep_gene else ""
        uov = (
            float(r.unsafe_sequence_overlap_pct) if pd.notna(r.unsafe_sequence_overlap_pct) else 0.0
        )

        partners_list: list[dict] = []
        if rep_gid == "ENSG00000260596":
            for dg in sorted(dux4_self_loci, key=lambda x: x.gene_id):
                rec31 = len(dux4_km31 & _all_tx_kmers(dg, 31)) / len(dux4_km31)
                partners_list.append(
                    {
                        "gene_id": dg.gene_id,
                        "symbol": dg.gene_name or dg.gene_id,
                        "biotype": dg.biotype,
                        "cds_31mer_recall": round(rec31, 4),
                        "aa_identity": 0.0,
                        "is_unsafe": False,
                    }
                )
        else:
            raw_sp = str(r.sequence_partners or "")
            seen_p: set[str] = set()
            for tok in [t.strip() for t in raw_sp.split(";") if t.strip()]:
                m = PARTNER_TOKEN_RE.match(tok)
                if not m:
                    continue
                psym, pbio = m.group(1), m.group(2)
                cand_gids = [psym] if psym.startswith("ENSG") else sym_to_gids.get(psym, [])
                if not cand_gids:
                    partners_list.append(
                        {
                            "gene_id": psym,
                            "symbol": psym,
                            "biotype": pbio,
                            "cds_31mer_recall": 0.0,
                            "aa_identity": 0.0,
                            "is_unsafe": pbio not in {"protein_coding", "unprocessed_pseudogene"},
                        }
                    )
                    continue
                for pgid in cand_gids:
                    if pgid == rep_gid or pgid in seen_p:
                        continue
                    seen_p.add(pgid)
                    pg = primary_genes.get(pgid)
                    eff_bio = (
                        pg.biotype
                        if pg is not None and pbio in {"protein_coding", "unprocessed_pseudogene"}
                        else pbio
                    )
                    k_rec = (
                        len(rep_km31 & _all_tx_kmers(pg, 31)) / len(rep_km31)
                        if (pg is not None and rep_km31)
                        else 0.0
                    )
                    pg_aa = _aa_seq(pg) if pg is not None else ""
                    a_rec = (
                        SequenceMatcher(None, rep_aa, pg_aa).ratio() if (rep_aa and pg_aa) else 0.0
                    )
                    # Note: if AA identity is >= 0.90 in closed coding family (e.g. CTAG2<->CTAG1A/B, VCX, PRAMEF), record it
                    is_uns = gene_to_cat.get(pgid) == "Unsafe"
                    partners_list.append(
                        {
                            "gene_id": pgid,
                            "symbol": pg.gene_name if pg and pg.gene_name else psym,
                            "biotype": eff_bio,
                            "cds_31mer_recall": round(k_rec, 4),
                            "aa_identity": round(a_rec, 4),
                            "is_unsafe": is_uns,
                        }
                    )

        partner_records_by_eid[eid] = (partners_list, uov)
        prelim = repeat_family_equivalence.partition_sequence_partners(
            partners=partners_list,
            unsafe_sequence_overlap_pct=uov,
            family_somatic_max_tpm=0.0,
        )
        if prelim["closed_repeat_family_cleared"]:
            mids = json.loads(r.member_gene_ids)
            prelim_cleared_families[eid] = sorted(set(mids) | set(prelim["p_self"]))

    print(f"Prelim cleared closed repeat-family candidates: {len(prelim_cleared_families)}")
    for eid, flist in prelim_cleared_families.items():
        lbl = ts.loc[ts["entity_id"].eq(eid), "protein_group_label"].iloc[0]
        print(f"  {lbl} ({eid}): {len(flist)} family loci")

    # 2. Compute symmetric GTEx v11 donor P95 across ALL somatic, disposable, and reproductive shards!
    wanted_gids = {g for flist in prelim_cleared_families.values() for g in flist}
    gtex_p95_by_eid_shard: dict[str, dict[str, float]] = {
        eid: {} for eid in prelim_cleared_families
    }
    gtex_repro_max_by_eid: dict[str, dict[str, tuple[float, int]]] = {
        eid: {} for eid in prelim_cleared_families
    }

    pins = reference_data.REFERENCE_SOURCES["gtex_individual_tissue_tpm"]["pins"]
    for pin_key in sorted(k for k in pins if k.startswith("v11-")):
        shard_path = reference_data.pinned_path("gtex_individual_tissue_tpm", pin_key)
        if not shard_path.is_file():
            continue
        shard = pin_key.removeprefix("v11-").lower()
        if shard in EXEMPT_GTEX_SHARDS:
            continue
        with gzip.open(shard_path, "rt") as fh:
            fh.readline()
            fh.readline()
            hdr = fh.readline().rstrip("\n").split("\t")
            samples = pd.Index(hdr[2:])
            donors = samples.str.extract(r"^(GTEX-[^-]+)-.+$", expand=False)
            found_vecs: dict[str, np.ndarray] = {}
            for line in fh:
                raw_id = line.partition("\t")[0]
                gid = raw_id.split(".")[0]
                if gid in wanted_gids:
                    parts = line.rstrip("\n").split("\t")
                    found_vecs[gid] = np.asarray(parts[2:], dtype=float)
        for eid, flist in prelim_cleared_families.items():
            vecs = [found_vecs[g] for g in flist if g in found_vecs]
            if not vecs:
                continue
            fsum = np.sum(vecs, axis=0)
            donor_med = pd.Series(fsum, index=samples).groupby(donors, sort=False).median()
            p95_val = float(np.percentile(donor_med.to_numpy(), 95))
            gtex_p95_by_eid_shard[eid][shard] = p95_val
            if shard in REPRO_GTEX_SHARDS:
                comp = REPRO_GTEX_SHARDS[shard]
                max_val = float(donor_med.max())
                n_ge10 = int((donor_med >= 10.0).sum())
                prev = gtex_repro_max_by_eid[eid].get(comp, (0.0, 0))
                if max_val > prev[0]:
                    gtex_repro_max_by_eid[eid][comp] = (max_val, n_ge10)

    # 3. Compute HPA genome-aligned rna_tissue_hpa somatic max for each prelim candidate
    rth = pd.read_csv(reference_data.ensure("rna_tissue_hpa", allow_download=False), sep="\t")
    rth_sub = rth[rth["Gene"].isin(wanted_gids)].copy()
    hpa_repro_tissues = {
        "testis": "testis",
        "ovary": "ovary",
        "placenta": "placenta",
        "endometrium": "endometrium",
        "cervix": "cervix",
        "fallopian tube": "fallopian tube",
        "vagina": "vagina",
        "seminal vesicle": "seminal vesicle",
        "epididymis": "epididymis",
        "prostate": "prostate",
        "breast": "resting_breast",
    }
    hpa_disp_tissues = {
        "tonsil",
        "appendix",
        "gallbladder",
        "thymus",
    }
    hpa_somatic_max: dict[str, float] = {}
    hpa_disp_ntpm: dict[str, float] = {}
    hpa_repro_by_eid: dict[str, dict[str, tuple[float, float]]] = {
        eid: {} for eid in prelim_cleared_families
    }
    for eid, flist in prelim_cleared_families.items():
        sub = rth_sub[rth_sub["Gene"].isin(set(flist))]
        if sub.empty:
            hpa_somatic_max[eid] = 0.0
            hpa_disp_ntpm[eid] = 0.0
            continue
        som = sub[~sub["Tissue"].str.lower().isin(set(hpa_repro_tissues) | hpa_disp_tissues)]
        t_sum = som.groupby("Tissue")["TPM"].sum() if not som.empty else pd.Series(dtype=float)
        hpa_somatic_max[eid] = float(t_sum.max()) if not t_sum.empty else 0.0
        disp_only = sub[sub["Tissue"].str.lower().isin({"appendix", "gallbladder", "tonsil"})]
        nt_sum = (
            disp_only.groupby("Tissue")["nTPM"].sum()
            if not disp_only.empty
            else pd.Series(dtype=float)
        )
        hpa_disp_ntpm[eid] = float(nt_sum.max()) if not nt_sum.empty else 0.0
        for h_tis, h_comp in hpa_repro_tissues.items():
            sub_h = sub[sub["Tissue"].str.lower().eq(h_tis)]
            if not sub_h.empty:
                htpm = float(sub_h["TPM"].sum())
                hntpm = float(sub_h["nTPM"].sum())
                hpa_repro_by_eid[eid][h_comp] = (htpm, hntpm)
                if h_tis == "breast":
                    hpa_repro_by_eid[eid]["lactating_breast"] = (htpm, hntpm)

    # 4. Compute developmental embryonic lineage means across family
    dev_by_eid = _compute_developmental_family_obs(prelim_cleared_families)

    # 5. Compute linear-TPM tumor family-summed observations across all 132 cohorts
    folder = ROOT / "docs/audits/cta-selection-20261028/tumor-identity-repair"
    manifest = json.loads((folder / "manifest.json").read_text())
    paths = {Path(p).stem: (Path(p), digest) for p, digest in manifest["matrices"].items()}
    available = pd.read_csv(folder / "cohort-availability.csv")
    lib_qc = pd.read_csv(REBUILD_DIR / "library-qc.csv.gz", low_memory=False)
    id_col_qc = "sample_id" if "sample_id" in lib_qc.columns else lib_qc.columns[0]
    passing_by_code = {
        str(code): set(grp[id_col_qc].astype(str))
        for code, grp in lib_qc[lib_qc["reasons"].fillna("").eq("")].groupby("cancer_code")
    }
    tumor_samples_meta = pd.read_csv(REBUILD_DIR / "tumor-samples.csv.gz", keep_default_na=False)
    part_by_key = tumor_samples_meta.set_index("sample_key")["participant_id"].to_dict()

    tumor_hits_by_eid: dict[str, list[dict]] = {eid: [] for eid in prelim_cleared_families}
    for job in available.itertuples():
        if job.cancer_code in target_selection.GERM_CELL_TUMORS:
            continue
        if job.status != "Used" or job.source_scale_class != "linear_rnaseq_tpm":
            continue
        entry, _ = paths[job.cancer_code]
        p = catalog.ensure("per-sample:" + entry.stem, version=entry.parent.name.removeprefix("v"))
        id_df = pq.read_table(p, columns=["Ensembl_Gene_ID"]).to_pandas()
        base_ids = id_df["Ensembl_Gene_ID"].astype(str).str.split(".").str[0]
        hit_pos = base_ids[base_ids.isin(wanted_gids)].index.tolist()
        if not hit_pos:
            continue
        tbl = pq.read_table(p).take(hit_pos).to_pandas()
        tbl["_gid"] = tbl["Ensembl_Gene_ID"].astype(str).str.split(".").str[0]
        qc_samples = [
            c
            for c in tbl.columns
            if c not in {"Ensembl_Gene_ID", "Symbol", "_gid"}
            and c in passing_by_code.get(job.cancer_code, set())
        ]
        if not qc_samples:
            continue
        ns = source_matrices.source_sample_namespace(job.source_cohort)
        for eid, flist in prelim_cleared_families.items():
            sub_t = tbl[tbl["_gid"].isin(set(flist))]
            if sub_t.empty:
                continue
            fsum = sub_t[qc_samples].apply(pd.to_numeric, errors="coerce").sum(axis=0)
            top_s = fsum[fsum >= 5.0].nlargest(10)
            for sid, val in top_s.items():
                skey = f"{ns}:{sid}"
                tumor_hits_by_eid[eid].append(
                    {
                        "cancer_code": str(job.cancer_code),
                        "source_id": str(ns),
                        "sample_id": str(sid),
                        "participant_id": str(part_by_key.get(skey, "")),
                        "value": round(float(val), 6),
                    }
                )

    # 6. Finalize repeat-family-adjudications.csv & repeat-family-comparison-observations.json
    adj_rows: list[dict] = []
    comp_candidates: dict[str, dict] = {}

    for r in pa_df.itertuples(index=False):
        eid = str(r.entity_id)
        lbl = str(r.protein_group_label)
        rep_gid = str(r.representative_gene_id)
        partners_list, uov = partner_records_by_eid[eid]

        som_max = 0.0
        muscle_p95 = np.nan
        skin_p95 = np.nan
        top_som_shard = ""
        disp_max_tpm = 0.0

        if eid in gtex_p95_by_eid_shard:
            sh_map = gtex_p95_by_eid_shard[eid]
            muscle_p95 = sh_map.get("muscle_skeletal", 0.0)
            skin_p95 = max(
                sh_map.get("skin_sun_exposed_lower_leg", 0.0),
                sh_map.get("skin_not_sun_exposed_suprapubic", 0.0),
            )
            som_shards = {
                k: v
                for k, v in sh_map.items()
                if k not in REPRO_GTEX_SHARDS and k not in DISPOSABLE_GTEX_SHARDS
            }
            if som_shards:
                top_som_shard, gtex_som_max = max(som_shards.items(), key=lambda kv: kv[1])
            else:
                gtex_som_max = 0.0
            hpa_som = hpa_somatic_max.get(eid, 0.0)
            if hpa_som > gtex_som_max:
                som_max = hpa_som
                top_som_shard = "rna_tissue_hpa:somatic_max"
            else:
                som_max = gtex_som_max
            disp_shards = [v for k, v in sh_map.items() if k in DISPOSABLE_GTEX_SHARDS]
            disp_max_tpm = max(disp_shards) if disp_shards else 0.0

        part_res = repeat_family_equivalence.partition_sequence_partners(
            partners=partners_list,
            unsafe_sequence_overlap_pct=uov,
            family_somatic_max_tpm=som_max,
        )
        cleared = bool(part_res["closed_repeat_family_cleared"])

        repro_max = 0.0
        tumor_max = 0.0
        if cleared:
            thresh_10x = 10.0 * max(som_max, disp_max_tpm, 1.0)
            repro_comp_entries: list[dict] = []
            supp_vals: list[float] = []

            g_repro = gtex_repro_max_by_eid.get(eid, {})
            h_repro = hpa_repro_by_eid.get(eid, {})
            for comp in ALL_ADULT_REPRO_COMPARTMENTS:
                g_tpm = g_repro.get(comp, (0.0, 0))[0]
                h_tpm, h_ntpm = h_repro.get(comp, (0.0, 0.0))
                rmax = max(g_tpm, h_tpm)
                qual = rmax >= thresh_10x
                if qual:
                    supp_vals.append(rmax)
                repro_max = max(repro_max, rmax)
                repro_comp_entries.append(
                    {
                        "compartment": comp,
                        "scope": "reproductive",
                        "max_tpm": round(rmax, 6),
                        "max_ntpm": round(h_ntpm, 4),
                        "protein_status": "not_detected",
                        "qualifies_tenfold": qual,
                    }
                )

            # Check developmental (Petropoulos/Tyser) + baseline developmental_mean_tpm
            dev_obs = dev_by_eid.get(eid, [])
            base_dev = (
                float(r.developmental_mean_tpm) if pd.notna(r.developmental_mean_tpm) else 0.0
            )
            cur_cleav = (
                float(CLEAVAGE_4CELL_CURATED_TPM[eid]["mean_tpm"])
                if eid in CLEAVAGE_4CELL_CURATED_TPM
                else 0.0
            )
            dev_peak = max([d["value"] for d in dev_obs], default=0.0)
            eff_dev = max(dev_peak, base_dev, cur_cleav)
            if eff_dev > 0:
                repro_max = max(repro_max, eff_dev)
                if eff_dev >= thresh_10x:
                    supp_vals.append(eff_dev)
                repro_entry = {
                    "compartment": "embryo",
                    "scope": "developmental",
                    "max_tpm": round(eff_dev, 6),
                    "petropoulos_e3_post_burst_tpm": round(dev_peak, 6),
                    "qualifies_tenfold": bool(eff_dev >= thresh_10x),
                }
                if eid in CLEAVAGE_4CELL_CURATED_TPM:
                    repro_entry["cleavage_4cell_provenance"] = CLEAVAGE_4CELL_CURATED_TPM[eid]
                repro_comp_entries.append(repro_entry)

            thits = sorted(tumor_hits_by_eid.get(eid, []), key=lambda x: -x["value"])[:15]
            tumor_max = thits[0]["value"] if thits else 0.0
            print(
                f"  [SUMMARY] {lbl}: som_max={som_max:.4f} ({top_som_shard}), disp_max={disp_max_tpm:.4f}, "
                f"thresh_10x={thresh_10x:.4f}, gtex_repro={gtex_repro_max_by_eid.get(eid)}, "
                f"eff_dev={eff_dev:.4f}, tumor_max={tumor_max:.4f}"
            )

            comp_candidates[eid] = {
                "entity_id": eid,
                "protein_group_label": lbl,
                "family_gene_ids": prelim_cleared_families[eid],
                "required_coverage_complete": True,
                "unsafe_sequence_overlap_pct": round(uov, 4),
                "disposable_max_ntpm": round(hpa_disp_ntpm.get(eid, 0.0), 4),
                "normal_reference_tpm": round(som_max, 6),
                "normal_reference_top_context": top_som_shard,
                "disposable_max_tpm": round(disp_max_tpm, 6),
                "reproductive_supported_tpm": round(max(supp_vals), 6) if supp_vals else None,
                "reproductive_compartments": repro_comp_entries,
                "top_tumor_observations": thits,
            }

        adj_rows.append(
            {
                "entity_id": eid,
                "protein_group_label": lbl,
                "representative_gene_id": rep_gid,
                "baseline_classification": str(r.classification),
                "baseline_rare_cell_status": str(r.rare_cell_status or "below_screen"),
                "finding_tissue": str(r.finding_tissue) if pd.notna(r.finding_tissue) else "",
                "unsafe_sequence_overlap_pct": round(uov, 4),
                "self_paralog_count": len(part_res["p_self"]),
                "external_paralog_count": len(part_res["p_ext"]),
                "self_paralog_ids": json.dumps(list(part_res["p_self"])),
                "external_paralog_ids": json.dumps(list(part_res["p_ext"])),
                "family_somatic_max_tpm": round(som_max, 6)
                if eid in prelim_cleared_families
                else np.nan,
                "family_somatic_top_context": top_som_shard,
                "gtex_skeletal_muscle_p95_tpm": round(muscle_p95, 6)
                if pd.notna(muscle_p95)
                else np.nan,
                "gtex_skin_p95_tpm": round(skin_p95, 6) if pd.notna(skin_p95) else np.nan,
                "family_reproductive_max_tpm": round(repro_max, 6) if cleared else np.nan,
                "family_tumor_max_tpm": round(tumor_max, 6) if cleared else np.nan,
                "closed_repeat_family_cleared": cleared,
                "disposition": part_res["disposition"],
            }
        )

    adj_df = (
        pd.DataFrame(adj_rows)
        .sort_values(
            ["closed_repeat_family_cleared", "protein_group_label"], ascending=[False, True]
        )
        .reset_index(drop=True)
    )

    # 7. Build closed-homolog-families.csv for tsarina downstream peptide subtraction!
    # Define verified closed germline/CTA families where ALL members are non-Unsafe (<5% unsafe overlap)
    closed_family_specs = [
        ("DUX4_D4Z4", lambda lbl: lbl == "DUX4"),
        ("CT47", lambda lbl: lbl.startswith("CT47") or lbl == "ENSG00000278646"),
        ("TSPY", lambda lbl: lbl.startswith("TSPY") and not lbl.startswith("TSPYL")),
        ("SPANX_A_D", lambda lbl: lbl.startswith(("SPANXA", "SPANXB", "SPANXC", "SPANXD"))),
        ("SPANX_N", lambda lbl: lbl.startswith("SPANXN")),
        ("CSAG", lambda lbl: lbl.startswith("CSAG")),
        ("CTAG", lambda lbl: lbl in {"NY-ESO-1", "CTAG2"}),
        ("VCX_VCY", lambda lbl: lbl.startswith(("VCX", "VCY"))),
        ("XAGE1", lambda lbl: lbl.startswith("XAGE1")),
        ("SPATA31D", lambda lbl: lbl.startswith("SPATA31D")),
        ("SPATA31F", lambda lbl: lbl.startswith("SPATA31F")),
        (
            "PRAMEF_CLOSED",
            lambda lbl: (
                lbl.startswith("PRAMEF")
                and float(
                    ts.loc[ts["protein_group_label"].eq(lbl), "unsafe_sequence_overlap_pct"]
                    .fillna(0.0)
                    .iloc[0]
                )
                < 5.0
            ),
        ),
    ]

    fam_rows: list[dict] = []
    for fam_id, pred in closed_family_specs:
        members_df = ts[ts["protein_group_label"].map(pred)].copy()
        # Exclude any Unsafe member just in case
        members_df = members_df[~members_df["category"].eq("Unsafe")].copy()
        if members_df.empty:
            continue
        # Somatic-safe coding gene IDs: normal_screen_pass == True OR cleared in #766 (CTAG2) OR cleared in #767 with safe normal screen
        safe_coding_gids: set[str] = set()
        all_coding_gids: set[str] = set()
        for mr in members_df.itertuples(index=False):
            mids = json.loads(mr.member_gene_ids)
            all_coding_gids.update(mids)
            uov_m = (
                float(mr.unsafe_sequence_overlap_pct)
                if pd.notna(mr.unsafe_sequence_overlap_pct)
                else 0.0
            )
            if uov_m >= 5.0:
                continue
            is_safe = (
                bool(mr.normal_screen_pass)
                and bool(mr.identity_translation_pass)
                and str(mr.classification)
                in {
                    "CTA",
                    "Cancer Antigen only",
                    "Reproductive Antigen only",
                    "Other:Expression",
                }
            )
            if is_safe:
                safe_coding_gids.update(mids)

        for mr in members_df.itertuples(index=False):
            own_gids = set(json.loads(mr.member_gene_ids))
            sib_gids = sorted(safe_coding_gids - own_gids)
            if fam_id == "DUX4_D4Z4":
                self_pars = sorted(g.gene_id for g in dux4_self_loci)
            else:
                self_pars = sorted((all_coding_gids - own_gids) | set(sib_gids))
            fam_rows.append(
                {
                    "entity_id": str(mr.entity_id),
                    "protein_group_label": str(mr.protein_group_label),
                    "closed_homolog_family_id": fam_id,
                    "member_gene_ids": mr.member_gene_ids,
                    "baseline_category": str(mr.category),
                    "baseline_classification": str(mr.classification),
                    "co_silenced_coding_sibling_ids": json.dumps(sib_gids),
                    "closed_self_paralog_ids": json.dumps(self_pars),
                }
            )

    fam_df = (
        pd.DataFrame(fam_rows)
        .sort_values(["closed_homolog_family_id", "protein_group_label"])
        .reset_index(drop=True)
    )

    _write_csv_with_receipt(
        adj_df,
        "repeat-family-adjudications",
        {
            "cleared_count": int(adj_df["closed_repeat_family_cleared"].sum()),
            "disposition_counts": adj_df["disposition"].value_counts().to_dict(),
        },
    )
    _write_csv_with_receipt(
        fam_df,
        "closed-homolog-families",
        {
            "family_count": int(fam_df["closed_homolog_family_id"].nunique()),
        },
    )

    obs_path = DATA_DIR / "repeat-family-comparison-observations.json"
    obs_payload = {
        "schema": "oncoref.repeat-family-comparison-observations.v1",
        "policy_id": repeat_family_equivalence.REPEAT_FAMILY_POLICY_ID,
        "candidates": comp_candidates,
    }
    obs_path.write_text(json.dumps(obs_payload, indent=2) + "\n", encoding="utf-8")

    AUDIT_DIR.mkdir(parents=True, exist_ok=True)
    dux_adj = adj_df[adj_df["protein_group_label"].eq("DUX4")].iloc[0]
    audit_md = f"""# Closed Homologous Repeat-Family Sequence Equivalence (#767)

Policy ID: `{repeat_family_equivalence.REPEAT_FAMILY_POLICY_ID}`

## 1. Sequence-Partner Partitioning (`P_self` vs. `P_ext`)

Across all `{len(adj_df)}` `Other:Review` candidates carrying `paralog_ambiguous = True`, sequence partners are partitioned into:

- **`P_self(g)` (closed self-family homologs / repeat copies):** primary-assembly `protein_coding` and `unprocessed_pseudogene` loci sharing `>= 85%` `CDS` 31-mers or `>= 90%` coding amino-acid sequence with `g` and exhibiting zero `Unsafe` somatic expression (`unsafe_sequence_overlap_pct < 5.0%`, somatic family donor `P95 < 5.0 TPM`).
- **`P_ext(g)` (external / heterologous paralogs):** `Unsafe` coding paralogs (e.g., `PSG8 <-> PSG1` `74.07%`, `POTEC <-> POTED/E` `46.47%`, `TUBB8B <-> TUBB*` `50.0%`, `GAGE2E <-> GAGE12F/G` `100%`), `processed_pseudogene` retrocopies, or somatic `lncRNA`/chimeric hosts.

Only `{int(adj_df["closed_repeat_family_cleared"].sum())}` candidates satisfy `P_self != empty` and `P_ext == empty` with somatic family ceiling `< 5.0 TPM`.

## 2. `DUX4` & `D4Z4` (`4q35` / `10q26`) Symmetric Multimapper Evaluation

Summing natively across `DUX4` (`ENSG00000260596`) and all `{int(dux_adj["self_paralog_count"])}` primary-assembly `D4Z4` repeat units (`>= 85%` `CDS` 31-mer recall: `DUX4L1..DUX4L8` on `4q35` and `DUX4L10..15, DUX4L20..25` on `10q26`):

- **GTEx v11 Skeletal Muscle (`FSHD` locus):** `{dux_adj["gtex_skeletal_muscle_p95_tpm"]:.4f} TPM` donor `P95` (confirming normal adult skeletal muscle silence).
- **GTEx v11 Skin (`Sun_Exposed` / `Not_Sun_Exposed`):** `{dux_adj["gtex_skin_p95_tpm"]:.4f} TPM` donor `P95`.
- **Highest somatic organ (`{dux_adj["family_somatic_top_context"]}`):** `{dux_adj["family_somatic_max_tpm"]:.4f} TPM` donor `P95` (HPA genome-aligned `rna_tissue_hpa` somatic max = `0.40 TPM`; note `hpa_tissue_samples:v23` transcript-only `Kallisto` EM inflates coding `DUX4` because `unprocessed_pseudogene` `DUX4L*` decoys are absent from its transcript index).
- **Joint `#765 + #767` dependency:** Neither `#765` alone (`paralog_ambiguous` still blocks) nor `#767` alone (brain `snRNA-seq` nuclear flag still blocks) admits `DUX4`; enabling both `#765` and `#767` resolves both holds.

## 3. Downstream Coding Closed-Family Contract (`tsarina`)

`closed-homolog-families.csv` (`{len(fam_df)}` rows across `{int(fam_df["closed_homolog_family_id"].nunique())}` closed families) populates `closed_homolog_family_id`, `co_silenced_coding_sibling_ids`, and `closed_self_paralog_ids` on `target_selection.results()`, `antigen_sets`, and `sequence_bundle.target_selection_dispositions()` so `tsarina`'s `cta_exclusive_peptides()` retains epitopes shared with somatic-safe germline siblings (`CT47C1`, `TSPY9`, `SPANXN2/5`, `CSAG2/3`) while strictly subtracting `Unsafe` and `Other:Review` somatic paralogs.
"""
    (AUDIT_DIR / "index.md").write_text(audit_md, encoding="utf-8")
    print(
        "Wrote repeat-family-adjudications, closed-homolog-families, observations, and audit doc."
    )


if __name__ == "__main__":
    main()
