#!/usr/bin/env python3
"""Build holistic GTEx v11 deconvolution & sample-swap/splash audit tables (#778).

Scans all 52 pinned GTEx v11 per-sample tissue shards, infers splash-robust donor
chromosomal sex (``XX`` vs ``XY``), integrates ``SMNABTCH``/``SMGEBTCH`` cross-donor
batch matches from ``docs/audits/gtex-cross-tissue-signals/signals.csv.gz``, selects
lineage-disjoint anchor genes across 20+ organ archetypes, solves a joint 52-shard
non-negative least squares (NNLS) mixture for native retention, and writes
SHA-256-verified packaged tables and ``docs/audits/gtex-holistic-deconvolution-20261124/index.md``.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from oncoref import gtex_deconvolution, reference_data  # noqa: E402

_SIGNAL_PANEL_TO_SITES: dict[str, tuple[str, ...]] = {
    "Adipose": ("adipose_subcutaneous", "adipose_visceral_omentum"),
    "Adrenal": ("adrenal_gland",),
    "Artery": ("artery_aorta", "artery_coronary", "artery_tibial"),
    "Brain": (
        "brain_amygdala",
        "brain_anterior_cingulate_cortex_ba24",
        "brain_caudate_basal_ganglia",
        "brain_cerebellar_hemisphere",
        "brain_cerebellum",
        "brain_cortex",
        "brain_frontal_cortex_ba9",
        "brain_hippocampus",
        "brain_hypothalamus",
        "brain_nucleus_accumbens_basal_ganglia",
        "brain_putamen_basal_ganglia",
        "brain_spinal_cord_cervical_c-1",
        "brain_substantia_nigra",
    ),
    "Breast": ("breast_mammary_tissue",),
    "Colon": ("colon_sigmoid", "colon_transverse"),
    "Esophagus": (
        "esophagus_gastroesophageal_junction",
        "esophagus_mucosa",
        "esophagus_muscularis",
    ),
    "Heart": ("heart_atrial_appendage", "heart_left_ventricle"),
    "Kidney": ("kidney_cortex", "kidney_medulla"),
    "Liver": ("liver",),
    "Lung": ("lung",),
    "Muscle": ("muscle_skeletal",),
    "Nerve": ("nerve_tibial",),
    "Ovary": ("ovary",),
    "Pancreas": ("pancreas",),
    "Pituitary": ("pituitary",),
    "Prostate": ("prostate",),
    "Salivary gland": ("minor_salivary_gland",),
    "Skin": ("skin_not_sun_exposed_suprapubic", "skin_sun_exposed_lower_leg"),
    "Small intestine": ("small_intestine_terminal_ileum",),
    "Spleen": ("spleen",),
    "Stomach": ("stomach",),
    "Testis": ("testis",),
    "Thyroid": ("thyroid",),
    "Whole blood": ("whole_blood",),
}


def _sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _load_cross_donor_batch_pairs() -> set[tuple[str, str]]:
    """Load verified SMNABTCH / SMGEBTCH cross-donor batch matches mapped to exact shard sites."""
    sig_path = ROOT / "docs/audits/gtex-cross-tissue-signals/signals.csv.gz"
    if not sig_path.is_file():
        return set()
    sig = pd.read_csv(sig_path)
    matched = sig[
        sig["different_donor"].fillna(False).astype(bool)
        & (
            sig["shared_extraction"].fillna(False).astype(bool)
            | sig["shared_expression"].fillna(False).astype(bool)
        )
    ]
    pairs: set[tuple[str, str]] = set()
    for r in matched.itertuples(index=False):
        sites = _SIGNAL_PANEL_TO_SITES.get(str(r.source_panel), (str(r.source_panel).lower(),))
        for st in sites:
            pairs.add((str(r.sample_id), st))
    return pairs


def main() -> None:
    med_table = pd.read_csv(
        reference_data.pinned_path("gtex_tissue_rna", "v11"), sep="\t", skiprows=2
    )
    med_table["gene_id"] = med_table["Name"].astype(str).str.split(".").str[0]
    med_table = med_table.drop_duplicates("gene_id").set_index("gene_id")
    sym_map = med_table["Description"].astype(str).to_dict()
    tcols = [c for c in med_table.columns if c not in ("Name", "Description")]

    # Include candidate genes >= 45 TPM in at least one GTEx column + Y-chromosome sentinels
    cand_mask = med_table[tcols].astype(float).max(axis=1) >= 45.0
    cand_gids = sorted(
        set(med_table.index[cand_mask]) | set(gtex_deconvolution.Y_CHROMOSOME_SENTINEL_GENES)
    )
    cand_idx = {g.encode("ascii"): i for i, g in enumerate(cand_gids)}

    pins = sorted(reference_data.REFERENCE_SOURCES["gtex_individual_tissue_tpm"]["pins"])
    sites = [pk.removeprefix("v11-") for pk in pins]

    shard_med = np.zeros((len(cand_gids), len(sites)), dtype=np.float64)
    shard_p95 = np.zeros((len(cand_gids), len(sites)), dtype=np.float64)
    shard_samples: dict[str, list[str]] = {}
    shard_cand_mats: dict[str, np.ndarray] = {}

    for k, (pk, site) in enumerate(zip(pins, sites)):
        p = reference_data.pinned_path("gtex_individual_tissue_tpm", pk)
        with gzip.open(p, "rb") as fh:
            fh.readline()
            fh.readline()
            samps = fh.readline().decode("utf-8").rstrip("\n").split("\t")[2:]
            mat = np.zeros((len(cand_gids), len(samps)), dtype=np.float64)
            for line in fh:
                rn, _, rest = line.partition(b"\t")
                idx = cand_idx.get(rn.split(b".", 1)[0])
                if idx is not None:
                    mat[idx] = np.asarray(rest.partition(b"\t")[2].split(), dtype=np.float64)
        shard_samples[site] = samps
        shard_cand_mats[site] = mat
        shard_med[:, k] = np.median(mat, axis=1)
        shard_p95[:, k] = np.percentile(mat, 95.0, axis=1)

    med_df = pd.DataFrame(shard_med, index=cand_gids, columns=sites)
    p95_df = pd.DataFrame(shard_p95, index=cand_gids, columns=sites)
    sample_counts = {site: len(shard_samples[site]) for site in sites}
    anchors = gtex_deconvolution.select_lineage_disjoint_anchors(
        med_df,
        p95_df,
        shard_sample_counts=sample_counts,
        symbols=sym_map,
        min_source_median_tpm=75.0,
        min_episodic_source_p95_tpm=250.0,
        max_off_lineage_p95_tpm=5.0,
        max_off_lineage_ratio=0.01,
        max_off_lineage_shard_fraction=0.12,
        hard_max_off_lineage_ratio=0.10,
        max_anchors_per_source=10,
        min_anchors_per_source=3,
    )

    basis_gids = gtex_deconvolution.select_nnls_basis_genes(med_df, anchors, top_per_shard=8)
    anchor_gids = sorted(anchors["gene_id"].unique())
    gid_pos = {g: i for i, g in enumerate(cand_gids)}
    y_gids = [g for g in gtex_deconvolution.Y_CHROMOSOME_SENTINEL_GENES if g in gid_pos]
    anchor_row_pos = [gid_pos[g] for g in anchor_gids]
    basis_row_pos = [gid_pos[g] for g in basis_gids]
    y_row_pos = [gid_pos[g] for g in y_gids]

    all_sample_ids: list[str] = []
    all_sample_sites: dict[str, str] = {}
    anchor_blocks: list[np.ndarray] = []
    basis_blocks: list[np.ndarray] = []
    y_blocks: list[np.ndarray] = []
    for site in sites:
        samps = shard_samples[site]
        all_sample_ids.extend(samps)
        for s in samps:
            all_sample_sites[s] = site
        anchor_blocks.append(shard_cand_mats[site][anchor_row_pos, :])
        basis_blocks.append(shard_cand_mats[site][basis_row_pos, :])
        y_blocks.append(shard_cand_mats[site][y_row_pos, :])

    anchor_tpm = pd.DataFrame(np.hstack(anchor_blocks), index=anchor_gids, columns=all_sample_ids)
    basis_sample_tpm = pd.DataFrame(
        np.hstack(basis_blocks), index=basis_gids, columns=all_sample_ids
    )
    y_tpm = pd.DataFrame(np.hstack(y_blocks), index=y_gids, columns=all_sample_ids)
    sample_sex, donor_sex = gtex_deconvolution.infer_chromosomal_sex(y_tpm)
    batch_pairs = _load_cross_donor_batch_pairs()

    events_df, sample_summary = gtex_deconvolution.deconvolve_gtex_samples(
        anchor_tpm,
        pd.Series(all_sample_sites),
        anchors,
        basis_medians=med_df.loc[basis_gids],
        basis_sample_tpm=basis_sample_tpm,
        donor_sex=donor_sex,
        sample_sex=sample_sex,
        cross_donor_batch_pairs=batch_pairs,
    )

    splash_events = (
        events_df[events_df["regime"].eq("tube_swap") | events_df["alpha_hat"].ge(5e-3)]
        .sort_values(["regime", "alpha_hat", "sample_id"], ascending=[True, False, True])
        .reset_index(drop=True)
    )
    swaps_df = (
        events_df[events_df["regime"].eq("tube_swap")]
        .sort_values(["alpha_hat", "sample_id"], ascending=[False, True])
        .reset_index(drop=True)
    )

    tissue_rows: list[dict] = []
    for site in sites:
        s_sum = sample_summary[sample_summary["site"].eq(site)]
        ev_s = events_df[events_df["target_site"].eq(site)]
        n_samps = len(s_sum)
        n_donors = int(s_sum["donor_id"].nunique())
        n_swaps = int(s_sum["quarantined_swap"].sum())
        splash_samps = set(ev_s.loc[ev_s["regime"].eq("plate_pipetting_splash"), "sample_id"])
        carry_samps = set(ev_s.loc[ev_s["regime"].eq("same_donor_carryover"), "sample_id"])
        hop_samps = set(ev_s.loc[ev_s["regime"].eq("flowcell_index_hopping"), "sample_id"]) | set(
            s_sum.loc[s_sum["plate_hop_beta"].gt(0.0), "sample_id"]
        )
        bound_samps = set(ev_s.loc[ev_s["regime"].eq("anatomical_boundary"), "sample_id"])
        top_src = (
            ev_s.loc[ev_s["regime"].ne("anatomical_boundary"), "source_site"]
            .value_counts()
            .head(3)
            .index.tolist()
        )
        tissue_rows.append(
            {
                "site": site,
                "lineage_family": gtex_deconvolution.organ_lineage_family(site),
                "total_samples": n_samps,
                "total_donors": n_donors,
                "quarantined_swap_samples": n_swaps,
                "plate_splash_samples": len(splash_samps),
                "same_donor_carryover_samples": len(carry_samps),
                "index_hopping_samples": len(hop_samps),
                "anatomical_boundary_samples": len(bound_samps),
                "median_native_retention": round(float(s_sum["native_retention"].median()), 6),
                "median_purity_weight": round(float(s_sum["purity_weight"].median()), 6),
                "p05_purity_weight": round(float(np.percentile(s_sum["purity_weight"], 5.0)), 6),
                "top_foreign_sources": ";".join(top_src),
            }
        )
    tissue_summary_df = pd.DataFrame(tissue_rows)

    data_dir = ROOT / "oncoref/data"
    anchors_path = data_dir / "gtex-deconvolution-anchors.csv"
    tissue_path = data_dir / "gtex-deconvolution-tissue-summary.csv"
    swaps_path = data_dir / "gtex-quarantined-swaps.csv"
    splash_path = data_dir / "gtex-sample-splash-events.csv.gz"
    dsex_path = data_dir / "gtex-donor-sex.csv"
    receipt_path = data_dir / "gtex-deconvolution-receipt.json"

    y_sum_by_sample = y_tpm.astype(float).sum(axis=0)
    d_ids_all = [gtex_deconvolution.parse_gtex_sample_barcode(s).donor_id for s in all_sample_ids]
    donor_y_med = y_sum_by_sample.groupby(d_ids_all, sort=True).median()
    donor_shard_n = y_sum_by_sample.groupby(d_ids_all, sort=True).size()
    dsex_df = pd.DataFrame(
        [
            {
                "donor_id": str(did),
                "chromosomal_sex": str(donor_sex[did]),
                "median_y_sentinel_tpm": round(float(donor_y_med[did]), 4),
                "sequenced_shards": int(donor_shard_n[did]),
            }
            for did in donor_y_med.index
        ]
    )
    dsex_df.to_csv(dsex_path, index=False, lineterminator="\n")

    anchors_fmt = anchors.copy()
    for fcol in ("source_median_tpm", "source_p95_tpm", "max_off_lineage_p95_tpm"):
        anchors_fmt[fcol] = anchors_fmt[fcol].round(4)
    anchors_fmt.to_csv(anchors_path, index=False, lineterminator="\n")
    tissue_summary_df.to_csv(tissue_path, index=False, lineterminator="\n")

    for df_out in (swaps_df, splash_events):
        for fcol in (
            "alpha_hat",
            "nnls_weight",
            "alpha_lower_q25",
            "log10_iqr",
            "native_retention",
        ):
            if fcol in df_out.columns:
                df_out[fcol] = df_out[fcol].round(6)
    swaps_df.to_csv(swaps_path, index=False, lineterminator="\n")
    splash_events.to_csv(
        splash_path,
        index=False,
        compression={"method": "gzip", "mtime": 0},
        lineterminator="\n",
    )

    files_sha = {
        p.name: _sha256_file(p)
        for p in (anchors_path, tissue_path, swaps_path, splash_path, dsex_path)
    }
    receipt = {
        "policy_id": gtex_deconvolution.GTEX_DECONVOLUTION_POLICY_ID,
        "gtex_release": "v11",
        "shards_evaluated": len(sites),
        "samples_evaluated": len(sample_summary),
        "unique_donors": int(sample_summary["donor_id"].nunique()),
        "donor_sex_counts": {k: int(v) for k, v in sorted(donor_sex.value_counts().items())},
        "unique_sm_plate_prefixes": int(
            sample_summary.loc[
                sample_summary["sm_plate_prefix"].ne(""), "sm_plate_prefix"
            ].nunique()
        ),
        "nnls_basis_genes": len(basis_gids),
        "anchor_genes": len(anchors),
        "anchor_source_archetypes": int(anchors["source_site"].nunique()),
        "anchor_source_sites": sorted(anchors["source_site"].unique().tolist()),
        "regime_event_counts": {
            k: int(v) for k, v in sorted(events_df["regime"].value_counts().items())
        },
        "quarantined_swap_samples": int(sample_summary["quarantined_swap"].sum()),
        "files_sha256": files_sha,
    }
    receipt_path.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    audit_dir = ROOT / "docs/audits/gtex-holistic-deconvolution-20261124"
    audit_dir.mkdir(parents=True, exist_ok=True)
    swap_md_rows = [
        f"| `{r.sample_id}` | `{r.donor_sex}` | `{r.target_site}` | `{r.source_site}` | "
        f"`SM-{r.sm_plate_prefix}` | `{r.alpha_hat:.3f}` | `{r.native_retention:.3f}` | `{r.anchor_symbols}` |"
        for r in swaps_df.head(15).itertuples(index=False)
    ]
    splash_ex = splash_events[splash_events["regime"].eq("plate_pipetting_splash")].head(15)
    splash_md_rows = [
        f"| `{r.sample_id}` | `{r.donor_sex}` | `{r.target_site}` | `{r.source_site}` | "
        f"`SM-{r.sm_plate_prefix}` | `{r.co_plated_on_sm_batch}` | `{r.cross_donor_batch_matched}` | "
        f"`{100.0 * r.alpha_hat:.2f}%` | `{r.anchor_symbols}` |"
        for r in splash_ex.itertuples(index=False)
    ]

    audit_md = "\n".join(
        [
            "# Holistic GTEx Tissue Deconvolution & Donor Quantile Estimation (`#778`)",
            "",
            f"**Policy ID:** `{gtex_deconvolution.GTEX_DECONVOLUTION_POLICY_ID}`  ",
            f"**Evaluated Scope:** `{len(sites)}` GTEx v11 shards · `{len(sample_summary):,}` RNA-seq libraries · "
            f"`{receipt['unique_donors']}` `GTEX-XXXXX` donors · `{receipt['unique_sm_plate_prefixes']}` Broad `SM-XXXX` 96-well prep-plate batches · "
            f"`{receipt['anchor_source_archetypes']}` disjoint organ archetypes (`{receipt['anchor_genes']}` anchor genes, `{receipt['nnls_basis_genes']}` 52-shard NNLS basis genes)",
            "",
            "## Four Contamination Regimes Separated from Histological Lineage & Surgical Boundaries",
            "",
            "1. **Patterned-flowcell index hopping (`flowcell_index_hopping`, `2×10⁻⁵ – 5×10⁻³`):** Pervasive multiplexed pool bleed across co-sequenced tissues; conservatively subtracted via `alpha_lower_q25 * source_median_tpm` when bleed accounts for `≥25%` of observed TPM.",
            "2. **Cross-donor 96-well plate / extraction-batch pipetting splash (`plate_pipetting_splash`, `≥0.5%`):** Requires cross-donor `SM-XXXX` co-plating (without same-donor co-plating), verified cross-donor `SMNABTCH`/`SMGEBTCH` batch sharing, or chromosomal sex impossibility (`XX` recipient carrying `testis`/`prostate`, or `XY` recipient carrying `ovary`/`uterus`). Heavy splashes (`≥0.5%` from expressing sources) are dropped from the gene's donor denominator rather than zero-filled.",
            "3. **Same-donor prosector carryover vs. surgical/anatomical boundaries (`same_donor_carryover` vs. `anatomical_boundary`):** Separates non-adjacent autopsy bench carryover within `GTEX-XXXXX` from true contiguous surgical planes (`prostate → bladder`, `epidermis → subcutaneous adipose`, `gastrocnemius → tibial nerve`, `pectoralis → breast`, `cardia → stomach`), which are **never** subtracted from somatic safety `N`.",
            "4. **Gross biopsy tube swaps & cross-sex chimeras (`tube_swap`, `≥50%` foreign chimerism, `XX`/`XY` sample-donor sex discordance at `≥20%`, or `≥20%` foreign + joint 52-shard NNLS `native_retention < 0.25`):** Quarantined globally across all genes.",
            "",
            "## Quarantined Foreign-Organ Swaps",
            "",
            "| Sample ID | Donor Sex | Labeled Shard | Detected Organ | Plate Batch | Foreign $\\hat{\\alpha}$ | NNLS Native Retention | Detected Anchors |",
            "| --- | --- | --- | --- | --- | ---: | ---: | --- |",
            *(swap_md_rows or ["| *(none)* | - | - | - | - | - | - | - |"]),
            "",
            "## Top Cross-Donor Wet-Lab Pipetting Splashes (`plate_pipetting_splash`)",
            "",
            "| Sample ID | Donor Sex | Target Shard | Splashed Organ | Plate Batch | SM-Plate Peer | Batch-Matched | Splash % | Detected Anchors |",
            "| --- | --- | --- | --- | --- | --- | --- | ---: | --- |",
            *splash_md_rows,
            "",
            "## Programmatic Inspection",
            "",
            "```python",
            "from oncoref import gtex_deconvolution",
            "",
            "receipt = gtex_deconvolution.gtex_deconvolution_receipt()",
            "swaps = gtex_deconvolution.gtex_quarantined_swaps()",
            "summary = gtex_deconvolution.gtex_deconvolution_tissue_summary()",
            "print(receipt['anchor_source_archetypes'], receipt['quarantined_swap_samples'], len(summary))",
            "```",
            "",
        ]
    )
    (audit_dir / "index.md").write_text(audit_md, encoding="utf-8")
    print(json.dumps(receipt, indent=2))


if __name__ == "__main__":
    main()
