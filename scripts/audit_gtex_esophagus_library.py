"""Trace the GTEx muscularis outlier and tumor evidence for the 98 CTA to RA-only changes."""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from cta_report_common import sha256, write_csv
from gtex_sample_evidence import SOURCE, sample_matrix

from oncoref import reference_data

ROOT = Path(__file__).resolve().parents[1]
AUDITS = ROOT / "docs/audits"
REBUILD = AUDITS / "target-selection-20261110/rebuild"
OUT = AUDITS / "gtex-esophagus-muscularis"
SAMPLE = "GTEX-1I1GR-1426-SM-CY8H8"
DONOR = "GTEX-1I1GR"
MUSCULARIS = "Esophagus_Muscularis"
MUSCLE_GENES = ["ACTG2", "MYH11", "CNN1", "TAGLN", "ACTA2"]
QC_GENES = ["PRM1", "PRM2", "TNP1", "PHF7", "LELP1", "SMCP"]


def build(annotations, samples_dir):
    if sha256(annotations) != "98c246b3855edf62cb050982a7ca50f3069d010624baea47ef84a9465900ca4b":
        raise ValueError("GTEx v11 sample annotations checksum mismatch")
    OUT.mkdir(exist_ok=True)
    comparison = pd.read_csv(AUDITS / "gtex-reproductive-witnesses/catalog-comparison.csv.gz")
    lost = comparison.query(
        "prior_classification == 'CTA' and classification == 'Reproductive Antigen only'"
    )
    catalog = pd.read_csv(
        ROOT / "oncoref/data/target-selection.csv.gz", low_memory=False
    ).set_index("entity_id")
    cases = catalog.loc[lost.entity_id].copy()
    assert len(cases) == 98
    qc = pd.read_csv(REBUILD / "gtex-library-qc.csv.gz", keep_default_na=False).set_index(
        "sample_id"
    )
    attr = pd.read_csv(annotations, sep="\t", low_memory=False).set_index("SAMPID")
    subject = qc.loc[SAMPLE]
    assert (
        subject.source_tissue == MUSCULARIS
        and attr.loc[SAMPLE, "SMTSD"] == "Esophagus - Muscularis"
    )
    native_tissues = sorted(qc.source_tissue.unique())
    median_path = ROOT / "tmp/cta-source-consolidation/gtex-v11-median-tpm.gct.gz"
    med = pd.read_csv(median_path, sep="\t", skiprows=2).set_index("Name")
    selected = pd.read_csv(
        AUDITS / "cta-selection-20261028/alternative-evidence/gtex_tissue_rna-row-decisions.csv.gz"
    ).dropna(subset="row_position")
    canonical = pd.Series(
        selected.gene_id.to_numpy(), index=med.index[selected.row_position.astype(int)]
    )
    # Broad tissue identity on genes detected in at least one tissue median.
    expressed = med[native_tissues].max(axis=1).ge(1)
    # Select the independent testis panel using published tissue medians only.
    other_median = med[[s for s in native_tissues if s != "Testis"]].max(axis=1)
    enriched = (
        med.Testis.ge(10)
        & med.Testis.ge(20 * other_median.clip(lower=1))
        & ~med.Description.isin(QC_GENES)
    )
    wanted = med.index[expressed | enriched | med.Description.isin(QC_GENES + MUSCLE_GENES)].union(
        canonical.index[
            canonical.isin({g for ids in cases.member_gene_ids for g in json.loads(ids)})
        ]
    )
    pins = reference_data.REFERENCE_SOURCES[SOURCE]["pins"]

    def matrix(tissue, exact):
        slug = tissue.lower()
        path = samples_dir / (slug + ".gct.gz")
        if sha256(path) != pins["v11-" + slug]["sha256"]:
            raise ValueError(f"GTEx checksum mismatch for {slug}")
        return sample_matrix(path, {g.encode(): g for g in exact})[0]

    eso = matrix(MUSCULARIS, wanted)
    others = eso.drop(columns=SAMPLE)
    symbols = med.Description
    panel = med.index[enriched]
    independent = eso.loc[panel]
    shared = pd.concat([med[native_tissues], eso[SAMPLE]], axis=1).loc[expressed]
    shared[MUSCULARIS] = others.loc[shared.index].median(axis=1)
    correl = shared.corr(method="spearman")[SAMPLE].drop(SAMPLE).sort_values(ascending=False)
    write_csv(
        correl.rename("spearman").rename_axis("tissue").reset_index(),
        OUT / "tissue-correlations.csv",
    )
    testis = matrix("Testis", panel.union(med.index[symbols.isin(QC_GENES + MUSCLE_GENES)]))
    phenotype = pd.DataFrame(
        {
            "symbol": symbols,
            "suspect_tpm": eso[SAMPLE],
            "peer_median_tpm": others.median(axis=1),
            "peer_max_tpm": others.max(axis=1),
            "testis_median_tpm": testis.median(axis=1),
        }
    )
    write_csv(
        phenotype.loc[med.index[symbols.isin(QC_GENES + MUSCLE_GENES)]].reset_index(),
        OUT / "marker-comparison.csv",
    )
    # Estimate profile similarity within testis-enriched genes after subtracting the peer baseline.
    excess = (independent[SAMPLE] - others.loc[panel].median(axis=1)).clip(lower=0)
    # Divide each gene by its typical testis abundance before comparing sample-specific deviations.
    typical = testis.loc[panel].median(axis=1)
    ranked = (
        testis.loc[panel]
        .div(typical, axis=0)
        .corrwith(excess / typical, method="spearman")
        .sort_values(ascending=False)
    )
    cohort = attr.loc[ranked.index, ["SMNABTCH", "SMGEBTCH"]].copy()
    cohort["same_extraction_batch"] = cohort.SMNABTCH.eq(attr.loc[SAMPLE, "SMNABTCH"])
    cohort["same_expression_batch"] = cohort.SMGEBTCH.str.split(",").map(
        lambda batches: attr.loc[SAMPLE, "SMGEBTCH"] in batches
    )
    cohort["same_donor"] = qc.loc[cohort.index, "participant_id"].eq(DONOR)
    cohort["profile_correlation"] = ranked
    for offset in range(2):
        subset = sorted(panel)[offset::2]
        cohort[f"gene_half_{offset + 1}_rho"] = (
            testis.loc[subset]
            .div(typical[subset], axis=0)
            .corrwith(excess[subset] / typical[subset], method="spearman")
        )
    donor_tissues = qc.loc[
        qc.participant_id.eq(DONOR), ["source_tissue", "testis_tissue_content", "library_qc_notes"]
    ]
    write_csv(donor_tissues.reset_index(), OUT / "same-donor-tissues.csv")
    write_csv(cohort.rename_axis("sample_id").reset_index(), OUT / "testis-profile-comparison.csv")
    # Retain all other normal observations, including existing tissue summaries. Replace
    # just the reduced muscularis row with its exact leave-one-out within-sample group maximum.
    alternate = []
    for chunk in pd.read_csv(
        REBUILD / "normal-comparison.csv.gz",
        usecols=[
            "entity_id",
            "value",
            "scope",
            "qc_pass",
            "identity_resolved",
            "sample_id",
            "source_context",
            "source_id",
        ],
        chunksize=150000,
    ):
        valid = chunk[
            chunk.entity_id.isin(cases.index)
            & chunk.qc_pass
            & chunk.identity_resolved
            & chunk.scope.eq("normal")
            & ~chunk.sample_id.eq(SAMPLE)
        ]
        alternate.append(valid.loc[valid.groupby("entity_id").value.idxmax()])
    candidates = pd.concat(alternate, ignore_index=True)
    next_normal = candidates.loc[candidates.groupby("entity_id").value.idxmax()].set_index(
        "entity_id"
    )
    cache = ROOT / "tmp/target-catalog-cache"
    receipt = json.loads((cache / "receipt.json").read_text())
    tumor_path = cache / "tumor-comparison.parquet"
    if (
        sha256(tumor_path) != receipt["matrix_sha256"]
        or sha256(REBUILD / "tumor-samples.csv.gz") != receipt["samples_sha256"]
    ):
        raise ValueError("Tumor sample matrix or identity receipt mismatch")
    matrix_t = pd.read_parquet(tumor_path, columns=cases.index.tolist())
    meta = (
        pd.read_csv(REBUILD / "tumor-samples.csv.gz", keep_default_na=False)
        .set_index("sample_key")
        .loc[matrix_t.index]
    )
    prior = lost.set_index("entity_id")
    rows = []
    for entity, r in cases.iterrows():
        ms = json.loads(r.member_gene_ids)
        src = canonical[canonical.isin(ms)]
        peer = others.reindex(src.index).sum(min_count=len(ms))
        runner = next_normal.loc[entity]
        n, alternative_sample, alternative_context = (
            runner.value,
            runner.sample_id,
            runner.source_context,
        )
        if peer.notna().any() and peer.max() > n:
            n, alternative_sample, alternative_context = peer.max(), peer.idxmax(), MUSCULARIS
        tumor = matrix_t[entity].dropna()
        old = 10 * max(
            1,
            prior.loc[entity, "prior_normal_reference_tpm"],
            prior.loc[entity, "prior_disposable_max_tpm"],
        )
        hits = meta.loc[tumor[tumor.ge(old)].index]
        assert np.isclose(tumor.max(), r.tumor_max_tpm) and r.tumor_fold < 10 and len(hits) > 0
        top = json.loads(r.normal_reference_contexts)[0]
        rows.append(
            {
                "protein": r.protein_group_label,
                "top_cancer": r.tumor_top_cancer_code,
                "top_tumor_sample": r.tumor_top_sample,
                "tumor_tpm": tumor.max(),
                "measured_tumors": len(tumor),
                "tumors_ge1_tpm": int(tumor.ge(1).sum()),
                "prior_tenfold_libraries": len(hits),
                "prior_tenfold_known_patients": hits.loc[
                    hits.participant_id.ne(""), "participant_id"
                ].nunique(),
                "prior_tenfold_anonymous_libraries": int(hits.participant_id.eq("").sum()),
                "prior_tenfold_cancer_codes": "|".join(sorted(hits.cancer_code.unique())),
                "prior_threshold_tpm": old,
                "prior_tumor_fold": prior.loc[entity, "prior_tumor_fold"],
                "normal_tpm": r.normal_reference_tpm,
                "normal_tissue": top["tissue"],
                "normal_sample": top.get("sample_id", ""),
                "disposable_tpm": r.disposable_max_tpm,
                "threshold_tpm": r.single_tumor_threshold_tpm,
                "tumor_fold": r.tumor_fold,
                "reproductive_fold": r.reproductive_fold,
                "normal_tpm_without_sample": n,
                "next_normal_sample": alternative_sample,
                "next_normal_context": alternative_context,
                "tumor_fold_without_sample": tumor.max() / max(1, n, r.disposable_max_tpm),
            }
        )
    result = pd.DataFrame(rows).sort_values("protein")
    write_csv(result, OUT / "tumor-evidence-98.csv")
    summary = {
        "normal_sample": SAMPLE,
        "pathology": attr.loc[SAMPLE, "SMPTHNTS"],
        "rna_integrity_number": float(attr.loc[SAMPLE, "SMRIN"]),
        "extraction_batch": attr.loc[SAMPLE, "SMNABTCH"],
        "expression_batch": attr.loc[SAMPLE, "SMGEBTCH"],
        "esophagus_samples": eso.shape[1],
        "testis_samples": testis.shape[1],
        "genes_for_correlations": int(expressed.sum()),
        "independent_testis_enriched": len(panel),
        "above_every_esophagus_peer": int(
            independent[SAMPLE].gt(others.loc[panel].max(axis=1)).sum()
        ),
        "median_ratio_to_testis": float((excess / typical).median()),
        "proteins": len(result),
        "sample_sets_normal_maximum": int(result.normal_sample.eq(SAMPLE).sum()),
        "regain_tumor_tenfold_without_sample": int(result.tumor_fold_without_sample.ge(10).sum()),
        "subject_libraries": int(qc.participant_id.eq(DONOR).sum()),
        "annotation_sha256": sha256(annotations),
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("annotations", type=Path)
    parser.add_argument("samples_dir", type=Path)
    args = parser.parse_args()
    build(args.annotations, args.samples_dir)
