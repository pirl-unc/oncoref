"""Compare the frozen ten testis leads with the live symmetric catalog threshold."""

import json
from pathlib import Path

import pandas as pd

from cta_report_common import sha256, write_csv

from oncoref import target_selection

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "docs/audits/gtex-reproductive-witnesses"


def main():
    baseline = pd.read_csv(OUT / "testis-summary.csv")
    observed = pd.read_csv(OUT / "testis-samples.csv.gz")
    catalog = target_selection.results()
    source = json.loads((OUT / "baseline.json").read_text())
    path = OUT / source["projection"]
    if sha256(path) != source["projection_sha256"]:
        raise ValueError("Unexpected pre-GTEx catalog projection")
    previous = pd.read_csv(path).set_index("entity_id")
    current = catalog.set_index("entity_id")
    if set(previous.index) != set(current.index):
        raise ValueError("Protein identities changed across the GTEx comparison")
    compared = current[["protein_group_label", *previous.columns]].join(
        previous.add_prefix("prior_"), validate="one_to_one"
    )
    compared["changed"] = compared.classification.ne(compared.prior_classification)
    selected = (
        compared.normal_reference_tpm.notna()
        | compared.prior_normal_reference_tpm.notna()
        | compared.changed
    )
    write_csv(compared.loc[selected].reset_index(), OUT / "catalog-comparison.csv.gz")
    rows = catalog.set_index("protein_group_label").loc[baseline.protein]
    bound = (
        rows[["normal_reference_tpm", "disposable_max_tpm"]].max(axis=1, skipna=False).clip(lower=1)
    )
    if bound.isna().any():
        raise ValueError("Missing normal or Disposable measurement for a pilot protein")
    output = baseline.set_index("protein")[["threshold_tpm", "donors_at_tenfold"]].add_prefix(
        "prior_"
    )
    output["normal_tpm"] = rows.normal_reference_tpm
    output["disposable_tpm"] = rows.disposable_max_tpm
    output["threshold_tpm"] = bound * target_selection.ENRICHMENT_FOLD
    inventory = pd.read_csv(
        ROOT / "docs/audits/target-selection-20261110/rebuild/gtex-library-qc.csv.gz",
        keep_default_na=False,
    )
    good = inventory.loc[
        inventory.tissue.eq("testis") & inventory.library_qc_review.eq(""), "sample_id"
    ]
    valid = observed[observed.sample_id.isin(good)]
    meets = valid[valid.tpm.ge(valid.protein.map(output.threshold_tpm))]
    output["testis_donors_at_tenfold"] = (
        meets.groupby("protein").participant_id.nunique().reindex(output.index, fill_value=0)
    )
    output["testis_max_tpm"] = valid.groupby("protein").tpm.max()
    output["tumor_fold"] = rows.tumor_fold
    output["reproductive_fold"] = rows.reproductive_fold
    output["classification"] = rows.classification
    output.to_csv(OUT / "symmetric-comparison.csv")


if __name__ == "__main__":
    main()
