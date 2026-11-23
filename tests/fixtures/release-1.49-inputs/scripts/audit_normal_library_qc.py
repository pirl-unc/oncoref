"""Report advisory pancreatic and esophageal RNA signatures across normal sample sources."""

from pathlib import Path

import pandas as pd
from cta_report_common import write_csv

from oncoref import normal_library_qc as qc
from oncoref import reference_data

OUT = (
    Path(__file__).resolve().parents[1] / "docs/audits/target-selection-20261110/rebuild/normal-rna"
)


def immune_libraries():
    median = pd.read_csv(
        reference_data.pinned_path("gtex_tissue_rna", "v11"),
        sep="\t",
        skiprows=2,
        usecols=["Name", "Description", *(native for native, _ in qc.CARRYOVER_PANELS.values())],
    )
    median["gene_id"] = median.Name.str.split(".").str[0]
    reference = qc.carryover_reference(median)
    path = reference_data.pinned_path("hpa_immune_cell_samples", "v23")
    metadata, observations = [], []
    for chunk in pd.read_csv(
        path,
        sep="\t",
        chunksize=200_000,
        usecols=["Sample ID", "Donor", "Immune cell", "ENSG ID", "TPM"],
    ):
        metadata.append(chunk[["Sample ID", "Donor", "Immune cell"]].drop_duplicates())
        observations.append(chunk[chunk["ENSG ID"].isin(reference.index)])
    meta = pd.concat(metadata).drop_duplicates().set_index("Sample ID", verify_integrity=True)
    values = pd.concat(observations).pivot(index="ENSG ID", columns="Sample ID", values="TPM")
    panels = [
        qc.published_carryover(values.reindex(columns=group.index), cell, reference)
        for cell, group in meta.groupby("Immune cell")
    ]
    return (
        meta.join(pd.concat(panels))
        .rename_axis("sample_id")
        .reset_index()
        .rename(columns={"Donor": "participant_id", "Immune cell": "tissue"})
    )


def main():
    OUT.mkdir(exist_ok=True, parents=True)
    immune = immune_libraries()
    write_csv(immune, OUT / "immune-library-qc.csv")
    counts = []
    for label, data in [
        ("GTEx v11 tissue", pd.read_csv(OUT.parent / "gtex-library-qc.csv.gz")),
        ("HPA v23 tissue", pd.read_csv(OUT.parent / "normal-library-qc.csv")),
        ("HPA v23 sorted immune", immune),
    ]:
        marks = data[["carryover_P", "carryover_E"]].astype("boolean")
        counts.append(
            {
                "source": label,
                "libraries": len(data),
                "assessed": marks.notna().any(axis=1).sum(),
                "P_assessed": marks.carryover_P.count(),
                "P": marks.carryover_P.sum(),
                "E_assessed": marks.carryover_E.count(),
                "E": marks.carryover_E.sum(),
                "either": marks.any(axis=1).sum(),
            }
        )
    write_csv(pd.DataFrame(counts), OUT / "carryover-summary.csv")


if __name__ == "__main__":
    main()
