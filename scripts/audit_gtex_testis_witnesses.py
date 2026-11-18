"""Trace the ten #716 testis leads to samples in a pinned public GTEx v11 GCT."""

import argparse
import gzip
import json
from pathlib import Path

import numpy as np
import pandas as pd

from oncoref.curation import digest

ROOT = Path(__file__).resolve().parents[1]
REVIEW = ROOT / "docs/audits/antigen-evidence-review-20261118"
DECISIONS = (
    ROOT
    / "docs/audits/cta-selection-20261028/alternative-evidence/gtex_tissue_rna-row-decisions.csv.gz"
)
SHA256 = "d2f63d00c9cfd7151e108d6521b47142435eef6fc79371ec6a637e9ce23130fd"


def audit(matrix, output):
    if digest(matrix) != SHA256:
        raise ValueError("Unexpected GTEx testis matrix checksum")
    leads = pd.read_csv(REVIEW / "ca-only-78.csv").query("bulk_R_fold >= 10")
    leads = leads.merge(
        pd.read_csv(REVIEW / "requested-cohorts.csv"), validate="one_to_one", on="entity_id"
    )
    identities = pd.read_csv(DECISIONS).set_index("gene_id")
    wanted = {}
    for row in leads.itertuples():
        (gene,) = json.loads(row.member_gene_ids)
        decision = identities.loc[gene]
        (source,) = json.loads(decision.source_gene_ids)
        if decision.status != "canonical" or source.split(".")[0] != gene:
            raise ValueError(f"Unresolved GTEx identity for {gene}")
        wanted[source] = row
    summaries, observations = [], []
    with gzip.open(matrix, "rt") as handle:
        version, dimensions = handle.readline().strip(), handle.readline().split()
        columns = handle.readline().rstrip().split("\t")
        samples = pd.Index(columns[2:], name="sample_id")
        participants = samples.str.extract(r"^(GTEX-[^-]+)-.+$", expand=False)
        if (
            version != "#1.2"
            or len(samples) != int(dimensions[1])
            or samples.has_duplicates
            or participants.hasnans
            or columns[:2] != ["Name", "Description"]
        ):
            raise ValueError("Unexpected GTEx sample header")
        for line in handle:
            source, _, vector = line.rstrip().split("\t", 2)
            if source not in wanted:
                continue
            row = wanted[source]
            values = np.asarray(vector.split("\t"), dtype=float)
            if len(values) != len(samples) or not np.isfinite(values).all() or min(values) < 0:
                raise ValueError(f"Invalid TPM values for {source}")
            bounds = [row.N_reference_TPM, row.D_TPM]
            if not np.isfinite(bounds).all():
                raise ValueError(f"Incomplete normal reference for {source}")
            bar = 10 * max(*bounds, 1)
            hits = values >= bar
            observations.append(
                pd.DataFrame(
                    {
                        "protein": row.protein,
                        "source_gene": source,
                        "sample_id": samples,
                        "participant_id": participants,
                        "tpm": values,
                        "at_tenfold": hits,
                    }
                )
            )
            summaries.append(
                {
                    "protein": row.protein,
                    "source_gene": source,
                    "threshold_tpm": bar,
                    "median_tpm": np.median(values),
                    "max_tpm": max(values),
                    "top_sample": samples[values.argmax()],
                    "samples": len(samples),
                    "samples_at_tenfold": hits.sum(),
                    "donors_at_tenfold": participants[hits].nunique(),
                }
            )
    if sorted(r["source_gene"] for r in summaries) != sorted(wanted):
        raise ValueError("Missing or duplicate exact-version GTEx measurements")
    output.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(summaries).sort_values("protein").to_csv(
        output / "testis-summary.csv", index=False
    )
    pd.concat(observations).sort_values(["protein", "sample_id"]).to_csv(
        output / "testis-samples.csv.gz", index=False, compression={"method": "gzip", "mtime": 0}
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("matrix", type=Path)
    parser.add_argument(
        "--output", type=Path, default=ROOT / "docs/audits/gtex-reproductive-witnesses"
    )
    args = parser.parse_args()
    audit(args.matrix, args.output)
