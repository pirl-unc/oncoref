#!/usr/bin/env python3
"""Trace CDY1/B and OR2T7 through every registered cancer cohort, in native scales."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pandas as pd

from oncoref import expression, source_matrices
from oncoref.gene_matching import match_reference_genes

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "docs/audits/cta-other-literature-20261030"
BASE = ROOT / "docs/audits/cta-expanded-funnel-20261027/tumor-union"
REPAIR = ROOT / "docs/audits/cta-selection-20261028/tumor-identity-repair"


def build():
    genes = pd.read_csv(OUT / "genes.csv")
    target = genes[genes.gene_symbol.isin(["CDY1", "CDY1B", "OR2T7"])]
    inventory = pd.concat(
        [
            pd.read_csv(REPAIR / "cohort-availability.csv"),
            pd.read_csv(BASE / "cohort-availability.csv"),
        ]
    ).drop_duplicates("cancer_code")
    rows, inputs = [], {}
    for cohort in inventory.itertuples():
        if cohort.eligible_scale:
            path = REPAIR / "cohorts" / f"{cohort.cancer_code}.csv.gz"
            if not path.exists():
                path = BASE / "cohorts" / f"{cohort.cancer_code}.csv.gz"
            frame = pd.read_csv(path)
            inputs[str(path.relative_to(ROOT))] = hashlib.sha256(path.read_bytes()).hexdigest()
            for gene in target.itertuples():
                hit = frame[frame.gene_id.eq(gene.gene_id)]
                r = hit.iloc[0] if len(hit) == 1 else None
                rows.append(
                    {
                        "gene_id": gene.gene_id,
                        "gene_symbol": gene.gene_symbol,
                        "cancer_code": cohort.cancer_code,
                        "source_cohort": cohort.source_cohort,
                        "source_matrix_version": cohort.source_matrix_version,
                        "scale": cohort.source_scale_class,
                        "maximum": r.maximum_tpm if r is not None else None,
                        "sample": r.maximum_sample if r is not None else None,
                        "unit": "TPM",
                        "use": "cancer presence" if r is not None else "unmeasured",
                        "note": "Passing source-specific sample QC; no sums of maxima from different samples.",
                    }
                )
        else:
            path = (
                ROOT
                / f"tmp/cta-other-literature/cancer-matrices/v{cohort.source_matrix_version}/{cohort.cancer_code}.parquet"
            )
            if not path.exists():
                path = source_matrices.local_path(
                    cohort.cancer_code, version=str(cohort.source_matrix_version)
                )
            if not path.exists():
                raise FileNotFoundError(
                    f"Fetch the registered native source with source_matrices.ensure: {cohort.cancer_code}"
                )
            inputs[str(path)] = hashlib.sha256(path.read_bytes()).hexdigest()
            raw = pd.read_parquet(path)
            samples = expression.sample_columns(raw)
            # Identity matching does not reinterpret log/proxy measurements as
            # nonnegative abundance. Keep unique source rows; never sum them.
            identities = match_reference_genes(
                raw,
                id_column="Ensembl_Gene_ID",
                source=cohort.source_cohort,
                source_version=str(cohort.source_matrix_version),
                symbol_column="Symbol" if "Symbol" in raw else None,
            )
            matched = raw.reset_index(drop=True).copy()
            matched["Ensembl_Gene_ID"] = identities.canonical_gene_id.to_numpy()
            matched = matched[identities.identity_status.eq("matched").to_numpy()]
            for gene in target.itertuples():
                hit = matched[matched.Ensembl_Gene_ID.eq(gene.gene_id)]
                values = hit[samples].iloc[0] if len(hit) == 1 else pd.Series(dtype=float)
                values = pd.to_numeric(values, errors="coerce").dropna()
                rows.append(
                    {
                        "gene_id": gene.gene_id,
                        "gene_symbol": gene.gene_symbol,
                        "cancer_code": cohort.cancer_code,
                        "source_cohort": cohort.source_cohort,
                        "source_matrix_version": cohort.source_matrix_version,
                        "scale": cohort.source_scale_class,
                        "maximum": values.max() if len(values) else None,
                        "sample": values.idxmax() if len(values) else None,
                        "unit": cohort.source_scale_class,
                        "use": "descriptive only; incompatible scale",
                        "note": "Native source value, not absolute TPM; not used for admission, ceilings or cross-source ranking. Invalid legacy values remain invalid.",
                    }
                )
    frame = pd.DataFrame(rows).sort_values(
        ["gene_symbol", "use", "maximum"], ascending=[True, True, False]
    )
    frame.to_csv(OUT / "requested-cancer-sources.csv", index=False)
    (OUT / "requested-cancer-manifest.json").write_text(
        json.dumps(
            {
                "cohorts": len(inventory),
                "member_genes": target.gene_symbol.tolist(),
                "inputs": inputs,
                "output_sha256": hashlib.sha256(
                    (OUT / "requested-cancer-sources.csv").read_bytes()
                ).hexdigest(),
            },
            indent=2,
        )
        + "\n"
    )
    print(frame.groupby(["gene_symbol", "use"]).maximum.max().to_string())


if __name__ == "__main__":
    build()
