#!/usr/bin/env python3
"""Pin two distinct >=10 TPM tumor observations for normal-low proteins.

Uses the existing source-matrix cache, QC-selected samples, reference identity
resolver and pinned cohort receipts. It stops at two supporting sample identities:
counts are lower bounds, not prevalence. HPA pTPM is not substituted for TPM.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from cta_report_common import sha256, write_csv  # noqa: E402

from oncoref import expression, source_matrices  # noqa: E402
from oncoref.gene_matching import harmonize_reference_expression  # noqa: E402

BASELINE = ROOT / "docs/audits/cta-neither-20261030/baseline-1.24.0.csv.gz"
ORIGINAL = ROOT / "docs/audits/cta-expanded-funnel-20261027/tumor-union"
CURRENT = ROOT / "docs/audits/cta-selection-20261028"
REPAIRED = CURRENT / "tumor-identity-repair"


def sample_key(namespace: str, sample: str) -> str:
    """Deduplicate source views and TCGA aliquots; never count them as replication."""
    if sample.startswith("TCGA-"):
        return "TCGA:" + sample[:12]
    return namespace + ":" + sample


def build(*, download: bool = False) -> None:
    facts = pd.read_csv(BASELINE, low_memory=False)
    potential = facts[
        facts.required_coverage
        & ~facts.unresolved_evidence
        & facts.normal_maximum_floor_ratio.lt(1)
        & ~facts.any_detected_protein
    ].copy()
    members = {
        g: r.protein_equivalence_group_id
        for r in potential.itertuples()
        for g in json.loads(r.member_gene_ids)
    }
    expected_members = potential.set_index("protein_equivalence_group_id").member_count
    seen = {g: set() for g in potential.protein_equivalence_group_id}
    evidence, inventory, hashes = [], [], {}
    jobs = {}
    for folder in [ORIGINAL, REPAIRED]:
        receipt = json.loads((folder / "manifest.json").read_text())
        for name, digest in receipt["outputs"].items():
            if sha256(folder / name) != digest:
                raise ValueError(f"Changed pinned tumor evidence: {folder / name}")
        available = pd.read_csv(folder / "cohort-availability.csv")
        samples = pd.read_csv(folder / "sample-identity.csv.gz")
        paths = {
            Path(name).stem: (Path(name), digest) for name, digest in receipt["matrices"].items()
        }
        for row in available[available.status.eq("Used")].itertuples():
            if row.source_scale_class not in {"linear_rnaseq_tpm", "ncbi_gene_feature_count_tpm"}:
                raise ValueError(f"Tumor support requires native linear TPM: {row.cancer_code}")
            path = folder / "cohorts" / f"{row.cancer_code}.csv.gz"
            peaks = pd.read_csv(path, usecols=["gene_id", "maximum_tpm"])
            high = set(peaks.loc[peaks.gene_id.isin(members) & peaks.maximum_tpm.ge(10), "gene_id"])
            groups = folder / "protein-groups" / f"{row.cancer_code}.csv.gz"
            if groups.exists():
                group_peaks = pd.read_csv(groups)
                high_groups = set(
                    group_peaks.loc[group_peaks.maximum_tpm.ge(10), "protein_equivalence_group_id"]
                )
                high |= {g for g, entity in members.items() if entity in high_groups}
            jobs[row.cancer_code] = (
                row,
                folder,
                paths[row.cancer_code],
                samples[samples.cancer_code.eq(row.cancer_code)],
                high,
            )
        hashes[str(folder / "manifest.json")] = sha256(folder / "manifest.json")

    # The cache layout and fetcher are the same as other OncoRef expression data.
    cache = ROOT / "tmp/tumor-restricted-matrices"
    for code, (row, folder, (original_path, digest), samples, high) in sorted(
        jobs.items(), key=lambda item: (-len(item[1][-1]), item[0])
    ):
        needed_groups = {members[g] for g in high if len(seen[members[g]]) < 2}
        record = {"cancer_code": code, "source_cohort": row.source_cohort, "matrix_sha256": digest}
        if not needed_groups:
            inventory.append(record | {"status": "No additional qualifying support needed"})
            continue
        local = source_matrices.local_path(code, str(row.source_matrix_version))
        old_override = os.environ.get(source_matrices.CACHE_DIR_ENV_VAR)
        os.environ[source_matrices.CACHE_DIR_ENV_VAR] = str(cache)
        try:
            cached = source_matrices.local_path(code, str(row.source_matrix_version))
        finally:
            if old_override is None:
                os.environ.pop(source_matrices.CACHE_DIR_ENV_VAR, None)
            else:
                os.environ[source_matrices.CACHE_DIR_ENV_VAR] = old_override
        path = next((p for p in [original_path, local, cached] if p.exists()), local)
        if not path.exists() and download:
            old_override = os.environ.get(source_matrices.CACHE_DIR_ENV_VAR)
            os.environ[source_matrices.CACHE_DIR_ENV_VAR] = str(cache)
            try:
                path = source_matrices.ensure(code, str(row.source_matrix_version))
            finally:
                if old_override is None:
                    os.environ.pop(source_matrices.CACHE_DIR_ENV_VAR, None)
                else:
                    os.environ[source_matrices.CACHE_DIR_ENV_VAR] = old_override
        if not path.exists():
            inventory.append(record | {"status": "Pinned matrix unavailable"})
            continue
        if sha256(path) != digest:
            raise ValueError(f"Pinned matrix checksum mismatch: {code}")
        identity_path = folder / "identity" / f"{code}.csv.gz"
        identity = pd.read_csv(identity_path, low_memory=False)
        need_genes = {g for g, entity in members.items() if entity in needed_groups}
        source_ids = (
            identity.loc[
                identity.identity_status.eq("matched")
                & identity.canonical_gene_id.isin(need_genes),
                "source_gene_id",
            ]
            .astype(str)
            .tolist()
        )
        raw = pq.read_table(path, filters=[("Ensembl_Gene_ID", "in", source_ids)]).to_pandas()
        sample_ids = samples.sample_id.astype(str).tolist()
        if not set(sample_ids) <= set(expression.sample_columns(raw)):
            raise ValueError(f"QC sample identifiers absent from matrix: {code}")
        matrix, decisions = harmonize_reference_expression(
            raw[[*expression.id_columns(raw), *sample_ids]],
            sample_cols=sample_ids,
            source=row.source_cohort,
            source_version=str(row.source_matrix_version),
            symbol_column="Symbol" if "Symbol" in raw else None,
            duplicate_policy="quarantine",
        )
        matrix = matrix[matrix.Ensembl_Gene_ID.isin(need_genes)].copy()
        matrix["entity_id"] = matrix.Ensembl_Gene_ID.map(members)
        summed = matrix.groupby("entity_id")[sample_ids].sum(min_count=1)
        measured = matrix.groupby("entity_id")[sample_ids].count()
        namespace = source_matrices.source_sample_namespace(row.source_cohort)
        for entity, values in summed.iterrows():
            if np.isinf(values).any() or values.lt(0).any():
                raise ValueError("Expected finite nonnegative native TPM or missing")
            for sample in values[values.ge(10)].sort_values(ascending=False).index:
                key = sample_key(namespace, sample)
                # Require all loci measured; do not silently supply a missing member as zero.
                if (
                    key in seen[entity]
                    or measured.loc[entity, sample] != expected_members.loc[entity]
                ):
                    continue
                if len(seen[entity]) >= 2:
                    break
                seen[entity].add(key)
                evidence.append(
                    {
                        "protein_equivalence_group_id": entity,
                        "cancer_code": code,
                        "source_cohort": row.source_cohort,
                        "source_version": row.source_matrix_version,
                        "sample_id": sample,
                        "sample_key": key,
                        "unit": "TPM",
                        "value": values[sample],
                        "matrix_sha256": digest,
                        "identity_resolved": True,
                    }
                )
        inventory.append(record | {"status": "Assessed", "path": str(path)})
        write_csv(decisions, CURRENT / "tumor-restricted-identities" / f"{code}.csv.gz")
        print(
            f"{code}: {sum(len(x) >= 2 for x in seen.values())} proteins with two samples",
            flush=True,
        )

    support = facts[["protein_equivalence_group_id"]].copy()
    support["tumor_samples_ge10"] = support.protein_equivalence_group_id.map(
        lambda g: len(seen[g]) if g in seen else np.nan
    )
    support["tumor_identity_resolved"] = support.protein_equivalence_group_id.map(
        lambda g: bool(seen.get(g))
    )
    # Aggregate cancer-IHC counts lack reagent validation; they are not silently
    # promoted to reliable protein evidence. Reviewed records use a separate ledger.
    reviewed = pd.read_csv(ROOT / "oncoref/data/tumor-protein-evidence.csv", keep_default_na=False)
    accepted = reviewed[reviewed.status.eq("Accepted")]
    all_members = {g for ids in facts.member_gene_ids for g in json.loads(ids)}
    for row in accepted.itertuples():
        if (
            not row.source_url.startswith("https://")
            or not row.protein_evidence
            or row.identity_scope not in {"Locus", "Protein sequence"}
            or row.gene_id not in all_members
        ):
            raise ValueError(
                "Accepted tumor protein evidence needs source, assay and resolved identity"
            )
    supported_ids = set(accepted.gene_id)
    support["reliable_tumor_protein"] = facts.member_gene_ids.map(
        lambda ids: bool(set(json.loads(ids)) & supported_ids)
    )
    write_csv(support, CURRENT / "tumor-restricted-support.csv.gz")
    write_csv(pd.DataFrame(evidence), CURRENT / "tumor-restricted-observations.csv.gz")
    write_csv(pd.DataFrame(inventory), CURRENT / "tumor-restricted-source-coverage.csv")
    (CURRENT / "tumor-restricted-manifest.json").write_text(
        json.dumps(
            {
                "inputs": {
                    str(BASELINE): sha256(BASELINE),
                    "oncoref/data/tumor-protein-evidence.csv": sha256(
                        ROOT / "oncoref/data/tumor-protein-evidence.csv"
                    ),
                    **hashes,
                },
                "script_sha256": sha256(Path(__file__)),
                "threshold_tpm": 10,
                "distinct_samples_required": 2,
                "sample_counts": "lower bounds capped at two; source views and TCGA participants deduplicated",
                "protein_route": "explicit reviewed evidence; aggregate unvalidated cancer IHC is insufficient",
                "supported_proteins": sum(len(x) >= 2 for x in seen.values()),
                "outputs": {
                    name: sha256(CURRENT / name)
                    for name in [
                        "tumor-restricted-support.csv.gz",
                        "tumor-restricted-observations.csv.gz",
                        "tumor-restricted-source-coverage.csv",
                    ]
                },
            },
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--download", action="store_true")
    build(download=parser.parse_args().download)
