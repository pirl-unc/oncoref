#!/usr/bin/env python3
"""Rebuild classification inputs after changing representative-protein membership.

Historical measurements are pinned at the audit commit. Normal measurements are
regrouped from source data; tumor sums use pinned native samples. No split protein
inherits its former partners' expression or safety classifications.
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
from cta_report_common import sha256, write_csv

from oncoref.gene_matching import gene_identity_catalog
from oncoref.protein_equivalence import equivalence_groups, group_identifier

ROOT = Path(__file__).resolve().parents[1]
BASELINE_COMMIT = "ab30edd5565714433beb7245f8916bd7dbbdde96"
STAGE = ROOT / "tmp/transcript-selection-rebuild"
INPUTS = ROOT / "docs/audits/cta-transcript-selection-20261108/inputs"
FROZEN = INPUTS / "frozen"
BASE = INPUTS / "regrouped"
CURRENT = ROOT / "docs/audits/cta-selection-20261028"
EXPANDED = ROOT / "docs/audits/cta-expanded-funnel-20261027"
HISTORICAL = ROOT / "docs/audits/cta-release-1.21.0"


def baseline(relative: str) -> Path:
    """Materialize exact historical bytes from the declared local Git commit."""
    target = STAGE / "baseline" / relative
    content = subprocess.check_output(["git", "show", f"{BASELINE_COMMIT}:{relative}"], cwd=ROOT)
    if not target.exists() or target.read_bytes() != content:
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(content)
    return target


def prepare():
    FROZEN.mkdir(parents=True, exist_ok=True)
    old = pd.read_csv(HISTORICAL / "protein-group-classifications.csv.gz", low_memory=False)
    genes = {g for ids in old.member_gene_ids.map(json.loads) for g in ids}
    catalog = gene_identity_catalog()
    known = set(catalog.canonical_gene_id)
    groups = equivalence_groups(genes & known)
    retained = {g for ids in groups.member_gene_ids.map(json.loads) for g in ids}
    extra = []
    historical_identities = {
        m["canonical_gene_id"]: m for value in old.member_identities.map(json.loads) for m in value
    }
    for gene_id in sorted(genes - retained):
        identity = historical_identities[gene_id]
        extra.append(
            {
                "protein_equivalence_group_id": group_identifier([gene_id]),
                "protein_group_label": identity["canonical_symbol"] or gene_id,
                "proteoform_key": gene_id,
                "representative_gene_id": gene_id,
                "member_gene_ids": json.dumps([gene_id]),
                "member_count": 1,
                "member_identities": json.dumps([identity]),
                "group_registry_sha256": groups.group_registry_sha256.iloc[0],
                "identity_catalog_sha256": groups.identity_catalog_sha256.iloc[0],
            }
        )
    groups = pd.concat([groups, pd.DataFrame(extra)], ignore_index=True)
    old_by_id = old.set_index("protein_equivalence_group_id")
    gene = pd.read_csv(EXPANDED / "gene-classifications.csv.gz", low_memory=False).set_index(
        "gene_id"
    )
    old_results = pd.read_csv(baseline("oncoref/data/cta-selection-human.csv.gz"))
    previous = old_results.set_index("entity_id").category
    rows = []
    boolean_columns = [c for c in old if old[c].dtype == bool]
    changed = []
    for group in groups.to_dict("records"):
        entity = group["protein_equivalence_group_id"]
        members = json.loads(group["member_gene_ids"])
        if entity in old_by_id.index:
            row = old_by_id.loc[entity].to_dict()
        else:
            part = gene.reindex(members)
            row = dict.fromkeys(old.columns, np.nan)
            for col in boolean_columns:
                row[col] = bool(part[col].fillna(False).any()) if col in part else False
            for col in set(part.columns) & set(old.columns):
                if col in boolean_columns:
                    continue
                row[col] = part.iloc[0][col] if len(members) == 1 else np.nan
            row["baseline_member_gene_categories"] = json.dumps(
                part.category.fillna("unassessed").to_dict()
            )
            row["member_named_controls"] = ""
            # A measured member is a valid presence lower bound, never a sum of maxima.
            row["hpa_tumor_group_lower_bound_ptpm"] = part.tumor_maximum_max_ptpm.max()
            row["category"] = "Other"
            changed.append(group)
        row.update(group)
        row["category"] = previous.get(entity, "Other")
        rows.append(row)
    facts = pd.DataFrame(rows)
    write_csv(facts, FROZEN / "protein-group-classifications.csv.gz")
    write_csv(pd.DataFrame(changed), STAGE / "changed-identities.csv")
    member_map = {
        g: r.protein_equivalence_group_id
        for r in groups.itertuples()
        for g in json.loads(r.member_gene_ids)
    }
    for name in ["known-cta-decisions.csv.gz", "gene-to-protein-group-decisions.csv.gz"]:
        frame = pd.read_csv(HISTORICAL / name, low_memory=False)
        frame["protein_equivalence_group_id"] = frame.gene_id.map(member_map)
        if frame.protein_equivalence_group_id.isna().any():
            raise ValueError("Historical locus lost in regrouping")
        write_csv(frame, FROZEN / name)
    shutil.copyfile(
        HISTORICAL / "imported-protein-review-decisions.csv",
        FROZEN / "imported-protein-review-decisions.csv",
    )
    for name in [
        "protein-group-classifications.csv.gz",
        "summary.json",
        "manifest.json",
        "tumor-restricted-support.csv.gz",
        "tumor-restricted-manifest.json",
    ]:
        baseline(f"docs/audits/cta-selection-20261028/{name}")
    print(
        f"Prepared {len(facts)} protein identities, including {len(changed)} new membership sets",
        flush=True,
    )


def normal():
    import cta_compartment_report as compartments
    import cta_protein_group_normal as group_normal

    sources = STAGE / "bulk-sources"
    sources.mkdir(parents=True, exist_ok=True)
    for name in ["rna_tissue_hpa", "rna_tissue_gtex"]:
        shutil.copyfile(ROOT / f"tmp/cta-source-consolidation/{name}.tsv.zip", sources / name)
    group_normal.build(BASE / "group-normal", sources, EXPANDED)
    facts = pd.read_csv(FROZEN / "protein-group-classifications.csv.gz", low_memory=False)
    compartments.BASE = BASE
    ihc = compartments.group_ihc_evidence(facts, FROZEN)
    compartments.bulk_assessments(facts, FROZEN, ihc)
    compartments.other_assay_assessments(facts, FROZEN)
    compartments.uncertain_group_assays(facts, FROZEN)
    compartments.apply_hpa_protein_levels(facts, FROZEN)
    write_csv(facts, FROZEN / "protein-group-classifications.csv.gz")
    print("Rebuilt normal bulk, cell, immune, DICE and cornea inputs", flush=True)


def alternative():
    import build_cta_alternative_evidence as alt

    alt.FACTS = FROZEN / "protein-group-classifications.csv.gz"
    alt.build()


def tumor():
    """Recompute same-sample sums using exactly the earlier QC-selected samples."""
    import pyarrow.parquet as pq
    from build_tumor_restricted_support import sample_key
    from cta_tumor_union import summarize_matrix, winning_observations

    from oncoref import expression
    from oncoref.gene_matching import harmonize_reference_expression

    original = EXPANDED / "tumor-union"
    repaired = STAGE / "baseline-repair"
    if not repaired.exists():
        relative = "docs/audits/cta-selection-20261028/tumor-identity-repair"
        receipt_path = baseline(f"{relative}/manifest.json")
        receipt = json.loads(receipt_path.read_text())
        repaired.mkdir(parents=True)
        for name in ["manifest.json", *receipt["outputs"]]:
            target = repaired / name
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(baseline(f"{relative}/{name}"), target)
    jobs = {}
    for folder in [original, repaired]:
        receipt = json.loads((folder / "manifest.json").read_text())
        for name, expected in receipt["outputs"].items():
            if sha256(folder / name) != expected:
                raise ValueError(f"Changed pinned tumor evidence: {folder / name}")
        available = pd.read_csv(folder / "cohort-availability.csv")
        samples = pd.read_csv(folder / "sample-identity.csv.gz")
        paths = {Path(p).stem: (Path(p), digest) for p, digest in receipt["matrices"].items()}
        for r in available.loc[available.status.eq("Used")].itertuples():
            jobs[r.cancer_code] = (
                folder,
                r,
                paths[r.cancer_code],
                samples.loc[samples.cancer_code.eq(r.cancer_code)],
            )
    out = CURRENT / "tumor-identity-repair"
    for name in ["cohorts", "protein-groups", "identity"]:
        (out / name).mkdir(parents=True, exist_ok=True)
    facts = pd.read_csv(FROZEN / "protein-group-classifications.csv.gz", low_memory=False)
    changed = pd.read_csv(STAGE / "changed-identities.csv")
    wanted = facts.loc[
        facts.member_count.gt(1)
        | facts.protein_equivalence_group_id.isin(changed.protein_equivalence_group_id)
    ]
    group_columns = [
        "protein_equivalence_group_id",
        "protein_group_label",
        "member_gene_ids",
        "member_count",
        "representative_gene_id",
    ]
    groups = wanted[group_columns]
    members = {
        g: r.protein_equivalence_group_id
        for r in groups.itertuples()
        for g in json.loads(r.member_gene_ids)
    }
    counts = groups.set_index("protein_equivalence_group_id").member_count
    changed_ids = set(changed.protein_equivalence_group_id)
    seen = {g: set() for g in changed_ids}
    observations, available_rows, sample_rows, native_hashes = [], [], [], {}
    for code, (folder, r, (source, expected), samples) in sorted(jobs.items()):
        paths = [
            source,
            *[
                base / source.parent.name / source.name
                for base in [
                    ROOT / "tmp/tumor-restricted-matrices",
                    Path.home() / ".cache/oncoref/source-matrices",
                    STAGE / "matrices",
                ]
            ],
        ]
        path = next((p for p in paths if p.exists()), None)
        if path is None or sha256(path) != expected:
            raise ValueError(f"Missing or changed pinned native matrix: {code}")
        native_hashes[str(path)] = expected
        identity = pd.read_csv(folder / "identity" / f"{code}.csv.gz", low_memory=False)
        source_ids = (
            identity.loc[
                identity.identity_status.eq("matched") & identity.canonical_gene_id.isin(members),
                "source_gene_id",
            ]
            .astype(str)
            .tolist()
        )
        raw = pq.read_table(path, filters=[("Ensembl_Gene_ID", "in", source_ids)]).to_pandas()
        sample_ids = sorted(samples.sample_id.astype(str))
        if not set(sample_ids) <= set(expression.sample_columns(raw)):
            raise ValueError(f"Pinned samples absent: {code}")
        matrix, decisions = harmonize_reference_expression(
            raw[[*expression.id_columns(raw), *sample_ids]],
            sample_cols=sample_ids,
            source=r.source_cohort,
            source_version=str(r.source_matrix_version),
            symbol_column="Symbol" if "Symbol" in raw else None,
            duplicate_policy="quarantine",
        )
        matrix = matrix.loc[matrix.Ensembl_Gene_ID.isin(members)].copy()
        matrix["group_id"] = matrix.Ensembl_Gene_ID.map(members)
        summed = matrix.groupby("group_id")[sample_ids].sum(min_count=1)
        measured = matrix.groupby("group_id")[sample_ids].count()
        peaks = summarize_matrix(
            summed.reset_index().rename(columns={"group_id": "Ensembl_Gene_ID"}), sample_ids
        ).rename(columns={"gene_id": "protein_equivalence_group_id"})
        peaks["measured_members_at_maximum"] = [
            int(measured.loc[v.protein_equivalence_group_id, v.maximum_sample])
            if v.maximum_sample
            else 0
            for v in peaks.itertuples()
        ]
        peaks["cancer_code"], peaks["source_cohort"], peaks["unit"] = code, r.source_cohort, "TPM"
        peaks = peaks.merge(groups, on="protein_equivalence_group_id", validate="one_to_one")
        write_csv(peaks.loc[peaks.member_count.gt(1)], out / "protein-groups" / f"{code}.csv.gz")
        write_csv(decisions, out / "identity" / f"{code}.csv.gz")
        # Reuse the published per-gene tables; only changed multi-locus sums are new.
        if folder == repaired:
            shutil.copyfile(
                folder / "cohorts" / f"{code}.csv.gz", out / "cohorts" / f"{code}.csv.gz"
            )
        else:
            (out / "cohorts" / f"{code}.csv.gz").unlink(missing_ok=True)
        namespace = samples.namespace.iloc[0]
        for entity in changed_ids & set(summed.index):
            for sample in sample_ids:
                if len(seen[entity]) >= 2:
                    break
                key = sample_key(namespace, sample)
                if (
                    key in seen[entity]
                    or measured.loc[entity, sample] != counts[entity]
                    or summed.loc[entity, sample] < 10
                    or pd.isna(summed.loc[entity, sample])
                ):
                    continue
                seen[entity].add(key)
                observations.append(
                    {
                        "protein_equivalence_group_id": entity,
                        "cancer_code": code,
                        "source_cohort": r.source_cohort,
                        "source_version": r.source_matrix_version,
                        "sample_id": sample,
                        "sample_key": key,
                        "unit": "TPM",
                        "value": summed.loc[entity, sample],
                        "matrix_sha256": expected,
                        "identity_resolved": True,
                    }
                )
        available_rows.append(r._asdict())
        sample_rows.append(samples)
        print(f"Rebuilt {code}: {len(sample_ids)} pinned samples", flush=True)
    for dirname, entity in [
        ("cohorts", "gene_id"),
        ("protein-groups", "protein_equivalence_group_id"),
    ]:
        name = (
            "tumor-maxima.csv.gz" if dirname == "cohorts" else "protein-group-tumor-maxima.csv.gz"
        )
        write_csv(
            winning_observations(
                (
                    (jobs[code][0] if dirname == "cohorts" else out) / dirname / f"{code}.csv.gz"
                    for code in sorted(jobs)
                ),
                entity_column=entity,
            ),
            out / name,
        )
    write_csv(pd.DataFrame(available_rows), out / "cohort-availability.csv")
    write_csv(pd.concat(sample_rows, ignore_index=True), out / "sample-identity.csv.gz")
    write_csv(
        pd.DataFrame(
            [
                {
                    "cancer_code": code,
                    "gene_table": str(
                        (
                            (out if values[0] == repaired else values[0])
                            / "cohorts"
                            / f"{code}.csv.gz"
                        ).relative_to(ROOT)
                    ),
                    "sha256": sha256(values[0] / "cohorts" / f"{code}.csv.gz"),
                }
                for code, values in sorted(jobs.items())
            ]
        ),
        out / "gene-table-inputs.csv",
    )
    (out / "summary.json").write_text(
        json.dumps(
            {
                "scope": "complete_pinned_tumor_union",
                "cohorts": len(jobs),
                "sample_occurrences": sum(len(x) for x in sample_rows),
                "protein_membership_records_rebuilt": len(groups),
                "baseline_commit": BASELINE_COMMIT,
            },
            indent=2,
        )
        + "\n"
    )
    (out / "manifest.json").write_text(
        json.dumps(
            {
                "scope": "complete_pinned_tumor_union",
                "baseline_commit": BASELINE_COMMIT,
                "inputs": {
                    str(Path(__file__).relative_to(ROOT)): sha256(Path(__file__)),
                    "oncoref/data/proteoform-groups-genome.csv": sha256(
                        ROOT / "oncoref/data/proteoform-groups-genome.csv"
                    ),
                },
                "matrices": native_hashes,
                "outputs": {
                    str(p.relative_to(out)): sha256(p)
                    for p in sorted(out.rglob("*"))
                    if p.is_file() and p.name != "manifest.json"
                },
            },
            indent=2,
        )
        + "\n"
    )
    old = pd.read_csv(
        baseline("docs/audits/cta-selection-20261028/tumor-restricted-support.csv.gz")
    ).set_index("protein_equivalence_group_id")
    support = old.reindex(facts.protein_equivalence_group_id).reset_index()
    for row in support.index[support.protein_equivalence_group_id.isin(changed_ids)]:
        entity = support.loc[row, "protein_equivalence_group_id"]
        support.loc[row, "tumor_samples_ge10"] = len(seen[entity])
        support.loc[row, "tumor_identity_resolved"] = bool(seen[entity])
    protein = pd.read_csv(ROOT / "oncoref/data/tumor-protein-evidence.csv", keep_default_na=False)
    accepted = set(protein.loc[protein.status.eq("Accepted"), "gene_id"])
    support["reliable_tumor_protein"] = facts.member_gene_ids.map(
        lambda value: bool(set(json.loads(value)) & accepted)
    ).to_numpy()
    write_csv(support, CURRENT / "tumor-restricted-support.csv.gz")
    write_csv(pd.DataFrame(observations), CURRENT / "regrouped-tumor-support-observations.csv.gz")
    old_receipt = json.loads(
        baseline("docs/audits/cta-selection-20261028/tumor-restricted-manifest.json").read_text()
    )
    old_receipt["inputs"]["oncoref/data/tumor-protein-evidence.csv"] = sha256(
        ROOT / "oncoref/data/tumor-protein-evidence.csv"
    )
    old_receipt["inputs"]["transcript_regrouping"] = sha256(out / "manifest.json")
    for name in ["tumor-restricted-support.csv.gz", "regrouped-tumor-support-observations.csv.gz"]:
        old_receipt["outputs"][name] = sha256(CURRENT / name)
    (CURRENT / "tumor-restricted-manifest.json").write_text(
        json.dumps(old_receipt, indent=2) + "\n"
    )


def replay():
    import replay_cta_selection as selection

    selection.FROZEN, selection.BASE = FROZEN, BASE
    selection.build()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=["prepare", "normal", "alternative", "tumor", "replay"])
    args = parser.parse_args()
    globals()[args.stage]()
