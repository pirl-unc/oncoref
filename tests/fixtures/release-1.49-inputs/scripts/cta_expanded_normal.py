#!/usr/bin/env python3
"""Evaluate every reference gene against the linked two-axis normal-RNA policy."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from cta_cell_policy_report import summarize  # noqa: E402
from cta_quantitative_threshold_report import positive_profile_controls  # noqa: E402
from cta_report_common import input_hashes  # noqa: E402
from cta_report_common import sha256 as digest  # noqa: E402
from cta_report_common import write_csv as write  # noqa: E402

from oncoref import cell_rna_policy as policy  # noqa: E402
from oncoref import (  # noqa: E402
    gene_ids,
    hpa,
    normal_cell_data,
    normal_cell_screen,
    reference_data,
)
from oncoref.cta_normal_screen import _routine_by_tissue, _tissue  # noqa: E402
from oncoref.normal_tissues import normal_tissue_screening_tier  # noqa: E402


def bulk_assessment(bulk: pd.DataFrame, ihc: pd.DataFrame, controls, conflicts) -> pd.DataFrame:
    """Vectorized equivalent of the existing per-gene complete-tissue IHC match."""
    annotations = {t: _tissue(t, "hpa_rna_consensus", "v23") for t in bulk.Tissue.unique()}
    d = bulk.rename(
        columns={"Gene": "Ensembl_Gene_ID", "Tissue": "tissue", "nTPM": "rna_ntpm"}
    ).copy()
    d["tissue_concept"] = d.tissue.map({t: a.concept_id for t, a in annotations.items()})
    d["screening_tier"] = d.tissue.map(
        {
            t: normal_tissue_screening_tier(a, policy="ihc_adaptive_cta")
            for t, a in annotations.items()
        }
    )
    routine = _routine_by_tissue(
        tuple(ihc[["Tissue", "Cell type"]].drop_duplicates().itertuples(index=False, name=None)),
        "hpa_normal_tissue",
        "v23",
    )
    h = ihc.copy()
    concepts = {t: _tissue(t, "hpa_normal_tissue", "v23").concept_id for t in h.Tissue.unique()}
    h["tissue_concept"] = h.Tissue.map(concepts)
    h["negative"] = h.Level.eq("Not detected")
    h["grade"] = h.Reliability.map(policy.GRADES).fillna(0)
    h["required_pair"] = [
        (t, c) in routine.get(concept, set())
        for t, c, concept in zip(h.Tissue, h["Cell type"], h.tissue_concept, strict=True)
    ]
    if h.duplicated(["Gene", "Tissue", "Cell type"]).any():
        raise ValueError("Duplicate IHC context")
    groups = h.groupby(["Gene", "tissue_concept"]).agg(
        all_negative=("negative", "all"),
        grade=("grade", "min"),
        observed_required=("required_pair", "sum"),
    )
    groups["required"] = (
        groups.index.get_level_values("tissue_concept")
        .map({k: len(v) for k, v in routine.items()})
        .fillna(0)
    )
    groups["usable"] = (
        groups.all_negative
        & groups.grade.gt(0)
        & groups.required.gt(0)
        & groups.observed_required.eq(groups.required)
    )
    groups["ihc_reliability_used"] = groups.grade.map(
        dict(enumerate(["Missing", *policy.GRADES]))
    ).where(groups.usable, "Missing")
    d = d.merge(
        groups[["ihc_reliability_used"]].reset_index().rename(columns={"Gene": "Ensembl_Gene_ID"}),
        on=["Ensembl_Gene_ID", "tissue_concept"],
        how="left",
        validate="many_to_one",
    )
    d["ihc_reliability_used"] = d.ihc_reliability_used.fillna("Missing")
    d.loc[d.screening_tier.isin(["permitted", "unresolved"]), "ihc_reliability_used"] = "Missing"
    d.loc[d.Ensembl_Gene_ID.isin(conflicts), "ihc_reliability_used"] = "Missing"
    d["profile_positive_control"] = d.Ensembl_Gene_ID.isin(controls)
    return policy.apply_bulk_ladder(d, mode="two_axis_doubling")


def build(out: Path, baseline: Path) -> None:
    """Build the versioned audit tables, report and source provenance."""
    out.mkdir(parents=True, exist_ok=True)
    facts = pd.read_csv(baseline / "gene-comparison.csv.gz", low_memory=False)
    ids = set(facts.gene_id)
    ihc = hpa._read_hpa("hpa_normal_tissue", "v23").dropna(subset=["Tissue", "Cell type", "Level"])
    controls = positive_profile_controls(ihc)
    conflicts = set(facts.loc[facts["Antibody concern"], "gene_id"])
    receipts, aggregate, peaks = [], [], []

    def receipt(source):
        r = reference_data.provenance(source, verify_content=True)
        receipts.append({k: r[k] for k in ["name", "version", "url", "sha256"]})
        return r

    def collect(d, assay):
        key = "Ensembl_Gene_ID" if assay == "Normal bulk" else "gene_id"
        value = "rna_ntpm" if assay == "Normal bulk" else "value"
        cutoff = "effective_ceiling_ntpm" if assay == "Normal bulk" else "cta_cutoff"
        d = d.copy()
        d["ratio"] = d[value] / d[cutoff]
        if assay == "Normal bulk":
            d["unit"], d["cell_type"] = "nTPM", "bulk"
        a = d.groupby(key).agg(
            candidate_exclusion=("local_exclusion", "any"),
            unsafe_evidence=("reference_review", "any"),
            observed_contexts=(value, "size"),
        )
        a.index.name = "gene_id"
        a = a.reset_index()
        a["assay"] = assay
        aggregate.append(a)
        peak = d.sort_values("ratio", ascending=False).drop_duplicates(key)
        peak = peak[
            [
                key,
                "tissue",
                "cell_type",
                "unit",
                value,
                cutoff,
                "ratio",
                "negative_ihc_grade",
                "local_exclusion",
                "reference_review",
            ]
        ].rename(columns={key: "gene_id", value: "value", cutoff: "candidate_ceiling"})
        peak["assay"] = assay
        peaks.append(peak)

    receipt("hpa_rna_consensus")
    receipt("hpa_normal_tissue")
    original = hpa._read_hpa("hpa_rna_consensus", "v23")
    raw = original[original.Gene.isin(ids)].copy()
    # Only registered stable-ID aliases may fill missing gene/tissue pairs.
    extras = original[~original.Gene.isin(ids)].copy()
    extras["source_gene_id"] = extras.Gene
    extras["Gene"] = extras.Gene.map(gene_ids.resolve_ensembl_id)
    extras = extras[extras.Gene.isin(ids)]
    extras = extras[
        ~pd.MultiIndex.from_frame(extras[["Gene", "Tissue"]]).isin(
            pd.MultiIndex.from_frame(raw[["Gene", "Tissue"]])
        )
    ]
    if extras.duplicated(["Gene", "Tissue"]).any():
        raise ValueError("Competing aliases for a missing bulk observation")
    write(extras, out / "bulk-alias-filled-observations.csv")
    raw = pd.concat([raw, extras.drop(columns="source_gene_id")], ignore_index=True)
    b = bulk_assessment(raw, ihc, controls, conflicts)
    previous = pd.read_csv(baseline / "bulk-contexts.csv.gz", low_memory=False)
    replay = previous.merge(
        b, on=["Ensembl_Gene_ID", "tissue"], suffixes=("_old", "_new"), validate="one_to_one"
    )
    for col in [
        "effective_ceiling_ntpm",
        "local_exclusion",
        "reference_review",
        "negative_ihc_grade",
    ]:
        left, right = replay[col + "_old"], replay[col + "_new"]
        if not (left.eq(right) | (left.isna() & right.isna())).all():
            raise ValueError(f"Bulk vectorized assessment does not replay {col}")
    collect(b, "Normal bulk")
    write(b, out / "bulk-contexts.csv.gz")
    missing = facts[facts.normal_bulk_tissues.lt(50)]
    rows = []
    all_tissues = set(original.Tissue.unique())
    symbols = original.groupby("Gene")["Gene name"].first()
    for r in missing.itertuples():
        same = original[original["Gene name"].eq(r.symbol)]
        native = original[original.Gene.eq(r.gene_id)]
        aliases = [
            a for a, canonical in gene_ids.ensembl_id_aliases().items() if canonical == r.gene_id
        ]
        rows.append(
            {
                "gene_id": r.gene_id,
                "symbol": r.symbol,
                "source_symbol": symbols.get(r.gene_id),
                "observed_tissues": native.Tissue.nunique(),
                "missing_tissues": " | ".join(sorted(all_tissues - set(native.Tissue))),
                "registered_aliases": " | ".join(aliases),
                "other_same_symbol_source_ids": " | ".join(sorted(set(same.Gene) - {r.gene_id})),
                "alias_rows_filled": int(extras.Gene.eq(r.gene_id).sum()),
            }
        )
    write(pd.DataFrame(rows), out / "bulk-missing-identity-audit.csv")
    print(f"Bulk replay verified: {len(replay)} contexts; {len(extras)} alias fills", flush=True)

    def assess(d, assay):
        if d.empty:
            return
        d = normal_cell_screen.contexts(
            d.drop(columns=["category", "mapping_status", "concept_id"], errors="ignore")
        )
        d = policy.reproductive_immune_contexts(d)
        d = policy.protected_immune_contexts(d)
        d = policy.match_negative_ihc(
            d,
            ihc[ihc.Gene.isin(d.gene_id.unique())],
            positive_controls=controls,
            conflicting_genes=conflicts,
        )
        d = policy.assess(
            d,
            policy.linked_thresholds(mode="two_axis_doubling"),
            ihc_factors=dict.fromkeys(["Missing/uncertain", *policy.GRADES], 1),
            inclusive_ceiling=True,
        )
        d["local_exclusion"] = d.cta_excluded | ~d.mapping_status.eq("mapped")
        d.loc[d.category.eq(policy.NEUTRAL_REPRODUCTIVE_IMMUNE), "local_exclusion"] = False
        collect(d, assay)

    r = receipt("hpa_single_cell_tissue")
    for number, chunk in enumerate(pd.read_csv(r["path"], sep="\t", chunksize=250000)):
        d = (
            chunk[chunk.Gene.isin(ids)]
            .rename(
                columns={
                    "Gene": "gene_id",
                    "Tissue": "tissue",
                    "Cell type": "cell_type",
                    "nTPM": "value",
                }
            )
            .copy()
        )
        d["source"], d["version"], d["taxonomy_id"] = "hpa_single_cell_tissue", "v23", 9606
        d["state"], d["stage"], d["unit"] = "unspecified", "adult", "nTPM"
        d["value"] = d.value.where(d["Read count"].ge(20))
        assess(d, "HPA tissue clusters")
        if number % 5 == 0:
            print(f"Clusters: chunk {number}", flush=True)
    for source, label in [
        ("hpa_immune_cell", "HPA immune"),
        ("hpa_immune_monaco", "Monaco immune"),
    ]:
        receipt(source)
        d = normal_cell_screen.contexts(normal_cell_data.observations(source, genes=ids))
        assess(summarize(d, support=d.value.notna()), label)
        print(label, flush=True)
    ambiguous = []
    for source in normal_cell_data.sources().source:
        if source.startswith("dice_"):
            receipt(source)
            d = normal_cell_data.dice_observations(source, genes=ids)
            duplicated = set(
                d.loc[d.duplicated(["gene_id", "replicate_slot"], keep=False), "gene_id"]
            )
            ambiguous.extend(
                {
                    "gene_id": gid,
                    "source": source,
                    "reason": "ambiguous duplicate replicate identity; not pooled",
                }
                for gid in duplicated
            )
            d = normal_cell_screen.contexts(d[~d.gene_id.isin(duplicated)])
            assess(
                summarize(d, support=d.value.notna(), rank=3, unique="replicate_slot"),
                "DICE immune",
            )
            print(source, flush=True)
    write(
        pd.DataFrame(ambiguous, columns=["gene_id", "source", "reason"]),
        out / "optional-ambiguous-observations.csv",
    )
    for source in ["cornea_collin_adult_counts", "cornea_collin_adultcornea_cells"]:
        receipt(source)
    d = normal_cell_screen.contexts(normal_cell_data.cornea_summary(genes=ids))
    assess(
        summarize(d, support=d.positive_cells.ge(5) & d.read_count.ge(20), rank=2, unique="sample"),
        "Adult cornea",
    )
    totals = (
        pd.concat(aggregate)
        .groupby(["gene_id", "assay"])
        .agg(
            candidate_exclusion=("candidate_exclusion", "any"),
            unsafe_evidence=("unsafe_evidence", "any"),
            observed_contexts=("observed_contexts", "sum"),
        )
        .reset_index()
    )
    write(totals, out / "normal-assessments.csv.gz")
    peak_table = (
        pd.concat(peaks).sort_values("ratio", ascending=False).drop_duplicates(["gene_id", "assay"])
    )
    write(peak_table, out / "normal-peak-evidence.csv.gz")
    (out / "manifest.json").write_text(
        json.dumps(
            {
                "sources": receipts,
                "inputs": input_hashes(
                    [
                        Path(__file__),
                        baseline / "gene-comparison.csv.gz",
                        baseline / "bulk-contexts.csv.gz",
                        ROOT / "oncoref/cell_rna_policy.py",
                    ]
                ),
                "outputs": {
                    p.name: digest(p)
                    for p in out.iterdir()
                    if p.is_file() and p.name != "manifest.json"
                },
            },
            indent=2,
        )
        + "\n"
    )
    print(f"Completed {len(totals)} gene/assay assessments", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path)
    parser.add_argument(
        "--baseline", type=Path, default=ROOT / "docs/audits/cta-geometric-thresholds-20261027"
    )
    args = parser.parse_args()
    build(args.output, args.baseline)
