#!/usr/bin/env python3
"""Recompute same-context normal expression for identical-protein locus groups.

Bulk components are summed before taking the HPA/GTEx maximum. Cell-cluster
counts are summed before support filtering. Unpaired DICE distributions and
unavailable cell-level detection unions retain bounds rather than invented pairs.
"""

from __future__ import annotations

import argparse
import json
import sys
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from cta_cell_policy_report import summarize  # noqa: E402
from cta_expanded_normal import bulk_assessment  # noqa: E402
from cta_quantitative_threshold_report import positive_profile_controls  # noqa: E402
from cta_report_common import (  # noqa: E402
    input_hashes,
    sha256,
    write_csv,
)

from oncoref import cell_rna_policy as policy  # noqa: E402
from oncoref import hpa, normal_cell_data, normal_cell_screen, reference_data  # noqa: E402
from oncoref.gene_matching import match_reference_genes  # noqa: E402
from oncoref.protein_equivalence import (  # noqa: E402
    equivalence_groups,
    sum_context_observations,
    unpaired_rank_bounds,
)

# HPA v23 groups GTEx sub-tissues by taking their maximum. Keep the source
# subregion, sum loci there first, then aggregate to the consensus anatomy.
# https://v23.proteinatlas.org/about/download (RNA consensus tissue gene data)
BULK_REGION_GROUPS = {
    "caudate": "basal ganglia",
    "putamen": "basal ganglia",
    "nucleus accumbens": "basal ganglia",
    "hippocampus": "hippocampal formation",
    "substantia nigra": "midbrain",
}


def consensus_bulk_regions(components: pd.DataFrame) -> pd.DataFrame:
    """Maximum of same-subregion locus sums, with original anatomy retained."""
    result = components.copy()
    result["source_tissue"] = result.Tissue
    result["Tissue"] = result.Tissue.replace(BULK_REGION_GROUPS)
    return result.sort_values(
        ["nTPM", "cohort", "source_tissue"], ascending=[False, True, True]
    ).drop_duplicates(["gene_id", "Tissue"])


def combine_ihc(ihc: pd.DataFrame, groups: pd.DataFrame) -> pd.DataFrame:
    """Any positive persists; negative certainty requires all expected members."""
    members = {
        gid: r.protein_equivalence_group_id
        for r in groups.itertuples()
        for gid in json.loads(r.member_gene_ids)
    }
    meta = groups.set_index("protein_equivalence_group_id")
    d = ihc[ihc.Gene.isin(members)].copy()
    d["group_id"] = d.Gene.map(members)
    rows = []
    for (gid, tissue, cell), part in d.groupby(["group_id", "Tissue", "Cell type"]):
        positive = part.Level.isin(policy.POSITIVE)
        complete = part.Gene.nunique() == meta.loc[gid, "member_count"]
        all_negative = complete and part.Level.eq("Not detected").all()
        grade = part.Reliability.map(policy.GRADES).fillna(0).min()
        rows.append(
            {
                "Gene": meta.loc[gid, "representative_gene_id"],
                "Tissue": tissue,
                "Cell type": cell,
                "Level": "Low" if positive.any() else "Not detected" if all_negative else "Missing",
                "Reliability": {v: k for k, v in policy.GRADES.items()}.get(grade, "Uncertain")
                if complete
                else "Uncertain",
            }
        )
    return pd.DataFrame(rows)


def build(output: Path, source_dir: Path, baseline: Path) -> None:
    """Write group measurements, coverage, identity audits and normal decisions."""
    output.mkdir(parents=True, exist_ok=True)
    facts = pd.read_csv(baseline / "gene-classifications.csv.gz", low_memory=False)
    groups = equivalence_groups(facts.gene_id)
    groups = groups[groups.member_count.gt(1)].copy()
    member_to_rep = {
        gid: r.representative_gene_id
        for r in groups.itertuples()
        for gid in json.loads(r.member_gene_ids)
    }
    member_ids = set(member_to_rep)
    write_csv(groups, output / "protein-equivalence-groups.csv")
    receipts, identity_audits, summaries, peaks = [], [], [], []

    def native(source: str) -> Path:
        receipt = reference_data.provenance(source, verify_content=True)
        receipts.append({k: receipt[k] for k in ["name", "version", "url", "sha256"]})
        return Path(receipt["path"])

    def match(
        d: pd.DataFrame,
        source: str,
        id_column: str = "Gene",
        symbol_column: str | None = "Gene name",
    ) -> pd.DataFrame:
        pairs = d[[id_column, *([symbol_column] if symbol_column else [])]].drop_duplicates()
        a = match_reference_genes(
            pairs,
            id_column=id_column,
            symbol_column=symbol_column,
            source=source,
            source_version="pinned-source",
        )
        identity_audits.append(a)
        pairs["gene_id"] = a.canonical_gene_id.to_numpy()
        pairs["matched"] = a.identity_status.eq("matched").to_numpy()
        joined = d.merge(
            pairs,
            on=[id_column, *([symbol_column] if symbol_column else [])],
            validate="many_to_one",
        )
        return joined[joined.matched & joined.gene_id.isin(member_ids)].drop(columns="matched")

    def sum_rows(
        d: pd.DataFrame, contexts: list[str], values: list[str], unit: str
    ) -> pd.DataFrame:
        result = sum_context_observations(
            d, groups, context_columns=contexts, value_columns=values, unit=unit
        )
        result["gene_id"] = result.representative_gene_id
        return result

    raw_ihc = hpa._read_hpa("hpa_normal_tissue", "v23").dropna(
        subset=["Tissue", "Cell type", "Level"]
    )
    ihc = combine_ihc(raw_ihc, groups)
    controls = {member_to_rep[g] for g in positive_profile_controls(raw_ihc) if g in member_to_rep}
    conflicts = {
        member_to_rep[g]
        for g in facts.loc[facts["Antibody concern"], "gene_id"]
        if g in member_to_rep
    }
    native("hpa_normal_tissue")
    write_csv(ihc, output / "group-ihc.csv.gz")

    def collect(d: pd.DataFrame, assay: str, *, bulk: bool = False) -> None:
        key, value, cutoff = (
            ("Ensembl_Gene_ID", "rna_ntpm", "effective_ceiling_ntpm")
            if bulk
            else ("gene_id", "value", "cta_cutoff")
        )
        d = d.copy()
        d["ratio"] = d[value] / d[cutoff]
        d["assay"] = assay
        if bulk:
            d["cell_type"], d["unit"] = "bulk", "nTPM"
        agg = (
            d.groupby(key)
            .agg(
                candidate_exclusion=("local_exclusion", "any"),
                unsafe_evidence=("reference_review", "any"),
                observed_contexts=(value, "size"),
            )
            .reset_index()
            .rename(columns={key: "gene_id"})
        )
        agg["assay"] = assay
        summaries.append(agg)
        peak = d.sort_values("ratio", ascending=False).drop_duplicates(key)
        peaks.append(
            peak[
                [
                    key,
                    "assay",
                    "tissue",
                    "cell_type",
                    "unit",
                    value,
                    cutoff,
                    "ratio",
                    "local_exclusion",
                    "reference_review",
                ]
            ].rename(columns={key: "gene_id", value: "value", cutoff: "candidate_ceiling"})
        )

    def assess(d: pd.DataFrame, assay: str) -> pd.DataFrame:
        d = normal_cell_screen.contexts(
            d.drop(columns=["category", "mapping_status", "concept_id"], errors="ignore")
        )
        d = policy.protected_immune_contexts(policy.reproductive_immune_contexts(d))
        d = policy.match_negative_ihc(
            d, ihc, positive_controls=controls, conflicting_genes=conflicts
        )
        d = policy.assess(
            d,
            policy.linked_thresholds(mode="two_axis_doubling"),
            ihc_factors=dict.fromkeys(["Missing/uncertain", *policy.GRADES], 1),
            inclusive_ceiling=True,
        )
        d["local_exclusion"] = d.cta_excluded | ~d.mapping_status.eq("mapped")
        d.loc[d.category.eq(policy.NEUTRAL_REPRODUCTIVE_IMMUNE), "local_exclusion"] = False
        return d

    # Never add per-gene consensus maxima from potentially different cohorts.
    bulk_parts = []
    for source in ["rna_tissue_hpa", "rna_tissue_gtex"]:
        path = source_dir / source
        with zipfile.ZipFile(path) as archive:
            raw = pd.read_csv(archive.open(archive.namelist()[0]), sep="\t")
        d = match(raw, source)
        d["cohort"] = source
        bulk_parts.append(sum_rows(d, ["Tissue", "cohort"], ["nTPM"], "nTPM"))
        receipts.append(
            {
                "name": source,
                "version": "v23",
                "url": f"https://v23.proteinatlas.org/download/{source}.tsv.zip",
                "sha256": sha256(path),
            }
        )
    components = pd.concat(bulk_parts, ignore_index=True)
    write_csv(components, output / "bulk-source-sums.csv.gz")
    bulk = consensus_bulk_regions(components)
    b = bulk_assessment(bulk.rename(columns={"gene_id": "Gene"}), ihc, controls, conflicts)
    write_csv(b, output / "bulk-contexts.csv.gz")
    collect(b, "Normal bulk", bulk=True)

    source = "hpa_single_cell_tissue"
    chunks = []
    for chunk in pd.read_csv(native(source), sep="\t", chunksize=250000):
        # Match every source identity before selecting the group loci.
        selected = match(chunk, source)
        if not selected.empty:
            chunks.append(selected)
    raw = pd.concat(chunks, ignore_index=True)
    d = sum_rows(raw, ["Tissue", "Cluster", "Cell type"], ["nTPM", "Read count"], "nTPM")
    d = d.rename(columns={"Tissue": "tissue", "Cell type": "cell_type", "nTPM": "value"})
    d["source"], d["version"], d["taxonomy_id"] = source, "v23", 9606
    d["state"], d["stage"] = "unspecified", "adult"
    d["value"] = d.value.where(d["Read count"].ge(20))
    assessed = assess(d, "HPA tissue clusters")
    write_csv(assessed, output / "cluster-contexts.csv.gz")
    collect(assessed, "HPA tissue clusters")
    print("Group clusters complete", flush=True)

    keys = ["source", "version", "taxonomy_id", "tissue", "cell_type", "state", "stage"]
    for source, label in [
        ("hpa_immune_cell", "HPA immune"),
        ("hpa_immune_monaco", "Monaco immune"),
    ]:
        native(source)
        d = normal_cell_data.observations(source, genes=member_ids)
        # Source observations already have stable IDs; retain an explicit audit.
        d = match(
            d.rename(columns={"gene_id": "source_gene_id"}), source, "source_gene_id", "gene_name"
        )
        d = sum_rows(d, keys, ["value"], "nTPM" if source == "hpa_immune_cell" else "TPM")
        assessed = assess(d, label)
        write_csv(assessed, output / f"{source}-contexts.csv.gz")
        collect(assessed, label)

    dice_rows = []
    for source in normal_cell_data.sources().source:
        if not source.startswith("dice_"):
            continue
        native(source)
        d = normal_cell_data.dice_observations(source, genes=member_ids)
        if d.empty:
            continue
        d = match(
            d.rename(columns={"gene_id": "source_gene_id"}), source, "source_gene_id", "gene_name"
        )
        for group in groups.itertuples():
            part = d[d.gene_id.isin(json.loads(group.member_gene_ids))]
            if part.empty:
                continue
            vectors = [
                p.value.to_numpy()
                for _, p in part.groupby("gene_id")
                if not p.replicate_slot.duplicated().any()
            ]
            if not vectors:
                continue
            low, high = unpaired_rank_bounds(vectors)
            if len(vectors) < group.member_count:
                high = np.inf
            row = part.iloc[0][keys].to_dict()
            row.update(
                gene_id=group.representative_gene_id,
                value=low,
                value_upper=high,
                unit="TPM",
                measured_members=len(vectors),
                member_count=group.member_count,
                bound_reason="unpaired gene distributions; no shared donor IDs",
            )
            dice_rows.append(row)
    dice = assess(pd.DataFrame(dice_rows), "DICE immune")
    # Positive lower bounds may exclude. A possible upper exceedance can only
    # cause a review hold, not manufacture strong evidence.
    dice["unpaired_sum_review"] = (
        dice.measured_members.eq(dice.member_count)
        & dice.value_upper.gt(dice.cta_cutoff)
        & ~dice.local_exclusion
    )
    dice["local_exclusion"] |= dice.unpaired_sum_review
    write_csv(dice, output / "dice-rank-bounds.csv.gz")
    collect(dice, "DICE immune")

    native("cornea_collin_adult_counts")
    native("cornea_collin_adultcornea_cells")
    cornea = normal_cell_data.cornea_summary(genes=member_ids)
    cornea = match(
        cornea.rename(columns={"gene_id": "source_gene_id"}),
        "cornea_collin_adult_counts",
        "source_gene_id",
        "gene_name",
    )
    summed = sum_rows(cornea, [*keys, "sample"], ["value", "read_count"], "CP10K")
    cornea["gene_id"] = cornea.gene_id.map(member_to_rep)
    detection = (
        cornea.groupby(["gene_id", *keys, "sample"], dropna=False)
        .agg(
            positive_cells_lower=("positive_cells", "max"),
            positive_cells_upper=("positive_cells", "sum"),
            cells=("cells", "min"),
        )
        .reset_index()
    )
    summed = summed.merge(detection, on=["gene_id", *keys, "sample"], validate="one_to_one")
    summed["positive_cells_upper"] = summed[["positive_cells_upper", "cells"]].min(axis=1)
    write_csv(summed, output / "cornea-detection-bounds.csv.gz")
    normalized = normal_cell_screen.contexts(summed)
    lower = summarize(
        normalized,
        support=normalized.positive_cells_lower.ge(5) & normalized.read_count.ge(20),
        rank=2,
        unique="sample",
    )
    upper = summarize(
        normalized,
        support=normalized.positive_cells_upper.ge(5) & normalized.read_count.ge(20),
        rank=2,
        unique="sample",
    )
    lower = assess(lower, "Adult cornea")
    upper = assess(upper, "Adult cornea")
    join_keys = ["gene_id", *keys]
    lower = lower.merge(
        upper[[*join_keys, "local_exclusion"]].rename(
            columns={"local_exclusion": "possible_supported_exclusion"}
        ),
        on=join_keys,
        validate="one_to_one",
    )
    lower["local_exclusion"] |= lower.possible_supported_exclusion
    write_csv(lower, output / "cornea-contexts.csv.gz")
    collect(lower, "Adult cornea")

    # Shared lineage means are additive even when individual-cell detections
    # cannot be recovered. Admission requires the conservative detection union.
    dev = pd.read_csv(ROOT / "docs/audits/cta-quantitative-20261026/developmental-lineages.csv.gz")
    dev = match(dev, "E-MTAB-9388", "Gene ID", None)
    dev["representative"] = dev.gene_id.map(member_to_rep)
    development = (
        dev.groupby(["representative", "Lineage"])
        .agg(
            mean_tpm=("Mean TPM", "sum"),
            cells=("Cells", "min"),
            positive_cells_lower=("Positive cells TPM>=1", "max"),
            positive_cells_upper=("Positive cells TPM>=1", "sum"),
            observed_members=("gene_id", "nunique"),
        )
        .reset_index()
        .rename(columns={"representative": "gene_id"})
    )
    development["positive_cells_upper"] = development[["positive_cells_upper", "cells"]].min(axis=1)
    development["supported"] = (
        development.cells.ge(20)
        & development.positive_cells_lower.ge(5)
        & (development.positive_cells_lower / development.cells).ge(0.1)
    )
    development["admitted"] = development.supported & development.mean_tpm.ge(5)
    write_csv(development, output / "developmental-lineage-bounds.csv.gz")
    write_csv(pd.concat(summaries, ignore_index=True), output / "normal-assessments.csv.gz")
    write_csv(pd.concat(peaks, ignore_index=True), output / "normal-peak-evidence.csv.gz")
    identities = pd.concat(identity_audits, ignore_index=True).drop_duplicates(
        ["source", "source_gene_id", "source_gene_name"]
    )
    write_csv(identities, output / "source-gene-identities.csv.gz")
    (output / "manifest.json").write_text(
        json.dumps(
            {
                "sources": receipts,
                "inputs": input_hashes(
                    [
                        Path(__file__),
                        ROOT / "oncoref/protein_equivalence.py",
                        ROOT / "oncoref/gene_matching.py",
                        ROOT / "oncoref/cell_rna_policy.py",
                        baseline / "gene-classifications.csv.gz",
                        ROOT
                        / "docs/audits/cta-quantitative-20261026/developmental-lineages.csv.gz",
                    ]
                ),
                "outputs": {
                    p.name: sha256(p)
                    for p in output.iterdir()
                    if p.is_file() and p.name != "manifest.json"
                },
            },
            indent=2,
        )
        + "\n"
    )
    print(f"Completed {len(groups)} protein equivalence groups", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument(
        "--source-dir", type=Path, default=Path("/tmp/oncoref-cta-protein-group-sources")
    )
    parser.add_argument(
        "--baseline", type=Path, default=ROOT / "docs/audits/cta-expanded-funnel-20261027"
    )
    args = parser.parse_args()
    build(args.output, args.source_dir, args.baseline)
