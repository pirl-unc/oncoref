#!/usr/bin/env python3
"""Grade rare normal-cell RNA per adult donor from the pinned Census count build.

Replaces the pooled 1% detection screen for every protein that the count build
measured. The panel keeps the earlier scope (brain, spinal cord, heart and
selected lung populations; immune, stromal, generic and developmental labels
excluded) so the change is in the evidence standard, not in coverage. Adult age,
primary tissue and donor identity now come from metadata, not cell-type names.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from build_cta_rare_cell_screen import (  # noqa: E402
    EXCLUDED_CELLS,
    GENERIC,
    LUNG_VITAL,
    MAX_SHARE,
    VITAL_TISSUES,
)

from oncoref import target_selection  # noqa: E402
from oncoref.cta_cell_tiers import cell_tiers, differentiated_airway  # noqa: E402
from oncoref.tissue_hazards import registry as hazard_registry  # noqa: E402

DATA = ROOT / "oncoref/data"
AUDIT = ROOT / "docs/audits/normal-cell-atlas-20261110/rare-cell"
CUTOFFS = (1, 2, 5)
SILETTI = ROOT / "docs/audits/normal-cell-atlas-20261110/siletti-granule-gene-detection.csv.gz"
#: Census collection of the same study, so its donors are not a second study.
SILETTI_SOURCE = "283d65eb-dd53-496d-adb7-7570c7caa443"
RISK = DATA / "gene-mapping-risk.csv.gz"
LITERATURE = DATA / "normal-expression-literature.json"
BRAIN_BULK = {
    "amygdala", "basal ganglia", "caudate", "cerebellum", "cerebral cortex", "choroid plexus",
    "hippocampal formation", "hippocampus", "hypothalamus", "medulla oblongata", "midbrain",
    "pons", "substantia nigra", "thalamus", "white matter",
}  # fmt: skip
ORGAN_BULK = {"heart muscle": "heart", "lung": "lung", "spinal cord": "spinal cord"}


def digest(path):
    with Path(path).open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def write_csv(frame, path):
    compression = {"method": "gzip", "mtime": 0} if str(path).endswith(".gz") else None
    frame.to_csv(path, index=False, compression=compression)


def panel(populations: pd.DataFrame) -> pd.DataFrame:
    """Adult, donor-resolved UMI populations inside the reviewed rare-cell panel."""
    d = populations[
        populations.age_scope.eq("adult")
        & populations.donor_identity_known
        & populations.count_unit.eq("UMI_protocol")
    ].copy()
    d["tissue"] = d.tissue_general.astype(str)
    d["cell_type"] = d.cell_type.astype(str)
    tissue_cells = d.groupby("tissue").n_cells.sum()
    type_cells = d.groupby(["tissue", "cell_type"]).n_cells.transform("sum")
    airway = differentiated_airway(d.tissue, d.cell_type)
    in_tissue = d.tissue.isin(VITAL_TISSUES) | (
        d.tissue.eq("lung") & (d.cell_type.str.contains(LUNG_VITAL) | airway)
    )
    keep = (
        in_tissue
        & ~d.cell_type.isin(GENERIC)
        & ~d.cell_type.str.contains(EXCLUDED_CELLS)
        & type_cells.le(MAX_SHARE * d.tissue.map(tissue_cells))
    )
    d = d[keep].copy()
    d["tier"] = cell_tiers(d.tissue, d.cell_type)
    d = d[d.tier.isin(["vital", "protected"])]
    d["source_id"] = d.collection_id.astype(str)
    d["donor_id"] = d.donor_id.astype(str)
    return d


def donor_rows(counts: pd.DataFrame, scoped: pd.DataFrame) -> pd.DataFrame:
    """One row per gene, donor and named population; donors are never pooled."""
    keys = ["source_id", "donor_id", "tissue", "cell_type", "tier"]
    sizes = scoped.groupby(keys, as_index=False).n_cells.sum()
    joined = counts.merge(scoped[["population_id", *keys]], on="population_id")
    rows = joined.groupby(["gene_id", *keys], as_index=False).agg(
        n_positive=("n_positive", "sum"),
        n_multi=("n_ge2", "sum"),
        count_sum=("count_sum", "sum"),
        count_max=("count_max", "max"),
    )
    return rows.merge(sizes, on=keys, validate="many_to_one")


def siletti_rows() -> pd.DataFrame:
    """Granule nuclei from the directly analysed dissections.

    Census labels the primary Siletti observations with generic terms, which
    the panel excludes. The dissection files carry the specific granule label.
    """
    d = pd.read_csv(SILETTI)
    d = d[d.scope.ne("all")]
    rows = pd.DataFrame(
        {
            "gene_id": d.gene_id,
            "source_id": SILETTI_SOURCE,
            "donor_id": d.scope,
            "tissue": "brain",
            "cell_type": "cerebellar granule cell",
            "n_positive": d.n_positive,
            "n_multi": d.n_ge2,
            "count_sum": d.count_sum,
            "count_max": d.max_count,
            "n_cells": d.n_cells,
        }
    )
    rows["tier"] = cell_tiers(rows.tissue, rows.cell_type)
    return rows[rows.n_positive.gt(0)]


def bulk_confirmation(proteins: pd.DataFrame) -> pd.DataFrame:
    """Bulk TPM per protein and panel organ reached in both HPA and GTEx samples."""
    import rebuild_target_catalog as rebuild

    base = pd.read_csv(rebuild.BASE / "protein-group-classifications.csv.gz", low_memory=False)
    facts = base.rename(columns={"protein_equivalence_group_id": "entity_id"})
    facts = facts[facts.entity_id.isin(proteins.entity_id)]
    normal, _ = rebuild.normal_panel(facts, write=False)
    normal = normal[normal.qc_pass & normal.identity_resolved & normal.value.notna()].copy()
    normal["organ"] = normal.tissue.map(lambda t: "brain" if t in BRAIN_BULK else ORGAN_BULK.get(t))
    normal = normal[normal.organ.notna() & normal.cell_type.eq("bulk")]
    # One high subregion among hundreds is not confirmation. Require both
    # independent sample sets; report the weaker of their two maxima.
    normal["dataset"] = normal.source_id.map(lambda v: "HPA" if "_hpa:" in v else "GTEx")
    top = normal.sort_values("value", ascending=False).drop_duplicates(
        ["entity_id", "organ", "dataset"]
    )
    both = top.groupby(["entity_id", "organ"]).dataset.transform("nunique").eq(2)
    weaker = top[both].sort_values("value").drop_duplicates(["entity_id", "organ"])
    return weaker[["entity_id", "organ", "value", "source_id", "source_context"]].rename(
        columns={"value": "bulk_tpm", "source_id": "bulk_source", "source_context": "bulk_context"}
    )


def published_organs() -> dict[str, set[str]]:
    """Panel organs where reviewed literature reports the gene by a specific assay."""
    entries = json.loads(LITERATURE.read_text())["entries"]
    return {
        e["entity_id"]: {e["decision_tissue"]}
        for e in entries
        if e["effect"] == "confirm_detection"
    }


def assess_proteins(proteins, rows, measured_genes, assessable, risk, bulk):
    published = published_organs()
    by_gene = dict(iter(rows.groupby("gene_id")))
    floor = target_selection.DETECTION_FLOOR
    confirmed = bulk[bulk.bulk_tpm.ge(floor)].groupby("entity_id").organ.agg(set).to_dict()
    bulk_by = bulk.set_index(["entity_id", "organ"])
    records, findings = [], []
    for protein in proteins.itertuples():
        members = json.loads(protein.member_gene_ids)
        receipts = []
        for gene in members:
            if gene not in measured_genes:
                continue
            flags = risk.loc[gene]
            receipt = target_selection.assess_rare_cell_rna(
                by_gene.get(gene, rows.iloc[:0]).drop(columns=["count_sum", "count_max"]),
                mapping_risk=bool(flags.mapping_risk),
                confirmed_tissues=confirmed.get(protein.entity_id, set()),
                literature_tissues=published.get(protein.entity_id, set()),
            )
            if receipt["status"] == "unassessed":
                # Measured in the panel, with no positive donor population to grade.
                receipt["status"] = "below_screen"
            receipt["gene_id"] = gene
            receipts.append(receipt)
        rank = {"established": 0, "dismissed": 1, "below_screen": 2}
        receipts.sort(key=lambda r: (rank[r["status"]], -(r["max_multi_molecule_pct"] or 0)))
        top = receipts[0] if receipts else None
        finding = top["findings"][0] if top and top["findings"] else {}
        flags = risk.loc[top["gene_id"]] if top else None
        key = (protein.entity_id, finding.get("tissue"))
        support = bulk_by.loc[key] if key in bulk_by.index else None
        notes = []
        if finding and flags.paralog_ambiguous:
            notes.append("sequence shared with another locus")
        if finding and flags.host_overlapped:
            notes.append("exon inside another same-strand gene")
        records.append(
            {
                "entity_id": protein.entity_id,
                "protein_group_label": protein.protein_group_label,
                "member_gene_ids": protein.member_gene_ids,
                "rare_cell_status": top["status"] if top else "unassessed",
                "rare_cell_action": top["action"] if top else "none",
                "rare_cell_tier": top["tier"] if top else "",
                "rare_cell_gene_id": top["gene_id"] if top else "",
                "rare_cell_basis": finding.get("basis", ""),
                "rare_cell_mapping_notes": " | ".join(notes),
                "measured_member_gene_ids": json.dumps([g for g in members if g in measured_genes]),
                "assessable_donor_populations": assessable if top else 0,
                "max_detection_pct": top["max_detection_pct"] if top else None,
                "max_multi_molecule_pct": top["max_multi_molecule_pct"] if top else None,
                "paralog_ambiguous": bool(flags.paralog_ambiguous) if top else None,
                "host_overlapped": bool(flags.host_overlapped) if top else None,
                "shared_kmer_fraction": flags.shared_kmer_fraction if top else None,
                "sequence_partners": flags.sequence_partners if top else None,
                "same_strand_hosts": flags.same_strand_hosts if top else None,
                **{
                    f"finding_{k}": v
                    for k, v in finding.items()
                    if k not in {"status", "tier", "basis"}
                },
                "finding_bulk_tpm": support.bulk_tpm if support is not None else None,
                "finding_bulk_context": (
                    f"{support.bulk_source}: {support.bulk_context}"
                    if support is not None
                    else None
                ),
                "n_findings": sum(len(r["findings"]) for r in receipts),
            }
        )
        for receipt in receipts:
            for item in receipt["findings"]:
                findings.append(
                    {
                        "entity_id": protein.entity_id,
                        "protein_group_label": protein.protein_group_label,
                        "gene_id": receipt["gene_id"],
                        "mapping_risk": receipt["mapping_risk"],
                        **item,
                    }
                )
    return pd.DataFrame(records), pd.DataFrame(findings)


def sensitivity(rows, proteins):
    """Reported comparisons at 1%, 2% and 5%; these are not category rules."""
    gene_to_entity = {
        g: p.entity_id for p in proteins.itertuples() for g in json.loads(p.member_gene_ids)
    }
    d = rows[rows.n_positive.ge(target_selection.RARE_CELL_MIN_CELLS)].copy()
    d["entity_id"] = d.gene_id.map(gene_to_entity)
    d["pct"] = 100 * d.n_positive / d.n_cells
    d["pct_multi"] = 100 * d.n_multi / d.n_cells
    d["donor"] = d.source_id + ":" + d.donor_id
    out = []
    for tier, cutoff in [(t, c) for t in ("vital", "protected") for c in CUTOFFS]:
        hit = d[d.tier.eq(tier) & d.pct.ge(cutoff)]
        studies = hit.groupby(["entity_id", "tissue"]).source_id.nunique()
        multi = d[d.tier.eq(tier) & d.pct_multi.ge(cutoff) & d.n_multi.ge(10)]
        multi_studies = multi.groupby(["entity_id", "tissue"]).source_id.nunique()
        out.append(
            {
                "tier": tier,
                "detection_cutoff_pct": cutoff,
                "proteins_any_donor": hit.entity_id.nunique(),
                "proteins_two_or_more_studies": studies[studies.ge(2)]
                .index.get_level_values(0)
                .nunique(),
                "proteins_multi_molecule_any_donor": multi.entity_id.nunique(),
                "proteins_multi_molecule_two_or_more_studies": multi_studies[multi_studies.ge(2)]
                .index.get_level_values(0)
                .nunique(),
            }
        )
    return pd.DataFrame(out)


def render(evidence, findings, table, provenance):
    """Evidence page; category consequences are reported by the catalog rebuild."""
    rules = provenance["rules"]
    scope = provenance["panel"]
    grades = evidence.groupby(["rare_cell_status", "rare_cell_tier", "rare_cell_action"]).size()
    grade_rows = "\n".join(
        f"| {status} | {tier or '—'} | {action} | {n:,} |"
        for (status, tier, action), n in grades.items()
    )
    sens = "\n".join(
        f"| {r.tier} | {r.detection_cutoff_pct}% | {r.proteins_any_donor:,} | "
        f"{r.proteins_two_or_more_studies:,} | {r.proteins_multi_molecule_any_donor:,} | "
        f"{r.proteins_multi_molecule_two_or_more_studies:,} |"
        for r in table.itertuples()
    )

    def rows(frame):
        return "\n".join(
            f"| {r.protein_group_label} | {r.finding_tissue}: {r.finding_top_cell_type} | "
            f"{r.finding_studies_multi_molecule:.0f} / {r.finding_donors_multi_molecule:.0f} | "
            f"{r.finding_top_n_multi:,.0f}/{r.finding_top_n_cells:,.0f} "
            f"({r.finding_top_multi_molecule_pct:.1f}%) | "
            f"{'—' if pd.isna(r.finding_bulk_tpm) else f'{r.finding_bulk_tpm:.2g}'} | "
            f"{r.rare_cell_basis}"
            + (
                f" ({r.rare_cell_mapping_notes})"
                if isinstance(r.rare_cell_mapping_notes, str) and r.rare_cell_mapping_notes
                else ""
            )
            + " |"
            for r in frame.itertuples()
        )

    head = (
        "| Protein | Organ: strongest population | Studies / donors with multi-molecule cells | "
        "Strongest donor: cells with ≥2 molecules | Bulk TPM reached in both HPA and GTEx | Basis |\n"
        "|---|---|---:|---:|---:|---|\n"
    )
    unsafe = evidence[evidence.rare_cell_status.eq("established")].sort_values(
        "finding_top_multi_molecule_pct", ascending=False
    )
    dismissed = evidence[evidence.rare_cell_status.eq("dismissed")].sort_values(
        "finding_top_multi_molecule_pct", ascending=False
    )
    risk = provenance["mapping_risk"]
    page = f"""# Rare normal-cell RNA by adult donor

**{len(evidence):,} proteins were graded from {scope["cells"]:,} adult cells and nuclei of
{scope["donors"]:,} donors in {scope["collections"]} collections, plus the directly analysed
Siletti cerebellar granule nuclei. {len(unsafe)} have established normal-cell RNA and are
Unsafe; {len(dismissed)} have a signal that is dismissed with its reason.** Category
consequences are in the [catalog rebuild](../../target-selection-20261110/rebuild/index.md).

A count matrix cannot show where a read aligned. One mis-assigned read per cell, from a
paralog, a pseudogene or an overlapping gene, looks exactly like sparse expression. The
grade therefore ignores single-molecule cells and asks five questions. Every outcome is
terminal; nothing waits for review.

| Question | Requirement |
|---|---|
| Is there a signal? | {rules["signal"]} {rules["protected"]} |
| Is it replicated? | {rules["replication"]} {rules["single_study"]} |
| Could the reads belong elsewhere? | {rules["mapping_risk"]} |
| Does another assay agree? | {rules["confirmation"]} |
| Is it published? | {rules["literature"]} |

| Outcome | Rule |
|---|---|
| Established | {rules["established"]} |
| Dismissed | {rules["dismissed"]} |
| Below the screen | No replicated or frequent multi-molecule signal. No finding. |

Below the screen is not a demonstration of absence or safety. The
[literature review](../../target-selection-20261110/rebuild/literature.md) covers every
candidate whose signal is dismissed or who left the pooled flag.
Donors are never pooled, frequencies are never added across paralogs, and unmeasured
genes stay unassessed.

## Grades

| Grade | Tier | Catalog action | Proteins |
|---|---|---|---:|
{grade_rows}

## Mapping risk

Of {risk["genes"]:,} extracted genes, {risk["paralog_ambiguous"]:,} share at least
{risk["shared_fraction_flag"]:.0%} of their transcript {risk["k"]}-mers, or of their 3′ kilobase,
with another locus, and {risk["host_overlapped"]:,} have an exon inside another same-strand
gene; {risk["mapping_risk"]:,} carry either flag. Ensembl {risk["ensembl_release"]}, sense strand,
primary chromosomes. [Flags, shared fractions and partner loci for every gene](gene-mapping-risk.csv.gz).
A flag marks a possible
confounder, not a demonstrated mis-assignment, and an unflagged gene can still collect
intronic reads.

## Sensitivity to the detection cutoff

Reported comparisons only, with at least ten positive cells in the donor population.
The first two columns count any molecule; they show why single-molecule detection is
not graded.

| Tier | Cutoff | Proteins above it in any donor | In two or more studies | With ≥2 molecules per cell, any donor | With ≥2 molecules per cell, two or more studies |
|---|---:|---:|---:|---:|---:|
{sens}

## Established

{head}{rows(unsafe)}

## Dismissed

{head}{rows(dismissed)}

## Scope and limits

The panel keeps the scope of the earlier pooled screen: brain, spinal cord, heart and
selected lung populations, without immune, stromal, generic, unknown or developmental
labels, and without populations above {scope["max_tissue_share"]:.0%} of a tissue's sampled
cells. Adult age, primary tissue and donor identity come from metadata. Only UMI protocols
are graded. Census release {provenance["census_release"]}; filter `{provenance["filter"]}`.

- {provenance["limitations"][0]}.
- {provenance["limitations"][1]}.
- {provenance["limitations"][2]}.
- {provenance["limitations"][3]}; their other exclusions stand.

[Every finding](findings.csv) · [Donor populations with at least ten positive cells](donor-populations.csv.gz) ·
[Panel](panel.csv) · [Cutoff comparison](threshold-sensitivity.csv) ·
[Proteins selected for extraction](selected-proteins.csv.gz)

```sh
python scripts/extract_rare_cell_panel_counts.py --counts tmp/normal-atlas/counts-2025-11-08
python scripts/build_gene_mapping_risk.py --proteins tmp/normal-atlas/counts-2025-11-08/selected-proteins.csv.gz
python scripts/build_rare_cell_evidence.py --counts tmp/normal-atlas/counts-2025-11-08
```
"""
    (AUDIT / "index.md").write_text(page)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--counts", type=Path, required=True)
    args = parser.parse_args()
    manifest = json.loads((args.counts / "manifest.json").read_text())
    features = pd.read_csv(args.counts / "selected-features.csv")
    extract = json.loads((args.counts / "panel-counts/manifest.json").read_text())
    blocks = sorted((args.counts / "panel-counts").glob("part-*.parquet"))
    if len(blocks) != extract["chunks"]:
        raise ValueError(f"Count build incomplete: {len(blocks)}/{extract['chunks']} chunks")
    counts = (
        pd.concat(pd.read_parquet(p) for p in blocks)
        .groupby(["population_id", "gene_id"], as_index=False)
        .agg(
            n_positive=("n_positive", "sum"),
            n_ge2=("n_ge2", "sum"),
            count_sum=("count_sum", "sum"),
            count_max=("count_max", "max"),
            n_fractional=("n_fractional", "sum"),
        )
    )
    populations = pd.read_parquet(args.counts / "populations.parquet")
    scoped = panel(populations)
    if extract["n_cells"] != int(scoped.n_cells.sum()) or extract["n_genes"] != len(features):
        raise ValueError("Extracted cells or genes differ from the current panel")
    if not set(counts.population_id) <= set(scoped.population_id):
        raise ValueError("Counts outside the panel")
    # Extracted with the panel, but an unlabelled population is not a named one.
    scoped = scoped[scoped.cell_type.ne("unknown")]
    floor = target_selection.RARE_CELL_MIN_CELLS * 100 / target_selection.RARE_CELL_DETECTION_PCT
    keys = ["source_id", "donor_id", "tissue", "cell_type"]
    sizes = scoped.groupby(keys).n_cells.sum()
    assessable = int(sizes.ge(floor).sum())
    proteins = pd.read_csv(args.counts / "selected-proteins.csv.gz")[
        ["entity_id", "protein", "member_gene_ids", "category", "Unsafe_reasons"]
    ].rename(
        columns={
            "protein": "protein_group_label",
            "category": "draft_category",
            "Unsafe_reasons": "draft_unsafe_reasons",
        }
    )
    wanted = {g for ids in proteins.member_gene_ids for g in json.loads(ids)}
    extra = siletti_rows()
    extra = extra[extra.gene_id.isin(wanted)]
    assessable += int(pd.read_csv(SILETTI).query("scope != 'all'").scope.nunique())
    rows = pd.concat([donor_rows(counts, scoped), extra], ignore_index=True)
    if rows.duplicated(["gene_id", *keys]).any():
        raise ValueError("A donor population is supplied by two sources")
    risk = pd.read_csv(RISK).set_index("gene_id")
    risk_provenance = json.loads((DATA / "gene-mapping-risk.json").read_text())
    if risk_provenance["sha256"] != digest(RISK):
        raise ValueError("Mapping-risk table does not match its provenance")
    AUDIT.mkdir(parents=True, exist_ok=True)
    write_csv(risk.reset_index(), AUDIT / "gene-mapping-risk.csv.gz")
    bulk = bulk_confirmation(proteins)
    write_csv(bulk, AUDIT / "bulk-confirmation.csv.gz")
    evidence, findings = assess_proteins(
        proteins, rows, set(features.feature_id) | set(extra.gene_id), assessable, risk, bulk
    )
    AUDIT.mkdir(parents=True, exist_ok=True)
    write_csv(evidence, DATA / "rare-cell-evidence.csv.gz")
    write_csv(findings, AUDIT / "findings.csv")
    detected = rows[rows.n_positive.ge(target_selection.RARE_CELL_MIN_CELLS)]
    write_csv(detected, AUDIT / "donor-populations.csv.gz")
    table = sensitivity(rows, proteins)
    write_csv(table, AUDIT / "threshold-sensitivity.csv")
    write_csv(proteins, AUDIT / "selected-proteins.csv.gz")
    summary = (
        scoped.groupby(["tissue", "tier"])
        .agg(
            cells=("n_cells", "sum"),
            cell_types=("cell_type", "nunique"),
            collections=("source_id", "nunique"),
        )
        .reset_index()
    )
    write_csv(summary, AUDIT / "panel.csv")
    provenance = {
        "schema": "oncoref.rare-cell-evidence",
        "policy_id": target_selection.POLICY_ID,
        "census_release": manifest["census_release"],
        "filter": manifest["filter"],
        "count_manifest_sha256": digest(args.counts / "manifest.json"),
        "panel_extract": extract,
        "count_parts_sha256": hashlib.sha256(
            "".join(digest(p) for p in blocks).encode()
        ).hexdigest(),
        "fractional_count_cells": int(counts.n_fractional.sum()),
        "supplement": {
            "source_id": SILETTI_SOURCE,
            "path": str(SILETTI.relative_to(ROOT)),
            "sha256": digest(SILETTI),
            "reason": "Specific granule labels exist only in the dissection files",
        },
        "hazard_registry_id": hazard_registry()["registry_id"],
        "rules": target_selection.definition()["rare_cell_rna"],
        "mapping_risk": risk_provenance,
        "confirmation_sources": sorted(bulk.bulk_source.unique()),
        "literature_sha256": digest(LITERATURE),
        "panel": {
            "scope": "adult, known donor, UMI protocol, primary normal tissue",
            "vital_tissues": sorted(VITAL_TISSUES),
            "lung_selector": LUNG_VITAL.pattern,
            "excluded_cell_types": EXCLUDED_CELLS.pattern,
            "excluded_generic_terms": sorted(GENERIC),
            "max_tissue_share": MAX_SHARE,
            "cells": int(scoped.n_cells.sum()),
            "donors": int(scoped.groupby(["source_id", "donor_id"]).ngroups),
            "collections": int(scoped.source_id.nunique()),
            "assessable_donor_populations": assessable,
        },
        "proteins": len(evidence),
        "unavailable_gene_ids": manifest["unavailable_genes"],
        "status_counts": evidence.rare_cell_status.value_counts().to_dict(),
        "action_counts": evidence.rare_cell_action.value_counts().to_dict(),
        "limitations": [
            "A gene-count matrix cannot check alignment, strand or exonic support; the "
            "mapping-risk flags are reference properties, not read-level checks",
            "No empty-droplet or ambient-RNA correction is applied to the pinned counts",
            "Donor overlap across collections has not been reviewed",
            "Proteins that were Unsafe for another reason were not extracted",
        ],
        "sha256": digest(DATA / "rare-cell-evidence.csv.gz"),
    }
    (DATA / "rare-cell-evidence.json").write_text(json.dumps(provenance, indent=2) + "\n")
    render(evidence, findings, table, provenance)
    print(json.dumps({k: provenance[k] for k in ["panel", "status_counts", "action_counts"]}))


if __name__ == "__main__":
    main()
