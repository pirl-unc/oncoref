#!/usr/bin/env python3
"""Replay current selection from pinned measurements without modifying releases."""

from __future__ import annotations

import argparse
import json
import sys
import zipfile
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from cta_compartment_report import group_ihc_evidence  # noqa: E402
from cta_quantitative_threshold_report import positive_profile_controls  # noqa: E402
from cta_report_common import input_hashes, sha256, write_csv  # noqa: E402

from oncoref import cell_rna_policy, hpa, normal_cell_screen, reference_data  # noqa: E402
from oncoref import cta_selection_policy as policy  # noqa: E402
from oncoref.cta_funnel import EvidenceLayer  # noqa: E402
from oncoref.gene_matching import (  # noqa: E402
    match_reference_genes,
    resolve_long_measurement_ids,
    resolve_measurement_rows,
)

FROZEN = ROOT / "docs/audits/cta-release-1.21.0"
BASE = ROOT / "docs/audits/cta-harmonized-20261027"
OUTPUT = ROOT / "docs/audits/cta-selection-20261028"
CROSSWALK = ROOT / "oncoref/data/cta-ihc-compartments.csv"
ALTERNATIVE = OUTPUT / "alternative-evidence"
ASSAYS = [
    "Normal bulk",
    "HPA tissue clusters",
    "HPA immune",
    "Monaco immune",
    "DICE immune",
    "Adult cornea",
]


def read(path: Path) -> pd.DataFrame:
    """Read a pinned evidence artifact, preserving missing numeric values."""
    return pd.read_csv(path, low_memory=False)


class ContextCollector:
    """Reduce full context streams without losing breadth, units or exclusions."""

    def __init__(self, facts: pd.DataFrame, output: Path):
        self.facts, self.output = facts, output
        self.summaries: list[pd.DataFrame] = []
        self.peaks: list[pd.DataFrame] = []
        self.breadth: list[pd.DataFrame] = []
        self.comparators: list[pd.DataFrame] = []
        self.mapped: list[pd.DataFrame] = []

    def collect(self, frame: pd.DataFrame, assay: str) -> None:
        d = policy.assess_contexts(frame)
        d["assay"] = assay
        if "value_upper" in d:
            complete = d.measured_members.eq(d.member_count)
            d["candidate_exclusion"] |= complete & d.value_upper.gt(d.candidate_ceiling)
            d["unresolved_context"] |= (
                complete & d.value_upper.gt(d.candidate_ceiling) & ~d.value.gt(d.candidate_ceiling)
            )
        if "upper_supported_value" in d:
            possible = d.upper_supported_value.gt(d.candidate_ceiling)
            d["candidate_exclusion"] |= possible
            d["unresolved_context"] |= possible & ~d.value.gt(d.candidate_ceiling)
        normal = d.selection_tier.isin(policy.TIER_BASES)
        raw = d.get("raw_maximum", d.value).fillna(d.value)
        d["normal_maximum_floor_ratio"] = raw / (
            policy.APPRECIABLE_RNA * d.unit.map(policy.UNIT_SCALE)
        )
        d["healthy_maximum_ntpm"] = d.value.where(normal & d.unit.eq("nTPM"))
        d["healthy_maximum_tpm"] = d.value.where(normal & d.unit.eq("TPM"))
        d["healthy_maximum_cp10k"] = d.value.where(normal & d.unit.eq("CP10K"))
        d["ratio"] = d.value / d.candidate_ceiling
        self.summaries.append(
            d.groupby("gene_id")
            .agg(
                candidate_exclusion=("candidate_exclusion", "any"),
                unsafe_evidence=("unsafe_evidence", "any"),
                unresolved_context=("unresolved_context", "any"),
                normal_maximum_floor_ratio=("normal_maximum_floor_ratio", "max"),
                healthy_maximum_ntpm=("healthy_maximum_ntpm", "max"),
                healthy_maximum_tpm=("healthy_maximum_tpm", "max"),
                healthy_maximum_cp10k=("healthy_maximum_cp10k", "max"),
            )
            .reset_index()
            .assign(assay=assay)
        )
        # Keep separate worst Unsafe and worst non-Unsafe contexts. The maximum
        # ratio alone can hide a second, Vital/Protected exclusion.
        chosen = d[d.value.notna() & d.candidate_ceiling.notna()].sort_values(
            "ratio", ascending=False
        )
        self.peaks.append(chosen.drop_duplicates(["gene_id", "unsafe_evidence"]))
        manageable = d.selection_tier.eq("manageable") & d.value.ge(
            policy.APPRECIABLE_RNA * d.unit.map(policy.UNIT_SCALE)
        )
        part = d.loc[manageable, ["gene_id", "tissue"]].copy()
        part["tissue"] = part.tissue.replace({"rectum": "large bowel", "colon": "large bowel"})
        self.breadth.append(part.drop_duplicates())
        names = set(
            self.facts.loc[
                self.facts.protein_group_label.isin(
                    [
                        "PRAME",
                        "CTAG2",
                        "MAGEA1",
                        "MAGEA3",
                        "MAGEA4",
                        "MAGEA12",
                        "NY-ESO-1",
                        "CT83",
                        "TBXT",
                        "XAGE1A/B",
                        "SSX2/B",
                        "SMC1B",
                        "PAGE1",
                    ]
                ),
                "representative_gene_id",
            ]
        )
        self.comparators.append(d[d.gene_id.isin(names)])
        if "ihc_match_scope" in d:
            self.mapped.append(d[d.ihc_match_scope.eq("reviewed_compartment") & d.value.notna()])

    def finish(self) -> pd.DataFrame:
        summary = (
            pd.concat(self.summaries)
            .groupby(["gene_id", "assay"])
            .agg(
                {
                    "candidate_exclusion": "any",
                    "unsafe_evidence": "any",
                    "unresolved_context": "any",
                    "normal_maximum_floor_ratio": "max",
                    "healthy_maximum_ntpm": "max",
                    "healthy_maximum_tpm": "max",
                    "healthy_maximum_cp10k": "max",
                }
            )
            .reset_index()
        )
        peaks = (
            pd.concat(self.peaks, ignore_index=True)
            .sort_values("ratio", ascending=False)
            .drop_duplicates(["gene_id", "assay", "unsafe_evidence"])
        )
        write_csv(summary, self.output / "assay-assessments.csv.gz")
        write_csv(peaks, self.output / "normal-peak-evidence.csv.gz")
        write_csv(
            pd.concat(self.comparators, ignore_index=True),
            self.output / "widely-studied-contexts.csv.gz",
        )
        write_csv(
            pd.concat(self.mapped, ignore_index=True),
            self.output / "compartment-ihc-matches.csv.gz",
        )
        breadth = pd.concat(self.breadth).drop_duplicates().groupby("gene_id").size()
        self.facts["manageable_tissue_count"] = (
            self.facts.representative_gene_id.map(breadth).fillna(0).astype(int)
        )
        return summary


def replay_contexts(facts: pd.DataFrame, output: Path) -> pd.DataFrame:
    """Stream all HPA clusters; reassess all pinned bulk/immune/corneal rows."""
    collector = ContextCollector(facts, output)
    # Bulk IHC requires complete routine tissue coverage, as in the frozen audit.
    old_bulk = read(FROZEN / "bulk-contexts.csv.gz")
    member_to_rep = {
        g: r.representative_gene_id
        for r in facts.itertuples()
        for g in json.loads(r.member_gene_ids)
    }
    raw_bulk = hpa._read_hpa("hpa_rna_consensus", "v23")
    lookup = (
        read(output / "hpa_rna_consensus-accepted-identities.csv.gz")
        .set_index("source_gene_id")
        .canonical_gene_id
    )
    raw_bulk["member"] = raw_bulk.Gene.map(lookup)
    raw_bulk = raw_bulk[raw_bulk.member.isin(member_to_rep)].copy()
    if raw_bulk.duplicated(["member", "Tissue"]).any():
        raise ValueError("Duplicate bulk measurement vectors require resolution")
    raw_bulk["gene_id"] = raw_bulk.member.map(member_to_rep)
    bulk = (
        raw_bulk.groupby(["gene_id", "Tissue"])
        .nTPM.sum(min_count=1)
        .reset_index()
        .rename(columns={"Tissue": "tissue", "nTPM": "value"})
    )
    bulk = normal_cell_screen.contexts(
        bulk.assign(
            source="hpa_rna_consensus",
            version="v23",
            unit="nTPM",
            taxonomy_id=9606,
            stage="adult",
            cell_type="bulk",
        )
    )
    bulk = bulk.merge(
        old_bulk[["gene_id", "tissue", "negative_ihc_grade"]],
        on=["gene_id", "tissue"],
        how="left",
        validate="one_to_one",
    )
    bulk["negative_ihc_grade"] = bulk.negative_ihc_grade.fillna("Missing")
    collector.collect(bulk, "Normal bulk")
    alternative = read(ALTERNATIVE / "normal-bulk.csv.gz")
    # Exact whole-organ correspondence permits the same routine bulk IHC grade
    # across RNA assays. Regional brain peaks retain Missing: a parent organ's
    # IHC cannot establish negative staining in every anatomical subregion.
    grades = old_bulk.set_index(["gene_id", "tissue"]).negative_ihc_grade
    equivalent_bulk = alternative.source.eq("gtex_tissue_rna")
    keys = pd.MultiIndex.from_frame(alternative[["gene_id", "tissue"]])
    alternative.loc[equivalent_bulk, "negative_ihc_grade"] = (
        grades.reindex(keys).fillna("Missing").to_numpy()[equivalent_bulk]
    )
    collector.collect(alternative, "Normal bulk")
    grouped_ihc, grouped_controls, grouped_conflicts = group_ihc_evidence(facts, output)
    ihc = hpa._read_hpa("hpa_normal_tissue", "v23").dropna(subset=["Tissue", "Cell type", "Level"])
    controls = positive_profile_controls(ihc)
    conflicts = set(facts.loc[facts["Antibody concern"], "representative_gene_id"])
    single = set(facts.loc[facts.member_count.eq(1), "representative_gene_id"])
    crosswalk = read(CROSSWALK)
    receipt = reference_data.provenance("hpa_single_cell_tissue", verify_content=True)
    lookup = (
        read(output / "hpa_single_cell_tissue-accepted-identities.csv.gz")
        .set_index("source_gene_id")
        .canonical_gene_id
    )
    for number, chunk in enumerate(pd.read_csv(receipt["path"], sep="\t", chunksize=250000)):
        chunk["source_gene_id"] = chunk.Gene
        chunk["Gene"] = chunk.Gene.map(lookup)
        d = (
            chunk[chunk.Gene.isin(single)]
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
        if d.empty:
            continue
        d = d.assign(
            source="hpa_single_cell_tissue",
            version="v23",
            taxonomy_id=9606,
            state="unspecified",
            stage="adult",
            unit="nTPM",
        )
        d["raw_maximum"] = d.value
        d["value"] = d.value.where(d["Read count"].ge(20))
        d = cell_rna_policy.protected_immune_contexts(
            cell_rna_policy.reproductive_immune_contexts(normal_cell_screen.contexts(d))
        )
        collector.collect(
            policy.match_ihc(
                d,
                ihc[ihc.Gene.isin(d.gene_id.unique())],
                positive_controls=controls,
                conflicting_genes=conflicts,
                crosswalk=crosswalk,
            ),
            "HPA tissue clusters",
        )
        if number % 10 == 0:
            print(f"Replayed cluster chunk {number}", flush=True)
    d = read(BASE / "group-normal/cluster-contexts.csv.gz")
    collector.collect(
        policy.match_ihc(
            d,
            grouped_ihc,
            positive_controls=grouped_controls,
            conflicting_genes=grouped_conflicts,
            crosswalk=crosswalk,
        ),
        "HPA tissue clusters",
    )
    for assay, d in read(FROZEN / "immune-cornea-contexts.csv.gz").groupby("assay"):
        if assay in {"HPA immune", "Monaco immune"}:
            continue
        collector.collect(d, assay)
    member_to_rep = {
        g: r.representative_gene_id
        for r in facts.itertuples()
        for g in json.loads(r.member_gene_ids)
    }
    for source, assay, unit in [
        ("hpa_immune_cell", "HPA immune", "nTPM"),
        ("hpa_immune_monaco", "Monaco immune", "TPM"),
    ]:
        raw = hpa._read_hpa(source, "v23")
        lookup = (
            read(output / f"{source}-accepted-identities.csv.gz")
            .set_index("source_gene_id")
            .canonical_gene_id
        )
        raw["member"] = raw.Gene.map(lookup)
        raw = raw[raw.member.isin(member_to_rep)].copy()
        # Duplicate aliases of a locus cannot count twice in a group sum.
        if raw.duplicated(["member", "Immune cell"]).any():
            raise ValueError(f"Duplicate immune measurement vectors require resolution: {source}")
        raw["gene_id"] = raw.member.map(member_to_rep)
        d = (
            raw.groupby(["gene_id", "Immune cell"])[unit]
            .sum(min_count=1)
            .reset_index()
            .rename(columns={"Immune cell": "cell_type", unit: "value"})
        )
        d = d.assign(
            source=source,
            version="v23",
            unit=unit,
            taxonomy_id=9606,
            stage="adult",
            tissue="blood",
            mapping_status="mapped",
            category="protected",
            negative_ihc_grade="Missing",
        )
        collector.collect(d, assay)
    collector.collect(read(FROZEN / "dice-group-contexts.csv.gz"), "DICE immune")
    collector.collect(read(FROZEN / "cornea-group-contexts.csv.gz"), "Adult cornea")
    (output / "cluster-source.json").write_text(
        json.dumps({k: receipt[k] for k in ["name", "version", "url", "sha256"]}, indent=2) + "\n"
    )
    return collector.finish()


def protein_layers(facts: pd.DataFrame, output: Path) -> list[EvidenceLayer]:
    """Retain reviewed protein eligibility; assign risk using measured anatomy."""
    raw = hpa._read_hpa("hpa_normal_tissue", "v23").dropna(subset=["Tissue", "Cell type", "Level"])
    d = raw.rename(
        columns={"Gene": "gene_id", "Tissue": "tissue", "Cell type": "cell_type"}
    ).assign(source="hpa_normal_tissue", version="v23", taxonomy_id=9606)
    d = normal_cell_screen.contexts(d)
    d["selection_tier"] = policy.compartment_tiers(d)
    normal = d.category.isin(policy.TIER_BASES)
    positive = d.Level.isin(cell_rna_policy.POSITIVE)
    nonlow = positive & ~d.Level.eq("Low") & normal
    d["candidate_exclusion"] = nonlow
    d["unsafe_evidence"] = (
        nonlow
        & d.Reliability.isin(cell_rna_policy.GRADES)
        & d.selection_tier.isin(["vital", "protected"])
    )
    d["unresolved_evidence"] = nonlow & ~d.Reliability.isin(cell_rna_policy.GRADES)
    # Preserve explicit earlier adjudications of unreliable staining. This mask
    # concerns protein evidence eligibility, not a previous final category.
    allowed = facts.set_index("representative_gene_id").ihc_detected
    member_to_rep = {
        g: r.representative_gene_id
        for r in facts.itertuples()
        for g in json.loads(r.member_gene_ids)
    }
    d["representative_gene_id"] = d.gene_id.map(member_to_rep)
    eligible = d.representative_gene_id.map(allowed).fillna(False).astype(bool)
    for col in ["candidate_exclusion", "unsafe_evidence", "unresolved_evidence"]:
        d[col] &= eligible
    summary = d.groupby("representative_gene_id")[
        ["candidate_exclusion", "unsafe_evidence", "unresolved_evidence"]
    ].any()

    def lookup(column):
        return facts.representative_gene_id.map(summary[column]).fillna(False).astype(bool)

    facts["unresolved_evidence"] |= lookup("unresolved_evidence")
    detected = set(d.loc[positive, "representative_gene_id"].dropna())
    facts["any_detected_protein"] = facts.representative_gene_id.isin(detected)
    write_csv(d[positive | d.candidate_exclusion], output / "hpa-protein-evidence.csv.gz")
    layers = [
        EvidenceLayer(
            "HPA normal protein", lookup("candidate_exclusion"), lookup("unsafe_evidence")
        )
    ]
    imported = read(FROZEN / "imported-protein-review-decisions.csv")
    imported["unsafe_evidence"] = imported.retained_protein_exclusion & imported[
        "Tissue category"
    ].isin(["vital", "protected"])
    members = facts.member_gene_ids.map(json.loads)

    def imported_mask(column: str) -> pd.Series:
        ids = set(imported.loc[imported[column], "Gene ID"])
        return members.map(lambda group: bool(set(group) & ids))

    layers.append(
        EvidenceLayer(
            "Accepted normal-protein evidence",
            imported_mask("retained_protein_exclusion"),
            imported_mask("unsafe_evidence"),
        )
    )
    detected_ids = set(imported["Gene ID"])
    facts["any_detected_protein"] |= members.map(lambda group: bool(set(group) & detected_ids))
    write_csv(imported, output / "imported-protein-review-decisions.csv")
    # Eligible external observations already require exact cell coverage and a
    # validated reagent. Resolve their anatomical tier independently of intensity.
    external = read(ROOT / "docs/audits/external-ihc-20261027/assessed-evidence.csv.gz")
    external = external[external.eligible & external.result.eq("detected")].copy()
    if external.empty:
        empty = pd.Series(False, index=facts.index)
        layers.append(EvidenceLayer("External validated normal protein", empty, empty))
    else:
        external = normal_cell_screen.contexts(
            external.assign(source="hpa_normal_tissue", version="v23")
        )
        members = facts.member_gene_ids.map(json.loads)
        hold_ids = set(external.loc[external.category.isin(policy.TIER_BASES), "gene_id"])
        unsafe_ids = set(external.loc[external.category.isin(["vital", "protected"]), "gene_id"])
        layers.append(
            EvidenceLayer(
                "External validated normal protein",
                members.map(lambda ids: bool(set(ids) & hold_ids)),
                members.map(lambda ids: bool(set(ids) & unsafe_ids)),
            )
        )
        facts["any_detected_protein"] |= members.map(
            lambda ids: bool(set(ids) & set(external.gene_id))
        )
    write_csv(external, output / "external-protein-evidence.csv.gz")
    return layers


def developmental_coverage(facts: pd.DataFrame, output: Path) -> None:
    """Recover verified missing IDs and retain conflicting alternate-contig rows.

    Canonical measurements take precedence. A missing canonical gene can use a
    unique registered alias, or identical alias vectors; distinct alias vectors
    remain unresolved. Protein-identical loci are still separate before sums.
    """
    from scipy.io import mmread

    d = read(ROOT / "docs/audits/cta-quantitative-20261026/developmental-lineages.csv.gz")
    directory = ROOT / "tmp/cta-quantitative"
    receipt = json.loads(
        (ROOT / "docs/audits/cta-quantitative-20261026/developmental-source.json").read_text()
    )
    for filename in ["atlas_tpm", "atlas_design"]:
        if sha256(directory / filename) != receipt["files"][filename]["sha256"]:
            raise ValueError(f"Developmental source checksum mismatch: {filename}")
    (output / "developmental-source.json").write_text(json.dumps(receipt, indent=2) + "\n")
    with zipfile.ZipFile(directory / "atlas_tpm") as archive:
        name = next(n for n in archive.namelist() if n.endswith(".mtx"))
        gene_rows = pd.read_csv(
            archive.open(name + "_rows"), sep="\t", header=None, names=["gene_id", "label"]
        )
        cells = pd.read_csv(archive.open(name + "_cols"), sep="\t", header=None)[0]
        identities = match_reference_genes(
            gene_rows,
            id_column="gene_id",
            source="tyser_E-MTAB-9388",
            source_version="Ensembl106_filtered_TPM",
        )
        member_ids = {g for ids in facts.member_gene_ids.map(json.loads) for g in ids}
        candidates = identities[
            identities.identity_status.eq("matched")
            & identities.canonical_gene_id.isin(member_ids - set(d["Gene ID"]))
        ]
        matrix = mmread(archive.open(name)).tocsr()
    meta = pd.read_csv(directory / "atlas_design", sep="\t").set_index("Assay").loc[cells]
    lineages = meta.groupby("Sample Characteristic[inferred cell type - authors labels]").indices
    additions = []
    decisions = resolve_measurement_rows(
        candidates, equivalent=lambda a, b: (matrix[a] != matrix[b]).nnz == 0
    )
    for row in decisions.itertuples():
        gid, chosen = row.gene_id, row.row_position
        if pd.isna(chosen):
            continue
        chosen = int(chosen)
        for lineage, indices in lineages.items():
            values = matrix[chosen, indices]
            count = int((values >= 1).sum())
            additions.append(
                {
                    "Gene ID": gid,
                    "Lineage": lineage,
                    "Cells": len(indices),
                    "Positive cells TPM>=1": count,
                    "Fraction TPM>=1": count / len(indices),
                    "Mean TPM": float(values.mean()),
                    "Embryos": 1,
                }
            )
    write_csv(decisions, output / "developmental-identity-recoveries.csv")
    d = pd.concat([d, pd.DataFrame(additions)], ignore_index=True)
    d = d.assign(source="tyser_E-MTAB-9388", stage="Carnegie stage 7", unit="TPM")
    d = pd.concat([d, read(ALTERNATIVE / "developmental-lineages.csv.gz")], ignore_index=True)
    write_csv(d, output / "developmental-lineages.csv.gz")
    member_to_rep = {
        g: r.representative_gene_id
        for r in facts.itertuples()
        for g in json.loads(r.member_gene_ids)
    }
    d["representative"] = d["Gene ID"].map(member_to_rep)
    groups = (
        d.groupby(["representative", "source", "stage", "Lineage"])
        .agg(
            mean_tpm=("Mean TPM", "sum"),
            cells=("Cells", "min"),
            positive_cells_lower=("Positive cells TPM>=1", "max"),
        )
        .reset_index()
    )
    groups["supported"] = (
        groups.cells.ge(20)
        & groups.positive_cells_lower.ge(5)
        & (groups.positive_cells_lower / groups.cells).ge(0.1)
    )
    best = groups[groups.supported].groupby("representative").mean_tpm.max()
    facts["developmental_mean_tpm"] = facts.representative_gene_id.map(best)
    facts["developmental_route"] = facts.developmental_mean_tpm.ge(5)
    facts["positive_context"] = facts.reproductive_route | facts.developmental_route
    facts["developmental_all_lineage_max_tpm"] = facts.representative_gene_id.map(
        groups.groupby("representative").mean_tpm.max()
    )
    members = facts.member_gene_ids.map(json.loads)
    measured_ids = set(d["Gene ID"])
    facts["developmental_measured"] = members.map(lambda ids: set(ids) <= measured_ids)


def required_rna_coverage(facts: pd.DataFrame, output: Path) -> None:
    """Require anatomical coverage across compatible sources, not one export.

    Reproductive tissues do not supply normal-organ exclusion coverage. Every
    group member must have finite measurements for every required normal tissue;
    all native subregions of an alternative tissue must be present. Tissue-cell
    and immune panels retain their distinct, non-interchangeable coverage roles.
    """
    from audit_cta_source_coverage import PANELS

    complete_by_source, rows = {}, []
    bulk_seen = {}
    for source, _title, columns, unit, expected in PANELS:
        receipt = reference_data.provenance(source, "v23", verify_content=True)
        values = [unit, "Read count"] if source == "hpa_single_cell_tissue" else [unit]
        d = pd.read_csv(receipt["path"], sep="\t", usecols=["Gene", "Gene name", *columns, *values])
        accepted, decisions = resolve_long_measurement_ids(
            d,
            id_column="Gene",
            symbol_column="Gene name",
            context_columns=columns,
            value_columns=values,
            source=source,
            source_version="v23",
        )
        write_csv(accepted, output / f"{source}-accepted-identities.csv.gz")
        write_csv(decisions, output / f"{source}-identity-decisions.csv.gz")
        lookup = accepted.set_index("source_gene_id").canonical_gene_id
        d["gene_id"] = d.Gene.map(lookup)
        if source in {"hpa_single_cell_tissue", "hpa_rna_consensus"}:
            tissue = (
                d[["Tissue", *(["Cell type"] if "Cell type" in d else [])]]
                .drop_duplicates()
                .rename(columns={"Tissue": "tissue", "Cell type": "cell_type"})
            )
            if "cell_type" not in tissue:
                tissue["cell_type"] = "bulk"
            tissue = normal_cell_screen.contexts(
                tissue.assign(source=source, version="v23", taxonomy_id=9606, stage="adult")
            )
            allowed = tissue[~tissue.category.isin(cell_rna_policy.PERMITTED)]
            keys = ["Tissue", "Cell type"] if "Cell type" in d else ["Tissue"]
            d = d.merge(
                allowed.rename(columns={"tissue": "Tissue", "cell_type": "Cell type"})[keys],
                on=keys,
                validate="many_to_one",
            )
        d["context"] = d[columns[0]].astype(str)
        for column in columns[1:]:
            d["context"] = d.context + " / " + d[column].astype(str)
        universe = set(d.context)
        if source != "hpa_rna_consensus" and len(universe) != expected:
            raise ValueError(f"Changed required source contexts: {source}")
        seen = d[d[unit].notna()].groupby("gene_id").context.agg(set).to_dict()
        complete_by_source[source] = {g for g, contexts in seen.items() if contexts == universe}
        if source == "hpa_rna_consensus":
            bulk_seen, bulk_universe = seen, universe
            bulk_before = complete_by_source[source].copy()
        print(f"Resolved required coverage: {source}", flush=True)
    for source in ["hpa_brain", "gtex_tissue_rna"]:
        coverage = read(ALTERNATIVE / f"{source}-member-coverage.csv.gz")
        for row in coverage.itertuples():
            bulk_seen.setdefault(row.gene_id, set()).update(
                set(json.loads(row.available_tissues)) & bulk_universe
            )
    complete_by_source["hpa_rna_consensus"] = {
        g for g, seen in bulk_seen.items() if bulk_universe <= seen
    }
    for row in facts.itertuples():
        ids = set(json.loads(row.member_gene_ids))
        missing = sorted(set().union(*(bulk_universe - bulk_seen.get(g, set()) for g in ids)))
        rows.append(
            {
                "entity_id": row.protein_equivalence_group_id,
                "label": row.protein_group_label,
                "required_bulk_tissues": len(bulk_universe),
                "missing_bulk_tissues": json.dumps(missing),
                "bulk_complete_before": ids <= bulk_before,
                "bulk_complete": not missing,
                "required_coverage": all(
                    ids <= measured for measured in complete_by_source.values()
                ),
            }
        )
    coverage = pd.DataFrame(rows)
    facts["required_coverage"] = coverage.required_coverage.to_numpy()
    write_csv(coverage, output / "combined-rna-coverage.csv.gz")


def build(output: Path = OUTPUT, *, reuse_contexts: bool = False) -> None:
    """Regenerate decisions, complete profile evidence and all counting exports."""
    output.mkdir(parents=True, exist_ok=True)
    facts = read(FROZEN / "protein-group-classifications.csv.gz")
    # Unnamed canonical loci retain their stable ID as the display label.
    facts["protein_group_label"] = facts.protein_group_label.fillna(facts.representative_gene_id)
    facts["previous_category"] = facts.category
    alternative_receipt = json.loads((ALTERNATIVE / "manifest.json").read_text())
    for name, digest in alternative_receipt["outputs"].items():
        if sha256(ALTERNATIVE / name) != digest:
            raise ValueError(f"Changed alternative evidence: {name}")
    required_rna_coverage(facts, output)
    repaired = output / "tumor-identity-repair"
    receipt = json.loads((repaired / "manifest.json").read_text())
    for name, expected in receipt["outputs"].items():
        if sha256(repaired / name) != expected:
            raise ValueError(f"Changed repaired tumor evidence: {name}")
    for multi, filename, key in [
        (False, "tumor-maxima.csv.gz", "gene_id"),
        (True, "protein-group-tumor-maxima.csv.gz", "protein_equivalence_group_id"),
    ]:
        revised = read(repaired / filename).set_index(key)
        identifiers = facts.protein_equivalence_group_id if multi else facts.representative_gene_id
        maximum = identifiers.map(revised.maximum_tpm)
        better = (
            facts.member_count.gt(1).eq(multi)
            & maximum.notna()
            & (facts.oncoref_maximum_tpm.isna() | maximum.gt(facts.oncoref_maximum_tpm))
        )
        for target, source in [
            ("oncoref_maximum_tpm", "maximum_tpm"),
            ("oncoref_maximum_sample", "maximum_sample"),
            ("oncoref_cancer_code", "cancer_code"),
        ]:
            facts.loc[better, target] = identifiers.map(revised[source])[better]
    facts["tumor_union_ge2"] = facts.oncoref_maximum_tpm.ge(
        2
    ) | facts.hpa_tumor_group_lower_bound_ptpm.ge(2)
    developmental_coverage(facts, output)
    inputs = [
        Path(__file__),
        ROOT / "oncoref/cta_selection_policy.py",
        ROOT / "oncoref/gene_matching.py",
        CROSSWALK,
        FROZEN / "protein-group-classifications.csv.gz",
        FROZEN / "bulk-contexts.csv.gz",
        FROZEN / "immune-cornea-contexts.csv.gz",
        FROZEN / "dice-group-contexts.csv.gz",
        FROZEN / "cornea-group-contexts.csv.gz",
        BASE / "group-normal/cluster-contexts.csv.gz",
        ROOT / "oncoref/normal_cell_screen.py",
        ROOT / "oncoref/cell_rna_policy.py",
        ROOT / "oncoref/normal_tissues.py",
        ROOT / "oncoref/compartment_rna_policy.py",
        ROOT / "oncoref/data/normal-tissue-taxonomy.yaml",
        ROOT / "oncoref/data/normal-tissue-mappings.csv",
        ROOT / "scripts/cta_compartment_report.py",
        ROOT / "scripts/cta_quantitative_threshold_report.py",
        ALTERNATIVE / "manifest.json",
    ]
    sources = []
    for name in ["hpa_single_cell_tissue", "hpa_normal_tissue"]:
        source = reference_data.provenance(name, verify_content=True)
        sources.append({key: source[key] for key in ["name", "version", "url", "sha256"]})
    cache = output / "context-replay.json"
    if reuse_contexts:
        receipt = json.loads(cache.read_text())
        if receipt["inputs"] != input_hashes(inputs) or receipt["sources"] != sources:
            raise ValueError("Context inputs changed; replay without --reuse-contexts")
        for name, digest in receipt["outputs"].items():
            if sha256(output / name) != digest:
                raise ValueError("Changed cached context output")
        summary = read(output / "assay-assessments.csv.gz")
        breadth = (
            read(output / "manageable-breadth.csv.gz").set_index("gene_id").manageable_tissue_count
        )
        facts["manageable_tissue_count"] = facts.representative_gene_id.map(breadth)
    else:
        replay_contexts(facts, output)
        # Both modes use the same persisted numeric representation.
        summary = read(output / "assay-assessments.csv.gz")
        write_csv(
            facts[["representative_gene_id", "manageable_tissue_count"]].rename(
                columns={"representative_gene_id": "gene_id"}
            ),
            output / "manageable-breadth.csv.gz",
        )
        names = [
            "assay-assessments.csv.gz",
            "normal-peak-evidence.csv.gz",
            "widely-studied-contexts.csv.gz",
            "compartment-ihc-matches.csv.gz",
            "manageable-breadth.csv.gz",
        ]
        cache.write_text(
            json.dumps(
                {
                    "inputs": input_hashes(inputs),
                    "sources": sources,
                    "outputs": {n: sha256(output / n) for n in names},
                },
                indent=2,
            )
            + "\n"
        )
    facts["unresolved_evidence"] = (
        facts.representative_gene_id.map(summary.groupby("gene_id").unresolved_context.any())
        .fillna(False)
        .astype(bool)
    )
    bulk_maximum = (
        summary[summary.assay.eq("Normal bulk")].set_index("gene_id").healthy_maximum_ntpm
    )
    facts["adult_bulk_max_ntpm"] = facts.representative_gene_id.map(bulk_maximum)
    facts["reproductive_route"] = facts.reproductive_ntpm.ge(5) & facts.reproductive_ntpm.gt(
        facts.adult_bulk_max_ntpm
    )
    facts["positive_context"] = facts.reproductive_route | facts.developmental_route
    for col in [
        "normal_maximum_floor_ratio",
        "healthy_maximum_ntpm",
        "healthy_maximum_tpm",
        "healthy_maximum_cp10k",
    ]:
        facts[col] = facts.representative_gene_id.map(summary.groupby("gene_id")[col].max())
    empty = pd.Series(False, index=facts.index)
    layers = [
        EvidenceLayer("Reproductive / embryonic association", ~facts.positive_context, empty),
        EvidenceLayer("Tumor sample ≥2 TPM/pTPM", ~facts.tumor_union_ge2, empty),
        EvidenceLayer("Required RNA coverage", ~facts.required_coverage, empty),
    ]
    proteins = protein_layers(facts, output)
    layers.insert(2, proteins[0])
    for assay in ASSAYS:
        part = summary[summary.assay.eq(assay)].set_index("gene_id")
        masks = [
            facts.representative_gene_id.map(part[col]).fillna(False).astype(bool)
            for col in ["candidate_exclusion", "unsafe_evidence"]
        ]
        layers.append(EvidenceLayer(assay, *masks))
    layers.extend(
        [
            proteins[1],
            EvidenceLayer("Clinical safety exclusion", facts["Safety veto"], facts["Safety veto"]),
            proteins[2],
        ]
    )
    facts["vital_protected_exclusion"] = pd.DataFrame(
        {name: unsafe for name, _, unsafe in layers}
    ).any(axis=1)
    tumor = facts[["oncoref_maximum_tpm", "hpa_tumor_group_lower_bound_ptpm"]]
    facts["tumor_measured"] = tumor.notna().any(axis=1)
    # A measured group lower bound is enough for admission but cannot prove low
    # expression of an incompletely measured group. Low profiles require singletons.
    facts["tumor_profile_complete"] = facts.tumor_measured & facts.member_count.eq(1)
    facts["tumor_maximum_floor_ratio"] = tumor.max(axis=1) / policy.APPRECIABLE_RNA
    facts = policy.neither_profiles(facts)
    result, funnel, independent = policy.classify(
        facts, layers, entity_column="protein_equivalence_group_id"
    )
    reverse = policy.classify(
        facts, list(reversed(layers)), entity_column="protein_equivalence_group_id"
    ).classifications
    if not result.category.equals(reverse.category):
        raise ValueError("Final assignments depend on layer order")
    write_csv(result, output / "protein-group-classifications.csv.gz")
    write_csv(funnel, output / "filtering-funnel.csv")
    write_csv(independent, output / "independent-layer-effects.csv")
    write_csv(
        result[result.category.fillna("Unassigned").ne(result.previous_category)],
        output / "classification-changes.csv.gz",
    )
    for name in ["known-cta-decisions.csv.gz", "gene-to-protein-group-decisions.csv.gz"]:
        d = read(FROZEN / name)
        lookup = result.set_index("protein_equivalence_group_id")
        for target, col in [
            ("protein_group_category", "category"),
            ("protein_group_exclusion_reasons", "exclusion_reasons"),
        ]:
            d[target] = d.protein_equivalence_group_id.map(lookup[col])
        d["policy"] = policy.POLICY_ID
        write_csv(d, output / name)
    counts = result.category.fillna("Unassigned").value_counts().to_dict()
    meta = {
        "policy": policy.POLICY_ID,
        "categories": counts,
        "selected_member_loci": int(result.loc[result.category.eq("CTA"), "member_count"].sum()),
        "neither_profiles": result.loc[result.category.eq("Neither"), "neither_profile"]
        .value_counts()
        .to_dict(),
        "order_independent": True,
    }
    (output / "summary.json").write_text(json.dumps(meta, indent=2) + "\n")
    print(json.dumps(meta), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reuse-contexts", action="store_true")
    args = parser.parse_args()
    build(reuse_contexts=args.reuse_contexts)
