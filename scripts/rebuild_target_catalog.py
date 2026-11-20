#!/usr/bin/env python3
"""Rebuild antigen memberships from pinned normal and per-sample tumor evidence."""

from __future__ import annotations

import gzip
import hashlib
import json
import re
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from antigen_sample_evidence import (  # noqa: E402
    HPA_SOURCE,
    LACTATION_MARKER_TPM,
    LACTATION_MARKERS,
    LACTATION_MIN_MARKERS,
    developmental_samples,
    tissue_library_qc,
    tissue_sample_genes,
    tissue_samples,
)
from build_cta_alternative_evidence import GTEX_TISSUES, anatomy  # noqa: E402
from build_cta_candidate_evidence import replay_updates  # noqa: E402
from cta_report_common import sha256, write_csv  # noqa: E402
from gtex_sample_evidence import SOURCE as GTEX_SOURCE  # noqa: E402
from gtex_sample_evidence import individual_samples  # noqa: E402

from oncoref import (  # noqa: E402
    cta_literature,
    cta_sequence_overlap,
    library_qc,
    normal_library_qc,
    reference_data,
    source_matrices,
    target_selection,
    tissue_tiers,
)
from oncoref import (  # noqa: E402
    samples as sample_api,
)
from oncoref.gene_matching import (  # noqa: E402
    harmonize_reference_expression,
    match_reference_genes,
)

BASE = ROOT / "docs/audits/cta-selection-20261028"
FROZEN = ROOT / "docs/audits/cta-transcript-selection-20261108/inputs/frozen"
OUT = ROOT / "docs/audits/target-selection-20261110/rebuild"
CACHE = ROOT / "tmp/target-catalog-cache"
COMPARISON_ID = "oncoref.native-gene-tpm-normal-reference.2026-11-20.1"
INPUTS = {}
#: Pooled 1% detection layers, replaced by donor-resolved rare-cell evidence.
LEGACY_RARE = {"Rare vital-tissue cell expression", "Rare Protected-cell expression review"}
RARE_UNSAFE = "Established rare normal-cell RNA"
LITERATURE_UNSAFE = "Published normal-tissue expression"
#: Checks that must pass before tumor evidence and expression bounds are assessed.
PRIOR_CHECKS = [c for c in target_selection.SHARED_CHECKS if c != "expression_bounds_pass"]
DISPOSABLE_NTPM_TISSUES = ["appendix", "gallbladder", "tonsil"]


def read(path, **kwargs):
    path = Path(path)
    INPUTS[str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else str(path.resolve())] = (
        sha256(path)
    )
    return pd.read_csv(path, low_memory=False, **kwargs)


def read_reference(name, version, **kwargs):
    path = reference_data.ensure(name, version)
    INPUTS[f"reference:{name}:{version}"] = sha256(path)
    return pd.read_csv(path, low_memory=False, **kwargs)


def shared_facts():
    for name in [
        "gene-knowledge.json",
        "gene-knowledge-updates.json",
        "normal-tissue-tiers.yaml",
        "normal-tissue-tier-updates.yaml",
        "readthrough-loci.csv",
        "nose-skin-reference.csv",
        "transcript-selection.csv.gz",
        "translation-qualification.csv.gz",
    ]:
        path = ROOT / "oncoref/data" / name
        INPUTS[str(path.relative_to(ROOT))] = sha256(path)
    baseline = read(BASE / "protein-group-classifications.csv.gz")
    f = replay_updates(baseline).classifications
    f["released_category"] = baseline.category
    f["released_classification"] = baseline.classification
    f["legacy_exclusion_reasons"] = f.exclusion_reasons
    f["curation_replay_category"] = f.category
    f["entity_id"] = f.protein_equivalence_group_id
    layers = [c[7:] for c in f if c.startswith("unsafe:")]
    if not f.category.eq("Unsafe").eq(f[["unsafe:" + n for n in layers]].any(axis=1)).all():
        raise ValueError("Unsafe calls do not decompose into their evidence layers")
    gate_receipt = OUT / "revised-tissue-gate-inputs.json"
    INPUTS[str(gate_receipt.relative_to(ROOT))] = sha256(gate_receipt)
    revised = read(OUT / "revised-tissue-gates.csv")
    for assay, group in revised.groupby("assay"):
        if (
            group.gene_id.duplicated().any()
            or not group.gene_id.isin(f.representative_gene_id).all()
        ):
            raise ValueError(f"Invalid tissue gate identities for {assay}")
        for prefix, field in [("hold:", "candidate_exclusion"), ("unsafe:", "unsafe_evidence")]:
            values = f.representative_gene_id.map(group.set_index("gene_id")[field])
            f.loc[values.notna(), prefix + assay] = values.dropna().astype(bool)
    for name in ["rare-cell-evidence.csv.gz", "rare-cell-evidence.json"]:
        INPUTS[f"oncoref/data/{name}"] = sha256(ROOT / "oncoref/data" / name)
    provenance = json.loads((ROOT / "oncoref/data/rare-cell-evidence.json").read_text())
    if provenance["sha256"] != INPUTS["oncoref/data/rare-cell-evidence.csv.gz"]:
        raise ValueError("Rare-cell evidence does not match its provenance")
    if provenance["policy_id"] != target_selection.POLICY_ID:
        raise ValueError("Rare-cell evidence was graded under another policy")
    rare = pd.read_csv(ROOT / "oncoref/data/rare-cell-evidence.csv.gz").set_index("entity_id")
    kept = [n for n in layers if n not in LEGACY_RARE]
    other_unsafe = f[["unsafe:" + n for n in kept]].any(axis=1)
    missing_rare = ~other_unsafe & ~f.entity_id.isin(rare.index)
    if not f.loc[missing_rare, "representative_gene_id"].isin(revised.gene_id).all():
        raise ValueError("A protein without another exclusion lacks rare-cell evidence")
    rare_columns = [c for c in rare if c not in {"protein_group_label", "member_gene_ids"}]
    f = f.drop(columns=[c for c in rare_columns if c in f]).join(rare[rare_columns], on="entity_id")
    f["rare_cell_status"] = f.rare_cell_status.fillna("not_extracted")
    f["rare_cell_action"] = f.rare_cell_action.fillna("none")
    f.loc[missing_rare, "rare_cell_action"] = "review"
    # A pooled flag that the donor build could not measure stays unresolved.
    f["rare_cell_legacy_only"] = (
        f["hold:Rare vital-tissue cell expression"]
        & f.rare_cell_status.eq("unassessed")
        & ~other_unsafe
    )
    path = ROOT / "oncoref/data/normal-expression-literature.json"
    INPUTS["oncoref/data/normal-expression-literature.json"] = sha256(path)
    if provenance["literature_sha256"] != sha256(path):
        raise ValueError("Rare-cell evidence was graded with another literature review")
    document = json.loads(path.read_text())
    literature = pd.DataFrame(document["entries"]).set_index("entity_id")
    for field in [
        "verdict",
        "effect",
        "rule",
        "decision_tissue",
        "decision_reason",
        "recorded_tissue",
        "recorded_tier",
    ]:
        f["literature_" + field] = f.entity_id.map(literature[field])
    f["literature_effect"] = f.literature_effect.fillna("not_reviewed")
    # A family rule names proteins by label and records the tissue the family shares.
    f["literature_family"], f["literature_family_tier"] = "", ""
    family_tissue = pd.Series("", index=f.index)
    for rule in document["family_rules"]:
        member = f.protein_group_label.str.contains(rule["pattern"]) & f.literature_family.eq("")
        f.loc[member, "literature_family"] = rule["id"]
        f.loc[member, "literature_family_tier"] = rule["tier"]
        # Every family names its home tissue; filter_acts governs only the filters' tier.
        family_tissue[member] = rule["tissue"]
    # Nose and skin references record a tissue on the row.
    reference = pd.read_csv(
        ROOT / "oncoref/data/nose-skin-reference.csv", keep_default_na=False
    ).set_index("entity_id")
    reference = reference[reference.recorded_tissue.ne("")]
    f["reference_recorded_tissue"] = f.entity_id.map(reference.recorded_tissue).fillna("")
    f["reference_recorded_basis"] = f.entity_id.map(reference.recorded_basis).fillna("")
    f["reference_recorded_tier"] = f.reference_recorded_tissue.map(
        lambda text: next(
            (
                tier
                for tier in tissue_tiers.TIERS
                if tier in {tissue_tiers.assess(t)["tier"] for t in text.split("; ") if t}
            ),
            "",
        )
    )
    # A reviewed adult home has the same coverage consequence regardless of publication.
    homes = defaultdict(set)
    for source in cta_literature.publication_registry()["sources"]:
        for claim in source["gene_claims"]:
            if "unscreened_adult_home" in claim:
                for gene_id in claim["gene_ids"]:
                    homes[gene_id].add(claim["unscreened_adult_home"])
    f["published_home_tissue"] = f.member_gene_ids.map(
        lambda ids: "; ".join(
            sorted({home for gene in json.loads(ids) for home in homes.get(gene, ())})
        )
    )
    # The recorded home, when it is a tissue no panel resolves (see classify).
    recorded = (
        family_tissue
        + "; "
        + f.literature_recorded_tissue.fillna("")
        + "; "
        + f.reference_recorded_tissue
        + "; "
        + f.published_home_tissue
    )
    f["home_tissue"] = recorded.map(
        lambda text: "; ".join(
            dict.fromkeys(
                t
                for t in (part.split(":")[0].strip() for part in text.split(";"))
                if t in target_selection.UNMEASURED_HOME_TISSUES
            )
        )
    )
    f["home_tissue_unmeasured"] = f.home_tissue.ne("")
    readthrough = (
        pd.read_csv(ROOT / "oncoref/data/readthrough-loci.csv")
        .fillna({"parent_genes": ""})
        .groupby("entity_id")
        .agg(parents=("parent_genes", "; ".join), basis=("basis", "first"))
    )
    f["readthrough_locus"] = f.entity_id.isin(readthrough.index)
    f["readthrough_parents"] = f.entity_id.map(readthrough.parents).fillna("")
    literature_unsafe = f.literature_effect.eq("exclude")
    rare_unsafe = f.rare_cell_action.eq("unsafe")
    rare_review = f.rare_cell_action.eq("review") | f.rare_cell_legacy_only
    f["unsafe"] = other_unsafe | rare_unsafe | literature_unsafe
    f["Unsafe_reasons"] = [
        " | ".join(
            [n for n, on in zip(kept, row) if on]
            + ([RARE_UNSAFE] if rare_hit else [])
            + ([LITERATURE_UNSAFE] if published else [])
        )
        for row, rare_hit, published in zip(
            f[["unsafe:" + n for n in kept]].to_numpy(), rare_unsafe, literature_unsafe
        )
    ]
    ninemers = cta_sequence_overlap.group_ninemers(f.member_gene_ids)
    overlap = cta_sequence_overlap.unsafe_overlap(ninemers, f.index)
    f["released_unsafe_sequence_overlap_pct"] = f.unsafe_sequence_overlap_pct
    f["unsafe_sequence_overlap_pct"] = (100 * overlap(f.unsafe)).round(2)
    admission = {
        "Reproductive / embryonic association",
        "Cancer expression",
        "Reproductive concentration",
        "Translation qualification",
        "Required RNA coverage",
    }
    normal_holds = [c for c in f if c.startswith("hold:") and c[5:] not in admission | LEGACY_RARE]
    f["normal_screen_pass"] = ~f[normal_holds].any(axis=1) & ~rare_review
    f["required_coverage_complete"] = f.required_coverage
    f["identity_translation_pass"] = f[
        [
            "translation_support_eligible",
            "translation_qualified",
            "representative_selection_approved",
        ]
    ].all(axis=1)
    f["evidence_resolved"] = ~f.unresolved_evidence
    f["sequence_pass"] = (
        ~f.nested_microprotein & f.unsafe_sequence_overlap_pct.lt(95) & ~f.readthrough_locus
    )
    f["shared_eligible"] = ~f.unsafe & f[PRIOR_CHECKS].all(axis=1)
    # Retain measured counterfactuals for newly held rare-cell signals. A review
    # hold still blocks membership through normal_screen_pass in classify().
    f["comparison_eligible"] = (
        ~f.unsafe
        & f[[c for c in PRIOR_CHECKS if c != "normal_screen_pass"]].all(axis=1)
        & ~f[normal_holds].any(axis=1)
        & ~f.rare_cell_legacy_only
    )
    branches = {
        "reproductive_bulk": f.reproductive_ntpm.ge(1),
        "reproductive_cells": f.reproductive_cell_ntpm.ge(1),
        "localized_reproductive_site": (
            f.reproductive_cell_indirect_ntpm.ge(1)
            & f.outside_cluster_complete
            & f.outside_cluster_max_ntpm.lt(1)
            & ~f.normal_induction_evidence
        ),
        "reproductive_ihc": f.reproductive_protein_association,
        "reviewed_reproductive_protein": f.knowledge_reproductive_support,
        "developmental_rna": f.developmental_mean_tpm.ge(1),
    }
    routes = pd.DataFrame(branches, index=f.index)
    f["reproductive_association_routes"] = routes.apply(lambda r: " | ".join(r.index[r]), axis=1)
    measured = (
        f.reproductive_ntpm.notna() | f.reproductive_cell_ntpm.notna() | f.developmental_measured
    )
    f["reproductive_association_evidence"] = np.select(
        [routes.any(axis=1), measured], ["qualifies", "does_not_qualify"], default="unassessed"
    )
    return f


def group_normal_rows(raw, facts, *, source, statistic, gene_col="gene_id"):
    """Sum member loci within one measured context, never sum their maxima."""
    gene_to_entity = {
        g: r.entity_id for r in facts.itertuples() for g in json.loads(r.member_gene_ids)
    }
    members = facts.set_index("entity_id").member_count
    d = raw[raw[gene_col].isin(gene_to_entity)].copy()
    d["entity_id"] = d[gene_col].map(gene_to_entity)
    keys = ["entity_id", "tissue", "source_context"]
    if d.duplicated([gene_col, "tissue", "source_context"]).any():
        raise ValueError(f"Duplicate normal locus/context in {source}")
    grouped = d.groupby(keys, sort=True).value.agg(["sum", "count"]).reset_index()
    grouped["member_count"] = grouped.entity_id.map(members)
    grouped["value"] = grouped["sum"].where(grouped["count"].eq(grouped.member_count))
    grouped = grouped.drop(columns="sum").rename(columns={"count": "measured_members"})
    return grouped.assign(source_id=source, statistic=statistic, unit="TPM", cell_type="bulk")


def write_hpa_context_inventory():
    """Keep the report rebuildable from the compact HPA v23 tissue and cell inventory."""
    rows = []
    for key in ["hpa_rna_consensus", "hpa_single_cell_tissue"]:
        fields = {"Tissue": "tissue"}
        if key == "hpa_single_cell_tissue":
            fields["Cell type"] = "cell_type"
        pieces = [
            x.drop_duplicates()
            for x in pd.read_csv(
                reference_data.ensure(key, "v23"),
                sep="\t",
                usecols=list(fields),
                chunksize=250000,
            )
        ]
        sites = pd.concat(pieces).drop_duplicates().rename(columns=fields)
        rows.append(
            sites.assign(
                family=key,
                source_id=key + ":v23",
                source_context=sites.tissue,
                role="Shared screen",
            )
        )
    write_csv(pd.concat(rows, ignore_index=True), OUT / "hpa-normal-source-contexts.csv")


def normal_panel(facts, write=True, *, include_individual_samples=True, normal_stat="p95"):
    rows = []
    hpa_ids = read(BASE / "hpa_rna_consensus-accepted-identities.csv.gz")
    lookup = hpa_ids.set_index("source_gene_id").canonical_gene_id
    for name, context in [
        ("rna_tissue_hpa", "Tissue"),
        ("rna_tissue_gtex", "Tissue"),
        ("rna_brain_gtex", "Brain region"),
    ]:
        raw = read_reference(name, "v23", sep="\t")
        raw["gene_id"] = raw.Gene.map(lookup)
        raw["tissue"] = raw[context].replace({"hippocampus": "hippocampal formation"})
        raw["source_context"] = raw[context]
        raw["value"] = raw.TPM
        rows.append(
            group_normal_rows(raw, facts, source=name + ":v23", statistic="tissue mean TPM")
        )
    raw = read_reference("rna_brain_hpa", "v23", sep="\t")
    crosswalk = read(BASE / "alternative-evidence/hpa-brain-anatomy.csv")
    if set(raw.Subregion) != set(crosswalk.Subregion):
        raise ValueError("Changed brain anatomy")
    decisions = read(BASE / "alternative-evidence/hpa_brain-row-decisions.csv.gz")
    brain_lookup = {
        g: r.gene_id
        for r in decisions[decisions.row_position.notna()].itertuples()
        for g in json.loads(r.source_gene_ids)
    }
    raw["gene_id"] = raw.Gene.map(brain_lookup)
    raw["tissue"] = raw.Subregion.map(crosswalk.set_index("Subregion")["Brain region"])
    raw["source_context"] = raw.Subregion
    raw["value"] = raw.TPM
    rows.append(
        group_normal_rows(raw, facts, source="rna_brain_hpa:v23", statistic="subregion mean TPM")
    )

    raw = read_reference("gtex_tissue_rna", "v11", sep="\t", skiprows=2)
    decisions = read(BASE / "alternative-evidence/gtex_tissue_rna-row-decisions.csv.gz")
    chosen = decisions[decisions.row_position.notna()]
    mapped = raw.iloc[chosen.row_position.astype(int)].copy()
    mapped["gene_id"] = chosen.gene_id.to_numpy()
    gtex_columns = {
        re.sub(r"[ ()-]+", "_", key.replace("c-1", "cMINUS1"))
        .strip("_")
        .replace("cMINUS1", "c-1"): value
        for key, value in GTEX_TISSUES.items()
    }
    # The historical crosswalk grouped uterine bulk with endometrium. Keep the
    # organ-level collection label when attributing individual witnesses.
    gtex_columns["Uterus"] = "uterus"
    long = mapped.melt(
        id_vars="gene_id",
        value_vars=list(gtex_columns),
        var_name="source_context",
        value_name="value",
    )
    long["tissue"] = long.source_context.map(gtex_columns)
    rows.append(
        group_normal_rows(
            long, facts, source="gtex_tissue_rna:v11", statistic="tissue donor median TPM"
        )
    )
    normal = pd.concat(rows, ignore_index=True)
    if include_individual_samples:
        from oncoref.curation import reference_input

        archive = reference_input("hpa_tissue_samples")
        genes = tissue_sample_genes(archive)
        carryover_reference = normal_library_qc.carryover_reference(mapped)
        audit = tissue_library_qc(genes, hpa_ids, carryover_reference)
        gtex, gtex_audit = individual_samples(
            facts, mapped, gtex_columns, carryover_reference, normal_stat=normal_stat
        )
        if write:
            write_csv(audit, OUT / "normal-library-qc.csv")
            write_csv(gtex_audit, OUT / "gtex-library-qc.csv.gz")
        individual = tissue_samples(facts, hpa_ids, genes=genes)
        individual = individual.merge(
            audit[["sample_id", "library_qc_review"]], on="sample_id", validate="many_to_one"
        )
        individual["member_count"] = individual.entity_id.map(
            facts.set_index("entity_id").member_count
        )
        individual["measured_members"] = individual.member_count.where(individual.identity_resolved)
        normal = pd.concat([normal, individual, gtex], ignore_index=True)
        INPUTS["reference:hpa_individual_tissue_tpm:v23"] = HPA_SOURCE["sha256"]
        INPUTS.update(
            {
                f"reference:{GTEX_SOURCE}:v11-{native.lower()}": reference_data.REFERENCE_SOURCES[
                    GTEX_SOURCE
                ]["pins"]["v11-" + native.lower()]["sha256"]
                for native in gtex_columns
            }
        )
    normal["tissue"] = normal.tissue.replace(
        {"nucleus accumbens": "basal ganglia", "putamen": "basal ganglia"}
    )
    normal = anatomy(normal.rename(columns={"source_id": "source"}))
    normal["category"] = tissue_tiers.context_tiers(normal)
    normal = normal.rename(columns={"source": "source_id"})
    normal["scope"] = normal.category.where(
        normal.category.isin(["reproductive", "disposable"]), "normal"
    )
    normal.loc[normal.tissue.eq("thymus"), "scope"] = "thymic"
    normal["context_id"] = normal.tissue
    normal["qc_pass"] = normal.value.notna()
    normal["identity_resolved"] = normal.measured_members.eq(normal.member_count)
    # Same combined anatomical coverage rule, now requiring native TPM support.
    consensus = read_reference("hpa_rna_consensus", "v23", sep="\t", usecols=["Tissue"])
    concepts = anatomy(pd.DataFrame({"tissue": sorted(consensus.Tissue.unique())}))
    concepts = concepts.assign(source="hpa_rna_consensus", cell_type="bulk")
    concepts["category"] = tissue_tiers.context_tiers(concepts)
    required = set(
        concepts.loc[~concepts.category.isin(["reproductive", "disposable"]), "tissue"]
    ) - {"thymus"}

    from bulk_normal_percentiles import dice_percentiles, hpa_immune_percentiles, tissue_percentiles

    normal = tissue_percentiles(normal, method=normal_stat)
    immune = pd.concat(
        [
            dice_percentiles(facts, INPUTS, method=normal_stat),
            hpa_immune_percentiles(facts, hpa_ids, INPUTS, method=normal_stat),
        ]
    )
    monaco = read_reference("hpa_immune_monaco", "v23", sep="\t")
    monaco_ids = read(BASE / "hpa_immune_monaco-accepted-identities.csv.gz")
    monaco["gene_id"] = monaco.Gene.map(monaco_ids.set_index("source_gene_id").canonical_gene_id)
    monaco["value"] = monaco.TPM
    monaco["tissue"] = "blood"
    monaco["source_context"] = monaco["Immune cell"]
    monaco = group_normal_rows(
        monaco, facts, source="hpa_immune_monaco:v23", statistic="immune cell mean TPM"
    )
    monaco["cell_type"] = monaco.source_context
    monaco["context_id"] = "blood"
    monaco["scope"] = "normal"
    monaco["qc_pass"] = monaco.value.notna()
    monaco["identity_resolved"] = monaco.measured_members.eq(monaco.member_count)
    normal = pd.concat([normal, immune, monaco], ignore_index=True)
    normal["comparison_id"] = COMPARISON_ID
    for column in (
        "compartment",
        "mammary_state_inference",
        "sample_id",
        "participant_id",
        "observation_kind",
        "library_qc_review",
    ):
        if column not in normal:
            normal[column] = ""
        normal[column] = normal[column].fillna("")
    normal["compartment"] = normal.compartment.where(normal.compartment.ne(""), normal.tissue)
    cols = [
        "entity_id",
        "source_id",
        "comparison_id",
        "tissue",
        "compartment",
        "mammary_state_inference",
        "cell_type",
        "source_context",
        "context_id",
        "scope",
        "statistic",
        "n_observations",
        "unit",
        "value",
        "qc_pass",
        "identity_resolved",
        "sample_id",
        "participant_id",
        "observation_kind",
        "library_qc_review",
    ]
    normal = (
        normal[cols]
        .sort_values(["entity_id", "source_id", "source_context"])
        .reset_index(drop=True)
    )
    if write:
        write_csv(normal, OUT / "normal-comparison.csv.gz")
    return normal, required


def disposable_ntpm(candidates):
    """HPA consensus nTPM in Disposable bulk tissues; a separate scale from TPM."""
    identities = read(BASE / "hpa_rna_consensus-accepted-identities.csv.gz")
    lookup = {
        g: r.entity_id for r in candidates.itertuples() for g in json.loads(r.member_gene_ids)
    }
    raw = read_reference("hpa_rna_consensus", "v23", sep="\t")
    raw = raw[raw.Tissue.isin(DISPOSABLE_NTPM_TISSUES)].copy()
    raw["entity_id"] = raw.Gene.map(identities.set_index("source_gene_id").canonical_gene_id).map(
        lookup
    )
    rows = raw.groupby(["entity_id", "Tissue"]).nTPM.agg(["sum", "count"]).reset_index()
    expected = candidates.set_index("entity_id").member_count
    rows["value"] = rows["sum"].where(rows["count"].eq(rows.entity_id.map(expected)))
    return rows.groupby("entity_id").value.max()


def expression_bounds(selected, evidence, normal):
    """Target maximum and Disposable caps for proteins passing the prior checks."""
    rows = normal[normal.qc_pass & normal.identity_resolved & normal.value.notna()]
    by_scope = rows.groupby(["scope", "entity_id"]).value.max()
    b = selected[["entity_id"]].merge(
        evidence[
            ["entity_id", "normal_reference_tpm", "tumor_max_tpm", "reproductive_max_sample_tpm"]
        ],
        on="entity_id",
    )
    b["disposable_max_tpm"] = b.entity_id.map(by_scope.get("disposable", pd.Series(dtype=float)))
    b["reproductive_max_tpm"] = b.entity_id.map(
        by_scope.get("reproductive", pd.Series(dtype=float))
    )
    b["disposable_max_ntpm"] = b.entity_id.map(disposable_ntpm(selected))
    # Use the same qualified observations as the two enrichment decisions.
    # Historical developmental averages remain annotations, not admission evidence.
    b["target_peak_tpm"] = b[["tumor_max_tpm", "reproductive_max_sample_tpm"]].max(axis=1)
    measured = b[["normal_reference_tpm", "disposable_max_tpm", "disposable_max_ntpm"]].notna()
    ceiling = b[["normal_reference_tpm", "disposable_max_tpm"]].max(axis=1)
    b["target_maximum_pass"] = b.target_peak_tpm.ge(target_selection.DETECTION_FLOOR) & b[
        "target_peak_tpm"
    ].gt(ceiling)
    b["disposable_cap_pass"] = b.disposable_max_tpm.le(
        target_selection.DISPOSABLE_CAP
    ) & b.disposable_max_ntpm.le(target_selection.DISPOSABLE_CAP)
    b["expression_bounds_status"] = np.where(
        ~measured.all(axis=1),
        "unassessed",
        np.where(b.target_maximum_pass & b.disposable_cap_pass, "pass", "fail"),
    )
    return b.drop(columns=["normal_reference_tpm", "tumor_max_tpm", "reproductive_max_sample_tpm"])


def participant_id(namespace, sample):
    """Compatibility wrapper; identity policy is owned by the public sample API."""
    rows = pd.DataFrame([{"source_cohort": namespace, "sample_id": sample}])
    return sample_api.sample_identities(rows).iloc[0].participant_id


def annotate_papers(decisions):
    """Retain the full curated occurrence ledger, joined by exact member loci."""
    links = cta_literature.current_attestations(catalog=decisions)
    for name in [
        "cta-literature.json.gz",
        "cta-literature-followup.json",
        "cta-selection-audit.json",
        "cta-publication-membership.csv",
        "cta-publication-sources.csv",
    ]:
        path = ROOT / "oncoref/data" / name
        INPUTS[str(path.relative_to(ROOT))] = sha256(path)
    rows = links.to_dict("records")
    with gzip.GzipFile(
        filename=str(OUT / "paper-attestations.json.gz"), mode="wb", mtime=0
    ) as handle:
        handle.write(json.dumps(rows, separators=(",", ":")).encode())
    by_gene, nominations = defaultdict(set), defaultdict(set)
    for i, row in enumerate(rows):
        for gene in row["gene_ids"]:
            by_gene[gene].add(i)
        for gene in row["nomination_gene_ids"]:
            nominations[gene].add(row["source_id"])
    result = decisions.copy()
    annotations = []
    for row in result.itertuples():
        genes = json.loads(row.member_gene_ids)
        indices = set().union(*(by_gene[g] for g in genes))
        source_ids = sorted({rows[i]["source_id"] for i in indices})
        nomination_sources = sorted(set().union(*(nominations[g] for g in genes)))
        covered_members = [g for g in genes if nominations[g]]
        annotations.append(
            {
                "paper_source_ids": json.dumps(source_ids),
                "paper_attestation_indices": json.dumps(sorted(indices)),
                "paper_source_count": len(source_ids),
                "published_nomination": bool(nomination_sources),
                "nomination_source_ids": json.dumps(nomination_sources),
                "nomination_member_ids": json.dumps(covered_members),
                "nomination_scope": "all_members"
                if len(covered_members) == len(genes)
                else "some_members"
                if covered_members
                else "none",
            }
        )
    return pd.concat([result.reset_index(drop=True), pd.DataFrame(annotations)], axis=1)


def pinned_tumor_matrix(path, digest):
    from oncoref import catalog

    entry = Path(path)
    resolved = catalog.ensure(
        "per-sample:" + entry.stem, version=entry.parent.name.removeprefix("v")
    )
    if sha256(resolved) != digest:
        raise ValueError(f"Changed tumor matrix: {resolved}")
    return resolved


def tumors(facts):
    """Read all pinned matrices; retain per-protein, per-physical-sample values."""
    folder = BASE / "tumor-identity-repair"
    manifest_path = folder / "manifest.json"
    INPUTS[str(manifest_path.relative_to(ROOT))] = sha256(manifest_path)
    manifest = json.loads(manifest_path.read_text())
    paths = {Path(p).stem: (Path(p), digest) for p, digest in manifest["matrices"].items()}
    available = read(folder / "cohort-availability.csv")
    samples = read(folder / "sample-identity.csv.gz")
    gene_to_entity = {
        g: r.entity_id for r in facts.itertuples() for g in json.loads(r.member_gene_ids)
    }
    expected = facts.set_index("entity_id").member_count
    pieces, coverage, qc_rows = [], [], []
    identity_cache = {}
    for job in available.itertuples():
        path, digest = paths[job.cancer_code]
        path = pinned_tumor_matrix(path, digest)
        if job.cancer_code in target_selection.GERM_CELL_TUMORS:
            coverage.append(
                {
                    "cancer_code": job.cancer_code,
                    "source_cohort": job.source_cohort,
                    "status": "germ cell tumor: not cancer evidence for germline genes",
                    "matrix_sha256": digest,
                }
            )
            continue
        if job.status != "Used" or job.source_scale_class != "linear_rnaseq_tpm":
            coverage.append(
                {
                    "cancer_code": job.cancer_code,
                    "status": "not comparable to native transcript TPM",
                    "matrix_sha256": digest,
                }
            )
            continue
        id_rows = pq.read_table(path, columns=["Ensembl_Gene_ID", "Symbol"]).to_pandas()
        key = hashlib.sha256(
            pd.util.hash_pandas_object(id_rows, index=False).to_numpy().tobytes()
        ).hexdigest()
        if key not in identity_cache:
            identity_cache[key] = match_reference_genes(
                id_rows,
                id_column="Ensembl_Gene_ID",
                symbol_column="Symbol",
                source=job.source_cohort,
                source_version=str(job.source_matrix_version),
            )
        identity = identity_cache[key]
        identity = identity[
            identity.identity_status.eq("matched") & identity.canonical_gene_id.isin(gene_to_entity)
        ]
        if identity.empty:
            coverage.append(
                {
                    "cancer_code": job.cancer_code,
                    "source_cohort": job.source_cohort,
                    "status": "no resolved eligible protein rows",
                    "matrix_sha256": digest,
                }
            )
            continue
        raw = pq.read_table(path).to_pandas()
        sample_ids = (
            samples.loc[samples.cancer_code.eq(job.cancer_code), "sample_id"].astype(str).tolist()
        )
        signatures = library_qc.signatures(raw, sample_ids)
        signatures["reasons"] = signatures.apply(lambda r: ";".join(library_qc.reasons(r)), axis=1)
        signatures["notes"] = signatures.apply(lambda r: ";".join(library_qc.notes(r)), axis=1)
        signatures["cancer_code"] = job.cancer_code
        signatures["source_cohort"] = job.source_cohort
        qc_rows.append(signatures.reset_index())
        sample_ids = signatures.index[signatures.reasons.eq("")].tolist()
        raw = raw[raw.Ensembl_Gene_ID.isin(identity.source_gene_id)]
        raw, decisions = harmonize_reference_expression(
            raw,
            sample_cols=sample_ids,
            source=job.source_cohort,
            source_version=str(job.source_matrix_version),
            symbol_column="Symbol",
            duplicate_policy="quarantine",
        )
        write_csv(decisions, OUT / "tumor-identities" / f"{job.cancer_code}.csv.gz")
        raw["gene_id"] = raw.Ensembl_Gene_ID
        raw["entity_id"] = raw.gene_id.map(gene_to_entity)
        values = raw.groupby("entity_id")[sample_ids].sum(min_count=1)
        counts = raw.groupby("entity_id")[sample_ids].count()
        values = values.where(counts.eq(expected.reindex(values.index), axis=0))
        if np.isinf(values.to_numpy()).any() or values.lt(0).any().any():
            raise ValueError("Invalid tumor expression in pinned QC-passing samples")
        namespace = source_matrices.source_sample_namespace(job.source_cohort)
        # A compact wide block preserves zeros and measurements for every candidate.
        values.columns = [namespace + ":" + s for s in sample_ids]
        pieces.append(values)
        coverage.append(
            {
                "cancer_code": job.cancer_code,
                "source_cohort": job.source_cohort,
                "namespace": namespace,
                "status": "assessed",
                "matrix_sha256": digest,
                "samples": len(sample_ids),
                "measured_proteins": int(values.notna().any(axis=1).sum()),
            }
        )
        print(
            f"Tumor {job.cancer_code}: {len(sample_ids)} samples, {len(values)} proteins",
            flush=True,
        )
    write_csv(pd.concat(qc_rows, ignore_index=True), OUT / "library-qc.csv.gz")
    INPUTS["oncoref/library_qc.py"] = sha256(ROOT / "oncoref/library_qc.py")
    INPUTS["oncoref/data/tissue-marker-panels.csv"] = sha256(
        ROOT / "oncoref/data/tissue-marker-panels.csv"
    )
    combined = pd.concat(pieces, axis=1).T
    if combined.index.duplicated().any():
        repeated = combined[combined.index.duplicated(keep=False)]
        if repeated.groupby(level=0).nunique(dropna=False).gt(1).any().any():
            raise ValueError("Overlapping cohort views disagree")
        combined = combined[~combined.index.duplicated()]
    sample_meta = samples.drop_duplicates(["namespace", "sample_id"]).copy()
    sample_meta = sample_api.sample_identities(sample_meta)
    for name in [
        "cancer-reference-expression-samples.csv.gz",
        "cancer-reference-sample-molecular-provenance.csv",
    ]:
        path = ROOT / "oncoref/data" / name
        INPUTS[str(path.relative_to(ROOT))] = sha256(path)
    INPUTS["oncoref/samples.py"] = sha256(ROOT / "oncoref/samples.py")
    sample_meta = sample_meta[sample_meta.sample_key.isin(combined.index)]
    combined.index.name = "sample_key"
    CACHE.mkdir(parents=True, exist_ok=True)
    combined.to_parquet(CACHE / "tumor-comparison.parquet", compression="zstd")
    write_csv(sample_meta, OUT / "tumor-samples.csv.gz")
    historical = read(
        ROOT / "docs/audits/cta-expanded-funnel-20261027/tumor-union/cohort-availability.csv"
    )
    coverage += historical[~historical.cancer_code.isin(available.cancer_code)].to_dict("records")
    write_csv(pd.DataFrame(coverage), OUT / "tumor-source-coverage.csv")
    cache_receipt = {
        "matrix_sha256": sha256(CACHE / "tumor-comparison.parquet"),
        "samples_sha256": sha256(OUT / "tumor-samples.csv.gz"),
        "identity_code_sha256": sha256(ROOT / "oncoref/gene_matching.py"),
        "members": dict(sorted(gene_to_entity.items())),
        "inputs": dict(INPUTS),
    }
    (CACHE / "receipt.json").write_text(json.dumps(cache_receipt, indent=2) + "\n")
    return combined, sample_meta


def cached_tumors(facts):
    """Reuse a verified comparison matrix after changes outside tumor extraction."""
    receipt = json.loads((CACHE / "receipt.json").read_text())
    members = {g: r.entity_id for r in facts.itertuples() for g in json.loads(r.member_gene_ids)}
    if members != receipt["members"]:
        raise ValueError("Cached tumor membership differs; rebuild the matrices")
    for name in ("oncoref/library_qc.py", "oncoref/data/tissue-marker-panels.csv"):
        if name not in receipt["inputs"]:
            raise ValueError("Cached tumor evidence predates library QC")
    # Non-tumor curation/normal inputs are freshly read by shared_facts and
    # normal_panel. Only tumor extraction inputs belong to this cache identity.
    tumor_inputs = {
        name: digest
        for name, digest in receipt["inputs"].items()
        if "tumor-identity-repair/" in name
        or name
        in {
            "oncoref/library_qc.py",
            "oncoref/data/tissue-marker-panels.csv",
            "oncoref/samples.py",
            "oncoref/data/cancer-reference-expression-samples.csv.gz",
            "oncoref/data/cancer-reference-sample-molecular-provenance.csv",
        }
    }
    for name, digest in tumor_inputs.items():
        if sha256(ROOT / name) != digest:
            raise ValueError(f"Changed cached comparison input: {name}")
    manifest = json.loads((BASE / "tumor-identity-repair/manifest.json").read_text())
    for name, digest in manifest["matrices"].items():
        pinned_tumor_matrix(name, digest)
    coverage = pd.read_csv(OUT / "tumor-source-coverage.csv")
    if coverage.cancer_code.isin(target_selection.GERM_CELL_TUMORS).any():
        raise ValueError("Excluded histologies require a fresh tumor extraction")
    if sha256(CACHE / "tumor-comparison.parquet") != receipt["matrix_sha256"]:
        raise ValueError("Cached comparison matrix checksum mismatch")
    if sha256(OUT / "tumor-samples.csv.gz") != receipt["samples_sha256"]:
        raise ValueError("Cached sample identity checksum mismatch")
    if sha256(ROOT / "oncoref/gene_matching.py") != receipt["identity_code_sha256"]:
        raise ValueError("Cached identity mapping changed")
    INPUTS.update(tumor_inputs)
    return pd.read_parquet(CACHE / "tumor-comparison.parquet"), pd.read_csv(
        OUT / "tumor-samples.csv.gz", keep_default_na=False
    )


def build(*, reuse_tumor_cache=False, normal_stat="p95"):
    OUT.mkdir(parents=True, exist_ok=True)
    facts = shared_facts()
    selected = facts[facts.comparison_eligible].copy()
    print(f"Common eligibility: {len(selected)} / {len(facts)}", flush=True)
    normal, required = normal_panel(selected, normal_stat=normal_stat)
    write_hpa_context_inventory()
    development, development_provenance = developmental_samples(selected)
    development["comparison_id"] = COMPARISON_ID
    reproductive = pd.concat(
        [normal[normal.scope.eq("reproductive") & normal.sample_id.ne("")], development],
        ignore_index=True,
    )
    if not set(reproductive.compartment) <= {*target_selection.REPRODUCTIVE_COMPARTMENTS, "breast"}:
        raise ValueError(
            f"Unexpected reproductive compartments: {sorted(set(reproductive.compartment))}"
        )
    reproductive["qc_pass"] &= reproductive.library_qc_review.fillna("").eq("")
    write_csv(reproductive, OUT / "reproductive-comparison.csv.gz")
    for source in development_provenance:
        INPUTS[f"reference:{source['source']}:{source['version']}"] = source["sha256"]
    matrix, meta = cached_tumors(selected) if reuse_tumor_cache else tumors(selected)
    meta = meta.set_index("sample_key").loc[matrix.index]
    summaries, witnesses, reproductive_witnesses = [], [], []
    grouped = dict(iter(normal.groupby("entity_id")))
    reproductive_groups = dict(iter(reproductive.groupby("entity_id")))
    for row in selected.itertuples():
        values = matrix.get(row.entity_id, pd.Series(np.nan, index=matrix.index))
        tumor = pd.DataFrame(
            {
                "entity_id": row.entity_id,
                "source_id": meta.namespace,
                "sample_id": meta.sample_id,
                "participant_id": meta.participant_id,
                "value": values,
                "comparison_id": COMPARISON_ID,
                "unit": "TPM",
                "qc_pass": values.notna(),
                "identity_resolved": values.notna(),
            }
        ).reset_index(drop=True)
        reference = grouped.get(row.entity_id, normal.iloc[:0])
        receipt = target_selection.assess_cancer_rna(
            tumor,
            reference,
            required_normal_contexts=required,
            comparison_id=COMPARISON_ID,
            comparability_reference="comparison.json",
            unit="TPM",
        )
        support = pd.DataFrame(receipt["supporting_tumors"])
        qualified = support[support.single_sample_pass] if not support.empty else support
        top = (
            support.sort_values("value", ascending=False).drop_duplicates("participant_id")
            if not support.empty
            else support
        )
        # Keep high anonymous samples as separate witnesses, not independent patients.
        anonymous = support[support.participant_id.eq("")] if not support.empty else support
        top = (
            pd.concat([top.head(3), anonymous.head(3)]).drop_duplicates(["source_id", "sample_id"])
            if not support.empty
            else top
        )
        # Preserve every tenfold witness; weaker examples are a bounded preview.
        top = (
            pd.concat([qualified, top]).drop_duplicates(["source_id", "sample_id"])
            if not support.empty
            else top
        )
        witnesses.extend(top.to_dict("records"))
        ra = reproductive_groups.get(row.entity_id, reproductive.iloc[:0])
        ra_receipt = target_selection.assess_reproductive_rna(
            ra,
            reference,
            required_normal_contexts=required,
            comparison_id=COMPARISON_ID,
            comparability_reference="comparison.json",
            unit="TPM",
        )
        ra_valid = ra[ra.qc_pass & ra.identity_resolved & ra.value.notna()]
        ra_top = ra_valid.sort_values("value", ascending=False).head(1)
        # Tissue loss is a separate axis from tenfold membership. Retain measured
        # reproductive expression even when it does not meet the RA threshold.
        observed_contexts = []
        observed = ra[ra.identity_resolved & ra.value.ge(target_selection.DETECTION_FLOOR)]
        for (tissue, compartment, cell, scope), contexts in observed.groupby(
            ["tissue", "compartment", "cell_type", "scope"], sort=True
        ):
            observed_contexts.append(
                {
                    "tissue": tissue,
                    "compartment": compartment,
                    "cell_type": cell,
                    "scope": scope,
                    "stage": "embryonic" if scope == "developmental" else "adult",
                    "native_tpm_max": float(contexts.value.max()),
                    "qc_supported": bool(contexts.qc_pass.any()),
                    "sources": sorted(contexts.source_id.unique()),
                }
            )
        comp_max = {
            comp: float(val)
            for comp, val in ra_valid.groupby("compartment", sort=True).value.max().items()
            if val > 0
        }
        ra_support = pd.DataFrame(ra_receipt["supporting_reproductive_samples"])
        ra_supported_top = ra_receipt["supported_maximum_observations"]
        ra_qualifying = (
            ra_support[ra_support.single_sample_pass] if not ra_support.empty else ra_support
        )
        qual_tpm = (
            {
                comp: float(val)
                for comp, val in ra_qualifying.groupby("compartment", sort=True)
                .support_ceiling.max()
                .items()
                if pd.notna(val) and val > 0
            }
            if not ra_qualifying.empty
            else {}
        )
        if not ra_support.empty:
            reproductive_witnesses.extend(ra_qualifying.to_dict("records"))
        cancer_keys = (
            set(qualified.source_id + ":" + qualified.sample_id) if not qualified.empty else set()
        )
        summaries.append(
            {
                "entity_id": row.entity_id,
                "cancer_evidence": receipt["status"],
                "cancer_routes": " | ".join(receipt["routes"]),
                "normal_reference_tpm": receipt["normal_maximum"],
                "normal_reference_contexts": json.dumps(
                    pd.DataFrame(receipt["normal_maximum_observations"][:10])
                    .astype(object)
                    .where(lambda d: d.notna(), None)
                    .to_dict("records"),
                    allow_nan=False,
                ),
                "normal_reference_tied_contexts": len(receipt["normal_maximum_observations"]),
                "missing_comparison_contexts": json.dumps(receipt["missing_normal_contexts"]),
                "independent_patients_above_normal": receipt[
                    "independent_participants_above_normal"
                ],
                "samples_above_normal": receipt["samples_above_normal"],
                "samples_at_tenfold": receipt["samples_at_tenfold"],
                "independent_patients_at_tenfold": receipt["independent_participants_at_tenfold"],
                "independent_patient_tumor_folds": json.dumps(
                    receipt["independent_participant_folds"]
                ),
                "anonymous_samples_at_tenfold": receipt["anonymous_samples_at_tenfold"],
                "patient_or_sample_witnesses_at_tenfold": receipt[
                    "participant_or_sample_witnesses_at_tenfold"
                ],
                "cancer_types_at_tenfold": int(
                    meta.loc[meta.index.isin(cancer_keys), "cancer_code"].nunique()
                ),
                "reproductive_evidence": ra_receipt["status"],
                "reproductive_routes": " | ".join(ra_receipt["routes"]),
                "reproductive_assessment_reason": ra_receipt["reason"],
                "reproductive_max_sample_tpm": None
                if ra_top.empty
                else float(ra_top.value.iloc[0]),
                "reproductive_fold": None
                if ra_receipt["supported_maximum"] is None
                or ra_receipt["single_sample_threshold"] is None
                else ra_receipt["supported_maximum"]
                / (ra_receipt["single_sample_threshold"] / target_selection.REPRODUCTIVE_FOLD),
                "reproductive_supported_tpm": ra_receipt["supported_maximum"],
                "reproductive_supported_context": ""
                if not ra_supported_top
                else ra_supported_top[0]["source_context"],
                "reproductive_supported_observations": json.dumps(
                    pd.DataFrame(ra_supported_top)
                    .astype(object)
                    .where(lambda d: d.notna(), None)
                    .to_dict("records"),
                    allow_nan=False,
                ),
                "reproductive_expression_contexts": json.dumps(observed_contexts),
                "reproductive_compartment_max_tpm": json.dumps(comp_max, allow_nan=False),
                "reproductive_qualifying_tpm": json.dumps(qual_tpm, allow_nan=False),
                "reproductive_qualifying_compartments": " | ".join(qual_tpm),
                "reproductive_top_sample": ""
                if ra_top.empty
                else ra_top.source_id.iloc[0] + ":" + ra_top.sample_id.iloc[0],
                "reproductive_top_context": "" if ra_top.empty else ra_top.source_context.iloc[0],
                "reproductive_samples_at_tenfold": ra_receipt["samples_at_tenfold"],
                "reproductive_donors_at_tenfold": ra_receipt["independent_participants_at_tenfold"],
                "reproductive_anonymous_samples_at_tenfold": ra_receipt[
                    "anonymous_samples_at_tenfold"
                ],
                "unresolved_patient_samples_above_normal": receipt[
                    "unresolved_participant_samples_above_normal"
                ],
                "tumor_samples_measured": int(values.notna().sum()),
                "tumor_max_tpm": values.max() if values.notna().any() else None,
                "tumor_top_sample": values.idxmax() if values.notna().any() else "",
                "tumor_top_cancer_code": meta.loc[values.idxmax(), "cancer_code"]
                if values.notna().any()
                else "",
                "single_tumor_threshold_tpm": receipt["single_sample_threshold"],
                # Highest tumor over max(N, Disposable, 1): the cancer rule needs >= 10.
                "tumor_fold": values.max()
                / (receipt["single_sample_threshold"] / target_selection.TUMOR_FOLD)
                if values.notna().any() and receipt["single_sample_threshold"]
                else None,
                "cancer_assessment_reason": receipt["reason"],
            }
        )
    evidence = pd.DataFrame(summaries)
    evidence = evidence.merge(
        expression_bounds(selected, evidence, normal), on="entity_id", validate="one_to_one"
    )
    facts = facts.merge(evidence, on="entity_id", how="left", validate="one_to_one")
    facts["expression_bounds_status"] = facts.expression_bounds_status.fillna("unassessed")
    # Unassessed bounds never pass; an eligible protein without them needs review.
    facts["expression_bounds_pass"] = facts.expression_bounds_status.eq("pass")
    facts["evidence_resolved"] &= ~(
        facts.shared_eligible & facts.expression_bounds_status.eq("unassessed")
    )
    facts["shared_eligible"] &= facts.expression_bounds_pass
    facts["cancer_evidence"] = facts.cancer_evidence.fillna("unassessed")
    facts["reproductive_evidence"] = facts.reproductive_evidence.fillna("unassessed")
    facts["reproductive_routes"] = facts.reproductive_routes.fillna("")
    facts["cancer_assessment_reason"] = facts.cancer_assessment_reason.fillna(
        "Shared eligibility fails; positive cancer membership cannot be assigned"
    )
    decisions = annotate_papers(target_selection.classify(facts))

    def _cta_tier(r) -> str:
        if not r.is_cta:
            return ""
        qual = (
            json.loads(r.reproductive_qualifying_tpm)
            if isinstance(r.reproductive_qualifying_tpm, str) and r.reproductive_qualifying_tpm
            else {}
        )
        cmax = (
            json.loads(r.reproductive_compartment_max_tpm)
            if isinstance(r.reproductive_compartment_max_tpm, str)
            and r.reproductive_compartment_max_tpm
            else {}
        )
        if not (set(qual) & set(target_selection.CORE_REPRODUCTIVE_COMPARTMENTS)):
            return "extended_only"
        h0 = max(float(r.normal_reference_tpm), float(r.disposable_max_tpm), 1.0)
        ext_max = max(
            (float(cmax.get(c, 0.0)) for c in target_selection.EXTENDED_REPRODUCTIVE_COMPARTMENTS),
            default=0.0,
        )
        if ext_max <= target_selection.DISPOSABLE_CAP and float(
            r.tumor_max_tpm
        ) >= target_selection.TUMOR_FOLD * max(h0, ext_max):
            return "core"
        return "core_with_extended_exposure"

    decisions["cta_tier"] = [_cta_tier(r) for r in decisions.itertuples()]
    reasons = {
        "normal_screen_pass": "Normal-expression evidence blocks admission",
        "required_coverage_complete": "Required normal coverage is incomplete",
        "identity_translation_pass": "Identity or translation requires review",
        "evidence_resolved": "Evidence is unresolved",
        "sequence_pass": "Nested, read-through or Unsafe-contained sequence",
    }
    bounds_failed = decisions.expression_bounds_status.eq("fail")

    decisions["decision_reason"] = [
        str(r.Unsafe_reasons)
        if r.category == "Unsafe"
        else " | ".join(
            [text for field, text in reasons.items() if not getattr(r, field)]
            + (
                ["Pooled rare-cell flag could not be assessed per donor"]
                if r.rare_cell_legacy_only
                else []
            )
            + (["Target maximum or Disposable cap fails"] if failed else [])
        )
        if not r.shared_eligible
        else f"Normal home is {r.home_tissue}, without sufficient locus-specific normal coverage for cancer admission"
        if r.subtype == "Unscreened tissue"
        else "Cancer and reproductive/developmental evidence qualify"
        if r.category == "CTA"
        else "Cancer tenfold evidence qualifies; reproductive tenfold evidence does not qualify"
        if r.category == "Cancer Antigen only"
        else f"Reproductive enrichment qualifies. Cancer admission lacks coverage of {r.home_tissue}."
        if r.category == "Reproductive Antigen only"
        and r.cancer_evidence == "qualifies"
        and r.home_tissue
        else "Reproductive/developmental evidence qualifies; " + r.cancer_assessment_reason
        if r.category == "Reproductive Antigen only"
        else "Neither tenfold compartment qualifies; "
        + r.cancer_assessment_reason
        + "; "
        + str(r.reproductive_assessment_reason)
        for r, failed in zip(decisions.itertuples(), bounds_failed)
    ]
    decisions["exclusion_reasons"] = decisions.decision_reason.where(
        decisions.category.isin(["Other", "Unsafe"]), ""
    )
    # Legacy logical flags must not masquerade as new-policy decisions.
    old_logic = [
        "policy",
        "decision_flow_id",
        "assessment_status",
        "positive_context",
        "reproductive_route",
        "developmental_route",
        "cancer_association",
        "reproductive_association",
        "developmental_association",
        "reproductive_rna_association",
        "reproductive_cell_association",
        "reproductive_cell_indirect_association",
        "reproductive_concentration_hold",
        "tumor_restricted_profile",
        "tumor_restricted_identity",
        "tumor_restricted_lead",
    ]
    decisions = decisions.rename(columns={c: "legacy_" + c for c in old_logic if c in decisions})
    decisions = decisions.rename(
        columns={
            c: "legacy_" + c
            for c in decisions
            if c.startswith(("hold:", "unsafe:"))
            or c in {"rare_cell_flag", "rare_cell_hold", "rare_cell_max_pct", "rare_cell_type"}
            or c in {"rare_cell_detail", "rare_cell_screened"}
        }
    )
    decisions["policy"] = target_selection.POLICY_ID
    decisions["assessment_status"] = np.where(
        decisions.classification.eq("Other:Review"), "review_required", "screened"
    )
    write_csv(decisions, OUT / "protein-decisions.csv.gz")
    write_csv(pd.DataFrame(witnesses), OUT / "tumor-witnesses.csv.gz")
    write_csv(pd.DataFrame(reproductive_witnesses), OUT / "reproductive-witnesses.csv.gz")
    write_csv(evidence, OUT / "cancer-evidence.csv.gz")
    summary = {
        "policy_id": target_selection.POLICY_ID,
        "proteins": len(decisions),
        "categories": decisions.category.value_counts().to_dict(),
        "classifications": decisions.classification.value_counts().to_dict(),
        "memberships": {
            c: int(decisions[c].sum())
            for c in [
                "is_cancer_antigen",
                "is_reproductive_antigen",
                "is_cta",
                "is_target_candidate",
            ]
        },
        "cta_tiers": decisions.loc[decisions.is_cta, "cta_tier"].value_counts().to_dict(),
        "changes": decisions.groupby(["released_category", "category"])
        .size()
        .reset_index(name="proteins")
        .to_dict("records"),
        "shared_eligible": int(decisions.eligible.sum()),
        "assessed_for_cancer_evidence": len(selected),
        "unsafe_changes": {
            "released_unsafe_now_other": int(
                (decisions.released_category.eq("Unsafe") & decisions.category.ne("Unsafe")).sum()
            ),
            "newly_unsafe": int(
                (decisions.released_category.ne("Unsafe") & decisions.category.eq("Unsafe")).sum()
            ),
        },
        "rare_cell": {
            "status": decisions.rare_cell_status.value_counts().to_dict(),
            "action": decisions.rare_cell_action.value_counts().to_dict(),
            "legacy_flag_unassessed": int(decisions.rare_cell_legacy_only.sum()),
        },
        "literature": decisions.literature_effect.value_counts().to_dict(),
        "expression_bounds": decisions.expression_bounds_status.value_counts().to_dict(),
        "comparison_required_contexts": sorted(required),
        "assessed_cancer_statuses": evidence.cancer_evidence.value_counts().to_dict(),
        "published_nomination_by_category": decisions.groupby("category")
        .published_nomination.sum()
        .astype(int)
        .to_dict(),
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    (OUT / "inputs.json").write_text(json.dumps(INPUTS, indent=2) + "\n")
    comparison = {
        "comparison_id": COMPARISON_ID,
        "unit": "TPM",
        "purpose": "Source-native gene-abundance screening; not calibrated cross-study differential expression",
        "compatibility_basis": [
            "Linear transcript-per-million gene abundance; reject microarray proxies, CP10K, CPM proxies, nTPM, pTPM and gene-feature-count TPM",
            "Current exact locus mapping; ambiguous duplicate vectors quarantined; protein-identical loci summed only in the same sample/context with all members measured",
            "Pinned tumor matrices and per-sample quality selection from the existing source-scale registry",
            f"Complete {len(required)}-context non-thymic, non-reproductive panel, plus native-TPM peripheral immune and brain measurements. Adult whole thymus is Disposable and retained in the separate thymic scope.",
        ],
        "expression_bounds": target_selection.definition()["expression_bounds"],
        "disposable_ntpm_tissues": DISPOSABLE_NTPM_TISSUES,
        "normal_statistic": "Maximum across GTEx v11 within-site linear donor P95, HPA v23 normal biopsy P95 and supported context statistics: tissue/subregion means, GTEx v11 donor medians, and immune summaries including DICE within-population P95 for single loci and a conservative bound for unpaired loci",
        "batch_corrected": False,
        "limitations": [
            "Quantifier, library preparation and study effects are retained; the fold rule is a screening convention",
            "A maximum of tissue summaries is not a maximum over every healthy donor",
            "Other assay scales continue to enforce the existing shared normal screen independently",
        ],
        "participant_identity": "Shared sample_identities API joins curated case/donor mappings and documented TCGA, TARGET, MMRF and Treehouse identifiers. Every sample retains identity status, basis and references. Unresolved or conflicting identities support only the single-sample route.",
        "participant_identity_policy": sample_api.PATIENT_IDENTITY_POLICY,
        "normal_sources": [
            {key: pin[key] for key in ("url", "sha256", "archive_sha256", "bytes")}
            | {"source": name, "version": "v23"}
            for name in ("rna_tissue_hpa", "rna_tissue_gtex", "rna_brain_gtex", "rna_brain_hpa")
            for pin in [reference_data.REFERENCE_SOURCES[name]["pins"]["v23"]]
        ],
        "individual_tissue_source": HPA_SOURCE,
        "gtex_individual_tissue_source": {
            "source": GTEX_SOURCE,
            "release": "v11",
            "unit": "TPM",
            "shards": sorted(
                v
                for v in reference_data.REFERENCE_SOURCES[GTEX_SOURCE]["pins"]
                if v.startswith("v11-")
            ),
            "registry": "oncoref/data/reference-releases.json",
            "sample_inventory": "gtex-library-qc.csv.gz",
            "scope": "All 52 primary tissue matrices; cultured cell lines excluded.",
            "reduction": "Non-reproductive sites use gtex_donor_p95:v11, the linearly interpolated 95th percentile across GTEX donors after complete within-sample protein-group sums and within-donor replicate medians. Reproductive sites retain individual sample witnesses.",
        },
        "mammary_lactation": {
            "state_is_inferred": True,
            "unknown_marker_coverage_is_resting": False,
            "resting_label_meaning": "no high secretory signature in sampled biopsies; not confirmed physiology",
            "markers": list(LACTATION_MARKERS),
            "marker_tpm_threshold": LACTATION_MARKER_TPM,
            "min_markers": LACTATION_MIN_MARKERS,
            "lactating_samples": sorted(
                reproductive.loc[
                    reproductive.compartment.eq("lactating_breast"), "sample_id"
                ].unique()
            ),
            "resting_samples": sorted(
                reproductive.loc[
                    reproductive.compartment.eq("resting_breast"), "sample_id"
                ].unique()
            ),
        },
        "normal_library_qc": {
            "path": "normal-library-qc.csv",
            "sha256": sha256(OUT / "normal-library-qc.csv"),
            "policy": "Signatures are evaluated within each tissue. Reproductive programs are advisory notes and never change sample eligibility. Other flags require adjudication; they never silently remove normal measurements from the denominator. Samples with blocking signatures cannot provide positive reproductive evidence.",
        },
        "developmental_sources": development_provenance,
        "positive_observations": "CA: individual passing tumor; RA: individual reproductive tissue or supported embryo-by-lineage mean. Qualitative protein and incompatible RNA scales cannot bypass the fold.",
        "witness_tables": "All tenfold witnesses are exported; tumor table also retains up to six weaker above-normal examples per protein. Filter single_sample_pass for qualifying counts. Unknown donor identity is never claimed independent.",
        "method_references": [
            "https://v23.proteinatlas.org/about/assays%2Bannotation",
            "https://pmc.ncbi.nlm.nih.gov/articles/PMC5903355/",
        ],
        "required_contexts": sorted(required),
        "tumor_matrix_cache_sha256": sha256(CACHE / "tumor-comparison.parquet"),
    }
    (OUT / "comparison.json").write_text(json.dumps(comparison, indent=2) + "\n")
    compact_columns = [
        "entity_id",
        "protein_group_label",
        "member_gene_ids",
        "representative_gene_id",
        "representative_transcripts",
        "representative_sequence_sha256s",
        "category",
        "subtype",
        "classification",
        "assessment_status",
        "decision_reason",
        "released_category",
        "released_classification",
        "policy",
        *target_selection.SHARED_CHECKS,
        "is_cancer_antigen",
        "is_reproductive_antigen",
        "is_cta",
        "cta_tier",
        "is_target_candidate",
        "cancer_evidence",
        "cancer_routes",
        "cancer_assessment_reason",
        "reproductive_evidence",
        "reproductive_routes",
        "reproductive_association_routes",
        "reproductive_association_evidence",
        "reproductive_assessment_reason",
        "reproductive_max_sample_tpm",
        "reproductive_fold",
        "reproductive_top_sample",
        "reproductive_top_context",
        "reproductive_samples_at_tenfold",
        "reproductive_donors_at_tenfold",
        "reproductive_anonymous_samples_at_tenfold",
        "reproductive_ntpm",
        "reproductive_cell_ntpm",
        "developmental_mean_tpm",
        "rare_cell_status",
        "rare_cell_action",
        "rare_cell_tier",
        "rare_cell_gene_id",
        "rare_cell_basis",
        "rare_cell_mapping_notes",
        "finding_tissue",
        "finding_top_cell_type",
        "finding_studies_multi_molecule",
        "finding_donors_multi_molecule",
        "finding_studies_detected",
        "finding_top_multi_molecule_pct",
        "finding_top_detection_pct",
        "finding_bulk_tpm",
        "finding_bulk_context",
        "paralog_ambiguous",
        "host_overlapped",
        "sequence_partners",
        "same_strand_hosts",
        "max_detection_pct",
        "max_multi_molecule_pct",
        "rare_cell_legacy_only",
        "literature_effect",
        "literature_verdict",
        "literature_rule",
        "literature_decision_tissue",
        "literature_decision_reason",
        "literature_recorded_tissue",
        "literature_recorded_tier",
        "literature_family",
        "literature_family_tier",
        "home_tissue",
        "published_home_tissue",
        "reference_recorded_tissue",
        "reference_recorded_tier",
        "reference_recorded_basis",
        "readthrough_locus",
        "readthrough_parents",
        "expression_bounds_status",
        "target_peak_tpm",
        "disposable_max_tpm",
        "disposable_max_ntpm",
        "normal_reference_tpm",
        "normal_reference_contexts",
        "normal_reference_tied_contexts",
        "single_tumor_threshold_tpm",
        "independent_patients_above_normal",
        "unresolved_patient_samples_above_normal",
        "samples_above_normal",
        "samples_at_tenfold",
        "independent_patients_at_tenfold",
        "independent_patient_tumor_folds",
        "reproductive_supported_tpm",
        "reproductive_supported_context",
        "reproductive_supported_observations",
        "reproductive_expression_contexts",
        "reproductive_compartment_max_tpm",
        "reproductive_qualifying_tpm",
        "reproductive_qualifying_compartments",
        "anonymous_samples_at_tenfold",
        "patient_or_sample_witnesses_at_tenfold",
        "cancer_types_at_tenfold",
        "tumor_samples_measured",
        "tumor_max_tpm",
        "tumor_fold",
        "tumor_top_sample",
        "tumor_top_cancer_code",
        "missing_comparison_contexts",
        "published_nomination",
        "nomination_scope",
        "nomination_source_ids",
        "nomination_member_ids",
        "paper_source_count",
        "paper_source_ids",
        "paper_attestation_indices",
        "knowledge_finding_ids",
        "unsafe_sequence_overlap_pct",
        "hpa_diffuseness",
        "healthy_maximum_ntpm",
        "healthy_maximum_tpm",
        "oncoref_maximum_tpm",
        "Unsafe_reasons",
        "legacy_exclusion_reasons",
    ]
    compact = decisions[compact_columns]
    write_csv(compact, ROOT / "oncoref/data/target-selection.csv.gz")
    gene_decisions = compact.assign(gene_id=compact.member_gene_ids.map(json.loads)).explode(
        "gene_id"
    )
    if not gene_decisions.gene_id.is_unique:
        raise ValueError("A member gene received multiple classifications")
    write_csv(gene_decisions, OUT / "gene-decisions.csv.gz")
    migrations = compact[compact.category.ne(compact.released_category)]
    write_csv(migrations, OUT / "category-changes.csv")
    for category in target_selection.CATEGORIES:
        name = category.lower().replace(" ", "-") + ".csv.gz"
        write_csv(compact[compact.category.eq(category)], OUT / name)
    metadata = {
        "schema": "oncoref.target-selection",
        "version": 1,
        "release_ready": False,
        "release_ready_meaning": "Scientific review completion, not software publication status. False while any listed research adjudication remains open; unresolved evidence stays excluded from positive sets.",
        "open_reviews": [
            "Tumor-source comparability and outlier support",
            "Read-level confirmation of the rare-cell mapping-risk flags",
            "Hair-follicle normal-protein replay (#687)",
        ],
        "summary": summary,
        "comparison": comparison,
        "data_sha256": sha256(ROOT / "oncoref/data/target-selection.csv.gz"),
        "report_url": "https://pirl-unc.github.io/oncoref/audits/target-selection-20261110/rebuild/report.html",
    }
    (ROOT / "oncoref/data/target-selection.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reuse-tumor-cache", action="store_true")
    parser.add_argument(
        "--normal-stat",
        choices=list(target_selection.NORMAL_TISSUE_STATISTICS),
        default="p95",
        help="Within-tissue healthy normal TPM summary (default: p95; optional: mean_plus_2std)",
    )
    args = parser.parse_args()
    build(reuse_tumor_cache=args.reuse_tumor_cache, normal_stat=args.normal_stat)
