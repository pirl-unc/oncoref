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

from build_cta_alternative_evidence import GTEX_TISSUES, anatomy  # noqa: E402
from build_cta_candidate_evidence import replay_updates  # noqa: E402
from cta_report_common import sha256, write_csv  # noqa: E402

from oncoref import (  # noqa: E402
    cta_literature,
    reference_data,
    source_matrices,
    target_selection,
    tissue_hazards,
)
from oncoref.gene_matching import (  # noqa: E402
    harmonize_reference_expression,
    match_reference_genes,
)

BASE = ROOT / "docs/audits/cta-selection-20261028"
FROZEN = ROOT / "docs/audits/cta-transcript-selection-20261108/inputs/frozen"
RAW = ROOT / "tmp/target-rebuild-inputs"
OUT = ROOT / "docs/audits/target-selection-20261109/rebuild"
CACHE = ROOT / "tmp/target-catalog-cache"
COMPARISON_ID = "oncoref.native-gene-tpm-normal-reference.2026-11-09.1"
INPUTS = {}


def read(path, **kwargs):
    path = Path(path)
    INPUTS[str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else str(path)] = sha256(path)
    return pd.read_csv(path, low_memory=False, **kwargs)


def shared_facts():
    for name in [
        "gene-knowledge.json",
        "gene-knowledge-updates.json",
        "normal-tissue-hazards.yaml",
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
    f["unsafe"] = f.category.eq("Unsafe")
    admission = {
        "Reproductive / embryonic association",
        "Cancer expression",
        "Reproductive concentration",
        "Translation qualification",
        "Required RNA coverage",
    }
    normal_holds = [c for c in f if c.startswith("hold:") and c[5:] not in admission]
    f["normal_screen_pass"] = ~f[normal_holds].any(axis=1)
    f["required_coverage_complete"] = f.required_coverage
    f["identity_translation_pass"] = f[
        [
            "translation_support_eligible",
            "translation_qualified",
            "representative_selection_approved",
        ]
    ].all(axis=1)
    f["evidence_resolved"] = ~f.unresolved_evidence
    f["sequence_pass"] = ~f.nested_microprotein & f.unsafe_sequence_overlap_pct.lt(95)
    f["shared_eligible"] = ~f.unsafe & f[list(target_selection.SHARED_CHECKS)].all(axis=1)
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
    f["reproductive_routes"] = routes.apply(lambda r: " | ".join(r.index[r]), axis=1)
    measured = (
        f.reproductive_ntpm.notna() | f.reproductive_cell_ntpm.notna() | f.developmental_measured
    )
    f["reproductive_evidence"] = np.select(
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


def normal_panel(facts):
    rows = []
    downloads = OUT / "downloaded-hpa-sources.json"
    INPUTS[str(downloads.relative_to(ROOT))] = sha256(downloads)
    for source in json.loads(downloads.read_text()):
        if sha256(ROOT / source["path"]) != source["sha256"]:
            raise ValueError(f"Changed normal source: {source['source']}")
    hpa_ids = read(BASE / "hpa_rna_consensus-accepted-identities.csv.gz")
    lookup = hpa_ids.set_index("source_gene_id").canonical_gene_id
    for name, context in [
        ("rna_tissue_hpa", "Tissue"),
        ("rna_tissue_gtex", "Tissue"),
        ("rna_brain_gtex", "Brain region"),
    ]:
        raw = read(RAW / f"{name}.tsv.zip", sep="\t")
        raw["gene_id"] = raw.Gene.map(lookup)
        raw["tissue"] = raw[context]
        raw["source_context"] = raw[context]
        raw["value"] = raw.TPM
        rows.append(
            group_normal_rows(raw, facts, source=name + ":v23", statistic="tissue mean TPM")
        )
    raw = read(ROOT / "tmp/cta-source-consolidation/rna_brain_hpa.tsv.zip", sep="\t")
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

    raw = read(
        ROOT / "tmp/cta-source-consolidation/gtex-v11-median-tpm.gct.gz", sep="\t", skiprows=2
    )
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
    normal["tissue"] = normal.tissue.replace(
        {"nucleus accumbens": "basal ganglia", "putamen": "basal ganglia"}
    )
    normal = anatomy(normal.rename(columns={"source_id": "source"}))
    normal["category"] = tissue_hazards.context_tiers(normal)
    normal = normal.rename(columns={"source": "source_id"})
    normal["scope"] = normal.category.where(
        normal.category.isin(["reproductive", "disposable"]), "normal"
    )
    normal["context_id"] = normal.tissue
    normal["qc_pass"] = normal.value.notna()
    normal["identity_resolved"] = normal.measured_members.eq(normal.member_count)
    # Same combined anatomical coverage rule, now requiring native TPM support.
    consensus = read(
        reference_data.local_path("hpa_rna_consensus", "v23"), sep="\t", usecols=["Tissue"]
    )
    concepts = anatomy(pd.DataFrame({"tissue": sorted(consensus.Tissue.unique())}))
    concepts = concepts.assign(source="hpa_rna_consensus", cell_type="bulk")
    concepts["category"] = tissue_hazards.context_tiers(concepts)
    required = set(concepts.loc[~concepts.category.isin(["reproductive", "disposable"]), "tissue"])

    # DICE is already a supported native-TPM statistic. Unpaired protein groups
    # retain the conservative upper bound rather than adding unrelated donors.
    immune = read(FROZEN / "immune-cornea-contexts.csv.gz")
    groups = read(FROZEN / "dice-group-contexts.csv.gz")
    immune = pd.concat(
        [immune[immune.unit.eq("TPM") & immune.source.str.startswith("dice_")], groups],
        ignore_index=True,
    )
    rep_to_entity = facts.set_index("representative_gene_id").entity_id
    immune = immune[immune.gene_id.isin(rep_to_entity.index)].copy()
    immune["entity_id"] = immune.gene_id.map(rep_to_entity)
    immune["value"] = immune.get("value_upper", pd.Series(np.nan, index=immune.index)).fillna(
        immune.value
    )
    immune["source_id"] = immune.source.astype(str) + ":" + immune.version.astype(str)
    immune["source_context"] = immune.cell_type
    immune["context_id"] = "blood"
    immune["scope"] = "normal"
    immune["statistic"] = "supported immune TPM; DICE third-highest, group upper bound if unpaired"
    immune["qc_pass"] = immune.value.notna()
    immune["identity_resolved"] = True
    monaco_path = reference_data.local_path("hpa_immune_monaco", "v23")
    monaco = read(monaco_path, sep="\t")
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
    cols = [
        "entity_id",
        "source_id",
        "comparison_id",
        "tissue",
        "cell_type",
        "source_context",
        "context_id",
        "scope",
        "statistic",
        "unit",
        "value",
        "qc_pass",
        "identity_resolved",
    ]
    normal = (
        normal[cols]
        .sort_values(["entity_id", "source_id", "source_context"])
        .reset_index(drop=True)
    )
    write_csv(normal, OUT / "normal-comparison.csv.gz")
    return normal, required


def participant_id(namespace, sample):
    """Only infer donor identity for documented project identifier schemes."""
    if re.match(r"^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}(?:-|$)", sample):
        return "TCGA:" + sample[:12]
    if re.match(r"^TARGET-\d{2}-[A-Z0-9]{6}(?:-|$)", sample):
        return "TARGET:" + "-".join(sample.split("-")[:3])
    if re.match(r"^MMRF[_-]\d+[_-]", sample):
        return "MMRF:" + re.split(r"[_-]", sample)[1]
    return ""


def annotate_papers(decisions):
    """Retain the full curated occurrence ledger, joined by exact member loci."""
    links = cta_literature.attestations()
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
    pieces, coverage = [], []
    identity_cache = {}
    for job in available.itertuples():
        path, digest = paths[job.cancer_code]
        if sha256(path) != digest:
            raise ValueError(f"Changed tumor matrix: {path}")
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
        raw = pq.read_table(
            path, filters=[("Ensembl_Gene_ID", "in", identity.source_gene_id.tolist())]
        ).to_pandas()
        sample_ids = (
            samples.loc[samples.cancer_code.eq(job.cancer_code), "sample_id"].astype(str).tolist()
        )
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
    combined = pd.concat(pieces, axis=1).T
    if combined.index.duplicated().any():
        repeated = combined[combined.index.duplicated(keep=False)]
        if repeated.groupby(level=0).nunique(dropna=False).gt(1).any().any():
            raise ValueError("Overlapping cohort views disagree")
        combined = combined[~combined.index.duplicated()]
    sample_meta = samples.drop_duplicates(["namespace", "sample_id"]).copy()
    sample_meta["sample_key"] = sample_meta.namespace + ":" + sample_meta.sample_id
    sample_meta["participant_id"] = [
        participant_id(r.namespace, r.sample_id) for r in sample_meta.itertuples()
    ]
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


def build():
    OUT.mkdir(parents=True, exist_ok=True)
    facts = shared_facts()
    selected = facts[facts.shared_eligible].copy()
    print(f"Common eligibility: {len(selected)} / {len(facts)}", flush=True)
    normal, required = normal_panel(selected)
    matrix, meta = tumors(selected)
    meta = meta.set_index("sample_key").loc[matrix.index]
    summaries, witnesses = [], []
    grouped = dict(iter(normal.groupby("entity_id")))
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
        witnesses.extend(top.to_dict("records"))
        summaries.append(
            {
                "entity_id": row.entity_id,
                "cancer_evidence": receipt["status"],
                "cancer_routes": " | ".join(receipt["routes"]),
                "normal_reference_tpm": receipt["normal_maximum"],
                "normal_reference_contexts": json.dumps(
                    receipt["normal_maximum_observations"][:10]
                ),
                "normal_reference_tied_contexts": len(receipt["normal_maximum_observations"]),
                "missing_comparison_contexts": json.dumps(receipt["missing_normal_contexts"]),
                "independent_patients_above_normal": receipt[
                    "independent_participants_above_normal"
                ],
                "samples_above_normal": receipt["samples_above_normal"],
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
                "cancer_assessment_reason": receipt["reason"],
            }
        )
    evidence = pd.DataFrame(summaries)
    facts = facts.merge(evidence, on="entity_id", how="left", validate="one_to_one")
    facts["cancer_evidence"] = facts.cancer_evidence.fillna("unassessed")
    facts["cancer_assessment_reason"] = facts.cancer_assessment_reason.fillna(
        "Shared eligibility fails; positive cancer membership cannot be assigned"
    )
    decisions = annotate_papers(target_selection.classify(facts))
    reasons = {
        "normal_screen_pass": "Normal-expression evidence blocks admission",
        "required_coverage_complete": "Required normal coverage is incomplete",
        "identity_translation_pass": "Identity or translation requires review",
        "evidence_resolved": "Evidence is unresolved",
        "sequence_pass": "Nested or Unsafe-contained sequence",
    }
    decisions["decision_reason"] = [
        str(r.Unsafe_reasons)
        if r.category == "Unsafe"
        else " | ".join(text for field, text in reasons.items() if not getattr(r, field))
        if not r.shared_eligible
        else "Cancer and reproductive/developmental evidence qualify"
        if r.category == "CTA"
        else "Cancer evidence qualifies; no qualifying reproductive/developmental association recorded"
        if r.category == "Cancer Antigen only"
        else "Reproductive/developmental evidence qualifies; " + r.cancer_assessment_reason
        if r.category == "Reproductive Antigen only"
        else "Neither positive association qualifies; " + r.cancer_assessment_reason
        for r in decisions.itertuples()
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
        columns={c: "legacy_" + c for c in decisions if c.startswith(("hold:", "unsafe:"))}
    )
    decisions["policy"] = target_selection.POLICY_ID
    decisions["assessment_status"] = np.where(
        decisions.classification.eq("Other:Review"), "review_required", "screened"
    )
    write_csv(decisions, OUT / "protein-decisions.csv.gz")
    write_csv(pd.DataFrame(witnesses), OUT / "tumor-witnesses.csv.gz")
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
        "changes": decisions.groupby(["released_category", "category"])
        .size()
        .reset_index(name="proteins")
        .to_dict("records"),
        "shared_eligible": len(selected),
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
            "Complete 36-context combined nonreproductive normal anatomy panel, plus observed native-TPM normal immune and brain measurements",
        ],
        "normal_statistic": "Maximum supported context statistic: tissue/subregion means, GTEx v11 donor medians, and supported immune summaries including DICE third-highest TPM",
        "batch_corrected": False,
        "limitations": [
            "Quantifier, library preparation and study effects are retained; the fold rule is a screening convention",
            "A maximum of tissue summaries is not a maximum over every healthy donor",
            "Other assay scales continue to enforce the existing shared normal screen independently",
        ],
        "participant_identity": "TCGA, TARGET and MMRF participant schemes are resolved across aliquots/cohort views. Other source IDs support the single-sample route; unverified patient IDs cannot establish replication.",
        "normal_sources": json.loads((OUT / "downloaded-hpa-sources.json").read_text()),
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
        "is_target_candidate",
        "cancer_evidence",
        "cancer_routes",
        "cancer_assessment_reason",
        "reproductive_evidence",
        "reproductive_routes",
        "reproductive_ntpm",
        "reproductive_cell_ntpm",
        "developmental_mean_tpm",
        "normal_reference_tpm",
        "normal_reference_contexts",
        "normal_reference_tied_contexts",
        "single_tumor_threshold_tpm",
        "independent_patients_above_normal",
        "unresolved_patient_samples_above_normal",
        "samples_above_normal",
        "tumor_samples_measured",
        "tumor_max_tpm",
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
        "open_reviews": [
            "Tumor-source comparability and outlier support",
            "Rare normal and developmental coverage",
        ],
        "summary": summary,
        "comparison": comparison,
        "data_sha256": sha256(ROOT / "oncoref/data/target-selection.csv.gz"),
        "report_url": "https://pirl-unc.github.io/oncoref/audits/target-selection-20261109/rebuild/report.html",
    }
    (ROOT / "oncoref/data/target-selection.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    build()
