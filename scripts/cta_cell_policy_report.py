#!/usr/bin/env python3
"""Audit proposed dual cell-RNA cutoffs without changing released membership."""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from oncoref import (  # noqa: E402
    cell_rna_policy,
    cta_provenance,
    hpa,
    normal_cell_data,
    normal_cell_screen,
    reference_data,
)

CELL = ROOT / "docs/audits/cta-normal-cells-20261023"
KEYS = [
    "gene_id",
    "source",
    "version",
    "taxonomy_id",
    "tissue",
    "cell_type",
    "state",
    "stage",
    "unit",
    "category",
    "mapping_status",
]


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write(frame, path):
    frame.to_csv(
        path,
        index=False,
        compression={"method": "gzip", "mtime": 0} if path.name.endswith(".gz") else None,
    )


def summarize(frame, *, support, rank=1, unique=None):
    """Nth supported value preserves recurrence cutoffs without donor invention."""
    if unique and frame.duplicated([*KEYS, unique]).any():
        raise ValueError("Repeated sample identities would inflate recurrence")
    d = frame.copy()
    d["supported_value"] = d.value.where(support)
    base = d.groupby(KEYS, dropna=False).agg(
        observations=("value", "size"),
        supported_observations=("supported_value", "count"),
        raw_maximum=("value", "max"),
    )
    eligible = d[d.supported_value.notna()].sort_values("supported_value", ascending=False)
    eligible["rank"] = eligible.groupby(KEYS, dropna=False).cumcount() + 1
    value = eligible.loc[eligible["rank"].eq(rank)].set_index(KEYS).supported_value
    return base.join(value.rename("value")).reset_index()


def scenario_thresholds():
    # Explicit exploratory choices, not estimated/validated clinical limits.
    # Existing 5-native-unit immune triggers and 1-CP10K corneal trigger anchor
    # parts of this grid. Other boundaries and IHC multipliers are hypotheses.
    rows = []
    for assay, unit in cell_rna_policy.EXPECTED_UNITS.items():
        tiers = {"vital": (2, 5), "protected": (5, 10), "manageable": (10, 20)}
        if assay == "cornea_collin_adult_counts":
            tiers = {"vital": (1, 2)}
        elif assay != "hpa_single_cell_tissue":
            tiers = {"vital": (2, 5)}
        for tier, (low, high) in tiers.items():
            rows.append(
                {
                    "assay": assay,
                    "unit": unit,
                    "category": tier,
                    "cta_cutoff": low,
                    "reference_cutoff": high,
                }
            )
    return pd.DataFrame(rows)


def build(out):
    if out.exists() and any(out.iterdir()):
        raise ValueError("Use an empty output directory")
    out.mkdir(parents=True, exist_ok=True)
    inputs = {}

    def read(path, **kwargs):
        inputs[str(path.relative_to(ROOT))] = digest(path)
        return pd.read_csv(path, low_memory=False, **kwargs)

    research_path = ROOT / "oncoref/data/cta-research-gene-sets.csv.gz"
    current = read(research_path)
    current = current[current.Nominated].copy()
    ids = set(current["Gene ID"])
    if len(ids) != 2540 or current.Target.sum() != 215:
        raise ValueError("Unexpected frozen research-policy baseline")
    atlas = read(ROOT / "docs/audits/cta-atlas-20261024/all-genes.csv")
    controls = set(atlas.loc[atlas["Reproductive positive control"], "Gene ID"])
    conflicts = set(current.loc[current["Antibody concern"], "Gene ID"])
    frames, coverage = [], []
    raw = read(CELL / "cell-rna.csv.gz")
    raw = raw[raw.gene_id.isin(ids)].copy()
    frames.append(summarize(raw, support=raw.read_count.ge(20)))
    coverage.append(
        {"source": "HPA tissue clusters", "genes": raw.gene_id.nunique(), "rows": len(raw)}
    )
    raw = read(CELL / "immune-rna.csv.gz")
    raw = raw[raw.gene_id.isin(ids)].copy()
    for source, group in raw.groupby("source"):
        coverage.append({"source": source, "genes": group.gene_id.nunique(), "rows": len(group)})
    # HPA DICE summary is explicitly excluded from independent decision evidence.
    raw = normal_cell_screen.contexts(raw[~raw.source.eq("hpa_immune_schmiedel")])
    frames.append(summarize(raw, support=raw.value.notna()))
    raw = normal_cell_screen.contexts(read(CELL / "corneal-rna.csv.gz"))
    raw = raw[raw.gene_id.isin(ids)].copy()
    frames.append(
        summarize(
            raw, support=raw.positive_cells.ge(5) & raw.read_count.ge(20), rank=2, unique="sample"
        )
    )
    coverage.append(
        {"source": "Adult corneal atlas", "genes": raw.gene_id.nunique(), "rows": len(raw)}
    )
    receipts = []
    for source in normal_cell_data.sources().source:
        if not source.startswith("dice_"):
            continue
        raw = normal_cell_data.dice_observations(source, genes=ids)
        receipt = reference_data.provenance(source, raw.version.iloc[0])
        if digest(receipt["path"]) != receipt["sha256"]:
            raise ValueError("DICE source checksum mismatch")
        receipts.append({k: receipt[k] for k in ["name", "version", "sha256", "url"]})
        raw = normal_cell_screen.contexts(raw)
        frames.append(summarize(raw, support=raw.value.notna(), rank=3, unique="replicate_slot"))
        coverage.append({"source": source, "genes": raw.gene_id.nunique(), "rows": len(raw)})
    contexts = pd.concat(frames, ignore_index=True)
    contexts = contexts[contexts.gene_id.isin(ids)].copy()
    receipt = reference_data.provenance("hpa_normal_tissue", "v23")
    if digest(receipt["path"]) != receipt["sha256"]:
        raise ValueError("IHC checksum mismatch")
    receipts.append({k: receipt[k] for k in ["name", "version", "sha256", "url"]})
    contexts = cell_rna_policy.match_negative_ihc(
        contexts,
        hpa._read_hpa("hpa_normal_tissue", "v23"),
        positive_controls=controls,
        conflicting_genes=conflicts,
    )
    thresholds = scenario_thresholds()
    factors = {"Missing/uncertain": 1, "Approved": 1.25, "Supported": 1.5, "Enhanced": 2}
    write(thresholds, out / "illustrative-thresholds.csv")
    result = cell_rna_policy.assess(contexts, thresholds, ihc_factors=factors)
    result["gene_name"] = result.gene_id.map(current.set_index("Gene ID").Gene)
    write(result.sort_values(KEYS), out / "context-assessments.csv.gz")
    flagged = result[result.cta_excluded].sort_values(KEYS)
    write(flagged, out / "exclusion-findings.csv.gz")
    summary = result.groupby("gene_id").agg(
        cta_excluded=("cta_excluded", "any"),
        reference_review=("reference_review", "any"),
        assessed_contexts=("value", "count"),
        total_contexts=("value", "size"),
        exact_ihc_contexts=("ihc_exact_match", "sum"),
    )
    reasons = flagged.groupby("gene_id").agg(
        flagged_contexts=("source", "size"),
        exclusion_sources=("source", lambda x: ";".join(sorted(set(x)))),
        exclusion_tissue_tiers=("category", lambda x: ";".join(sorted(set(x)))),
        exclusion_ihc_tiers=("negative_ihc_grade", lambda x: ";".join(sorted(set(x)))),
        exclusion_native_units=("unit", lambda x: ";".join(sorted(set(x)))),
    )
    summary = summary.join(reasons)
    genes = current[
        ["Gene ID", "Gene", "Category", "Detail", "Target", "Expression forbidden"]
    ].merge(summary, left_on="Gene ID", right_index=True, how="left", validate="one_to_one")
    genes["no_source_records"] = genes.total_contexts.isna()
    for col in ["cta_excluded", "reference_review"]:
        genes[col] = genes[col].fillna(False).astype(bool)
    genes["retained_after_overlay"] = genes.Target & ~genes.cta_excluded
    genes["middle_zone"] = genes.cta_excluded & ~genes.reference_review
    write(genes, out / "gene-assessments.csv")
    write(genes[genes.cta_excluded], out / "proposed-cta-exclusions.csv")
    write(genes[genes.reference_review], out / "proposed-reference-review.csv")
    write(
        genes[genes.middle_zone | genes.no_source_records], out / "uncertain-reference-review.csv"
    )
    write(genes[genes.Target & genes.cta_excluded], out / "selected-gene-changes.csv")
    sensitivity = []
    for lower_factor in [0.5, 1, 1.5]:
        for upper_factor in [1, 2, 3]:
            t = thresholds.copy()
            t.cta_cutoff *= lower_factor
            t.reference_cutoff *= upper_factor
            r = cell_rna_policy.assess(contexts, t, ihc_factors=factors)
            excluded = set(r.loc[r.cta_excluded, "gene_id"])
            reference = set(r.loc[r.reference_review, "gene_id"])
            sensitivity.append(
                {
                    "lower_factor": lower_factor,
                    "upper_factor": upper_factor,
                    "cta_excluded": len(excluded),
                    "reference_review": len(reference),
                    "retained_of_215": len(set(current.loc[current.Target, "Gene ID"]) - excluded),
                }
            )
    write(pd.DataFrame(sensitivity), out / "threshold-sensitivity.csv")
    no_relax = cell_rna_policy.assess(contexts, thresholds, ihc_factors=dict.fromkeys(factors, 1))
    write(
        no_relax.groupby("gene_id")[["cta_excluded", "reference_review"]].any().reset_index(),
        out / "without-ihc-relaxation.csv",
    )
    write(pd.DataFrame(coverage), out / "filter-source-coverage.csv")

    # Publication counts deduplicate by paper, never sum nested source tables.
    refs = cta_provenance.paper_membership()
    membership = refs[refs.protein_candidate_eligible].copy()
    chosen = set(current.loc[current.Target, "Gene ID"])
    sets = {
        doi: set(membership.loc[membership.doi.eq(doi), "Ensembl_Gene_ID"])
        for doi in sorted(refs.doi.unique())
    }
    paper_rows = []
    for doi, members in sets.items():
        other = set().union(*(v for k, v in sets.items() if k != doi))
        paper_rows.append(
            {
                "doi": doi,
                "source_tags": ";".join(sorted(set(refs.loc[refs.doi.eq(doi), "source_tag"]))),
                "genes": len(members),
                "selected": len(members & chosen),
                "unique_to_paper": len(members - other),
            }
        )
    write(
        pd.DataFrame(paper_rows).sort_values("genes", ascending=False),
        out / "nomination-sources.csv",
    )
    write(
        pd.DataFrame(
            [
                {"source_a": a, "source_b": b, "overlap": len(av & bv)}
                for a, av in sets.items()
                for b, bv in sets.items()
            ]
        ),
        out / "source-overlaps.csv",
    )
    write(
        membership[["doi", "source_tag", "Ensembl_Gene_ID", "Symbol"]].drop_duplicates(),
        out / "nomination-membership.csv.gz",
    )
    extras = current[~current["Gene ID"].isin(set().union(*sets.values()))]
    write(extras[["Gene ID", "Gene", "Category"]], out / "additional-trial-nominations.csv")
    # Disjoint stages of the current decision result, not a reconstructed history.
    remaining = current.copy()
    funnel = [{"stage": "Mapped nominations", "remaining": len(remaining), "removed": 0}]
    for label, mask in [
        ("After accepted normal exclusions", remaining.Category.ne("Unsafe")),
        ("After unresolved evidence and admission holds", remaining.Category.eq("CTA")),
    ]:
        n = len(remaining)
        remaining = remaining.loc[mask.reindex(remaining.index)]
        funnel.append({"stage": label, "remaining": len(remaining), "removed": n - len(remaining)})
    funnel.append(
        {
            "stage": "Illustrative cell-RNA overlay",
            "remaining": int(genes.retained_after_overlay.sum()),
            "removed": int((genes.Target & genes.cta_excluded).sum()),
        }
    )
    write(pd.DataFrame(funnel), out / "filtering-funnel.csv")
    write(
        current.groupby(["Category", "Detail"]).size().rename("genes").reset_index(),
        out / "current-decision-breakdown.csv",
    )
    plot = current.merge(
        atlas[["Gene ID", "Normal RNA maximum (nTPM)", "Normal RNA maximum (TPM)"]],
        on="Gene ID",
        validate="one_to_one",
    )
    scored = result[result.category.isin(cell_rna_policy.NORMAL)].copy()
    scored["cta_ratio"] = scored.value / scored.cta_cutoff
    plot["maximum_cell_cutoff_ratio"] = plot["Gene ID"].map(
        scored.groupby("gene_id").cta_ratio.max()
    )
    tissue = read(ROOT / "docs/audits/cta-screen-20261023/tissue-evidence.csv.gz")
    adult = tissue[tissue.screening_tier.isin(cell_rna_policy.NORMAL)]
    plot["adult_bulk_max_ntpm"] = plot["Gene ID"].map(
        adult.groupby("Ensembl_Gene_ID").rna_ntpm.max()
    )
    write(
        plot[
            [
                "Gene ID",
                "Gene",
                "Category",
                "Detail",
                "Target",
                "Reproductive RNA (nTPM)",
                "Normal RNA maximum (nTPM)",
                "Normal RNA maximum (TPM)",
                "maximum_cell_cutoff_ratio",
                "adult_bulk_max_ntpm",
            ]
        ],
        out / "plot-data.csv",
    )
    highlights = result[result.gene_name.isin(["PRAME", "MAGEA4", "CTAG1A", "CTAG1B", "TBXT"])]
    write(highlights, out / "highlight-contexts.csv.gz")
    summary_json = {
        "status": "Illustrative calibration only; not an activated classification policy",
        "current_selected": 215,
        "retained_after_overlay": int(genes.retained_after_overlay.sum()),
        "cta_excluded": int(genes.cta_excluded.sum()),
        "reference_review": int(genes.reference_review.sum()),
        "middle_zone": int(genes.middle_zone.sum()),
        "no_source_records": int(genes.no_source_records.sum()),
        "negative_ihc_factors": factors,
        "paper_union": len(set().union(*sets.values())),
        "papers": len(sets),
        "additional_trial_nominations": len(extras),
        "native_units": "Kept separate by assay; thresholds are not unit conversions",
        "ihc_limit": "Exact labels plus positive profile control; no matched donors or cell-specific sensitivity calibration",
        "overlay_limit": "Adds exclusions only. Does not rescue prior holds or treat incomplete coverage as absence.",
        "reference_limit": "RNA evidence lists only. No FASTA/reference membership is regenerated.",
    }
    (out / "summary.json").write_text(json.dumps(summary_json, indent=2) + "\n")
    for p in [
        Path(__file__),
        ROOT / "oncoref/cell_rna_policy.py",
        ROOT / "oncoref/data/cta-publication-membership.csv",
        ROOT / "oncoref/data/cta-publication-sources.csv",
    ]:
        if p.exists():
            inputs[str(p.relative_to(ROOT))] = digest(p)
    manifest = {
        "inputs": inputs,
        "retrieved_sources": receipts,
        "outputs": {p.name: digest(p) for p in sorted(out.iterdir())},
    }
    (out / "run-manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(summary_json, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    build(parser.parse_args().output)
