#!/usr/bin/env python3
"""Separate broad literature coverage from the optional tumor-paper filter."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from oncoref import cta_literature, curation, target_filters, target_selection  # noqa: E402

OUT = ROOT / "docs/audits/antigen-evidence-review-20261117"
ORIGINAL = ROOT / "docs/audits/cta-count-history-20261117/standard-paper-exclusions.csv"


def main():
    original = pd.read_csv(ORIGINAL)
    original = original[original.papers_counted.eq(0)]
    links = cta_literature.current_attestations()
    counts = target_filters.paper_counts()
    subset = links[links.entity_id.isin(original.entity_id)]
    rows = []
    for record in original.itertuples():
        papers = subset[subset.entity_id.eq(record.entity_id)]
        values = counts.loc[record.entity_id].to_dict() if record.entity_id in counts.index else {}
        rows.append(
            {
                "entity_id": record.entity_id,
                "protein_group_label": record.protein_group_label,
                "all_curated_sources": papers.source_id.nunique(),
                "coverage_eligible_sources": papers.loc[
                    papers.eligible_for_coverage, "source_id"
                ].nunique(),
                "counted_tumor_support_sources": values.get("papers_counted", 0),
                "reproductive_sources": values.get("papers_reproductive", 0),
                "source_ids": json.dumps(sorted(papers.source_id.unique())),
            }
        )
    broad = pd.DataFrame(rows)
    OUT.mkdir(parents=True, exist_ok=True)
    broad.to_csv(OUT / "original-49-literature-coverage.csv", index=False)
    columns = [
        "entity_id",
        "protein",
        "source_id",
        "source_url",
        "citation",
        "title",
        "source_roles",
        "gene_ids",
        "claim_scope",
        "eligible_for_coverage",
        "published_nomination",
        "recorded_entity_id",
    ]
    details = (
        subset[columns]
        .drop_duplicates(["entity_id", "source_id", "claim_scope", "eligible_for_coverage"])
        .copy()
    )
    for col in ("source_roles", "gene_ids"):
        details[col] = details[col].map(json.dumps)
    details.to_csv(OUT / "original-49-paper-links.csv", index=False)
    catalog = target_selection.results()
    filtered = target_filters.apply(
        target_filters.FILTER_MODES["standard"], catalog[catalog.is_cta]
    )
    alive = pd.Series(True, index=filtered.index)
    for layer in target_filters.LAYERS:
        if layer == "papers":
            exclusions = filtered[alive & ~filtered.pass_papers]
        alive &= filtered["pass_" + layer]
    summary = {
        "catalog_policy": target_selection.POLICY_ID,
        "filter_policy": target_filters.FILTER_POLICY_ID,
        "catalog_sha256": target_selection.metadata()["data_sha256"],
        "original_49_sha256": curation.digest(ORIGINAL),
        "original_zero_counted": len(broad),
        "with_any_literature": int(broad.all_curated_sources.gt(0).sum()),
        "with_coverage_eligible_literature": int(broad.coverage_eligible_sources.gt(0).sum()),
        "with_reproductive_literature": int(broad.reproductive_sources.gt(0).sum()),
        "no_curated_links": broad.loc[
            broad.all_curated_sources.eq(0), "protein_group_label"
        ].tolist(),
        "standard_sequential_paper_exclusions": len(exclusions),
        "excluded_by_count": {
            str(k): int(v) for k, v in exclusions.papers_counted.value_counts().items()
        },
        "cta_modes": {
            mode: int(target_filters.apply(filt, catalog[catalog.is_cta]).keep.sum())
            for mode, filt in target_filters.FILTER_MODES.items()
        },
        "limitations": [
            "A source count describes curated coverage, not an exhaustive literature search or independent validation.",
            "Normal counterevidence, reproductive biology, related-copy assays and nominations have distinct claims.",
            "No paper count establishes peptide presentation or clinical safety.",
        ],
    }
    (OUT / "literature-summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary))


if __name__ == "__main__":
    main()
