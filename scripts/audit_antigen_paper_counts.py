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


def publication_report():
    data = cta_literature.publication_registry()
    # A generated view of the central ledger, never a second editable curation store.
    (OUT / "reviewed-paper-claims.json").write_text(json.dumps(data, indent=2) + "\n")
    text = [
        "# Centrally reviewed papers and gene claims\n",
        "Generated from `oncoref/data/cta-literature-followup.json` (schema 2). "
        "Edit that central ledger, not this report. Original nomination tables are "
        "immutable source inputs. Reviewed claims override inherited scope for the "
        "same paper/gene in the current catalog, API and publication filters.\n",
        "Gene lists show reviewed claims, not every gene mentioned by a paper. "
        "Reproductive expression, translation and functional experiments are separate. "
        "A publication use is neither automatic tenfold qualification nor clinical clearance.\n",
    ]
    for paper in data["sources"]:
        text.append(f"## [{paper['citation']}]({paper['url']})\n\n{paper['title']}\n")
        text.append(
            "| Gene | Original name | Reported as | Context | Allowed uses |\n|---|---|---|---|---|"
        )
        for claim in paper["gene_claims"]:
            fields = [claim[k] for k in ("protein", "reported_name", "reported_as", "context")]
            fields.append(", ".join(claim["allowed_uses"]))
            text.append("| " + " | ".join(v.replace("|", " / ") for v in fields) + " |")
        text.append("")
        for claim in paper["gene_claims"]:
            text.append(
                f"**{claim['protein']}** (`{', '.join(claim['gene_ids'])}`): "
                f"{claim['assay']}. Source: {claim['source_location']}. "
                f"Limitations: {claim['limitations']}\n"
            )
    (OUT / "reviewed-paper-claims.md").write_text("\n".join(text) + "\n")


def version_comparison(catalog):
    history = pd.read_csv(ORIGINAL.parent / "catalog-counts.csv", dtype={"version": str})
    notes = {
        "1.36": "Earlier catalog definition",
        "1.37": "CA/RA set algebra: 74 leave, 93 enter",
        "1.38": "153 to RA-only; 11 to Other:Sequence; 9 to Unsafe",
        "1.39": "22 to RA-only after tumor QC; 6 CTAs restored",
        "1.40": "RA 10×, individual normal samples, review holds and restrictive embryo gate",
        "1.41": "Embryo gate repaired; 44 CTAs restored",
        "699": "Compartment and host views; membership unchanged",
        "700": "Explicit host choices and gene-specific paper uses; membership unchanged",
    }
    lines = [
        "# Catalog versions and optional selection losses\n",
        "Counts are identical selected-reference protein groups. CA-only and RA-only "
        "exclude the CTA intersection. The either-set count includes CTA once.\n",
        "| Version | CTA | CA-only | RA-only | Either CA or RA | Change in CTA | Main change |",
        "|---|---:|---:|---:|---:|---:|---|",
    ]
    previous = None
    for row in history.to_dict("records"):
        version = str(row["version"])
        cta = int(row["CTA"])
        ca, ra = row["Cancer Antigen only"], row["Reproductive Antigen only"]
        cells = [str(cta), *["—" if pd.isna(v) else str(int(v)) for v in (ca, ra)]]
        cells.append("—" if pd.isna(ca) else str(int(cta + ca + ra)))
        delta = "—" if previous is None else f"{cta - previous:+d}"
        label = "PR #" + version if version in {"699", "700"} else version
        lines.append("| " + " | ".join([label, *cells, delta, notes[version]]) + " |")
        previous = cta
    counts = catalog.category.value_counts()
    current_cta = int(counts["CTA"])
    lines.append(
        "| "
        + " | ".join(
            [
                f"Current ({target_selection.POLICY_ID})",
                str(current_cta),
                str(int(counts["Cancer Antigen only"])),
                str(int(counts["Reproductive Antigen only"])),
                str(int(catalog.is_target_candidate.sum())),
                f"{current_cta - previous:+d}",
                "GTEx individual reproductive witnesses and within-site normal donor P99; reproductive-program flags are advisory",
            ]
        )
        + " |"
    )
    lines += [
        "\nThe frozen 1.39 → 1.41 difference is 147: 78 CA-only, 35 Other:Review, "
        "15 Other:Tissue, 12 Other:Expression and 7 RA-only. None became Unsafe. "
        "These are observed destinations, not independent causal attributions for overlapping gates.\n",
        f"Each optional mode below starts from the same {int(catalog.is_cta.sum())} CTAs; "
        "losses are sequential within that mode. "
        "Default anatomy protects no reproductive tissue and removes zero.\n",
        "| Mode | Normal-RNA loss | Curated-home loss | Paper loss | MAGE-name loss | CTAs kept |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for mode, filt in target_filters.FILTER_MODES.items():
        assessed = target_filters.apply(filt, catalog[catalog.is_cta])
        alive = pd.Series(True, index=assessed.index)
        losses = {}
        for layer in target_filters.LAYERS:
            losses[layer] = int((alive & ~assessed["pass_" + layer]).sum())
            alive &= assessed["pass_" + layer]
        values = [losses[k] for k in ("off_target", "curation", "papers", "exclusions")]
        lines.append("| " + " | ".join([mode, *map(str, values), str(int(alive.sum()))]) + " |")
    lines += [
        "\nStandard requires two counted tumor-support papers; reproductive-only papers remain "
        "separate support. The [IHC counterfactual](index.md#how-much-would-reproductive-high-ihc-add) "
        "is not enabled.\n",
        f"Catalog `{target_selection.POLICY_ID}`; filter `{target_filters.FILTER_POLICY_ID}`. "
        "Historical snapshots and exact transitions: [count history](../cta-count-history-20261117/index.md).\n",
    ]
    (OUT / "version-comparison.md").write_text("\n".join(lines))


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
    publication_report()
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
        "allowed_uses",
    ]
    details = (
        subset[columns]
        .drop_duplicates(["entity_id", "source_id", "claim_scope", "eligible_for_coverage"])
        .copy()
    )
    for col in ("source_roles", "gene_ids", "allowed_uses"):
        details[col] = details[col].map(json.dumps)
    details.to_csv(OUT / "original-49-paper-links.csv", index=False)
    catalog = target_selection.results()
    version_comparison(catalog)
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
