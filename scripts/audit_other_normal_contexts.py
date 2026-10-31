#!/usr/bin/env python3
"""Trace appreciable normal measurements using the production source adapters."""

from __future__ import annotations

import shutil

import pandas as pd
from replay_cta_selection import OUTPUT, ROOT, ContextCollector, replay_contexts

from oncoref import cell_rna_policy
from oncoref import cta_selection_policy as policy

OUT = ROOT / "docs/audits/cta-other-literature-20261030"


class ReviewCollector(ContextCollector):
    def __init__(self, facts, output):
        super().__init__(facts, output)
        self.ids = set(facts.representative_gene_id)
        self.observations = []

    def collect(self, frame, assay):
        frame = frame[frame.gene_id.isin(self.ids)].copy()
        if frame.empty:
            return
        d = policy.assess_contexts(frame)
        raw = d.get("raw_maximum", d.value).fillna(d.value)
        floor = policy.APPRECIABLE_RNA * d.unit.map(policy.UNIT_SCALE)
        d["supported_at_least_floor"] = d.value.ge(floor)
        d["raw_at_least_floor"] = raw.ge(floor)
        d["raw_value"] = raw
        d["assay"] = assay
        d = d[d.supported_at_least_floor | d.raw_at_least_floor]
        columns = [
            "gene_id",
            "assay",
            "source",
            "version",
            "tissue",
            "cell_type",
            "category",
            "selection_tier",
            "value",
            "raw_value",
            "unit",
            "Read count",
            "negative_ihc_grade",
            "candidate_ceiling",
            "supported_at_least_floor",
            "raw_at_least_floor",
        ]
        self.observations.append(d.reindex(columns=columns))
        super().collect(frame, assay)


def main():
    facts = pd.read_csv(OUT / "baseline-proteins.csv.gz", low_memory=False)
    work = ROOT / "tmp/cta-other-literature/normal-replay"
    work.mkdir(parents=True, exist_ok=True)
    for source in [
        "hpa_rna_consensus",
        "hpa_single_cell_tissue",
        "hpa_immune_cell",
        "hpa_immune_monaco",
    ]:
        name = f"{source}-accepted-identities.csv.gz"
        shutil.copyfile(OUTPUT / name, work / name)
    collector = ReviewCollector(facts, work)
    replay_contexts(facts, work, collector=collector)
    frame = pd.concat(collector.observations, ignore_index=True).drop_duplicates()
    frame["protein_label"] = frame.gene_id.map(
        facts.set_index("representative_gene_id").protein_group_label
    )
    frame["baseline_subtype"] = frame.gene_id.map(facts.set_index("representative_gene_id").subtype)
    frame.to_csv(
        OUT / "normal-contexts.csv.gz", index=False, compression={"method": "gzip", "mtime": 0}
    )
    normal = frame[
        frame.baseline_subtype.eq("Normal")
        & frame.supported_at_least_floor
        & ~frame.category.isin(
            ["reproductive", "developmental", cell_rna_policy.NEUTRAL_REPRODUCTIVE_IMMUNE]
        )
    ]
    normal.to_csv(OUT / "normal-category-tissues.csv", index=False)
    summary = (
        normal.groupby(["tissue", "selection_tier", "unit"])
        .agg(proteins=("gene_id", "nunique"), maximum=("value", "max"))
        .reset_index()
        .sort_values("proteins", ascending=False)
    )
    summary.to_csv(OUT / "normal-tissue-summary.csv", index=False)
    print(summary.head(30).to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
