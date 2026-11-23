#!/usr/bin/env python3
"""Audit published targets and normal-expression peaks without changing categories."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import rebuild_target_catalog as rebuild  # noqa: E402
import replay_cta_selection as replay  # noqa: E402

from oncoref import cta_selection_policy, hpa  # noqa: E402

BASE = ROOT / "docs/audits/target-selection-20261109"
OUT = BASE / "category-profiles"
CACHE = ROOT / "tmp/target-category-profiles"
INPUTS = {}


def sha(path):
    with Path(path).open("rb") as f:
        return hashlib.file_digest(f, "sha256").hexdigest()


def read(path, **kwargs):
    path = Path(path)
    INPUTS[str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else str(path)] = sha(path)
    return pd.read_csv(path, low_memory=False, **kwargs)


def write(frame, path):
    frame.to_csv(
        path,
        index=False,
        compression={"method": "gzip", "mtime": 0} if path.suffix == ".gz" else None,
    )


class PeakCollector:
    """Keep the largest supported value per unit and tissue tier, not value/cutoff."""

    def __init__(self, facts):
        self.facts = facts.set_index("representative_gene_id")
        self.blocks = []
        self.gaps = []
        self.member_context_counts = {}

    def collect(self, frame, assay):
        d = frame[frame.gene_id.isin(self.facts.index)].copy()
        # The legacy replay's HPA bulk/immune adapters summed members without
        # carrying counts into the resulting frame. Recover them from raw rows.
        for source, context in [
            ("hpa_rna_consensus", "Tissue"),
            ("hpa_immune_cell", "Immune cell"),
            ("hpa_immune_monaco", "Immune cell"),
        ]:
            selected = d.source.eq(source) & d.gene_id.map(self.facts.member_count).gt(1)
            if not selected.any():
                continue
            if source not in self.member_context_counts:
                lookup = (
                    pd.read_csv(CACHE / f"{source}-accepted-identities.csv.gz")
                    .set_index("source_gene_id")
                    .canonical_gene_id
                )
                member_to_rep = {
                    g: rep
                    for rep, row in self.facts[self.facts.member_count.gt(1)].iterrows()
                    for g in json.loads(row.member_gene_ids)
                }
                raw = hpa._read_hpa(source, "v23").copy()
                raw["member"] = raw.Gene.map(lookup)
                raw = raw[raw.member.isin(member_to_rep)]
                raw["gene_id"] = raw.member.map(member_to_rep)
                self.member_context_counts[source] = (
                    raw.groupby(["gene_id", context]).member.nunique().to_dict()
                )
            context_column = "tissue" if context == "Tissue" else "cell_type"
            counts = self.member_context_counts[source]
            d.loc[selected, "measured_members"] = [
                counts.get((r.gene_id, getattr(r, context_column)), 0)
                for r in d[selected].itertuples()
            ]
            d.loc[selected, "member_count"] = d.loc[selected, "gene_id"].map(
                self.facts.member_count
            )
        incomplete = pd.Series(False, index=d.index)
        if {"measured_members", "member_count"} <= set(d):
            incomplete |= (
                d.measured_members.notna()
                & d.member_count.notna()
                & d.measured_members.lt(d.member_count)
            )
        if "complete_members" in d:
            incomplete |= d.complete_members.eq(False)
        if incomplete.any():
            self.gaps.append(
                d.loc[
                    incomplete,
                    [
                        c
                        for c in [
                            "gene_id",
                            "source",
                            "tissue",
                            "cell_type",
                            "unit",
                            "value",
                            "value_upper",
                            "measured_members",
                            "member_count",
                        ]
                        if c in d
                    ],
                ].assign(reason="incomplete_protein_members")
            )
            d = d[~incomplete].copy()
        for col in [
            "clinical_rna_credit",
            "clinical_product_ids",
            "clinical_cohort_ids",
            "clinical_credit_scope",
        ]:
            d[col] = d.gene_id.map(self.facts[col])
        d = cta_selection_policy.assess_contexts(d)
        d = d[
            d.selection_tier.isin(["vital", "protected", "manageable", "disposable"])
            & d.mapping_status.eq("mapped")
        ].copy()
        d["statistic"] = d.get("statistic", pd.Series("", index=d.index)).fillna("")
        for bound in ["value_upper", "upper_supported_value"]:
            if bound in d:
                replace = np.isfinite(d[bound]) & (d.value.isna() | d[bound].gt(d.value))
                d.loc[replace, "value"] = d.loc[replace, bound]
                d.loc[replace, "statistic"] = "unpaired member upper bound"
        if "complete_members" in d:
            d = d[~d.complete_members.eq(False)]  # missing means no supplied group flag
        d = d[d.value.notna()]
        d["assay"] = assay
        d["source_context"] = d.get("source_context", pd.Series("", index=d.index)).fillna("")
        d["scope"] = d.selection_tier.map(lambda t: "disposable" if t == "disposable" else "normal")
        d["entity_id"] = d.gene_id.map(self.facts.entity_id)
        fields = [
            "entity_id",
            "gene_id",
            "source",
            "version",
            "tissue",
            "cell_type",
            "source_context",
            "scope",
            "selection_tier",
            "unit",
            "value",
            "statistic",
            "assay",
            "candidate_ceiling",
            "negative_ihc_grade",
        ]
        d = d.sort_values(["value", "tissue", "cell_type"], ascending=[False, True, True])
        self.blocks.append(d.drop_duplicates(["entity_id", "unit", "scope"])[fields])

    def finish(self):
        d = pd.concat(self.blocks, ignore_index=True)
        return d.sort_values("value", ascending=False).drop_duplicates(
            ["entity_id", "unit", "scope"]
        )


def extract(facts):
    CACHE.mkdir(parents=True, exist_ok=True)
    for name in [
        "hpa_rna_consensus",
        "hpa_single_cell_tissue",
        "hpa_immune_cell",
        "hpa_immune_monaco",
    ]:
        source = replay.OUTPUT / f"{name}-accepted-identities.csv.gz"
        INPUTS[str(source.relative_to(ROOT))] = sha(source)
        shutil.copyfile(source, CACHE / source.name)
    collector = PeakCollector(facts)
    assessed = replay.replay_contexts(facts, CACHE, collector=collector)
    write(assessed, CACHE / "screen-peaks.csv.gz")
    write(pd.concat(collector.gaps, ignore_index=True), CACHE / "incomplete-normal-contexts.csv.gz")
    # The cancer comparison contains additional native-TPM HPA sources. Reuse
    # its identity/aggregation contract for the entire inventory, reducing each
    # source within anatomical tissue after member sums, before concatenation.
    shutil.copyfile(
        BASE / "rebuild/downloaded-hpa-sources.json", CACHE / "downloaded-hpa-sources.json"
    )
    original = rebuild.group_normal_rows
    duplicate_rows = []

    def reduced(*args, **kwargs):
        raw, selected = args
        gene_col = kwargs.get("gene_col", "gene_id")
        keys = [gene_col, "tissue", "source_context"]
        raw = raw.copy()
        duplicate = raw[gene_col].notna() & raw.duplicated(keys, keep=False)
        if duplicate.any():
            repeated = raw.loc[duplicate].copy()
            conflicts = repeated.groupby(keys, dropna=False).value.transform(
                lambda s: s.nunique(dropna=False) > 1
            )
            repeated["resolution"] = conflicts.map(
                {True: "quarantined_conflicting_vectors", False: "identical_replicate_collapsed"}
            )
            repeated["source"] = kwargs["source"]
            duplicate_rows.append(repeated[[*keys, "value", "resolution", "source"]])
            raw = raw.drop(index=repeated.index[conflicts]).drop_duplicates(keys)
        frame = original(raw, selected, **kwargs)
        return frame.sort_values("value", ascending=False).drop_duplicates(["entity_id", "tissue"])

    original_out = rebuild.OUT
    try:
        rebuild.OUT = CACHE
        rebuild.group_normal_rows = reduced
        normal, _ = rebuild.normal_panel(facts)
    finally:
        rebuild.OUT = original_out
        rebuild.group_normal_rows = original
    normal = normal[
        normal.scope.isin(["normal", "disposable"])
        & normal.qc_pass
        & normal.identity_resolved
        & np.isfinite(normal.value)
    ]
    normal = normal.sort_values("value", ascending=False).drop_duplicates(["entity_id", "scope"])
    write(normal, CACHE / "comparison-peaks.csv.gz")
    write(
        pd.concat(duplicate_rows, ignore_index=True)
        if duplicate_rows
        else pd.DataFrame(
            columns=["gene_id", "tissue", "source_context", "value", "resolution", "source"]
        ),
        CACHE / "duplicate-normal-observations.csv.gz",
    )
    INPUTS.update(rebuild.INPUTS)
    receipt = {
        "inputs": INPUTS,
        "script_sha256": sha(Path(__file__)),
        "outputs": {
            p.name: sha(p)
            for p in [
                CACHE / "screen-peaks.csv.gz",
                CACHE / "comparison-peaks.csv.gz",
                CACHE / "duplicate-normal-observations.csv.gz",
                CACHE / "incomplete-normal-contexts.csv.gz",
            ]
        },
    }
    (CACHE / "receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print("Normal peak extraction complete", flush=True)


def build(reuse=False):
    OUT.mkdir(parents=True, exist_ok=True)
    facts = read(BASE / "rebuild/protein-decisions.csv.gz")
    if not reuse:
        extract(facts)
    receipt = json.loads((CACHE / "receipt.json").read_text())
    for path, expected in receipt["inputs"].items():
        assert sha(ROOT / path) == expected, path
    for path, expected in receipt["outputs"].items():
        assert sha(CACHE / path) == expected, path
    # Rendering is implemented below; extraction remains independently cached.
    print(
        {
            "screen_peaks": len(read(CACHE / "screen-peaks.csv.gz")),
            "comparison_peaks": len(read(CACHE / "comparison-peaks.csv.gz")),
        },
        flush=True,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reuse-peaks", action="store_true")
    build(parser.parse_args().reuse_peaks)
