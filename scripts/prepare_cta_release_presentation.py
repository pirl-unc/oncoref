#!/usr/bin/env python3
"""Derive presentation data from the current audit without reclassifying genes."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
AUDIT = ROOT / "docs/audits/cta-release-1.21.0"


def records(frame: pd.DataFrame) -> list[dict]:
    """Encode pandas missing values as JSON nulls, never NaN literals."""
    return json.loads(frame.to_json(orient="records"))


def build() -> None:
    """Bind the deck's counts, paper unions and known-CTA tracking to the audit."""
    group = pd.read_csv(AUDIT / "protein-group-classifications.csv.gz", low_memory=False)
    known = pd.read_csv(AUDIT / "known-cta-decisions.csv.gz", low_memory=False)
    registry = pd.read_csv(ROOT / "oncoref/data/cta-publication-sources.csv")
    membership = known.dropna(subset=["membership_sources"]).copy()
    membership["source_tag"] = membership.membership_sources.str.split(";")
    membership = membership.explode("source_tag").merge(
        registry[["source_tag", "doi"]], on="source_tag", validate="many_to_one"
    )
    papers, selected_sets = [], {}
    for doi, source in registry.groupby("doi"):
        part = membership[membership.doi.eq(doi)].drop_duplicates("protein_equivalence_group_id")
        selected = set(
            part.loc[part.protein_group_category.eq("CTA"), "protein_equivalence_group_id"]
        )
        name = source.iloc[0].source_tag.split("_")[0].replace("daSilva", "da Silva ")
        selected_sets[name] = selected
        papers.append(
            {
                "paper": name,
                "doi": doi,
                "mapped_groups": len(part),
                **{
                    c: int(part.protein_group_category.eq(c).sum())
                    for c in ["CTA", "Unsafe", "Neither"]
                },
            }
        )
    papers.sort(key=lambda row: (-row["mapped_groups"], row["paper"]))
    top = [p["paper"] for p in papers[:6]]
    historical = known[known.historical_CTA].protein_group_category.value_counts().to_dict()
    association = {
        "reproductive_only": int((group.reproductive_route & ~group.developmental_route).sum()),
        "developmental_only": int((~group.reproductive_route & group.developmental_route).sum()),
        "both": int((group.reproductive_route & group.developmental_route).sum()),
    }
    data = {
        "summary": json.loads((AUDIT / "summary.json").read_text()),
        "report_sha256": hashlib.sha256((AUDIT / "manifest.json").read_bytes()).hexdigest(),
        "universe": len(group),
        "papers": papers,
        "paper_union": membership.protein_equivalence_group_id.nunique(),
        "overlap_order": top,
        "overlap": [[len(selected_sets[a] & selected_sets[b]) for b in top] for a in top],
        "historical_loci": historical,
        "association": association,
        **{
            key: records(pd.read_csv(AUDIT / file))
            for key, file in {
                "funnel": "filtering-funnel.csv",
                "components": "policy-component-effects.csv",
                "named": "named-antigen-decisions.csv",
                "normal": "named-normal-evidence.csv",
                "inventory": "species-source-inventory.csv",
            }.items()
        },
    }
    (AUDIT / "presentation-data.json").write_text(json.dumps(data, indent=2) + "\n")
    print(
        json.dumps(
            {
                "papers": len(papers),
                "paper_union": data["paper_union"],
                "association": association,
                "historical_loci": historical,
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    build()
