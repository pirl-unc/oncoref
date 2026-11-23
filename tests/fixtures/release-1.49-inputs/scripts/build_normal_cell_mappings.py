#!/usr/bin/env python3
"""Extend explicit source mappings from the pinned cell/antibody vocabularies."""

from __future__ import annotations

import re
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from oncoref import normal_cell_data as data  # noqa: E402
from oncoref import normal_tissues as anatomy  # noqa: E402
from oncoref import reference_data as references  # noqa: E402

MOBILE_IMMUNE = re.compile(
    r"t-cells?|b-cells?|nk-cells?|plasma cells?|macrophages?|dendritic|granulocytes?|"
    r"monocytes?|neutrophils?|eosinophils?|basophils?|platelets?|erythroid|megakaryocytes?|immune cells?|red blood cells?|lymphocytes?",
    re.I,
)


def concept(tissue, cell):
    if MOBILE_IMMUNE.search(cell) and "hofbauer" not in cell.lower():
        return (
            "blood",
            "Mobile hematopoietic lineage; does not inherit reproductive/disposable permission from sampling organ.",
        )
    if tissue == "cornea/conjunctiva":
        if cell.lower() in {"blood vessels", "lymphatic vessels"}:
            return (
                "blood_vessel" if cell.lower() == "blood vessels" else "cornea",
                "Vascular compartment in adult ocular surface; source cell label retained.",
            )
        cid = (
            "conjunctiva"
            if "conjunctiv" in cell.lower()
            else "corneal_limbus"
            if "limbal" in cell.lower()
            else "cornea"
        )
        return (
            cid,
            "Adult ocular surface; published cell label distinguishes cornea, limbus and conjunctiva.",
        )
    tissue = {"pbmc": "blood", "vascular": "blood vessel"}.get(tissue.lower(), tissue)
    a = anatomy.resolve_normal_tissue(tissue, species="human", cell_type=cell)
    if a.mapping_status != "mapped":
        a = anatomy.resolve_normal_tissue(tissue, species="human")
    return (
        a.concept_id or "",
        "Source organ context retained; unspecified cell anatomy is not narrowed.",
    )


def main():
    path = ROOT / "oncoref/data/normal-tissue-mappings.csv"
    existing = pd.read_csv(path, keep_default_na=False)
    managed = set(references._normal_manifest["artifacts"])
    existing = existing.loc[~existing.source_id.isin(managed)]
    rows = []

    def add(source, tissue, cell, cid=None, note=None):
        if cid is None:
            cid, note = concept(tissue, cell)
        spec = references.REFERENCE_SOURCES[source]
        rows.append(
            {
                "taxon_id": 9606,
                "source_id": source,
                "source_version": references.resolve_version(source),
                "source_tissue": tissue,
                "source_cell_type": cell,
                "concept_id": cid,
                "mapping_kind": "organ_context" if cid else "unresolved",
                "source_reference": spec.get(
                    "publication_url", spec["urls"][references.resolve_version(source)]
                ),
                "physiological_state": "",
                "notes": note,
            }
        )

    clusters = data.table("hpa_single_cell_clusters")
    for tissue, cell in (
        clusters[["Tissue", "Cell type"]].drop_duplicates().itertuples(index=False, name=None)
    ):
        add("hpa_single_cell_tissue", tissue, cell)
    primary = pd.read_parquet(
        data._antibody_cache(), columns=["Tissue", "Cell type"]
    ).drop_duplicates()
    for tissue, cell in primary.itertuples(index=False, name=None):
        a = anatomy.resolve_normal_tissue(
            tissue,
            species="human",
            cell_type=cell,
            source_id="hpa_normal_tissue",
            source_version="v23",
        )
        cid, note = concept(tissue, cell)
        if not cid and a.mapping_status == "mapped":
            cid = a.concept_id
        add("hpa_antibody_xml", tissue, cell, cid, note)
    for source in [
        "hpa_immune_cell",
        "hpa_immune_cell_samples",
        "hpa_immune_monaco",
        "hpa_immune_schmiedel",
    ]:
        p = references.ensure(source)
        cells = pd.read_csv(p, sep="\t", usecols=["Immune cell"])["Immune cell"].unique()
        for cell in cells:
            add(
                source,
                "blood",
                cell,
                "blood",
                "Sorted blood immune population; activation state remains separate metadata.",
            )
    for source, spec in references._normal_manifest["artifacts"].items():
        if source.startswith("dice_"):
            add(
                source,
                "blood",
                spec["cell_type"],
                "blood",
                "DICE immune population; anonymous per-gene specimen distribution, no paired donor identity.",
            )
    for cell in data.cornea_cells().cell_type.unique():
        add("cornea_collin_adult_counts", "cornea/conjunctiva", cell)
    new = pd.DataFrame(rows).drop_duplicates(
        ["source_id", "source_version", "source_tissue", "source_cell_type"]
    )
    pd.concat([existing, new], ignore_index=True).to_csv(path, index=False, lineterminator="\r\n")
    print(f"Added {len(new)} explicit mappings; unresolved: {int(new.concept_id.eq('').sum())}")


if __name__ == "__main__":
    main()
