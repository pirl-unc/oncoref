#!/usr/bin/env python3
"""Screen every gene for rare expression in vital-tissue cell types (CELLxGENE).

Bulk and cluster panels average away rare normal cells, so a gene can read near
zero everywhere while a small population of an essential tissue expresses it
(MAGEA12 in 1.5% of midbrain dopaminergic neurons). This screen records, per
gene, the vital-tissue cell types with the highest fraction of expressing cells.

Scope, set by policy
- Tissues: brain, spinal cord and heart. In lung, alveolar epithelium, endothelium
  and airway basal cells are Vital; differentiated airway epithelium (ciliated, club,
  goblet, secretory, ionocytes) is Manageable and recorded but never flagged
  (oncoref.cta_cell_tiers, shared with the RNA cascade).
- Ignored: immune lineages (including microglia), fibroblasts and stroma, and
  clearly developmental cell types, plus broad parent terms. Rare expression in
  spleen, adipose, immune subsets and similar tissues is allowed.
- A cell type needs at least 3,000 cells and 10 expressing cells; the filter fires at
  1% or more of its cells.

The service cannot restrict samples to adult donors, so developmental cell types are
removed by name only. Values are detection frequencies among cells, not abundance:
expressing cells typically carry one UMI, and ambient RNA, doublets or paralog
multimapping can create detections.

Two steps so rules can change without a second pull:
  collect   query the service for all genes and cache raw rows in --work
  build     apply the rules to the cache and write the packaged table
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
import time
import urllib.request
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from oncoref.cta_cell_tiers import DIFFERENTIATED_AIRWAY, differentiated_airway  # noqa: E402

CLASSIFICATIONS = ROOT / "docs/audits/cta-selection-20261028/protein-group-classifications.csv.gz"
OUTPUT = ROOT / "oncoref/data/cta-rare-cell-screen.csv.gz"
PROVENANCE = ROOT / "oncoref/data/cta-rare-cell-screen.json"
API = "https://api.cellxgene.cziscience.com/wmg/v2"
NORMAL, HUMAN = "PATO:0000461", "NCBITaxon:9606"
CACHED_TISSUES = frozenset({"brain", "spinal cord", "heart", "lung"})
VITAL_TISSUES = frozenset({"brain", "spinal cord", "heart"})
# Lung: the gas-exchange barrier (alveolar epithelium, endothelium) and airway basal
# cells (the stem pool) are Vital; differentiated airway epithelium is Manageable.
# The tier table is shared with the RNA cascade (oncoref.cta_cell_tiers).
LUNG_VITAL = re.compile(r"alveolar|pneumocyte|aerocyte|endothelial|basal", re.I)
EXCLUDED_CELLS = re.compile(
    r"macrophage|monocyte|microglia|\bt cell|\bb cell|\bnk cell|natural killer|lymphocyte|"
    r"leukocyte|dendritic|mast cell|neutrophil|plasma cell|myeloid|erythro|platelet|"
    r"hematopoietic|immune|fibroblast|mesenchymal|stromal|adipocyte|"
    r"neuroblast|precursor|progenitor|radial glia|immature|fetal|embryonic|primitive|"
    r"cycling|proliferating|germinal|ganglionic eminence|malignant|neoplastic",
    re.I,
)
GENERIC = frozenset(
    {
        "cell",
        "native cell",
        "animal cell",
        "eukaryotic cell",
        "somatic cell",
        "nucleate cell",
        "non-terminally differentiated cell",
        "epithelial cell",
        "stromal cell",
        "connective tissue cell",
        "secretory cell",
        "neural cell",
        "mononuclear cell",
        "electrically responsive cell",
        "contractile cell",
        "motile cell",
        "supporting cell",
        "stem cell",
        "neuron",
        "glial cell",
        "central nervous system neuron",
        "endothelial cell",
        "muscle cell",
        "cardiac muscle cell",
        "vascular associated smooth muscle cell",
    }
)
MIN_CELLS, MIN_EXPRESSING, FLAG_PCT, MAX_SHARE = 3000, 10, 1.0, 0.10
COLUMNS = [
    "gene_id", "tissue", "cell_type", "n_expressing", "n_cells", "pct_expressing",
    "mean_expr_log", "tier", "flagged", "review_status",
]  # fmt: skip


def post(path: str, body: dict) -> dict:
    request = urllib.request.Request(
        f"{API}/{path}", json.dumps(body).encode(), {"Content-Type": "application/json"}
    )
    for attempt in range(4):
        try:
            return json.load(urllib.request.urlopen(request, timeout=600))
        except Exception:  # transient service errors are retried
            if attempt == 3:
                raise
            time.sleep(10 * (attempt + 1))
    raise RuntimeError("unreachable")


def universe_genes(path: Path = CLASSIFICATIONS) -> list[str]:
    groups = pd.read_csv(path, low_memory=False)
    return sorted({g for ids in groups.member_gene_ids for g in json.loads(ids)})


def collect(
    work: Path,
    genes: list[str],
    batch: int = 8,
    shard: tuple[int, int] = (0, 1),
    also: tuple[Path, ...] = (),
) -> None:
    """Cache raw rows for the four tissues; resumable, shardable across workers.

    ``also`` lists other cache directories whose genes are already done. The remaining
    genes are split into ``shard[1]`` interleaved shards and this worker takes shard
    ``shard[0]``; start every worker from the same ``also`` state.
    """
    work.mkdir(parents=True, exist_ok=True)
    dims = json.load(urllib.request.urlopen(f"{API}/primary_filter_dimensions", timeout=60))
    supported = {k for g in dims["gene_terms"][HUMAN] for k in g}
    tissues = {k: v for t in dims["tissue_terms"][HUMAN] for k, v in t.items()}
    state = work / "done.json"
    done = set(json.loads(state.read_text())) if state.exists() else set()
    elsewhere = {
        g
        for d in also
        if (d / "done.json").exists()
        for g in json.loads((d / "done.json").read_text())
    }
    remaining = [g for g in genes if g in supported and g not in elsewhere and g not in done]
    todo = remaining[shard[0] :: shard[1]]
    print(f"{len(genes)} genes, {len(todo)} to query in this shard", flush=True)
    for start in range(0, len(todo), batch):
        chunk = todo[start : start + batch]
        reply = post(
            "query",
            {
                "filter": {
                    "gene_ontology_term_ids": chunk,
                    "organism_ontology_term_id": HUMAN,
                    "disease_ontology_term_ids": [NORMAL],
                }
            },
        )
        labels, rows = reply["term_id_labels"]["cell_types"], []
        for gene, by_tissue in reply["expression_summary"].items():
            for term, by_cell in by_tissue.items():
                tissue = tissues.get(term, term)
                if tissue not in CACHED_TISSUES:
                    continue
                whole = labels.get(term, {}).get("CL:0000000", {}).get("aggregated", {})
                for cell, value in by_cell.items():
                    summary = value.get("aggregated", {})
                    info = labels.get(term, {}).get(cell, {}).get("aggregated", {})
                    if summary.get("n", 0) >= 1:
                        rows.append(
                            (
                                gene,
                                tissue,
                                info.get("name", cell),
                                int(summary["n"]),
                                info.get("total_count"),
                                whole.get("total_count"),
                                summary.get("me"),
                            )
                        )
        pd.DataFrame(
            rows,
            columns=[
                "gene_id",
                "tissue",
                "cell_type",
                "n_expressing",
                "n_cells",
                "tissue_cells",
                "mean_expr_log",
            ],
        ).to_parquet(work / f"part-{start:06d}-{chunk[0]}.parquet")
        done.update(chunk)
        state.write_text(json.dumps(sorted(done)))
        if (start // batch) % 100 == 0:
            print(f"{len(done)} genes", flush=True)
    print("COLLECTED", flush=True)


def eligible(raw: pd.DataFrame) -> pd.DataFrame:
    """Apply the vital-tissue and cell-type rules to cached raw rows."""
    d = raw.copy()
    airway = differentiated_airway(d.tissue, d.cell_type)
    in_tissue = d.tissue.isin(VITAL_TISSUES) | (
        d.tissue.eq("lung") & (d.cell_type.str.contains(LUNG_VITAL) | airway)
    )
    d = d.assign(tier=np.where(airway, "manageable", "vital"))
    d = d[in_tissue & ~d.cell_type.isin(GENERIC) & ~d.cell_type.str.contains(EXCLUDED_CELLS)]
    d = d[
        d.n_cells.ge(MIN_CELLS)
        & d.n_expressing.ge(MIN_EXPRESSING)
        & d.n_cells.le(MAX_SHARE * d.tissue_cells)
    ]
    d["pct_expressing"] = 100 * d.n_expressing / d.n_cells
    return d


def build(works: list[Path]) -> pd.DataFrame:
    dims = json.load(urllib.request.urlopen(f"{API}/primary_filter_dimensions", timeout=60))
    supported = {k for g in dims["gene_terms"][HUMAN] for k in g}
    meta = {
        "snapshot_id": dims.get("snapshot_id"),
        "unsupported": sorted(set(universe_genes()) - supported),
    }
    parts = sorted(p for work in works for p in work.glob("part-*.parquet"))
    raw = pd.concat(pd.read_parquet(p) for p in parts)
    table = eligible(raw).sort_values(["gene_id", "pct_expressing"], ascending=[True, False])
    table = table.groupby(["gene_id", "tier"]).head(3).reset_index(drop=True)
    table["pct_expressing"] = table.pct_expressing.round(4)
    table["mean_expr_log"] = table.mean_expr_log.round(4)
    table["flagged"] = table.pct_expressing.ge(FLAG_PCT) & table.tier.eq("vital")
    table["review_status"] = "screen_unreviewed"
    table = table[COLUMNS]
    table.to_csv(OUTPUT, index=False)
    queried = sorted({g for w in works for g in json.loads((w / "done.json").read_text())})
    flagged = table[table.flagged].gene_id.nunique()
    PROVENANCE.write_text(
        json.dumps(
            {
                "schema": "oncoref.cta-rare-cell-screen",
                "source": "CZ CELLxGENE Discover where-is-my-gene summaries",
                "endpoint": f"{API}/query",
                "filters": {"organism": HUMAN, "disease": "normal (PATO:0000461)"},
                "queried": date.today().isoformat(),
                "snapshot_id": meta.get("snapshot_id"),
                "rules": {
                    "vital_tissues": sorted(VITAL_TISSUES),
                    "lung_vital": LUNG_VITAL.pattern,
                    "lung_manageable_recorded_not_flagged": DIFFERENTIATED_AIRWAY.pattern,
                    "excluded_cell_types": EXCLUDED_CELLS.pattern,
                    "excluded_generic_terms": sorted(GENERIC),
                    "min_cells": MIN_CELLS,
                    "min_expressing_cells": MIN_EXPRESSING,
                    "flag_pct": FLAG_PCT,
                    "max_tissue_share": MAX_SHARE,
                },
                "screened_genes": len(queried),
                "flagged_genes": flagged,
                "unsupported_genes": meta.get("unsupported", []),
                "screened_gene_ids": queried,
                "limitations": [
                    "Detection frequency among cells, not abundance; most expressing cells carry about one UMI",
                    "Samples cannot be restricted to adults; developmental cell types are removed by name only",
                    "Ambient RNA, doublets and multimapping among close paralogs can create detections",
                    "Automated screen; review status stays screen_unreviewed",
                ],
                "sha256": hashlib.sha256(OUTPUT.read_bytes()).hexdigest(),
            },
            indent=2,
        )
        + "\n"
    )
    return table


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("work", type=Path, help="cache directory for raw rows")
    parser.add_argument(
        "--also",
        type=Path,
        action="append",
        default=[],
        help="other cache directories to skip and to merge when building",
    )
    parser.add_argument("--shard", default="0/1", help="K/N: take shard K of N remaining genes")
    parser.add_argument("--build-only", action="store_true", help="apply rules to the caches")
    args = parser.parse_args()
    if not args.build_only:
        k, n = map(int, args.shard.split("/"))
        collect(args.work, universe_genes(), shard=(k, n), also=tuple(args.also))
    else:
        table = build([args.work, *args.also])
        print(len(table), "rows;", table[table.flagged].gene_id.nunique(), "flagged genes")
