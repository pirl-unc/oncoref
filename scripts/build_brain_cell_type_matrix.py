#!/usr/bin/env python3
"""Build a gene by brain-cell-type detection matrix from CELLxGENE summaries.

For every gene in the human protein universe this records, in normal samples,
the number of expressing cells and total cells for each brain and nervous-system
cell type, plus each gene's dominant neuron type. Values are detection
frequencies among cells, not abundance. The service returns ontology rollups
(leaf-only queries are rejected), so broad parent terms and any type covering
over 10% of its tissue are excluded from "dominant" calls.

Output is large and is not committed. Resumable: rerun to continue.
"""

from __future__ import annotations

import argparse
import json
import re
import time
import urllib.request
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
API = "https://api.cellxgene.cziscience.com/wmg/v2"
HUMAN, NORMAL = "NCBITaxon:9606", "PATO:0000461"
NERVOUS = {"brain", "spinal cord", "nervous system", "peripheral nervous system", "cortex"}
NEURON = re.compile(
    r"neuron|neuroblast|purkinje|granule|pyramidal|interneuron|medium spiny|"
    r"dopaminergic|serotonergic|motor neuron|cajal|chandelier|martinotti|basket",
    re.I,
)
MIN_CELLS, MIN_EXPRESSING, MAX_SHARE = 3000, 10, 0.10


def post(body: dict) -> dict:
    request = urllib.request.Request(
        f"{API}/query", json.dumps(body).encode(), {"Content-Type": "application/json"}
    )
    for attempt in range(4):
        try:
            return json.load(urllib.request.urlopen(request, timeout=600))
        except Exception:
            if attempt == 3:
                raise
            time.sleep(10 * (attempt + 1))
    raise RuntimeError("unreachable")


def universe_genes() -> list[str]:
    path = ROOT / "docs/audits/cta-selection-20261028/protein-group-classifications.csv.gz"
    groups = pd.read_csv(path, low_memory=False)
    return sorted({g for ids in groups.member_gene_ids for g in json.loads(ids)})


def collect(out: Path, batch: int = 8) -> None:
    dims = json.load(urllib.request.urlopen(f"{API}/primary_filter_dimensions", timeout=60))
    supported = {k for g in dims["gene_terms"][HUMAN] for k in g}
    tissues = {k: v for t in dims["tissue_terms"][HUMAN] for k, v in t.items()}
    genes = [g for g in universe_genes() if g in supported]
    state = out / "done.json"
    done = set(json.loads(state.read_text())) if state.exists() else set()
    parts = out / "parts"
    parts.mkdir(parents=True, exist_ok=True)
    todo = [g for g in genes if g not in done]
    print(f"{len(genes)} genes, {len(todo)} to do", flush=True)
    for start in range(0, len(todo), batch):
        chunk = todo[start : start + batch]
        reply = post(
            {
                "filter": {
                    "gene_ontology_term_ids": chunk,
                    "organism_ontology_term_id": HUMAN,
                    "disease_ontology_term_ids": [NORMAL],
                }
            }
        )
        labels = reply["term_id_labels"]["cell_types"]
        rows = []
        for gene, by_tissue in reply["expression_summary"].items():
            for term, by_cell in by_tissue.items():
                tissue = tissues.get(term, term)
                if tissue not in NERVOUS:
                    continue
                whole = labels.get(term, {}).get("CL:0000000", {}).get("aggregated", {})
                for cell, value in by_cell.items():
                    summary = value.get("aggregated", {})
                    info = labels.get(term, {}).get(cell, {}).get("aggregated", {})
                    n = summary.get("n", 0)
                    if n >= 1:
                        rows.append(
                            (
                                gene, tissue, info.get("name", cell), int(n),
                                info.get("total_count"), whole.get("total_count"),
                                summary.get("me"),
                            )
                        )
        pd.DataFrame(
            rows,
            columns=[
                "gene_id", "tissue", "cell_type", "n_expressing", "n_cells", "tissue_cells",
                "mean_expr_log",
            ],
        ).to_parquet(parts / f"{start:06d}-{chunk[0]}.parquet")
        done.update(chunk)
        state.write_text(json.dumps(sorted(done)))
        if (start // batch) % 50 == 0:
            print(f"{len(done)}/{len(genes)}", flush=True)
    print("COLLECTED", flush=True)


def summarize(out: Path) -> pd.DataFrame:
    long = pd.concat(pd.read_parquet(p) for p in sorted((out / "parts").glob("*.parquet")))
    long["pct_expressing"] = 100 * long.n_expressing / long.n_cells
    long.to_parquet(out / "brain-cell-type-long.parquet")
    neurons = long[
        long.cell_type.str.contains(NEURON)
        & long.n_cells.ge(MIN_CELLS)
        & long.n_expressing.ge(MIN_EXPRESSING)
        & long.n_cells.le(MAX_SHARE * long.tissue_cells)
    ]
    ranked = neurons.sort_values(["gene_id", "pct_expressing"], ascending=[True, False])
    top = ranked.groupby("gene_id").head(3)
    best = top.groupby("gene_id").first()
    summary = pd.DataFrame(
        {
            "dominant_neuron_type": best.cell_type,
            "dominant_tissue": best.tissue,
            "dominant_pct_expressing": best.pct_expressing.round(4),
            "dominant_n_expressing": best.n_expressing,
            "dominant_n_cells": best.n_cells,
            "top3_neuron_types": top.groupby("gene_id").apply(
                lambda d: "; ".join(f"{c} ({p:.2f}%)" for c, p in zip(d.cell_type, d.pct_expressing))
            ),
            "neuron_types_detected": neurons.groupby("gene_id").cell_type.nunique(),
        }
    )
    total = long[(long.tissue == "brain") & (long.cell_type == "cell")].set_index("gene_id")
    summary["brain_expressing_cells"] = total.n_expressing
    summary["brain_cells"] = total.n_cells
    summary.to_csv(out / "gene-dominant-neuron-type.csv.gz")
    return summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--summarize-only", action="store_true")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    if not args.summarize_only:
        collect(args.output)
    print(summarize(args.output).shape, flush=True)
