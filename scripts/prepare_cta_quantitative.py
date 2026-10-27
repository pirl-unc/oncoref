#!/usr/bin/env python3
"""Prepare genome-wide measured evidence for the experimental quantitative audit.

Requires SciPy for the Matrix Market reader. Source files are downloaded separately
from the pinned URLs in developmental-source.json; --sources names that directory.
No membership, control or publication columns participate in these summaries.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.io import mmread

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from oncoref import (  # noqa: E402
    cta_membership,
    hpa,
    normal_cell_data,
    normal_tissues,
    reference_data,
)


def prepare_development(source_dir, out):
    receipt = json.loads(
        (ROOT / "docs/audits/cta-quantitative-20261026/developmental-source.json").read_text()
    )
    for name, record in receipt["files"].items():
        if hashlib.sha256((source_dir / name).read_bytes()).hexdigest() != record["sha256"]:
            raise ValueError(f"Developmental source checksum mismatch: {name}")
    b = cta_membership.results()
    base = b[["Gene ID", "Gene"]]
    with zipfile.ZipFile(source_dir / "atlas_tpm") as z:
        print("TPM files", z.namelist(), flush=True)
        name = next(n for n in z.namelist() if n.endswith(".mtx"))
        matrix = mmread(z.open(name)).tocsr()
        genes = pd.read_csv(z.open(name + "_rows"), sep="\t", header=None)[0].astype(str)
        cells = pd.read_csv(z.open(name + "_cols"), sep="\t", header=None)[0].astype(str)
        print(
            "matrix",
            matrix.shape,
            "column sums",
            np.quantile(np.asarray(matrix.sum(axis=0)).ravel(), [0, 0.5, 1]),
            flush=True,
        )
    meta = pd.read_csv(source_dir / "atlas_design", sep="\t").set_index("Assay").loc[cells]
    assert cells.is_unique and genes.is_unique
    label = "Sample Characteristic[inferred cell type - authors labels]"
    print(meta[label].value_counts().to_string(), flush=True)
    frames = []
    for group, ids in meta.groupby(label).indices.items():
        a = matrix[:, ids]
        n = len(ids)
        count = np.asarray((a >= 1).sum(axis=1)).ravel()
        frames.append(
            pd.DataFrame(
                {
                    "Gene ID": genes,
                    "Lineage": group,
                    "Cells": n,
                    "Positive cells TPM>=1": count,
                    "Fraction TPM>=1": count / n,
                    "Mean TPM": np.asarray(a.mean(axis=1)).ravel(),
                    "Embryos": 1,
                }
            )
        )
    d = pd.concat(frames, ignore_index=True)
    d = d[d["Gene ID"].isin(base["Gene ID"])]
    d.to_csv(
        out / "developmental-lineages.csv.gz",
        index=False,
        compression={"method": "gzip", "mtime": 0},
    )
    f = d[(d.Cells >= 20) & (d["Positive cells TPM>=1"] >= 5) & (d["Fraction TPM>=1"] >= 0.1)]
    f = f.sort_values(["Mean TPM", "Gene ID"], ascending=[False, True]).drop_duplicates("Gene ID")

    print(
        "qualifiable development genes", len(f), "mean>=5", (f["Mean TPM"] >= 5).sum(), flush=True
    )
    f.to_csv(out / "developmental-maxima.csv", index=False)
    # The normal IHC table is already a gene-level curated HPA consensus table.
    ihc = hpa._read_hpa("hpa_normal_tissue", "v23")
    print("IHC", ihc.columns.tolist(), flush=True)
    contexts = ihc[["Tissue", "Cell type"]].drop_duplicates()
    mappings = []
    for tissue, cell in contexts.itertuples(index=False, name=None):
        a = normal_tissues.resolve_normal_tissue(
            tissue,
            cell_type=cell,
            species="human",
            source_id="hpa_normal_tissue",
            source_version="v23",
            stage="adult",
        )
        mappings.append(
            {
                "Tissue": tissue,
                "Cell type": cell,
                "Category": normal_tissues.normal_tissue_category(a),
                "Mapping": a.mapping_status,
            }
        )
    ihc = ihc.merge(
        pd.DataFrame(mappings), on=["Tissue", "Cell type"], how="left", validate="many_to_one"
    )
    ihc = ihc[ihc.Gene.isin(base["Gene ID"])]
    ihc.to_csv(out / "normal-ihc.csv.gz", index=False, compression={"method": "gzip", "mtime": 0})


def prepare_cells(out):
    b = cta_membership.results()[["Gene ID", "Gene"]]
    ids = set(b["Gene ID"])
    result = b[["Gene ID", "Gene"]].set_index("Gene ID")
    receipts = []

    def receipt(source):
        r = reference_data.provenance(source, verify_content=True)
        receipts.append({k: r[k] for k in ["name", "version", "url", "sha256"]})
        return r

    source = "hpa_single_cell_tissue"
    r = receipt(source)
    summaries = []
    maps = {}
    for i, d in enumerate(pd.read_csv(r["path"], sep="\t", chunksize=250000)):
        d = d[d.Gene.isin(ids)].copy()
        for tissue, cell in (
            d[["Tissue", "Cell type"]].drop_duplicates().itertuples(index=False, name=None)
        ):
            if (tissue, cell) not in maps:
                a = normal_tissues.resolve_normal_tissue(
                    tissue,
                    cell_type=cell,
                    species="human",
                    source_id=source,
                    source_version=r["version"],
                    stage="adult",
                )
                maps[(tissue, cell)] = (normal_tissues.normal_tissue_category(a), a.mapping_status)
        d["Category"] = [maps[(t, c)][0] for t, c in zip(d.Tissue, d["Cell type"])]
        d["Mapping"] = [maps[(t, c)][1] for t, c in zip(d.Tissue, d["Cell type"])]
        d = d[~d.Category.isin(["reproductive", "disposable"])].copy()
        d["Supported nTPM"] = d.nTPM.where(d["Read count"].ge(20))
        d["Unresolved nTPM"] = d.nTPM.where(d.Mapping.ne("mapped"))
        summaries.append(
            d.groupby("Gene").agg(
                cell_records=("nTPM", "size"),
                cell_supported_records=("Supported nTPM", "count"),
                cell_max_ntpm=("Supported nTPM", "max"),
                cell_raw_max_ntpm=("nTPM", "max"),
                cell_unresolved_max_ntpm=("Unresolved nTPM", "max"),
            )
        )

        if i % 10 == 0:
            print("cluster chunk", i, flush=True)
    s = pd.concat(summaries)
    s = s.groupby(level=0).agg(
        {
            "cell_records": "sum",
            "cell_supported_records": "sum",
            "cell_max_ntpm": "max",
            "cell_raw_max_ntpm": "max",
            "cell_unresolved_max_ntpm": "max",
        }
    )
    result = result.join(s)
    result.to_csv(out / "cluster-summary.csv")
    print("clusters done", flush=True)
    for source, unit in [("hpa_immune_cell", "nTPM"), ("hpa_immune_monaco", "TPM")]:
        r = receipt(source)
        d = pd.read_csv(r["path"], sep="\t")
        s = d.groupby("Gene")[unit].agg(["max", "count"])
        s.columns = [source + "_max_" + unit, source + "_records"]
        result = result.join(s)
    dice = []
    for source in normal_cell_data.sources().source:
        if not source.startswith("dice_"):
            continue
        r = receipt(source)
        d = pd.read_csv(r["path"])
        d = d[d.Feature_name.str.split(".").str[0].isin(ids)]
        gene = d.Feature_name.str.split(".").str[0]
        values = d[[c for c in d if c.isdigit()]].to_numpy(float)
        assert np.isfinite(values).all() and (values >= 0).all()
        v = np.partition(values, -3, axis=1)[:, -3]
        dice.append(
            pd.DataFrame(
                {
                    "Gene ID": gene,
                    "Source": source,
                    "Third largest TPM": v,
                    "Reported values": values.shape[1],
                }
            )
        )
        print(source, flush=True)
    d = pd.concat(dice)
    d = d[d["Gene ID"].isin(ids)]
    d["PAR identity duplicated"] = d.duplicated(["Gene ID", "Source"], keep=False)
    d.to_csv(out / "dice-summaries.csv.gz", index=False, compression={"method": "gzip", "mtime": 0})
    result = result.join(
        d.groupby("Gene ID").agg(
            dice_third_max_tpm=("Third largest TPM", "max"), dice_sources=("Source", "nunique")
        )
    )
    for source in ["cornea_collin_adult_counts", "cornea_collin_adultcornea_cells"]:
        receipt(source)
    d = normal_cell_data.cornea_summary(genes=ids)
    d = d[d.gene_id.ne("")].copy()
    assert not d.duplicated(["gene_id", "cell_type", "sample"]).any()
    d["supported"] = d.value.where(d.positive_cells.ge(5) & d.read_count.ge(20))
    v = d[d.supported.notna()].sort_values("supported", ascending=False)
    v["rank"] = v.groupby(["gene_id", "cell_type"]).cumcount() + 1
    s = v[v["rank"].eq(2)].groupby("gene_id").supported.max().rename("cornea_second_max_cp10k")
    c = d.groupby("gene_id").size().rename("cornea_records")
    result = result.join(s).join(c)
    print("cornea done", flush=True)
    result.to_csv(out / "adult-cell-summary.csv")
    (out / "cell-source-receipts.json").write_text(json.dumps(receipts, indent=2) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sources", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    prepare_development(args.sources, args.out)
    prepare_cells(args.out)


if __name__ == "__main__":
    main()
