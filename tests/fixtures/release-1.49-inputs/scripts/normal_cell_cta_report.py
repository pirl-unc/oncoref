#!/usr/bin/env python3
"""Reproduce the added normal-cell evidence screen from checksum-pinned sources."""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from oncoref import normal_cell_data as data  # noqa: E402
from oncoref import normal_cell_screen as screen  # noqa: E402
from oncoref.cta_regen import cross_reactive_ihc  # noqa: E402
from oncoref.gene_ids import resolve_symbol  # noqa: E402

BASE = ROOT / "docs/audits/target-sets-20261023"
OUT = ROOT / "docs/audits/cta-normal-cells-20261023"
BENCHMARK = [
    "XAGE1A",
    "XAGE1B",
    "CTAG1A",
    "CTAG1B",
    "CTAG2",
    "PRAME",
    "MAGEA4",
    "CT55",
    "FOXR2",
    "ABCB5",
    "IL5",
    "IL13",
    "PAPPA",
    "ACRBP",
    "SAGE1",
    "CTCFL",
    "GAGE1",
    "SSX4",
    "SSX5",
    "HORMAD1",
    "TDRD1",
    "TEX101",
    "BRDT",
    "DAZL",
    "TRIP13",
    "GAGE12F",
    "PAGE2B",
    "CTAGE1",
]


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def write(frame, path):
    frame.to_csv(
        path,
        index=False,
        compression={"method": "gzip", "mtime": 0} if path.suffix == ".gz" else None,
    )


def collect(out):
    baseline = pd.read_csv(BASE / "gene-results.csv.gz", keep_default_na=False)
    pinned = json.loads((BASE / "run-manifest.json").read_text())
    if digest(BASE / "gene-results.csv.gz") != pinned["outputs"]["gene-results.csv.gz"]:
        raise ValueError("Baseline hash changed")
    clinical = pd.read_csv(
        ROOT / "docs/audits/cta-screen-20261023/clinical-results.csv", keep_default_na=False
    )
    benchmark_names = set(BENCHMARK) | set(clinical.Gene)
    selected = baseline.loc[~baseline.Category.eq("Unsafe") | baseline.Gene.isin(benchmark_names)]
    ids = set(selected["Gene ID"])
    out.mkdir(parents=True, exist_ok=True)
    print(
        f"Reviewing {len(ids)} genes; all previous CTA/Neither plus the clinical benchmark",
        flush=True,
    )
    summary = pd.DataFrame({"Gene ID": sorted(ids)})

    def flag(frame, boolean, textcol, describe):
        flagged = frame.loc[frame.review]
        summary[boolean] = summary["Gene ID"].isin(flagged.gene_id)
        descriptions = {
            gid: "; ".join(describe(row) for row in g.itertuples(index=False))
            for gid, g in flagged.groupby("gene_id")
        }
        summary[textcol] = summary["Gene ID"].map(descriptions).fillna("")

    cluster = screen.cluster_findings(data.observations("hpa_single_cell_tissue", genes=ids))
    write(cluster, out / "cell-rna.csv.gz")
    flag(
        cluster,
        "Cell RNA review",
        "Cell RNA findings",
        lambda r: (
            f"{r.tissue}/{r.cell_type} cluster {r.cluster}: {r.value:g} nTPM, {r.read_count:g} reads; reproductive-cell max {r.reproductive_cell_max:g}"
        ),
    )
    summary["HPA cell RNA measured"] = summary["Gene ID"].isin(cluster.gene_id)
    print(f"HPA cluster review: {summary['Cell RNA review'].sum()}", flush=True)

    original = []
    for source in data.sources().source:
        if source.startswith("dice_"):
            original.append(screen.dice_findings(data.observations(source, genes=ids)))
    dice = pd.concat(original, ignore_index=True)
    write(dice, out / "dice-distributions.csv.gz")
    immune = pd.concat(
        [
            screen.immune_findings(data.observations(source, genes=ids))
            for source in ["hpa_immune_cell", "hpa_immune_monaco", "hpa_immune_schmiedel"]
        ],
        ignore_index=True,
    )
    write(immune, out / "immune-rna.csv.gz")
    flag(
        dice,
        "Immune RNA review",
        "Immune RNA findings",
        lambda r: (
            f"DICE {r.cell_type} ({r.state}): {r.high_measurements}/{r.measurements} values ≥5 TPM, max {r.maximum:g}"
        ),
    )
    extra = pd.Series(
        {
            gid: "; ".join(
                f"{r.source} {r.cell_type}: {r.value:g} {r.unit}" for r in g.itertuples()
            )
            for gid, g in immune.loc[immune.review].groupby("gene_id")
        },
        dtype=str,
    )
    summary["Immune RNA review"] |= summary["Gene ID"].isin(extra.index)
    summary["Immune RNA findings"] = (
        summary["Immune RNA findings"] + "; " + summary["Gene ID"].map(extra).fillna("")
    ).str.strip("; ")
    summary["DICE measured"] = summary["Gene ID"].isin(dice.gene_id)
    summary["HPA immune measured"] = summary["Gene ID"].isin(immune.gene_id)
    print(f"Immune review: {summary['Immune RNA review'].sum()}", flush=True)

    # Preserve original symbols and audit all mappings. Do not substitute paralogs.
    cornea = data.cornea_summary()
    symbol_to_id = baseline.loc[baseline.Gene.ne("")].set_index("Gene")["Gene ID"].to_dict()
    mapping = pd.DataFrame({"Source symbol": sorted(cornea.gene_name.unique())})
    mapping["Canonical symbol"] = mapping["Source symbol"].map(resolve_symbol)
    mapping["Gene ID"] = mapping["Canonical symbol"].map(symbol_to_id).fillna("")
    # Many-to-one aliases would require count-level reconciliation, not summing CP10K.
    duplicate = mapping["Gene ID"].ne("") & mapping["Gene ID"].duplicated(keep=False)
    mapping["Mapping"] = "Mapped"
    mapping.loc[mapping["Gene ID"].eq(""), "Mapping"] = "Unmapped"
    mapping.loc[duplicate, "Mapping"] = "Ambiguous"
    write(mapping, out / "corneal-gene-mapping.csv.gz")
    usable = mapping.loc[mapping.Mapping.eq("Mapped")].set_index("Source symbol")["Gene ID"]
    cornea["gene_id"] = cornea.gene_name.map(usable).fillna("")
    cornea = cornea.loc[cornea.gene_id.isin(ids)]
    write(cornea, out / "corneal-rna.csv.gz")
    corneal = screen.cornea_findings(cornea)
    corneal["gene_id"] = corneal.gene_name.map(usable)
    write(corneal, out / "corneal-findings.csv.gz")
    flag(
        corneal,
        "Corneal RNA review",
        "Corneal RNA findings",
        lambda r: (
            f"{r.cell_type}: max {r.maximum:g} mean CP10K, {r.supported_groups} supported sample groups"
        ),
    )
    summary["Corneal measured"] = summary["Gene ID"].isin(cornea.gene_id)
    print(f"Corneal review: {summary['Corneal RNA review'].sum()}", flush=True)

    primary = data.table("hpa_antibody_xml", genes=ids)
    primary = screen.antibody_findings(primary, ignored_genes=cross_reactive_ihc())
    # Raw positive antibody staining cannot overturn a curated negative profile.
    curated = data.table("hpa_normal_tissue", genes=ids)
    curated_genes = set(curated.Gene)
    primary["reproductive_protein"] &= ~primary.gene_id.isin(curated_genes)
    primary["review"] = primary.uncertain_normal | (
        primary.normal_signal & ~primary.gene_id.isin(curated_genes)
    )
    write(primary, out / "primary-antibodies.csv.gz")
    flag(
        primary,
        "Antibody review",
        "Antibody findings",
        lambda r: (
            f"{r.tissue}/{r.cell_type}: {r.level}, {r.Antibody}; conflicting or uncurated primary profile"
        ),
    )
    summary["Additional reproductive protein"] = summary["Gene ID"].isin(
        primary.loc[primary.reproductive_protein, "gene_id"]
    )
    summary["Primary IHC measured"] = summary["Gene ID"].isin(primary.gene_id)
    summary["Cell evidence assessed"] = summary[
        [
            "HPA cell RNA measured",
            "DICE measured",
            "HPA immune measured",
            "Corneal measured",
            "Primary IHC measured",
        ]
    ].any(axis=1)
    write(summary, out / "added-evidence.csv.gz")
    print(f"Antibody review: {summary['Antibody review'].sum()}", flush=True)


def classify(out):
    baseline = pd.read_csv(BASE / "gene-results.csv.gz", keep_default_na=False)
    summary = pd.read_csv(out / "added-evidence.csv.gz", keep_default_na=False)
    pubs = pd.read_csv(
        ROOT / "oncoref/data/normal-cell-publication-evidence.csv", keep_default_na=False
    )
    result = screen.classify(baseline, summary, pubs)
    write(result, out / "gene-results.csv.gz")
    changes = result.loc[result.Change.ne("Unchanged")]
    write(changes, out / "changed-genes.csv")
    nominated = result.loc[result["Previous CTA screen"].ne("Not nominated")]
    write(nominated, out / "cta-candidates.csv.gz")
    for category in ["CTA", "Unsafe", "Neither"]:
        write(result.loc[result.Category.eq(category)], out / f"{category.lower()}-genes.csv.gz")
    clinical = pd.read_csv(
        ROOT / "docs/audits/cta-screen-20261023/clinical-results.csv", keep_default_na=False
    )
    benchmark = result.loc[result.Gene.isin(set(BENCHMARK) | set(clinical.Gene))]
    write(benchmark, out / "important-ctas.csv")
    print(result.groupby(["Category", "Detail"]).size().to_string(), flush=True)
    print(
        changes[["Gene", "Previous category", "Category", "Detail"]].to_string(index=False),
        flush=True,
    )


def protein_assets(out, assets, previous):
    from target_sets_report import proteins

    genes = pd.read_csv(out / "gene-results.csv.gz", keep_default_na=False, low_memory=False)
    assets.mkdir(parents=True, exist_ok=True)
    proteins(assets, previous, genes)
    path = assets / "exclusion-basis.json"
    recipe = json.loads(path.read_text())
    recipe["policy"] = screen.POLICY_ID
    path.write_text(json.dumps(recipe, indent=2) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUT)
    parser.add_argument("--reuse-evidence", action="store_true")
    parser.add_argument(
        "--assets", type=Path, help="Regenerate the complete protein artifacts here"
    )
    parser.add_argument(
        "--previous-sequences",
        type=Path,
        help="Pinned Ensembl 115 sequence assets from the 21 October report",
    )
    args = parser.parse_args()
    if bool(args.assets) != bool(args.previous_sequences):
        parser.error("--assets and --previous-sequences must be supplied together")
    if not args.reuse_evidence:
        collect(args.output)
    classify(args.output)
    if args.assets:
        protein_assets(args.output, args.assets, args.previous_sequences)
