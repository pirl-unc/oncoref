#!/usr/bin/env python3
"""Audit every v1.17 CTA/Neither gene for high normal RNA and record follow-up."""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from normal_cell_cta_report import digest, write  # noqa: E402

from oncoref import normal_cell_data as data  # noqa: E402
from oncoref import normal_cell_screen as screen  # noqa: E402
from oncoref import normal_rna_review as review  # noqa: E402
from oncoref import reference_data  # noqa: E402

BASE = ROOT / "docs/audits/cta-normal-cells-20261023"
OUT = ROOT / "docs/audits/cta-neither-review-20261024"
SOURCES = [
    "hpa_rna_consensus",
    "hpa_single_cell",
    "hpa_single_cell_tissue",
    "hpa_immune_cell",
    "hpa_immune_cell_samples",
    "hpa_immune_monaco",
    "hpa_immune_schmiedel",
    *[s for s in reference_data._normal_manifest["artifacts"] if s.startswith("dice_")],
]


def load_base():
    manifest = json.loads((BASE / "run-manifest.json").read_text())
    if digest(BASE / "gene-results.csv.gz") != manifest["outputs"]["gene-results.csv.gz"]:
        raise ValueError("Frozen v1.17 decisions changed")
    return pd.read_csv(BASE / "gene-results.csv.gz", keep_default_na=False, low_memory=False)


def collect(out):
    base = load_base()
    important = pd.read_csv(BASE / "important-ctas.csv", keep_default_na=False)
    ids = set(base.loc[base.Category.ne("Unsafe"), "Gene ID"]) | set(important["Gene ID"])
    # Complete and incomplete reference translations include loci outside the
    # canonical protein-coding table (e.g. immune-receptor segments). They must
    # not evade the same ceiling by defaulting to an unassessed Neither protein.
    sequence_path = BASE / "unclassified-sequence-genes.csv.gz"
    manifest = json.loads((BASE / "run-manifest.json").read_text())
    if digest(sequence_path) != manifest["outputs"][sequence_path.name]:
        raise ValueError("Frozen sequence-locus inventory changed")
    sequence_ids = set(pd.read_csv(sequence_path)["Gene ID"])
    ids |= sequence_ids
    flags, maxima, provenance = [], [], []
    seen = set()

    def ingest(observations, label):
        d = review.findings(screen.contexts(observations))
        seen.update(d.gene_id)
        flags.append(d.loc[d.exclude].copy())
        # No pooling of native units, sources or anatomical contexts.
        keys = [
            "gene_id",
            "source",
            "version",
            "unit",
            "tissue",
            "cell_type",
            "category",
            "mapping_status",
            "permitted",
        ]
        grouped = d.groupby(keys, dropna=False, as_index=False).agg(
            measurements=("value", "size"),
            minimum=("value", "min"),
            maximum=("value", "max"),
            above_ceiling=("above_ceiling", "sum"),
        )
        maxima.append(grouped)
        print(
            f"{label}: {len(d):,} observations; {d.loc[d.exclude, 'gene_id'].nunique()} high genes",
            flush=True,
        )

    for source in SOURCES:
        provenance.append(reference_data.provenance(source, verify_content=True))
        if provenance[-1]["verification_state"] != "verified":
            raise ValueError(f"Unverified source {source}")
        observations = data.observations(source, genes=ids)
        ingest(observations, source)
        if source == "hpa_immune_cell_samples":
            # The native HPA sample table contains BOTH units. The standard
            # adapter selects nTPM; also audit its original TPM, without conversion.
            native = data.table(source, genes=ids)
            if len(native) != len(observations):
                raise ValueError("HPA sample RNA row alignment changed")
            tpm = observations.copy()
            tpm["unit"], tpm["value"] = "TPM", native.TPM.to_numpy()
            ingest(tpm, source + " native TPM")

    donors = pd.read_csv(ROOT / "oncoref/data/cta-donor-rna.csv", keep_default_na=False)
    d = donors.rename(
        columns={
            "Ensembl_Gene_ID": "gene_id",
            "Symbol": "gene_name",
            "rna_ntpm": "value",
            "source_version": "version",
            "donor_id": "donor",
            "sample_id": "sample",
        }
    ).copy()
    d = d.loc[d.gene_id.isin(ids)]
    d["source"], d["unit"], d["modality"], d["stage"] = "hpa_rna_consensus", "nTPM", "rna", "adult"
    d["cell_type"] = ""
    ingest(d, "reviewed vital-tissue donor RNA")

    high = pd.concat(flags, ignore_index=True)
    # Anonymous DICE replicate slots are retained; never fabricate donor IDs.
    high["finding_id"] = [
        hashlib.sha256(json.dumps(row, sort_keys=True, default=str).encode()).hexdigest()[:20]
        for row in high.fillna("").to_dict("records")
    ]
    if high.finding_id.duplicated().any():
        raise ValueError("Duplicate evidence observations")
    maximum = pd.concat(maxima, ignore_index=True)
    write(high, out / "high-normal-rna.csv.gz")
    write(maximum, out / "rna-context-maxima.csv.gz")
    nonpermitted = maximum.loc[~maximum.permitted]
    summary = pd.DataFrame({"Gene ID": sorted(ids)})
    summary["Normal RNA assessed"] = True
    summary["Normal RNA measured"] = summary["Gene ID"].isin(seen)
    summary["High normal RNA"] = summary["Gene ID"].isin(high.gene_id)
    for unit in ["TPM", "nTPM"]:
        m = nonpermitted.loc[nonpermitted.unit.eq(unit)].groupby("gene_id").maximum.max()
        summary[f"Normal RNA maximum ({unit})"] = summary["Gene ID"].map(m)
    descriptions = {}
    for gene, rows in high.groupby("gene_id"):
        groups = rows.groupby(["source", "unit", "tissue", "cell_type"], dropna=False)
        descriptions[gene] = "; ".join(
            f"{source}: {tissue}/{cell}, max {g.value.max():g} {unit}, {len(g)} observations >30"
            for (source, unit, tissue, cell), g in groups
        )
    summary["High normal RNA findings"] = summary["Gene ID"].map(descriptions).fillna("")
    write(summary, out / "normal-rna-summary.csv.gz")
    (out / "source-provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")


def classify(out):
    base = load_base()
    summary = pd.read_csv(out / "normal-rna-summary.csv.gz", keep_default_na=False)
    result = review.classify(base, summary)
    write(result, out / "gene-results.csv.gz")
    write(result.loc[result.Change.ne("Unchanged")], out / "changed-genes.csv")
    for category in ["CTA", "Unsafe", "Neither"]:
        write(result.loc[result.Category.eq(category)], out / f"{category.lower()}-genes.csv.gz")
    important = pd.read_csv(BASE / "important-ctas.csv", keep_default_na=False)
    write(result.loc[result["Gene ID"].isin(important["Gene ID"])], out / "important-ctas.csv")
    extra = summary.loc[~summary["Gene ID"].isin(base["Gene ID"])].copy()
    extra["Category"] = extra["High normal RNA"].map({True: "Unsafe", False: "Neither"})
    extra["Detail"] = extra["High normal RNA"].map({True: "High RNA", False: "Unassessed"})
    extra["Gene"] = ""
    extra["Target"], extra["Forbidden"] = False, extra["High normal RNA"]
    extra["Policy"] = review.POLICY_ID
    extra["Reason"] = extra["High normal RNA"].map(
        {
            True: "Exact reference gene ID has RNA >30 TPM/nTPM outside permitted contexts; excluded despite absence from the canonical gene table",
            False: "Reference sequence locus lacks a canonical gene classification; no promotion or inferred absence",
        }
    )
    write(extra, out / "sequence-gene-reviews.csv.gz")
    notes = pd.read_csv(ROOT / "oncoref/data/cta-followup-notes.csv", keep_default_na=False)
    if notes.Gene.duplicated().any() or not set(notes.Gene) <= set(result.Gene):
        raise ValueError("Follow-up notes require unique, present genes")
    queue = notes.merge(
        result.loc[result.Gene.isin(notes.Gene)], on="Gene", how="left", validate="one_to_one"
    )
    queue["Review date"] = "2026-10-24"
    queue["Evidence report"] = (
        "https://pirl-unc.github.io/oncoref/audits/cta-neither-review-20261024/"
    )
    queue["Previous report"] = (
        "https://pirl-unc.github.io/oncoref/audits/cta-normal-cells-20261023/"
    )
    write(queue, out / "follow-up.csv")
    # Auditable invariant across all measured non-permitted TPM/nTPM contexts.
    high = pd.read_csv(out / "high-normal-rna.csv.gz", keep_default_na=False)
    remaining = (
        set(result.loc[result.Category.ne("Unsafe"), "Gene ID"])
        | set(extra.loc[extra.Category.ne("Unsafe"), "Gene ID"])
    ) & set(high.gene_id)
    if remaining:
        raise AssertionError(f"High normal RNA remains outside forbidden: {remaining}")
    print(result.groupby(["Category", "Detail"]).size().to_string(), flush=True)
    print(
        result.loc[
            result.Change.ne("Unchanged"), ["Gene", "Previous category", "Category"]
        ].to_string(index=False),
        flush=True,
    )
    return result


def protein_assets(out, assets, previous):
    from target_sets_report import proteins

    assets.mkdir(parents=True, exist_ok=True)
    genes = pd.read_csv(out / "gene-results.csv.gz", keep_default_na=False, low_memory=False)
    extra = pd.read_csv(out / "sequence-gene-reviews.csv.gz", keep_default_na=False)
    genes = pd.concat([genes, extra], ignore_index=True)
    proteins(assets, previous, genes)
    # Keep the noncanonical identity coverage visible even after assigning RNA
    # exclusions; a decision does not resolve a locus into a canonical gene.
    write(extra, assets / "unclassified-sequence-genes.csv.gz")
    path = assets / "exclusion-basis.json"
    recipe = json.loads(path.read_text())
    recipe["policy"] = review.POLICY_ID
    path.write_text(json.dumps(recipe, indent=2) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUT)
    parser.add_argument("--reuse-evidence", action="store_true")
    parser.add_argument("--assets", type=Path)
    parser.add_argument("--previous-sequences", type=Path)
    args = parser.parse_args()
    if bool(args.assets) != bool(args.previous_sequences):
        parser.error("--assets and --previous-sequences must be supplied together")
    args.output.mkdir(parents=True, exist_ok=True)
    if not args.reuse_evidence:
        collect(args.output)
    classify(args.output)
    if args.assets:
        protein_assets(args.output, args.assets, args.previous_sequences)
