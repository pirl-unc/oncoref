#!/usr/bin/env python3
"""Export genome-wide Reproductive, Normal and Other expression categories."""

from __future__ import annotations

import argparse
import hashlib
import json
import multiprocessing
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from cta_permitted_report import annotate_archive  # noqa: E402
from cta_scope_report import annotate_measurements, sha256, write_csv  # noqa: E402

from oncoref import cta_regen, cta_sources, hpa, reference_data  # noqa: E402
from oncoref.gene_ids import canonical_gene_space  # noqa: E402
from oncoref.normal_expression_categories import POLICY_ID, classify_normal_expression  # noqa: E402
from oncoref.normal_tissues import normal_tissue_screening_policy  # noqa: E402


def classify(task):
    gid, symbol, rows, expected, unreliable = task
    return {
        "Ensembl_Gene_ID": gid,
        "Symbol": symbol,
        **classify_normal_expression(
            rows,
            routine_rna_contexts=expected["bulk_rna"],
            routine_ihc_contexts=expected["ihc"],
            reviewed_unreliable=unreliable,
        ),
    }


def build(out, workers):
    if out.exists() and any(out.iterdir()):
        raise ValueError("Use an empty output directory")
    out.mkdir(parents=True, exist_ok=True)
    space = canonical_gene_space()
    space = space.loc[space.biotype.eq("protein_coding")]
    print(f"Canonical coding universe: {len(space)}", flush=True)
    baseline = ROOT / "docs/audits/cta-recalibration-20261021"
    prior_manifest = json.loads((baseline / "run-manifest.json").read_text())
    expected = {k: set(map(tuple, v)) for k, v in prior_manifest["routine_contexts"].items()}
    receipts, parts = {}, []
    for modality, name in [("bulk_rna", "hpa_rna_consensus"), ("ihc", "hpa_normal_tissue")]:
        receipt = reference_data.provenance(name, "v23")
        assert sha256(receipt["path"]) == receipt["sha256"]
        receipts[name] = {k: receipt[k] for k in ["url", "sha256", "bytes", "version"]}
        parts.append(
            annotate_measurements(hpa._read_hpa(name, "v23"), modality, set(space.ensembl_gene_id))
        )
    observations = annotate_archive(pd.concat(parts, ignore_index=True))
    del parts
    groups = dict(iter(observations.groupby("Ensembl_Gene_ID", sort=False)))
    empty = observations.iloc[:0]
    unreliable = cta_regen.cross_reactive_ihc()
    tasks = (
        (
            r.ensembl_gene_id,
            r.symbol,
            groups.get(r.ensembl_gene_id, empty),
            expected,
            r.ensembl_gene_id in unreliable,
        )
        for r in space.itertuples()
    )
    rows = []
    with ProcessPoolExecutor(
        max_workers=workers, mp_context=multiprocessing.get_context("spawn")
    ) as pool:
        for i, result in enumerate(pool.map(classify, tasks, chunksize=20), 1):
            rows.append(result)
            if i % 1000 == 0:
                print(f"Classified {i}/{len(space)}", flush=True)
    frame = pd.DataFrame(rows)
    nominations = cta_sources.reviewed_publication_membership()
    frame["published_cta_nomination"] = frame.Ensembl_Gene_ID.isin(nominations.Ensembl_Gene_ID)
    path = baseline / "candidate-comparison.csv"
    assert sha256(path) == prior_manifest["outputs"][path.name]
    curated = pd.read_csv(path)
    frame["retained_recalibration_cta"] = frame.Ensembl_Gene_ID.isin(
        curated.loc[curated.retained, "Ensembl_Gene_ID"]
    )
    linked = frame.set_index("Ensembl_Gene_ID").loc[curated.Ensembl_Gene_ID]
    assert (linked.category.eq("Reproductive").to_numpy() == curated.evidence_pass.to_numpy()).all()
    assert frame.loc[frame.retained_recalibration_cta, "category"].eq("Reproductive").all()
    write_csv(frame, out / "gene-categories.csv.gz")
    for name, selected in {
        "reproductive": frame.category.eq("Reproductive"),
        "normal": frame.category.eq("Normal"),
        "other": frame.category.eq("Other"),
        "disposable": frame.other_subgroup.eq("Disposable"),
        "silent": frame.other_subgroup.eq("Silent"),
        "unassessed": frame.category.eq("Unassessed"),
        "conservative-normal-reference": frame.conservative_normal_reference,
    }.items():
        write_csv(frame.loc[selected], out / f"{name}.csv.gz")
    summary = (
        frame.groupby(["category", "other_subgroup"], dropna=False).size().reset_index(name="genes")
    )
    write_csv(summary, out / "summary.csv")
    summary_text = "\n".join(
        f"| {r.category} | {r.other_subgroup or '—'} | {r.genes:,} |" for r in summary.itertuples()
    )
    (out / "index.md").write_text(f"""# Genome-wide normal-expression categories

Policy `{POLICY_ID}`; {len(frame):,} canonical human protein-coding genes, pinned HPA v23. This universe includes genes outside published CTA nominations.

| Category | Other subgroup | Genes |
|---|---|---:|
{summary_text}

## Definitions

- **Reproductive:** positive reproductive RNA >0 or IHC, complete routine RNA, and all approved absolute RNA/selectivity/IHC gates. Publication nomination and manual CTA review are independent: this is broader than the retained curated CTA set.
- **Normal:** measured expression fails an absolute ceiling, the positive-IHC gate, or (when reproductive-positive) the pairwise selectivity gate. It means normal-expression counterevidence under this policy, not demonstrated clinical toxicity.
- **Other:** no positive reproductive evidence, complete routine RNA, no disqualifying positive IHC, and all approved absolute category/IHC ceilings pass. No reproductive ratio is required. **Disposable** has a measured positive in Disposable tissue; **Silent** does not. Both mean **below-threshold expression**, not literal absence from every normal tissue.
- **Unassessed:** no decisive normal-expression counterexample, but missing routine RNA, unresolved anatomy, unusable/reviewed-unreliable evidence, or an undefined required reproductive ratio prevents classification into Reproductive/Other. It is a separate evidence status and never part of Silent.

Missing/partial IHC retains baseline RNA limits and is explicitly recorded; this does not establish protein absence. The 50 RNA/80 IHC routine contexts do not cover all developmental stages or anatomical compartments. Numerical limits apply only to human HPA v23. [Exact category/IHC grid](../../normal-tissue-scopes.md).

## Preparing for subsequence comparisons

The three main categories are disjoint research screening groups. Sequence membership is a separate question: a reproductive gene can share a subsequence with a normal paralog, and normal RNA below a screening limit is still observed expression.

`normal_expression_observed` preserves any valid RNA >0 or positive IHC in Manageable/Protected/Vital tissue, even below a cutoff. The conservative normal-reference export includes those genes, all Normal genes, and Unassessed genes. It can overlap Reproductive and Other. A later subsequence comparison must examine exact protein/transcript occurrences and their reference versions across both sets; this report does not establish sequence specificity. Disposable expression remains separately visible.

## Downloads

- [All genes, categories, coverage and reasons](gene-categories.csv.gz)
- [Reproductive](reproductive.csv.gz) · [Normal](normal.csv.gz) · [Other union](other.csv.gz)
- [Disposable subgroup](disposable.csv.gz) · [Silent subgroup](silent.csv.gz) · [Unassessed](unassessed.csv.gz)
- [Conservative normal sequence-reference candidates](conservative-normal-reference.csv.gz)
- [Counts](summary.csv) · [Input and implementation checksums](run-manifest.json)

The [CTA recalibration report](../cta-recalibration-20261021/index.md) applies publication, family and nomination reviews in addition to the Reproductive expression screen. Its expression decisions reconcile exactly with this genome-wide classification. `retained_recalibration_cta` records that pre-external-review stage; use the [current clinical curation decisions](../cta-curation-20261021/index.md) for the final retained set and external safety holds.

Reproduce with `python scripts/normal_expression_categories_report.py --out <empty-directory> --workers 4`.
""")
    implementation = [
        "scripts/normal_expression_categories_report.py",
        "oncoref/normal_expression_categories.py",
        "oncoref/cta_normal_screen.py",
        "oncoref/normal_tissues.py",
        "scripts/cta_permitted_report.py",
        "scripts/cta_scope_report.py",
        "oncoref/data/normal-tissue-taxonomy.yaml",
        "oncoref/data/normal-tissue-mappings.csv",
        "oncoref/data/cta-ihc-unreliable.csv",
    ]
    manifest = {
        "schema": "oncoref.normal-expression-categories-report.v1",
        "policy_id": POLICY_ID,
        "screening_policy": normal_tissue_screening_policy("ihc_adaptive_cta"),
        "routine_contexts": prior_manifest["routine_contexts"],
        "sources": receipts,
        "canonical_sha256": hashlib.sha256(space.to_csv(index=False).encode()).hexdigest(),
        "inputs": {str(path.relative_to(ROOT)): sha256(path)},
        "implementation_sha256": {name: sha256(ROOT / name) for name in implementation},
        "outputs": {p.name: sha256(p) for p in sorted(out.iterdir()) if p.is_file()},
    }
    (out / "run-manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(summary.to_string(index=False), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()
    build(args.out, args.workers)
