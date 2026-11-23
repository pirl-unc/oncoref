#!/usr/bin/env python3
"""Build CTA/Unsafe/Neither genes and exact-protein backgrounds from pinned assets."""

from __future__ import annotations

import argparse
import json
import shutil
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from cta_clinical_curation_report import read_verified  # noqa: E402
from cta_permitted_report import annotate_archive  # noqa: E402
from cta_recalibration_report import markdown_table  # noqa: E402
from cta_scope_report import annotate_measurements, sha256, write_csv  # noqa: E402
from cta_screen_report import explorer as old_explorer  # noqa: E402
from normal_sequence_sets_report import fasta_records, write_fasta  # noqa: E402

from oncoref.cta_regen import cross_reactive_ihc  # noqa: E402
from oncoref.load_dataset import get_data  # noqa: E402
from oncoref.sequence_groups import protein_sequence_id  # noqa: E402
from oncoref.target_sets import (  # noqa: E402
    GROUPS,
    POLICY_ID,
    classify_genes,
    repartition_protein_table,
)

EXPRESSION = ROOT / "docs/audits/normal-expression-categories-20261021-r2"
SCREEN = ROOT / "docs/audits/cta-screen-20261023"
SEQUENCES = ROOT / "docs/audits/normal-sequence-sets-20261021"
ASSET_URL = "https://github.com/pirl-unc/oncoref/releases/download/v1.16.0/"


def explorer(out, result):
    # Keep the browser payload compact; full per-tissue evidence stays in CSV.
    columns = [
        "Gene",
        "Gene ID",
        "Category",
        "Detail",
        "Reason",
        "Reproductive RNA (nTPM)",
        "Maximum RNA (nTPM)",
        "IHC coverage",
        "Donor review",
        "Donor findings",
    ]
    old_explorer(out, result[columns])
    path = out / "explore.html"
    head, separator, body = path.read_text().partition('<script type="application/json"')
    head = head.replace('href="gene-results.csv"', 'href="gene-results.csv.gz"')
    head = head.replace(
        'href="variable.csv">Variable candidates', 'href="neither-genes.csv.gz">Neither genes'
    )
    start = head.index('<p class="note">')
    end = head.index("</p>", start) + 4
    head = (
        head[:start]
        + '<p class="note"><b>CTA = target candidates. Unsafe = forbidden background. Neither = neither.</b> These are research screening decisions. Low-expression and Variable genes are not silently added to the forbidden set. Donor RNA is reviewed for CTAG2 only.</p>'
        + head[end:]
    )
    head = head.replace(
        "<option>Retained</option><option>Variable</option><option>Excluded</option><option>Unassessed</option>",
        "<option>CTA</option><option>Unsafe</option><option>Neither</option>",
    )
    head = (
        head.replace(".Variable{", ".Neither{")
        .replace(".Excluded{", ".Unsafe{")
        .replace(".Retained{", ".CTA{")
    )
    path.write_text(head + separator + body)


def verify_assets(directory, names):
    manifest = json.loads((SEQUENCES / "run-manifest.json").read_text())
    for name in names:
        if sha256(directory / name) != manifest["outputs"][name]:
            raise ValueError(f"Sequence asset checksum mismatch: {name}")
    return manifest


def proteins(assets, previous, genes):
    names = [
        "proteins.csv.gz",
        "protein-occurrences.csv.gz",
        "translations-for-review.csv.gz",
        "translations-for-review.fasta.gz",
        "normal.fasta.gz",
        "unassessed.fasta.gz",
        "reproductive.fasta.gz",
        "other.fasta.gz",
    ]
    manifest = verify_assets(previous, names)
    reference_id = manifest["reference"]["reference_id"]
    old = pd.read_csv(previous / "proteins.csv.gz", keep_default_na=False)
    table = repartition_protein_table(old, genes, reference_id=reference_id)
    seqs = {}
    for name in names[-4:]:
        for pid, sequence in fasta_records(previous / name):
            if pid in seqs or protein_sequence_id(sequence) != pid:
                raise ValueError("Duplicate or incorrectly identified protein sequence")
            seqs[pid] = sequence
    if set(seqs) != set(table["Protein sequence"]):
        raise ValueError("Protein table and FASTA sequence coverage differ")
    write_csv(table, assets / "target-proteins.csv.gz")
    group_for_protein = table.set_index("Protein sequence").Group
    group_for_gene = genes.set_index("Gene ID").Category
    members = pd.read_csv(previous / "protein-occurrences.csv.gz", keep_default_na=False)
    members = members.rename(
        columns={
            "group": "previous_group",
            "expression_group": "previous_expression_group",
            "exclude_by_default": "previous_exclusion",
        }
    )
    members["group"] = members.protein_sequence_id.map(group_for_protein)
    members["gene_category"] = members.gene_id.map(group_for_gene).fillna("Neither")
    members["forbidden"] = members.group.eq("Unsafe")
    if members.group.isna().any():
        raise ValueError("Lost a complete protein occurrence")
    write_csv(members, assets / "target-protein-occurrences.csv.gz")
    for group in GROUPS:
        selected = table.loc[table.Group.eq(group)]
        write_csv(selected, assets / f"{group.lower()}-proteins.csv.gz")
        write_fasta(
            assets / f"{group.lower()}.fasta.gz",
            ((pid, seqs[pid]) for pid in selected["Protein sequence"]),
        )
    review = pd.read_csv(previous / "translations-for-review.csv.gz", keep_default_na=False)
    review = review.rename(
        columns={
            "expression_group": "previous_expression_group",
            "exclude_by_default": "previous_exclusion",
        }
    )
    review["gene_category"] = review.gene_id.map(group_for_gene).fillna("Neither")
    review["forbidden"] = review.gene_category.eq("Unsafe")
    review["target"] = False
    write_csv(review, assets / "translations-for-review.csv.gz")
    forbidden_records = set(review.loc[review.forbidden, "occurrence_id"])
    write_fasta(
        assets / "additional-forbidden.fasta.gz",
        (
            (oid, seq)
            for oid, seq in fasta_records(previous / "translations-for-review.fasta.gz")
            if oid in forbidden_records
        ),
    )
    shutil.copyfile(
        previous / "translations-for-review.fasta.gz", assets / "translations-for-review.fasta.gz"
    )
    missing = genes.loc[~genes["Gene ID"].isin(members.gene_id)]
    write_csv(missing, assets / "genes-without-complete-proteins.csv.gz")
    unknown = set(members.gene_id) | set(review.gene_id)
    unknown -= set(genes["Gene ID"])
    write_csv(
        pd.DataFrame(
            {
                "Gene ID": sorted(unknown),
                "Category": "Neither",
                "Reason": "No canonical expression classification",
            }
        ),
        assets / "unclassified-sequence-genes.csv.gz",
    )
    recipe = {
        "policy": POLICY_ID,
        "reference_id": reference_id,
        "target": ["cta.fasta.gz"],
        "forbidden": ["unsafe.fasta.gz", "additional-forbidden.fasta.gz"],
        "neither": ["neither.fasta.gz"],
        "neither_forbidden": False,
        "sequence_precedence": ["Unsafe", "Neither", "CTA"],
        "rule": "Enumerate contiguous standard-amino-acid k-mers within records; subtract the forbidden union from CTA k-mers. Never join across ambiguous residues, stops or record boundaries.",
        "limits": "Neither is neither a target nor proven absent. Shared k-mers can cross otherwise disjoint proteins; gene-level RNA does not establish isoform expression. Translation gaps remain explicit.",
    }
    (assets / "exclusion-basis.json").write_text(json.dumps(recipe, indent=2) + "\n")
    (assets / "reference.json").write_text(json.dumps(manifest["reference"], indent=2) + "\n")
    return table, manifest, names, len(review), len(missing)


def build(out, assets, previous, rna_path, ihc_path, dataset=None):
    for folder in [out, assets]:
        if folder.exists() and any(folder.iterdir()):
            raise ValueError("Use empty output and asset directories")
        folder.mkdir(parents=True, exist_ok=True)
    source_manifest = json.loads((EXPRESSION / "run-manifest.json").read_text())
    for key, path in [("hpa_rna_consensus", rna_path), ("hpa_normal_tissue", ihc_path)]:
        if sha256(path) != source_manifest["sources"][key]["sha256"]:
            raise ValueError(f"Unpinned HPA input: {key}")
    categories = read_verified(EXPRESSION, "gene-categories.csv.gz")
    screen = read_verified(SCREEN, "gene-results.csv")
    ids = set(categories.Ensembl_Gene_ID)
    parts = [
        annotate_measurements(pd.read_csv(path, sep="\t"), modality, ids)
        for modality, path in [("bulk_rna", rna_path), ("ihc", ihc_path)]
    ]
    observations = annotate_archive(pd.concat(parts, ignore_index=True))
    observations["taxonomy_id"] = 9606
    safety = get_data("cta-external-safety-reviews")
    genes = classify_genes(
        categories,
        screen,
        observations,
        safety_holds=safety.loc[safety.action.eq("hold"), "Ensembl_Gene_ID"],
        unreliable_genes=cross_reactive_ihc(),
        donor_results=read_verified(SCREEN, "donor-results.csv"),
    )
    # Current-candidate RNA exclusion must agree with the published ordered audit.
    tissue = read_verified(SCREEN, "tissue-evidence.csv.gz")
    rna = pd.to_numeric(tissue.rna_ntpm, errors="coerce")
    ceiling = pd.to_numeric(tissue.rna_ceiling_ntpm, errors="coerce")
    observed_high = set(tissue.loc[rna.gt(ceiling), "Ensembl_Gene_ID"])
    chosen = genes.loc[genes["Gene ID"].isin(screen["Gene ID"])]
    assert set(chosen.loc[chosen["Normal RNA above ceiling"], "Gene ID"]) == observed_high
    write_csv(genes, out / "gene-results.csv.gz")
    if dataset:
        shutil.copyfile(out / "gene-results.csv.gz", dataset)
    for group in GROUPS:
        write_csv(genes.loc[genes.Category.eq(group)], out / f"{group.lower()}-genes.csv.gz")
    write_csv(chosen, out / "cta-candidate-results.csv.gz")
    lost = chosen.loc[chosen["Previous CTA screen"].eq("Retained") & ~chosen.Target]
    write_csv(lost, out / "removed-cta-candidates.csv")
    write_csv(genes.loc[genes.Detail.eq("Variable")], out / "variable.csv")
    sensitivity = pd.DataFrame(
        [
            {
                "Minimum RNA (nTPM)": threshold,
                "Removed without IHC": int(
                    (
                        chosen["Previous CTA screen"].eq("Retained")
                        & chosen["IHC measurements"].eq(0)
                        & chosen["Reproductive RNA (nTPM)"].lt(threshold)
                    ).sum()
                ),
                "Removed without positive IHC": int(
                    (
                        chosen["Previous CTA screen"].eq("Retained")
                        & chosen["Positive IHC measurements"].eq(0)
                        & chosen["Reproductive RNA (nTPM)"].lt(threshold)
                    ).sum()
                ),
            }
            for threshold in [1, 2, 3, 5, 10]
        ]
    )
    write_csv(sensitivity, out / "admission-sensitivity.csv")
    del observations, parts
    table, sequence_manifest, asset_names, incomplete, missing = proteins(assets, previous, genes)
    shutil.copyfile(assets / "exclusion-basis.json", out / "exclusion-basis.json")
    counts = pd.DataFrame(
        {
            "Category": GROUPS,
            "Genes": [int(genes.Category.eq(g).sum()) for g in GROUPS],
            "Exact protein sequences": [int(table.Group.eq(g).sum()) for g in GROUPS],
        }
    )
    write_csv(counts, out / "counts.csv")
    subgroups = genes.groupby(["Category", "Detail"]).size().reset_index(name="Genes")
    write_csv(subgroups, out / "subgroups.csv")
    explorer(out, genes)
    named = genes.loc[
        genes.Gene.isin(
            [
                "CTAG1A",
                "CTAG1B",
                "CTAG2",
                "XAGE1A",
                "XAGE1B",
                "PRAME",
                "MAGEA4",
                "MAGEA3",
                "CT83",
                "ZNF730",
                "IL5",
                "MAGEA12",
            ]
        ),
        ["Gene", "Category", "Detail", "Reproductive RNA (nTPM)", "IHC coverage"],
    ]
    count = int(genes.Target.sum())
    (out / "index.md").write_text(f"""# CTA, Unsafe and Neither

**{count} CTA candidates**, down from 740. **{len(lost)} previous candidates move to Neither**, not the forbidden background. The rule is **at least 5 nTPM reproductive RNA when there is no positive IHC anywhere measured**. All 21 currently retained clinical/requested benchmark genes remain CTA candidates.

[Find any gene](explore.html) · [Removed candidates and reasons](removed-cta-candidates.csv) · [Gene sets](#the-three-sets) · [Protein downloads](#protein-sets-and-k-mer-exclusions) · [Threshold comparison](admission-sensitivity.csv)

## The three sets

| Set | Meaning | Target? | Forbidden background? |
|---|---|---|---|
| **CTA** | Curated candidate passing the existing RNA/IHC/selectivity rules and the new admission floor | Yes | No |
| **Unsafe** | Positive IHC in Vital/Protected/Manageable tissue, above-ceiling normal RNA, or an independent reviewed safety hold | No | Yes |
| **Neither** | Low-expression, Variable, unresolved, uncurated or selectivity-only failures without decisive normal-expression exclusion | No | **No** |

These are research labels. CTA means eligible under this screen, not proven safe treatment; Unsafe means operationally forbidden under the policy, not demonstrated clinical toxicity. Neither is not evidence of protein absence.

{markdown_table(counts)}

The gene universe is **{len(genes):,} canonical human protein-coding genes**, not just the 2,540 source-nominated CTA candidates. [All genes](gene-results.csv.gz) · [CTA](cta-genes.csv.gz) · [Unsafe](unsafe-genes.csv.gz) · [Neither](neither-genes.csv.gz) · [Candidate funnel](cta-candidate-results.csv.gz) · [Neither subgroups](subgroups.csv).

## What changed

- **150** previously retained genes have reproductive RNA <5 nTPM and no measured IHC at all. Their RNA is also <5 nTPM in every measured tissue.
- **IL5** adds one more removal: reproductive maximum 4.9 nTPM, measured IHC but no positive IHC anywhere. “No positive IHC” therefore removes **151**, rather than the 150 estimated for “no IHC.”
- **CTAG2 remains Neither / Variable**, with its reviewed donor RNA and negative local IHC preserved. Negative IHC and RNA donors are different people; only CTAG2 has imported individual-donor RNA evidence.
- **ZNF730 is Neither / Low**: it fails CTA selectivity and has low RNA without positive protein evidence, but does not exceed its normal RNA ceilings. It is not automatically a forbidden protein.
- Positive normal IHC, including Low, and independent safety holds still take precedence. **MAGEA12 stays Unsafe.**

{markdown_table(named)}

The 5 nTPM value is an **admission floor**, not a replacement normal-tissue safety ceiling. A gene below 5 nTPM everywhere can still be Unsafe if, for example, brain RNA exceeds its applicable 2 nTPM baseline. Missing measurements never become zero. Low requires usable routine RNA coverage; missing/unreliable evidence remains explicitly Unassessed within Neither.

IHC counts describe curated tissue/cell-type annotations, not independent donor counts. Missing IHC and measured-negative IHC are reported separately.

## The existing tissue rules remain

Each cell is **maximum nTPM / nominal minimum selectivity**. Relaxation requires reliable negative IHC across routine local cell types plus reliable positive reproductive IHC in the curated profile.

| Tissue | Baseline | Approved negative | Supported negative | Enhanced negative |
|---|---:|---:|---:|---:|
| Vital | 2 / 80% | 4 / 75% | 6 / 70% | 12 / 65% |
| Protected | 6 / 75% | 8 / 70% | 12 / 65% | 24 / 60% |
| Manageable | 10 / 70% | 20 / 65% | 30 / 60% | 60 / 55% |

CTA selection still requires `N ≤ C`, `N ≤ R × (1−s)/s + 0.1`, and `N < R`. A selectivity failure alone removes CTA eligibility; it does **not** create an Unsafe assignment. Reviewed Variable evidence is a separate reconsideration route, excluded from both target and forbidden sets unless competing protein/safety evidence wins. [Prior ordered rules and donor evidence](../cta-screen-20261023/index.md).

## Protein sets and k-mer exclusions

Complete identical amino-acid sequences are assigned once, with every encoding gene/transcript/protein occurrence retained. Priority is **Unsafe → Neither → CTA**. Thus an Unsafe encoding gene makes the shared protein forbidden; an unresolved/Neither encoding gene blocks target selection without making that protein forbidden. XAGE1A/B remain separate loci while their identical proteins are grouped.

| Set | Protein inventory | FASTA |
|---|---|---|
| CTA | [Proteins]({ASSET_URL}cta-proteins.csv.gz) | [Sequences]({ASSET_URL}cta.fasta.gz) |
| Unsafe | [Proteins]({ASSET_URL}unsafe-proteins.csv.gz) | [Sequences]({ASSET_URL}unsafe.fasta.gz) |
| Neither | [Proteins]({ASSET_URL}neither-proteins.csv.gz) | [Sequences]({ASSET_URL}neither.fasta.gz) |

**The forbidden union is `unsafe.fasta.gz` plus `additional-forbidden.fasta.gz`. Neither, Variable, Low and Unassessed are not added.** The additional file contains incomplete/ambiguous translation records only from Unsafe genes. The historical Normal + Unassessed exclusion recipe is not the current default.

Full proteins are disjoint; k-mers can still overlap. Subtract all forbidden k-mers from CTA k-mers and never join across record boundaries, ambiguous residues or stops. **{incomplete:,} translation records** require review; **{missing:,} classified genes** have no complete usable protein. Their gaps are preserved, not interpreted as absence. Gene RNA does not establish isoform-specific translation or peptide-HLA presentation.

**Gene and protein counts differ.** For example, MAGEA3 passes at the gene level, but its complete reference protein is also encoded by unresolved alternate-contig gene `ENSG00000292186`. That shared protein stays Neither pending identity review, as it did in the previous Unassessed sequence export. Matching the gene name alone does not resolve the locus. XAGE1A/B and CTAG1A/B shared proteins remain in CTA.

[Machine-readable recipe]({ASSET_URL}exclusion-basis.json) · [Additional forbidden records]({ASSET_URL}additional-forbidden.fasta.gz) · [All proteins]({ASSET_URL}target-proteins.csv.gz) · [All occurrences]({ASSET_URL}target-protein-occurrences.csv.gz) · [Translation review]({ASSET_URL}translations-for-review.csv.gz) · [Missing complete proteins]({ASSET_URL}genes-without-complete-proteins.csv.gz) · [Unclassified sequence genes]({ASSET_URL}unclassified-sequence-genes.csv.gz).

## Provenance and API

Pinned **human HPA v23** expression and **Ensembl 115 / GRCh38.p14** complete-sequence inventory. Existing source tissue/cell mappings and antibody confidence are preserved; human thresholds are not transferred to dogs or mice. Original dated reports remain available. [HPA sources and antibodies](../cta-curation-20261021-r2/index.md) · [Reference identity]({ASSET_URL}reference.json) · [Input/code/output checksums](run-manifest.json).

```python
from oncoref import target_sets
genes = target_sets.results()
targets = target_sets.gene_ids("CTA")
forbidden = target_sets.gene_ids("Unsafe")
neither = target_sets.gene_ids("Neither")
```

Use `target_sets` for the current three-way sets. Older `cta_screen` results and normal-sequence exports retain their explicitly dated policies.
""")
    inputs = [
        EXPRESSION / "run-manifest.json",
        EXPRESSION / "gene-categories.csv.gz",
        SCREEN / "run-manifest.json",
        SCREEN / "gene-results.csv",
        SCREEN / "tissue-evidence.csv.gz",
        SCREEN / "donor-results.csv",
        SEQUENCES / "run-manifest.json",
        ROOT / "oncoref/data/cta-external-safety-reviews.csv",
        ROOT / "oncoref/data/cta-ihc-unreliable.csv",
    ]
    implementation = [
        Path(__file__),
        ROOT / "oncoref/target_sets.py",
        ROOT / "oncoref/cta_screen.py",
        ROOT / "oncoref/normal_tissues.py",
        ROOT / "scripts/cta_scope_report.py",
        ROOT / "scripts/cta_permitted_report.py",
    ]
    manifest = {
        "policy": POLICY_ID,
        "admission_floor_ntpm": 5,
        "no_positive_ihc_floor": True,
        "forbidden_categories": ["Unsafe"],
        "variable_category": "Neither",
        "sources": source_manifest["sources"],
        "sequence_asset_base": sequence_manifest["asset_base_url"],
        "sequence_inputs": {name: sequence_manifest["outputs"][name] for name in asset_names},
        "inputs": {str(p.relative_to(ROOT)): sha256(p) for p in inputs},
        "implementation_sha256": {str(p.relative_to(ROOT)): sha256(p) for p in implementation},
        "outputs": {p.name: sha256(p) for p in sorted(out.iterdir())},
        "release_assets": {p.name: sha256(p) for p in sorted(assets.iterdir())},
        "published_dataset_sha256": sha256(dataset) if dataset else None,
    }
    (out / "run-manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(counts.to_string(index=False), flush=True)
    print(f"Removed CTA candidates: {len(lost)}; assets: {assets}", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ["out", "assets", "previous-sequences", "rna", "ihc"]:
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--dataset", type=Path)
    args = parser.parse_args()
    build(args.out, args.assets, args.previous_sequences, args.rna, args.ihc, args.dataset)
