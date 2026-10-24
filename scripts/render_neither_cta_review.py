#!/usr/bin/env python3
"""Render the high-normal-RNA audit without changing historical reports."""

from __future__ import annotations

import argparse
import json
import shutil
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from normal_cell_cta_report import digest, write  # noqa: E402
from review_neither_ctas import BASE, OUT  # noqa: E402

from oncoref import normal_rna_review as review  # noqa: E402
from oncoref.version import __version__  # noqa: E402


def number(value):
    return f"{float(value):.3f}".rstrip("0").rstrip(".") if value != "" else "Unmeasured"


def render(out, assets, publish=False):
    genes = pd.read_csv(out / "gene-results.csv.gz", keep_default_na=False, low_memory=False)
    queue = pd.read_csv(out / "follow-up.csv", keep_default_na=False)
    important = pd.read_csv(out / "important-ctas.csv", keep_default_na=False)
    high = pd.read_csv(out / "high-normal-rna.csv.gz", keep_default_na=False, low_memory=False)
    maxima = pd.read_csv(out / "rna-context-maxima.csv.gz", keep_default_na=False)
    sequence_genes = pd.read_csv(out / "sequence-gene-reviews.csv.gz", keep_default_na=False)
    counts = genes.Category.value_counts().to_dict()
    changes = genes.loc[genes.Change.ne("Unchanged")]
    prior_cta = genes["Previous category"].eq("CTA")
    prior_neither = genes["Previous category"].eq("Neither")
    original_important = pd.read_csv(BASE / "important-ctas.csv", keep_default_na=False)
    benchmark_ids = set(
        original_important.loc[
            original_important.Category.eq("CTA")
            & original_important["Previous category"].eq("CTA"),
            "Gene ID",
        ]
    )
    retained_important = important.loc[important["Gene ID"].isin(benchmark_ids)]
    lost_important = retained_important.loc[retained_important.Category.ne("CTA")]
    nominated = genes.loc[genes["Previous CTA screen"].ne("Not nominated")]
    neither_ids = set(genes.loc[genes.Category.eq("Neither"), "Gene ID"]) | set(
        sequence_genes.loc[sequence_genes.Category.eq("Neither"), "Gene ID"]
    )
    neither_maxima = maxima.loc[maxima.gene_id.isin(neither_ids) & ~maxima.permitted]
    assert not neither_maxima.maximum.gt(30).any()
    stats = {
        "policy": review.POLICY_ID,
        "genome": counts,
        "nominated": nominated.Category.value_counts().to_dict(),
        "previous_cta": int(prior_cta.sum()),
        "previous_neither": int(prior_neither.sum()),
        "cta_to_unsafe": int((prior_cta & genes.Category.eq("Unsafe")).sum()),
        "neither_to_unsafe": int((prior_neither & genes.Category.eq("Unsafe")).sum()),
        "remaining_neither_above_30": 0,
        "neither_maxima": neither_maxima.groupby("unit").maximum.max().to_dict(),
        "assessed_genes": int(genes["Normal RNA assessed"].sum()) + len(sequence_genes),
        "assessed_canonical_genes": int(genes["Normal RNA assessed"].sum()),
        "sequence_loci": len(sequence_genes),
        "sequence_loci_measured": int(sequence_genes["Normal RNA measured"].sum()),
        "sequence_loci_excluded": int(sequence_genes.Category.eq("Unsafe").sum()),
        "rna_measurements": int(maxima.measurements.sum()),
        "high_observations": len(high),
        "previous_benchmark_retained": int(retained_important.Category.eq("CTA").sum()),
        "previous_benchmark_total": len(retained_important),
        "benchmark_losses": lost_important.Gene.tolist(),
        "watchlist_genes": len(queue),
    }
    # Raw RNA observations from overlapping HPA summaries are not independent replication.
    newly_excluded_ids = set(changes["Gene ID"]) | set(
        sequence_genes.loc[sequence_genes.Category.eq("Unsafe"), "Gene ID"]
    )
    changed_high = high.loc[high.gene_id.isin(newly_excluded_ids)]
    per_source = changed_high.groupby(["source", "unit"], as_index=False).agg(
        genes=("gene_id", "nunique"),
        observations=("gene_id", "size"),
        maximum=("value", "max"),
    )
    write(per_source, out / "high-rna-by-source.csv")
    observed_normal = set(
        high.loc[high.category.isin(["vital", "protected", "manageable"]), "gene_id"]
    )
    stats["new_exclusions_only_unsited"] = len(set(changes["Gene ID"]) - observed_normal)

    proteins = pd.read_csv(assets / "target-proteins.csv.gz", keep_default_na=False)
    stats["proteins"] = proteins.Group.value_counts().to_dict()
    members = pd.read_csv(assets / "target-protein-occurrences.csv.gz", keep_default_na=False)
    by_gene = (
        members.loc[members.gene_id.isin(important["Gene ID"])]
        .groupby(["gene_id", "group"])
        .protein_sequence_id.nunique()
        .unstack(fill_value=0)
    )
    by_gene = by_gene.rename(columns={g: f"{g} complete proteins" for g in by_gene})
    write(
        important[["Gene", "Gene ID", "Category"]].merge(
            by_gene, left_on="Gene ID", right_index=True, how="left"
        ),
        out / "important-proteins.csv",
    )
    for name in [
        "exclusion-basis.json",
        "reference.json",
        "genes-without-complete-proteins.csv.gz",
        "unclassified-sequence-genes.csv.gz",
    ]:
        shutil.copyfile(assets / name, out / name)
    (out / "protein-assets.json").write_text(
        json.dumps(
            {
                "release": f"v{__version__}",
                "base_url": f"https://github.com/pirl-unc/oncoref/releases/download/v{__version__}/",
                "policy": review.POLICY_ID,
                "files": {
                    p.name: {"sha256": digest(p), "bytes": p.stat().st_size}
                    for p in sorted(assets.iterdir())
                    if p.is_file()
                },
            },
            indent=2,
        )
        + "\n"
    )

    watch_lines = []
    followup = [
        "# Candidate follow-up record",
        "",
        "Reviewed 24 October 2026 by Codex during the follow-up to PR #618. This is an implementation/evidence review, not an independent clinical adjudication. These records never override the general classification rules.",
        "",
        "[Report](index.md) · [Full follow-up CSV](follow-up.csv) · [Measured high RNA](high-normal-rna.csv.gz) · [All context maxima](rna-context-maxima.csv.gz) · [Pinned source URLs and checksums](source-provenance.json)",
        "",
    ]
    for row in queue.to_dict("records"):
        gene = row["Gene"]
        anchor = gene.lower()
        ntpm, tpm = (
            number(row["Normal RNA maximum (nTPM)"]),
            number(row["Normal RNA maximum (TPM)"]),
        )
        watch_lines.append(
            f"| [{gene}](follow-up.md#{anchor}) | {row['Category']} / {row['Detail']} | {ntpm} | {tpm} |"
        )
        followup += [
            f"## {gene}",
            "",
            f"**{row['Previous category']} → {row['Category']} / {row['Detail']}.** {row['Reason']}",
            "",
            f"Maximum outside explicitly permitted contexts: **{ntpm} nTPM** and **{tpm} TPM** (separate units). Reproductive bulk maximum: {row['Reproductive RNA (nTPM)']} nTPM.",
            "",
            f"**Open question:** {row['Review question']}",
            "",
            f"**Evidence needed:** {row['Reconsideration evidence']}",
            "",
            f"**Limits:** {row['Limitations']}",
            "",
        ]
        for field in [
            "High normal RNA findings",
            "Immune RNA findings",
            "Cell RNA findings",
            "Antibody findings",
            "Donor findings",
            "Publication sources",
        ]:
            if row[field]:
                followup += [f"**{field}:** {row[field]}", ""]
        followup += [
            f"[HPA v23 gene record](https://v23.proteinatlas.org/{row['Gene ID']}-{gene}) · [Previous decision and measured antibodies](../cta-normal-cells-20261023/explore.html)",
            "",
        ]
    (out / "follow-up.md").write_text("\n".join(followup))

    nc = stats["nominated"]
    pc = stats["proteins"]
    source_rows = "\n".join(
        f"| {r.source} | {r.unit} | {r.genes} | {r.observations:,} |"
        for r in per_source.itertuples()
    )
    (out / "index.md").write_text(f"""# CTA and Neither review: high normal RNA

**{counts["CTA"]:,} CTA genes remain. {stats["neither_to_unsafe"]} previous Neither genes and {stats["cta_to_unsafe"]} previous CTAs now enter the forbidden background. A further {stats["sequence_loci_excluded"]} noncanonical reference loci gain explicit RNA exclusions. No remaining Neither gene or reference locus has a measured non-permitted/unsited value above 30 TPM or nTPM in the audited sources.**

[Search genes](explore.html) · [Review findings and fixes](review.md) · [Candidate follow-up](follow-up.md) · [All changes](changed-genes.csv) · [Previous report](../cta-normal-cells-20261023/index.md)

## What changed?

The v1.17 rules put normal-cell RNA conflicts into Neither regardless of magnitude. The tissue-cluster rule also compared normal RNA with the reproductive maximum, allowing large normal signals to escape review when reproductive expression was still greater. The new absolute rule is independent of that ratio and of IHC confidence.

**Any observed normal RNA >30 TPM or >30 nTPM outside explicitly mapped Reproductive/Disposable contexts yields Unsafe / High RNA.** Exactly 30 does not trigger this additional rule. It applies to both former CTA and Neither genes, including Variable. Existing stricter protein, RNA and safety exclusions remain. No candidate is promoted and no clinical-target exception is made.

Unsafe here means a research exclusion and membership in the forbidden background, not demonstrated clinical toxicity. Unresolved anatomy does not establish permission: {stats["new_exclusions_only_unsited"]} newly excluded genes have only unsited high findings in this audit. Low-read clusters, ambiguous loci and possible contamination remain matters for adjudication; they do not silently waive the requested hard ceiling.

| Scope | CTA | Unsafe | Neither |
|---|---:|---:|---:|
| {len(nominated):,} mapped CTA nominations | {nc["CTA"]:,} | {nc["Unsafe"]:,} | {nc["Neither"]:,} |
| {len(genes):,} canonical human genes | {counts["CTA"]:,} | {counts["Unsafe"]:,} | {counts["Neither"]:,} |

Across remaining Neither genes the largest audited non-permitted/unsited observations are **{stats["neither_maxima"]["TPM"]:.4f} TPM** and **{stats["neither_maxima"]["nTPM"]:g} nTPM**. Missing observations remain missing. High reproductive/disposable values are permitted and therefore can still occur in Neither.

## Important targets

**MAGEA3 and ROPN1 leave the retained benchmark.** MAGEA3 has 53.1 nTPM and 251 reads in a colon undifferentiated-cell cluster. ROPN1 has 32.9 nTPM and 496 reads in salivary serous cells. Both cross the general rule; these are RNA-based research exclusions, not established target-mediated toxicity.

The other **{stats["previous_benchmark_retained"]} of the previous {stats["previous_benchmark_total"]}** benchmark targets remain, including XAGE1A/B, CTAG1A/B, PRAME, MAGEA4, CT83, SAGE1, CTCFL, GAGE1 and SSX4/5. CT55 remains retained. ACRBP, TEX101, GAGE12F, PAGE2B and TRIP13 now leave Neither because the expanded audit finds high normal RNA. HORMAD1, TDRD1, BRDT, DAZL, FOXR2 and CTAG2 remain Neither with explicit reconsideration questions.

| Watched gene | Current result | Maximum non-permitted/unsited nTPM | Maximum non-permitted/unsited TPM |
|---|---|---:|---:|
{chr(10).join(watch_lines)}

The [follow-up record](follow-up.md) contains measurements, uncertainties and evidence needed for reconsideration for each gene. It is a review queue, not an allowlist. [All important genes](important-ctas.csv) · [Important complete-protein decisions](important-proteins.csv).

## Coverage and units

The audit checks **{stats["rna_measurements"]:,} RNA observations across {stats["assessed_genes"]:,} genes**: all previous CTA and Neither genes plus the established important-target audit. Other Unsafe assignments remain inherited; they are not newly cleared. Sources are HPA bulk consensus, cell-type and tissue-cluster RNA, sorted immune summaries, individual HPA immune specimens, original DICE distributions, and the four reviewed CTAG2 heart RNA specimens. Both original TPM and nTPM columns from the HPA immune-sample table are checked separately.

The sequence inventory includes loci outside the canonical protein-coding table, including immune-receptor segments. Of {stats["sequence_loci"]:,} such exact reference IDs, {stats["sequence_loci_measured"]} have matching RNA observations and {stats["sequence_loci_excluded"]} cross the exclusion ceiling. Their complete proteins and incomplete translations now enter the forbidden background. The rest remain unassessed Neither; no absent measurement becomes zero and no paralog is substituted. [Sequence-locus decisions](sequence-gene-reviews.csv.gz).

- No averaging across donors, cells, tissues, studies or units can hide a value >30. A single reported high value suffices under this conservative policy.
- The older HPA cell-type table preserves its unsited contexts; it is not assigned an organ by guessing. Immune cells sampled in reproductive/disposable organs retain their existing hematopoietic mapping.
- HPA summaries and their original data overlap. They are not counted as independent biological replication.
- DICE slots are not stable donor identifiers across genes. RNA and IHC specimens are not assumed paired.
- Corneal CP10K is **not TPM** and does not enter the numerical ceiling. The existing corneal rule and uncertainty remain; developmental/culture data are not pooled into adult-normal evidence.
- This is measured-source coverage, not a claim that every normal tissue, physiological state or donor has been sampled. Human thresholds are not transferred to mouse or dog.

The following source counts overlap and must not be added as independent genes or cohorts:

| Source with newly excluding observations | Native unit | Changed genes | Observations >30 |
|---|---|---:|---:|
{source_rows}

## Files and API

[All genes](gene-results.csv.gz) · [CTA](cta-genes.csv.gz) · [Neither](neither-genes.csv.gz) · [Unsafe](unsafe-genes.csv.gz) · [All high observations](high-normal-rna.csv.gz) · [Context maxima including permitted tissues](rna-context-maxima.csv.gz) · [Source provenance](source-provenance.json) · [Run manifest](run-manifest.json)

```python
from oncoref import target_sets, normal_rna_review
target_sets.results()                 # current decisions and native-unit maxima
normal_rna_review.followup_reviews()   # watched genes, questions and evidence
normal_rna_review.sequence_gene_reviews()  # additional reference loci
normal_rna_review.forbidden_gene_ids()     # Unsafe canonical genes + reference loci
```

Complete protein sets are **{pc["CTA"]:,} CTA, {pc["Unsafe"]:,} Unsafe and {pc["Neither"]:,} Neither**. Shared exact sequences use Unsafe > Neither > CTA. The high-RNA exclusions also update incomplete translations in the additional forbidden background. [Protein assets and checksums](protein-assets.json) · [Exclusion recipe](exclusion-basis.json). Expression is gene-level and does not establish isoform-specific expression.

Reproduce from the same shared cache:

```bash
python scripts/review_neither_ctas.py --assets /path/to/new-assets --previous-sequences /path/to/pinned-ensembl115-assets
python scripts/render_neither_cta_review.py --assets /path/to/new-assets --publish-dataset
```

Policy: `{review.POLICY_ID}`. Baseline: `{review.BASE_POLICY_ID}`. The prior published report remains unchanged.
""")
    columns = [
        "Gene",
        "Gene ID",
        "Category",
        "Detail",
        "Reason",
        "Change",
        "Previous category",
        "Reproductive RNA (nTPM)",
        "Normal RNA maximum (TPM)",
        "Normal RNA maximum (nTPM)",
        "High normal RNA findings",
        "Immune RNA findings",
        "Cell RNA findings",
        "Antibody findings",
        "Donor findings",
    ]
    payload = genes[columns].copy()
    payload["Important"] = payload["Gene ID"].isin(important["Gene ID"])
    template = (ROOT / "scripts/templates/neither-review-explorer.html").read_text()
    for key, value in {
        "__CTA__": counts["CTA"],
        "__NEITHER__": counts["Neither"],
        "__CHANGED__": len(changes),
    }.items():
        template = template.replace(key, str(value))
    (out / "explore.html").write_text(
        template.replace("__DATA__", payload.to_json(orient="records").replace("<", "\\u003c"))
    )
    if publish:
        shutil.copyfile(out / "gene-results.csv.gz", ROOT / "oncoref/data/target-gene-sets.csv.gz")
        shutil.copyfile(out / "follow-up.csv", ROOT / "oncoref/data/cta-followup-reviews.csv")
        shutil.copyfile(
            out / "sequence-gene-reviews.csv.gz",
            ROOT / "oncoref/data/sequence-gene-rna-reviews.csv.gz",
        )
    (out / "summary.json").write_text(json.dumps(stats, indent=2) + "\n")
    inputs = [
        "oncoref/normal_rna_review.py",
        "scripts/review_neither_ctas.py",
        "scripts/render_neither_cta_review.py",
        "scripts/templates/neither-review-explorer.html",
        "oncoref/data/cta-followup-notes.csv",
        "oncoref/data/cta-donor-rna.csv",
        "oncoref/data/normal-tissue-taxonomy.yaml",
        "oncoref/data/normal-tissue-mappings.csv",
        "oncoref/normal_cell_data.py",
        "oncoref/normal_cell_screen.py",
        "oncoref/data/normal-cell-sources.json",
        "docs/audits/cta-normal-cells-20261023/gene-results.csv.gz",
        "docs/audits/cta-normal-cells-20261023/run-manifest.json",
        "docs/audits/cta-normal-cells-20261023/unclassified-sequence-genes.csv.gz",
    ]
    (out / "run-manifest.json").write_text(
        json.dumps(
            {
                "policy": review.POLICY_ID,
                "baseline_policy": review.BASE_POLICY_ID,
                "implementation_and_inputs": {p: digest(ROOT / p) for p in inputs},
                "outputs": {
                    p.name: digest(p)
                    for p in sorted(out.iterdir())
                    if p.is_file() and p.name != "run-manifest.json"
                },
                "forbidden_categories": ["Unsafe"],
                "summary": stats,
            },
            indent=2,
        )
        + "\n"
    )
    print(json.dumps(stats, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUT)
    parser.add_argument("--assets", type=Path, required=True)
    parser.add_argument("--publish-dataset", action="store_true")
    args = parser.parse_args()
    render(args.output, args.assets, args.publish_dataset)
