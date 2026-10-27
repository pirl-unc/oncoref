#!/usr/bin/env python3
"""Build a versioned research shortlist and an honest literature-review ledger."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from oncoref import cta_policy  # noqa: E402
from oncoref import cta_research_policy as policy  # noqa: E402


def digest(path):
    """Return the content SHA-256 used to pin an audit input or output."""
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write(frame: pd.DataFrame, path):
    """Write a data table with deterministic compression when requested."""
    frame.to_csv(
        path,
        index=False,
        compression={"method": "gzip", "mtime": 0} if path.name.endswith(".gz") else None,
    )


def build(out: Path, *, expansion, paper_receipts, publish=False) -> None:
    """Build the versioned audit tables, report and source provenance."""
    if out.exists() and any(out.iterdir()):
        raise ValueError("Use an empty output directory")
    out.mkdir(parents=True, exist_ok=True)
    base = ROOT / "docs/audits/cta-policy-20261025"
    old = ROOT / "docs/audits/cta-literature-20261025"
    manifest = json.loads((base / "run-manifest.json").read_text())
    for name in ["decision-facts.csv.gz", "reviewed-findings.csv", "candidate-results.csv.gz"]:
        if digest(base / name) != manifest["outputs"][name]:
            raise ValueError(f"Pinned input mismatch: {name}")
    facts = pd.read_csv(base / "decision-facts.csv.gz")
    findings = cta_policy.findings()
    clinical = policy.clinical_evidence()
    decisions = policy.classify(facts, findings, clinical)
    nominated = decisions.loc[decisions.Nominated]
    if len(nominated) != 2540:
        raise ValueError("Mapped nomination universe changed")
    changes = decisions.loc[decisions.Change.ne("Unchanged")]
    write(decisions, out / "gene-results.csv.gz")
    write(nominated, out / "candidate-results.csv.gz")
    write(decisions.loc[decisions.Target], out / "cta-genes.csv")
    write(changes, out / "changes.csv")
    write(clinical, out / "clinical-evidence.csv")
    identity_reviews = ROOT / "oncoref/data/cta-intake-identity-reviews.csv"
    shutil.copyfile(identity_reviews, out / "unmapped-identity-reviews.csv")
    # Preserve source-review and search depth; never label retrieved papers read.
    literature = pd.read_csv(old / "all-genes.csv", keep_default_na=False)
    columns = [
        "Gene ID",
        "Assessment basis",
        "New finding",
        "Next action",
        "New review PMIDs",
        "New review sources",
        "Search coverage",
        "Truncated searches",
        "Zero-hit topic searches",
    ]
    review = nominated[
        ["Gene ID", "Gene", "Category", "Expression category", "Selection basis"]
    ].merge(literature[columns], on="Gene ID", how="left", validate="one_to_one")
    review["Review status"] = "Unassessed primary papers"
    specific = review["Assessment basis"].ne(
        "Source flags and policy only; new papers not individually adjudicated"
    )
    review.loc[specific, "Review status"] = "Prior specific assessment; see recorded depth"
    for row in clinical.to_dict("records"):
        selected = review["Gene ID"].eq(row["Gene ID"])
        review.loc[selected, "Review status"] = "Clinical comparison reviewed 2026-10-26"
        review.loc[selected, "Assessment basis"] = row["Review depth"]
        review.loc[selected, "New finding"] = row["Clinical outcome"] + " " + row["Adverse events"]
        review.loc[selected, "New review sources"] = (
            row["Evidence sources"] + " | " + row["Normal selectivity sources"]
        )
        review.loc[selected, "Next action"] = row["Limitations"]
    cases = cta_policy.adjudications().set_index("Gene ID")
    review["Antibody adjudication"] = (
        review["Gene ID"].map(cases.Decision).fillna("No case in this adjudication queue")
    )
    review["Remaining work"] = [
        (
            "Primary papers remain individually unassessed. "
            if status == "Unassessed primary papers"
            else "Existing assessment is limited to the cited evidence; comprehensive coverage remains unverified. "
        )
        + (
            f"Expand {truncated} truncated searches and resolve gene/alias attribution. "
            if truncated
            else ""
        )
        + (
            "Both topic queries returned zero leads; follow aliases and citation chains before interpreting absence. "
            if no_hits == 2
            else ""
        )
        + "Preserve controls, contradictory evidence, full-text access and exact review portions."
        for status, truncated, no_hits in zip(
            review["Review status"], review["Truncated searches"], review["Zero-hit topic searches"]
        )
    ]
    write(review, out / "literature-review-status.csv")
    summary = {
        "policy": policy.POLICY_ID,
        "expression_policy": policy.BASE_POLICY_ID,
        "scope": "Evidence-weighted CTA target-development research selection",
        "nominated": nominated.Category.value_counts().to_dict(),
        "selected_genes": int(decisions.Target.sum()),
        "changes": changes[["Gene", "Gene ID", "Expression category", "Category"]].to_dict(
            "records"
        ),
        "numerical_rna_threshold_changes": 0,
        "new_breadth_admission_cutoff": None,
        "clinical_comparison_records": len(clinical),
        "review_status_counts": review["Review status"].value_counts().to_dict(),
        "issue_638_complete": False,
        "issue_638_reason": "Comprehensive individual primary-paper review, search expansion and full-text coverage remain incomplete.",
        "normal_protein_background_changed": False,
        "calibration": "User-requested PRAME/MAGE-A4 benchmarks; not independent validation",
    }
    extra_inputs = []
    if expansion is not None:
        summary_path = expansion / "expanded-summary.json"
        if not summary_path.exists():
            summary_path = expansion / "summary.json"
        retrieval = json.loads(summary_path.read_text())
        if retrieval["original_ledger_sha256"] != digest(old / "searches.csv"):
            raise ValueError("Expanded queries do not identify the pinned original ledger")
        for name, expected in retrieval["outputs"].items():
            source = expansion / name
            if not source.exists():
                source = expansion / ("expanded-" + name)
            if digest(source) != expected:
                raise ValueError("Expanded-search artifact checksum mismatch")
            shutil.copyfile(source, out / ("expanded-" + name))
            extra_inputs.append(out / ("expanded-" + name))
        summary["expanded_queries"] = retrieval["queries"]
        summary["expanded_queries_paginated_completely"] = retrieval["pagination_complete"]
        summary["expanded_genes"] = pd.read_csv(out / "expanded-searches.csv")["Gene ID"].nunique()
        shutil.copyfile(summary_path, out / "expanded-summary.json")
        extra_inputs.append(out / "expanded-summary.json")
    if paper_receipts is not None:
        shutil.copyfile(paper_receipts, out / "paper-receipts.json")
        extra_inputs.append(out / "paper-receipts.json")
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    lines = "\n".join(
        f"| {r['Gene']} | {r['Expression category']} | {r['Category']} |"
        for r in changes.to_dict("records")
    )
    (out / "index.md").write_text(f"""# Recalibrated CTA research selection

**{int(decisions.Target.sum())} CTA research candidates** under `{policy.POLICY_ID}`.
[Complete list](cta-genes.csv) · [All nominations](candidate-results.csv.gz) ·
[Every change](changes.csv) · [Reviewed clinical comparisons](clinical-evidence.csv).

The user requested a common evidence policy calibrated to include PRAME and
MAGE-A4 alongside genes of comparable quality. Both already pass the numerical
RNA gates. This revision changes the evidence threshold for research selection,
not the RNA measurements or cutoffs. NY-ESO-1 remains represented by CTAG1A/B.
This is benchmark-driven calibration, not independent prospective validation.

| Gene | Earlier expression category | Research selection |
|---|---|---|
{lines}

## Common evidence rule

The previous expression-qualified candidates remain selected. A new clinical
qualification route requires a prospective treated cohort, repeated objective
tumor responses (more than one; stable disease and isolated lesion shrinkage do
not count), target-specific evidence, and review of normal selectivity. Every
active normal finding must be linked and considered; adding a new finding
invalidates the review until reconsidered. An unresolved target-mediated organ
toxicity concern or substantial accepted normal protein prevents this route.
An unresolved target-mediated organ toxicity concern also holds an otherwise
expression-qualified candidate for review.

All measured RNA/IHC gates, reproductive qualification, coverage, hard >30 RNA
exclusions, donor variability and cell-RNA holds still apply. Neither a familiar
gene symbol nor trial registration qualifies a candidate. No breadth cutoff is
used. These criteria are a research evidence threshold, not a clinical safety
score or a claim that all encoded peptides can be targeted.

PRAME retains the accepted low adrenal-protein observation. Its selected
clinical product has prospective responses and tested peptide selectivity;
these support continued target development without asserting that adrenal
protein is absent. MAGE-A4 retains the unresolved HPA antibody findings;
independent assay-validation evidence and clinical responses support research
selection without declaring those stains experimentally resolved.

MAGE-A1 and MAGE-A10 were assessed as comparisons. MAGE-A1's IMA202 cohort had
no RECIST objective response; MAGE-A10's NSCLC study reported one response after
repeat treatment with subsequent severe complications. Neither establishes
the repeated-response criterion. MAGE-A3/12 remain subject to their existing
hard RNA/safety exclusions. All findings and qualification limits are exported.

## API and sequence scope

```python
from oncoref import cta_research_policy
genes = cta_research_policy.gene_ids()
decisions = cta_research_policy.results()
```

This is the new **research selection**. `target_sets` and `cta_policy` preserve
the released expression-restriction policy and its exact protein/background
assets. The new table explicitly retains `Expression category` and `Expression
forbidden`; it does not export a misleading new whole-protein forbidden set.
Clinical evidence for one peptide/receptor cannot remove every peptide of that
protein from a normal-expression background. Downstream target development
must carry the evidence scope and adjudicate the intended peptide/product.

The near-self policy boundary is tracked in
[#640](https://github.com/pirl-unc/oncoref/issues/640); the verified sequence,
occurrence and coverage-gap loader is tracked in
[#641](https://github.com/pirl-unc/oncoref/issues/641). This research selection
does not implement that loader or validate HLA/similarity thresholds. Consumers
must not subtract the new candidate list from the forbidden reference, reapply
an older Medium/High-IHC filter, or treat Neither as absent normal expression.

## Literature work (#638)

[Per-gene review status and remaining work](literature-review-status.csv) covers
all 2,540 mapped nominations. Search leads, prior specific assessments, completed
antibody dispositions and the new clinical comparisons remain distinct.
[Unmapped identity reviews](unmapped-identity-reviews.csv) retain BAGE and GAGE2B
outside that denominator, without substituting a paralog or inferring absence.
**#638 is not complete.** The per-gene table does not substitute for reading the
remaining primary papers or expanding truncated searches. No-hit searches are
not evidence of absence. The initial scripts/tables were already published in
PR #639; the older issue's statement that they remained local is stale.

The new retrieval mode supports cursor pagination with per-page byte receipts
and gene-context variants for ambiguous IMPACT, SPATIAL, FATE and MRP8 terms.
The original ledger reported 2,491,175 hits across 1,688 truncated queries;
these counts include substantial name ambiguity. The selected expansion covers
{summary["expanded_queries"]} queries for {summary["expanded_genes"]} clinical/ambiguity cases;
{summary["expanded_queries_paginated_completely"]} have complete pagination of those
query variants. [Search ledger](expanded-searches.csv) and [paper leads](expanded-papers.csv)
remain retrieval evidence, not paper adjudication. The variant searches can
lose ambiguous aliases and do not claim systematic completeness; the original
queries remain available for sensitivity review. [Primary-paper receipts](paper-receipts.json)
record the retrieved versions used for the targeted clinical review.

The [breadth report](../normal-breadth-20261026/index.md) supplies descriptive
profiles and missing-coverage bounds without changing memberships. Fetal and
embryonic targets need the separate [developmental-source assessment](../../developmental-targets.md).

## Reproduction

Run `.venv/bin/python scripts/build_cta_research_policy.py --out <empty-directory>
--expansion docs/audits/cta-research-20261026
--paper-receipts docs/audits/cta-research-20261026/paper-receipts.json`.
`--publish-dataset` copies the reviewed research selection to the package;
historical expression decisions and immutable reports remain unchanged.
""")
    inputs = [
        base / "decision-facts.csv.gz",
        base / "candidate-results.csv.gz",
        ROOT / "oncoref/data/cta-policy-findings.csv",
        ROOT / "oncoref/data/cta-research-clinical-evidence.csv",
        identity_reviews,
        old / "all-genes.csv",
        Path(__file__),
        ROOT / "oncoref/cta_research_policy.py",
    ]
    (out / "run-manifest.json").write_text(
        json.dumps(
            {
                "inputs": {str(p.resolve().relative_to(ROOT)): digest(p) for p in inputs},
                "retrieval_inputs": {p.name: digest(p) for p in extra_inputs},
                "outputs": {p.name: digest(p) for p in sorted(out.iterdir())},
            },
            indent=2,
        )
        + "\n"
    )
    if publish:
        shutil.copyfile(
            out / "gene-results.csv.gz", ROOT / "oncoref/data/cta-research-gene-sets.csv.gz"
        )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--publish-dataset", action="store_true")
    parser.add_argument("--expansion", type=Path, required=True)
    parser.add_argument("--paper-receipts", type=Path, required=True)
    args = parser.parse_args()
    build(
        args.out,
        publish=args.publish_dataset,
        expansion=args.expansion,
        paper_receipts=args.paper_receipts,
    )
