#!/usr/bin/env python3
"""Render completed antibody adjudications from pinned observations and findings."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "docs/audits/cta-adjudication-20261025"


def render():
    cases = pd.read_csv(ROOT / "oncoref/data/cta-antibody-adjudications.csv").fillna("")
    findings = pd.read_csv(ROOT / "oncoref/data/cta-policy-findings.csv").fillna("")
    current = pd.read_csv(ROOT / "oncoref/data/target-gene-sets.csv.gz", low_memory=False)
    assert len(cases) == 178 and cases["Review complete"].all()
    assert (
        cases.merge(current, on="Gene ID", validate="one_to_one").eval("Decision == Category").all()
    )
    lines = [
        """# Completed CTA adjudications — 25 October 2026

**All 178 antibody cases remaining in the local 391-gene draft have a completed decision: 177 Neither and 1 Unsafe.** No member of this discrepancy queue remains CTA. The released list has **213 CTA genes**.

[Current rules and funnel](../cta-policy-20261025/index.md) · [All case records](adjudications.csv) · [Raw positive observations](primary-observations.csv.gz) · [Reviewed findings](../cta-policy-20261025/reviewed-findings.csv)

## Decision boundary

- Accepted endogenous protein in non-permitted normal cells excludes unrestricted gene-level use: **Unsafe**.
- A material raw normal stain whose molecular specificity is unresolved gives **Neither**. Low intensity is not, by itself, proof of nonspecific staining. A generic off-target annotation does not identify which specimens are false positives.
- Curated negative IHC and passing RNA remain useful evidence, but do not silently erase a conflicting raw stain. Enhanced profile reliability does not validate every individual primary stain.
- A clinical product recognizes a particular peptide/HLA/receptor combination. Product-specific efficacy or tolerability does not clear every protein or sequence encoded by its target gene.

These are conservative research classifications. Neither is not proven toxicity; Unsafe is not a claim that an existing therapy caused injury. Future changes require source-specific evidence revisions, not an allowlist of familiar genes.

## What was actually reviewed

The 178-case review covers the pinned mapped positive primary-antibody rows, their named normal cells, intensity, primary reagents, HPA profile notes and staining summaries, with the RNA context preserved. This is a review of annotations and reagent evidence, **not an independent histopathology reread of every image** or new experimental validation. All 178 decisions are complete; molecular specificity remains unresolved in the raw-antibody findings. PRAME has additional accepted independent normal-protein evidence.

The eight initial cases were rechecked against their primary studies at the access depth recorded in the findings table. Full original methods were available for FTMT, YBX2, ACTRT1 and PLCZ1; other records explicitly retain abstract/excerpt and reagent-access limitations. A search result is not a full-paper adjudication. The broader 2,540-gene literature review and a continuous breadth metric remain open work.

## Eight initial cases and the important clinical targets

| Gene | Decision | Adjudication and limits |
|---|---|---|
| FTMT | Unsafe | Normal human nigral neurons with antibody specificity controls; small postmortem series. [Study](https://pubmed.ncbi.nlm.nih.gov/28386150/) |
| ABCC11 | Unsafe | Reported normal neural protein with multimodal support. Full reagent methods remain an access limitation. [Study](https://pubmed.ncbi.nlm.nih.gov/16359813/) |
| YBX2 | Unsafe | Vascular smooth-muscle/adrenal protein; not evidence of cardiomyocyte expression. Polyclonal and specimen limitations remain. [Study](https://pubmed.ncbi.nlm.nih.gov/16479255/) |
| CCDC158 | Unsafe | Reported normal renal tubular protein/function. Renal versus testicular isoforms remain unresolved; no isoform-specific clearance is claimed. [Study](https://pubmed.ncbi.nlm.nih.gov/39319391/) |
| ACTRT1 | Unsafe | Endogenous protein in differentiating normal human keratinocytes. Cultured-cell state is explicit; an RPE cell line is not normal retinal tissue evidence. [Study](https://pubmed.ncbi.nlm.nih.gov/33972689/) |
| MAGEC2 | Neither | Purkinje staining with TC-1 requires gene/family specificity resolution. [Study](https://pubmed.ncbi.nlm.nih.gov/12920247/) |
| LGALS13 | Neither | Historical extra-placental protein reactivity versus placental RNA restriction; galectin identity and anatomy are unresolved. [Early study](https://pubmed.ncbi.nlm.nih.gov/10527825/), [placental study](https://pubmed.ncbi.nlm.nih.gov/19497882/) |
| PLCZ1 | Neither | **Assay correction:** NYD-SP27 was measured by RT-PCR; the paper's pancreatic immunostaining measured CFTR. Reject the old PLCZ1 protein interpretation. Qualitative pancreatic RNA and isoform identity still require quantitative reconciliation. [Original paper](https://www.cuhk.edu.hk/ecbrc/cliparts/Publication%20list/Lancet2003.pdf) |
| PRAME | Unsafe | Independent monoclonal work reports low normal adrenal-cortex protein. A later negative panel using other reagents does not negate it. The initial positive study was reviewed at abstract level; clinical injury is not established. [Positive study](https://pubmed.ncbi.nlm.nih.gov/33890816/), [negative panel](https://pmc.ncbi.nlm.nih.gov/articles/PMC9588667/) |
| MAGEA4 | Neither | HPA normal neuronal/epidermal stains remain unresolved. An independent OTI1F9 assay has reassuring normal-panel data but does not directly resolve the older HPA reagents. Tumor diagnostic assay performance is not whole-gene safety. [Validation](https://pmc.ncbi.nlm.nih.gov/articles/PMC11170170/) |

**Preserved therapeutic evidence:** [PRAME/IMA203](https://pubmed.ncbi.nlm.nih.gov/40205198/) and [MAGE-A4/afami-cel](https://pubmed.ncbi.nlm.nih.gov/38554725/) remain product-specific clinical evidence. Neither publication is discarded, and neither supplies an unrestricted gene-level rescue. CTAG1A/B and XAGE1A/B remain CTA. CTAG2 remains Neither/Variable, with the donor evidence preserved.

## Every antibody case

Each heading is searchable by gene. The CSV carries exact IDs, source-native RNA maxima, reagent summaries and all contexts. “Evidence resolved = false” is a completed hold decision, not an unfinished review.
"""
    ]
    for _, row in cases.sort_values("Gene").iterrows():
        gid = row["Gene ID"]
        fs = findings.loc[findings["Gene ID"].eq(gid)]
        lines += [
            f"### {row['Gene']} — {row['Decision']}\n",
            f"`{gid}` · {int(row['Raw positive rows'])} positive raw rows; {int(row['Medium/High rows'])} Medium/High; {int(row['Vital Medium/High rows'])} Vital Medium/High. Reagents: {row['Antibodies']}.\n",
            f"**Contexts:** {row['Contexts']}.\n",
            f"**Decision:** {row['Rationale']}\n",
            f"**To reconsider:** {row['Reconsideration']}\n",
            f"[Pinned primary data]({row['Sources'].split(' | ')[0]}). Review records: {', '.join(fs['Review ID'])}.\n",
        ]
    lines += [
        "## Reproduce\n\nRun `python scripts/build_cta_adjudication_report.py` after rebuilding the current policy. This renders reviewed records; it does not perform new adjudication. [Checksums](run-manifest.json).\n"
    ]
    (OUT / "index.md").write_text("\n".join(lines))
    inputs = [
        ROOT / "oncoref/data/cta-antibody-adjudications.csv",
        ROOT / "oncoref/data/cta-policy-findings.csv",
        ROOT / "oncoref/data/target-gene-sets.csv.gz",
        Path(__file__),
    ]
    outputs = sorted(p for p in OUT.iterdir() if p.is_file() and p.name != "run-manifest.json")

    def sha(p):
        return hashlib.sha256(p.read_bytes()).hexdigest()

    (OUT / "run-manifest.json").write_text(
        json.dumps(
            {
                "inputs": {str(p.relative_to(ROOT)): sha(p) for p in inputs},
                "outputs": {p.name: sha(p) for p in outputs},
            },
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    render()
