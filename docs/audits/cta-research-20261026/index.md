# Recalibrated CTA research selection

**215 CTA research candidates** under `oncoref.cta-research-selection.2026-10-26.1`.
[Complete list](cta-genes.csv) · [All nominations](candidate-results.csv.gz) ·
[Every change](changes.csv) · [Reviewed clinical comparisons](clinical-evidence.csv).

The user requested a common evidence policy calibrated to include PRAME and
MAGE-A4 alongside genes of comparable quality. Both already pass the numerical
RNA gates. This revision changes the evidence threshold for research selection,
not the RNA measurements or cutoffs. NY-ESO-1 remains represented by CTAG1A/B.
This is benchmark-driven calibration, not independent prospective validation.

| Gene | Earlier expression category | Research selection |
|---|---|---|
| MAGEA4 | Neither | CTA |
| PRAME | Unsafe | CTA |

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
27 queries for 10 clinical/ambiguity cases;
27 have complete pagination of those
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
