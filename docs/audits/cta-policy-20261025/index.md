# CTA definitions — 25 October 2026

**213 CTA genes**, 2179 Neither and 17705 Unsafe across 20,097 canonical human genes.
Among 2,540 mapped source nominations: **213 CTA / 741 Neither / 1586 Unsafe**.

[Plain-language rules and complete threshold table](../../target-sets.md) · [All candidates](candidate-results.csv.gz) · [All genes](gene-results.csv.gz) · [CTA](cta-genes.csv.gz) · [Neither](neither-genes.csv.gz) · [Unsafe](unsafe-genes.csv.gz)

The updated policy applies reviewed literature as biological findings, not gene-name overrides. The unchanged measurement facts plus the previously accepted findings reproduce all 20,097 old categories exactly. New reviewed findings change 186 assignments from the released 399-gene CTA set (6 to Unsafe, 180 to Neither). RNA thresholds and tissue mappings are unchanged.

## Named targets — retained and excluded

| Gene | Category | Detail |
| --- | --- | --- |
| MAGEC2 | Neither | Conflict |
| CTCFL | Neither | Conflict |
| MAGEA10 | Neither | Conflict |
| AKAP4 | Neither | Conflict |
| MAGEA4 | Neither | Conflict |
| MAGEC1 | Neither | Conflict |
| SSX5 | CTA | Supported |
| CT55 | CTA | Supported |
| SAGE1 | Neither | Conflict |
| CTAG1B | CTA | Pass |
| PRAME | Unsafe | Normal |
| MAGEA1 | Neither | Conflict |
| CT83 | CTA | Pass |
| XAGE1A | CTA | Pass |
| XAGE1B | CTA | Pass |
| GAGE1 | CTA | Supported |
| SPANXB1 | CTA | Pass |
| SSX2 | CTA | Pass |
| SSX4 | CTA | Supported |
| CTAG1A | CTA | Pass |
| DAZL | Neither | Cellular |
| TDRD1 | Neither | Cellular |
| CTAG2 | Neither | Variable |
| BRDT | Neither | Cellular |
| HORMAD1 | Neither | Cellular |
| FOXR2 | Neither | Antibody |
| CTAGE1 | Neither | Antibody |
| ABCB5 | Unsafe | Normal |
| SPAG9 | Unsafe | Protein |
| DEPDC1 | Unsafe | Protein |
| SPA17 | Unsafe | Protein |
| ROPN1 | Unsafe | High RNA |
| TRIP13 | Unsafe | High RNA |
| BIRC5 | Unsafe | Protein |
| ACRBP | Unsafe | High RNA |
| TTK | Unsafe | Protein |
| KIF20A | Unsafe | Protein |
| IL5 | Unsafe | Normal |
| ROPN1B | Unsafe | High RNA |
| TEX101 | Unsafe | High RNA |
| IGF2BP3 | Unsafe | Protein |
| KIF20B | Unsafe | Protein |
| NUF2 | Unsafe | High RNA |
| LY6K | Unsafe | High RNA |
| PTTG1 | Unsafe | Protein |
| IL13 | Unsafe | Normal |
| PAPPA | Unsafe | Normal |
| CLDN6 | Unsafe | High RNA |
| MAGEA12 | Unsafe | Safety |
| MAGEA3 | Unsafe | High RNA |
| GAGE12F | Unsafe | High RNA |
| PAGE2B | Unsafe | High RNA |

Every changed gene has a [review record](reviewed-findings.csv) with source, species, reference, assay, identity/specificity judgement, review depth, limitations and linked issue. Abstract-level decisions are marked. Unsafe means exclusion from unrestricted gene-level CTA targeting, not demonstrated clinical injury. Neither conflicts have a completed curation decision; biological uncertainty remains. They need evidence that resolves the conflict before reconsideration.

## Evidence and follow-up

[Decision facts](decision-facts.csv.gz) · [Reviewed literature](reviewed-findings.csv) · [Changed rows](changes.csv) · [Follow-up questions](follow-up.csv) · [Named CTA checks](important-ctas.csv) · [Literature search and review coverage](../cta-literature-20261025/index.md)

CT55 remains CTA through affirmative reproductive and tumor-recognition evidence below the RNA floor; normal gates still apply. CTAG1A/B and XAGE1A/B remain CTA. CTAG2 remains Neither / Variable and FOXR2 Neither / Antibody. MAGEA4 is now Neither / Conflict; PRAME is Unsafe because an independent human adrenal-protein observation is accepted. This does not relabel their selected therapeutic products as clinically unsafe.

[All 178 completed antibody adjudications, evidence and reconsideration requirements](../cta-adjudication-20261025/index.md). PLCZ1's pancreatic paper measures its RNA, while its immunostaining detects CFTR; the erroneous protein interpretation is rejected and the qualitative RNA finding is held separately.

The literature audit is incomplete: 32 genes received new primary-paper assessment, with 46 specific notes including prior reviews. All **178** raw-antibody discrepancies remaining in the local 391-gene CTA draft received documented adjudications: **177 Neither and 1 Unsafe**. The two other genes in the original 180-gene queue were already addressed by the eight-gene review. **0 unresolved raw-antibody cases remain retained.** Review completion does not mean the molecular identity of every stain was resolved. The retained CTA set remains a research selection, not validated clinical safety. No new continuous breadth gate or calibrated RNA uncertainty model is claimed.

No remaining CTA or Neither has an audited non-permitted/unsited TPM or nTPM observation >30. Missing observations remain gaps; this is not proof of absence in unmeasured cells. Original source versions, native units and donor limitations are preserved in the [RNA audit](../cta-neither-review-20261024/index.md).

## Exact proteins

Complete exact-sequence groups: **359 CTA / 3,230 Neither / 126,396 Unsafe**. Priority for a protein encoded by multiple genes is Unsafe → Neither → CTA. The additional 1,721 reference loci retain their measured decisions, reissued under this policy. Neither is excluded from the default forbidden union.

[Sequence recipe](exclusion-basis.json) · [Reference](reference.json) · [Release asset inventory and SHA-256](protein-assets.json) · [Additional loci](sequence-gene-reviews.csv.gz)

Use the matching **[v1.20.0 release assets](https://github.com/pirl-unc/oncoref/releases/tag/v1.20.0)** and their recorded checksums. Earlier FASTAs retain their dated policies and must not be combined with the new gene list. The historical 624-gene CTA API is unchanged; downstream applications must explicitly choose the current `target_sets` policy.

## Reproduction

```sh
.venv/bin/python scripts/build_cta_policy.py --output /tmp/cta-policy-report --assets /tmp/cta-policy-proteins --sequence-input /path/to/pinned-original-sequence-assets
```

Use empty output directories. Add `--publish-dataset` only to update the repository's packaged current data after all invariants pass. The builder verifies pinned input and original sequence hashes, replays the baseline, evaluates every reason, and checks the >30 invariant before publishing.

[Policy](policy.json) · [Counts](summary.json) · [Input/code/output hashes](run-manifest.json). Historical reports remain frozen.
