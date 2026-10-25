# XAGE5: historical rescue and current antibody hold

Reviewed 24 October 2026 against OncoRef 1.19.1 and its pinned HPA v23 inputs.
XAGE5 (`ENSG00000171405`) remains **Neither / Antibody** in the current target
sets: **Target = false; Forbidden = false**. No gene-specific rescue or new
exclusion is introduced by this review.

## What happened after PR #59

[PR #59](https://github.com/pirl-unc/oncoref/pull/59) left the one-gene
`MANUALLY_EXPRESSED_CTA` exception as an open question.
[PR #60](https://github.com/pirl-unc/oncoref/pull/60), merged on 12 June 2026,
resolved that question by replacing the XAGE5 exception with a general
low-expression rule. The historical expression accessors retain candidates with
MODERATE restriction confidence and STRICT reproductive RNA, subject to their
other nomination and observed-expression requirements. XAGE5 qualifies at
1.1 nTPM in testis, as do eligible peers with the same evidence signature.

That historical membership is not the current CTA/Unsafe/Neither decision.
`cta.cta_gene_names()` still includes XAGE5; `target_sets.gene_names("CTA")`
does not. The general rescue does not bypass the current admission or
normal-cell review rules. Use `target_sets` for current target and forbidden
sets, as described in the [target-set guide](target-sets.md).

## Recorded expression

| Evidence | XAGE5 result | Interpretation |
|---|---|---|
| HPA v23 consensus bulk RNA | Testis 1.1 nTPM; all 49 other measured tissues 0 nTPM | Low reproductive RNA; reported zeros are source measurements |
| Maximum audited non-permitted RNA, TPM sources | 2.2666887504 TPM | Below the hard >30 TPM exclusion |
| Maximum audited non-permitted RNA, nTPM sources | 0.4 nTPM | Below the hard >30 nTPM exclusion |
| Curated HPA normal-tissue IHC | No records | Missing coverage, not negative staining |
| Primary antibody HPA052507 | 78 tissue/cell records; profile and antibody validation both Uncertain | Cannot supply reliable reproductive protein admission evidence |
| Primary staining in Vital/Protected/Manageable compartments | 12 positive records among 48 assessed records | Triggers the existing uncertain-antibody review rule |

The RNA maxima keep TPM and nTPM separate. They summarize the measured
contexts; they do not establish absence in every normal cell or donor.
The IHC counts are antibody/tissue/cell aggregates, not donor-positive fractions.

All 12 non-permitted positive primary records use HPA052507:

| Tissue | Cells | Staining |
|---|---|---|
| Caudate | Neuronal cells | Medium |
| Cerebral cortex | Neuronal cells | Medium |
| Hippocampus | Neuronal cells | Low |
| Duodenum | Glandular cells | Medium |
| Small intestine | Glandular cells | Medium |
| Liver | Hepatocytes | Low |
| Lung | Macrophages | Low |
| Lymph node | Non-germinal center cells | Medium |
| Pancreas | Exocrine glandular cells | Medium |
| Stomach 1 | Glandular cells | Low |
| Stomach 2 | Glandular cells | Low |
| Thyroid gland | Glandular cells | Medium |

The same uncertain antibody reports Low staining in testicular seminiferous
ducts and High staining in placental trophoblasts. Those findings do not
validate the antibody's specificity or resolve its normal-tissue staining.

Sources: the pinned [HPA v23 bulk RNA archive](https://v23.proteinatlas.org/download/rna_tissue_consensus.tsv.zip),
[curated IHC archive](https://v23.proteinatlas.org/download/normal_tissue.tsv.zip),
and [primary antibody XML archive](https://v23.proteinatlas.org/download/proteinatlas.xml.gz).
The reproducible extracted observations are in the
[primary-antibody table](audits/cta-normal-cells-20261023/primary-antibodies.csv.gz)
and [RNA context table](audits/cta-neither-review-20261024/rna-context-maxima.csv.gz).
The [current gene table](audits/cta-neither-review-20261024/gene-results.csv.gz)
records the decision and evidence flags.

## Why it remains Neither

The current screen has two unresolved requirements for XAGE5:

1. Its 1.1 nTPM reproductive bulk RNA is below the 5 nTPM admission floor.
   The reviewed evidence table currently supplies neither reliable positive
   reproductive protein nor the affirmative reproductive-plus-tumor evidence
   needed to bypass that floor. Published RNA nominations remain in intake;
   this is not a claim that XAGE5 lacks cancer literature.
2. Positive normal staining from an Uncertain antibody independently requires
   review. Even sufficient reproductive/tumor admission evidence would not
   remove that hold.

The current Cell RNA, Immune RNA and Corneal RNA review flags are false. No
audited non-permitted RNA value exceeds 30 TPM or nTPM. The unresolved staining
therefore supports **Neither / Antibody**, without asserting confirmed normal
protein production, clinical toxicity or inclusion in forbidden proteins.

## Evidence needed for reconsideration

The open question is whether HPA052507's neuronal and other normal staining is
specific to XAGE5. Reconsideration needs independently supported protein
identity in those compartments and reliable reproductive protein evidence or
reviewed reproductive-plus-tumor evidence. Resolving only the low RNA admission
floor would leave the normal-antibody hold in place. Any later evidence should
be recorded with its assay, tissue, source and limitations and evaluated by the
same rules used for other candidates.

```python
from oncoref import cta, target_sets

assert "XAGE5" in cta.cta_gene_names()  # Historical expression definition.
assert "XAGE5" not in target_sets.gene_names("CTA")
assert "XAGE5" in target_sets.gene_names("Neither")
assert "XAGE5" not in target_sets.gene_names("Unsafe")

row = target_sets.results().set_index("Gene").loc["XAGE5"]
assert row["Detail"] == "Antibody"
assert not row["Target"] and not row["Forbidden"]
```
