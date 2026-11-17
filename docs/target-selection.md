# Cancer and Reproductive Antigens

**CTA = CA ∩ RA. Candidate antigen proteins = CA ∪ RA.** Cancer-only proteins remain
Cancer Antigens; calling them CTAs would erase the reproductive distinction.

Policy `oncoref.target-selection.2026-11-17.1`: **290 CTA** · **106 Cancer Antigen only** · **113 Reproductive Antigen only** · **1,131 Other** · **18,253 Unsafe**.
These are research screening memberships, not demonstrated antigen presentation or safety.

[Browse every protein](audits/target-selection-20261110/rebuild/report.html) ·
[Funnel, Venn plots and counts](audits/target-selection-20261110/rebuild/summary.md) ·
[Machine-readable rules](audits/target-selection-20261110/definition.json)

![Antigen membership](audits/target-selection-20261110/rebuild/summary/antigen-sets.png)

## One threshold, two compartments

Let **H = max(N, D, 1 TPM)**: N is the highest supported non-reproductive normal
measurement and D the highest Disposable-tissue measurement. Both sets require
**at least 10 × H**, independently, after the same eligibility screen.

| Set | Required positive RNA evidence |
|---|---|
| CA | At least one QC-passing tumor sample ≥10 × H |
| RA | At least one reproductive tissue sample, or supported embryo-specific lineage mean, ≥10 × H |
| CTA | Both CA and RA |
| CA∪RA | Either set; the union does not imply tumor expression for RA-only proteins |

All comparisons use the declared native-TPM panel. Qualitative protein staining,
single-cell nTPM and UMI counts cannot be divided by TPM or bypass its numerical fold.
Reproductive IHC and older association routes remain visible as supporting annotations.

The normal reference includes **individual HPA v23 tissue samples from all measured
compartments**, plus tissue/subregion means, GTEx donor medians and supported immune
summaries. Reproductive observations are positive candidates; non-reproductive
observations tighten the denominator. It is still not a maximum over every healthy
donor, and cross-study batch effects are not corrected. Missing required normal
coverage never becomes zero. [Source receipt and limits](audits/target-selection-20261110/rebuild/comparison.json).

Developmental means are computed within one embryo and lineage. One well-sampled
population requires ≥20 cells, ≥5 cells with ≥1 TPM and ≥10% positive cells.
Smaller populations require ≥5 measured cells, ≥50% with ≥1 TPM, and two independent
embryos at the same age and lineage reaching 10 × H. The mean itself must reach
10 × H. Higher RA fold filters retain the replication requirement through
`reproductive_supported_tpm`. These are screening support conventions, not biological
limits on embryo size or lineage abundance. Embryos are
biological donors; multiple cells or lineages from one embryo are not independent donors.
The developmental scope is preimplantation and Carnegie-stage-7 embryos, not all fetal stages.

![Both enrichment comparisons](audits/target-selection-20261110/rebuild/summary/tumor-normal.png)

## Shared eligibility and Other

| Layer | Consequence |
|---|---|
| Accepted normal RNA, normal protein or applicable clinical exclusion | Unsafe |
| Missing coverage, uncertain identity/translation or unresolved normal evidence | Other:Review |
| Unsupported nested protein, ≥95% Unsafe 9-mer containment or read-through locus | Other:Sequence |
| Unmeasured normal home without qualifying reproductive evidence | Other:Tissue |
| Neither enrichment compartment qualifies, or expression bounds fail | Other:Expression |

Unsafe takes precedence. Other:Review then precedes the other Other subtypes.
Every failed check remains on the row, even when another determines its category.
**Other is neither an accepted hazard nor a claim of absent normal expression.**
Some Other assignments are technical; others reflect missing biological evidence.

Disposable expression retains separate ceilings of 32 native TPM and 32 nTPM.
The target maximum must exceed all measured non-reproductive contexts. Tissue
allowances apply to resolved cell compartments: immune or vascular cells inside a
reproductive organ do not automatically inherit that organ's allowance.

## Consequences of tissue loss

Reproductive membership and the consequences of losing the expressing tissue are
separate questions. Broad RA includes gonadal, developmental and accessory reproductive
contexts. Breast or prostate evidence is retained with its actual tissue context.
The [tissue registry](tissue-tiers.md) records **acute loss, lasting loss, repair and
management**. Whole testis and ovary retain Protected assessments; mammary tissue has
a conditional Manageable assessment. Bulk organ measurements do not establish which
cell population expresses a protein or whether only that population would be injured.

```python
from oncoref import tissue_tiers
loss = tissue_tiers.consequences("testis")
loss["persistent_loss"], loss["management"], loss["evidence"]
# Tissue loss alongside this protein's observed reproductive expression:
from oncoref import antigen_sets
antigen_sets.loss_consequences("CTAG2")
```

These records describe functional consequences, not receptor-specific probabilities.
Sex labels, presumed immune privilege, negative bulk RNA and negative IHC do not by
themselves demonstrate absence of an off-target compartment or dismiss a sparse signal.
Actual anatomy, prior removal and pregnancy/lactation state require explicit context.

## Rare normal-cell RNA

Single-cell and single-nucleus data are independent normal-expression layers.
They can identify a normal population diluted in bulk tissue. The donor-resolved
UMI rubric is distinct from HPA cell-cluster nTPM and bulk TPM.

| Check | Rule |
|---|---|
| Signal | >= 1% of a donor's cells and >= 10 cells carry >= 2 molecules. Protected populations need 2 times these frequencies. |
| Replication | The signal in >= 2 independent studies of the same organ, which agree on a cell type or include >= 3 donors. Donors of one study share its libraries and reference; one donor in each of two studies, in different cell types, agrees on nothing but the organ. |
| Frequent single-study signal | Without replication, a signal counts only when >= 5% of a donor's cells carry multiple molecules. |
| Orthogonal confirmation | Bulk RNA >= 1 native TPM in the same organ in both the HPA and GTEx sample sets, a different assay with its own read assignment. |
| Established | A replicated signal at a locus without mapping risk, or any signal with bulk confirmation. Unsafe. |
| Unresolved | A replicated signal with unresolved locus attribution, or a frequent single-study signal without independent confirmation. Other:Review pending adjudication. |
| Dismissed | A weak or inconsistent signal below replication and frequency requirements. Recorded with its reason; no category effect. |

A reference mapping-risk flag does not show where the observed reads aligned.
DCST1's replicated retinal signal and NKX2-6's frequent single-study cardiac signal
therefore require review, rather than being dismissed as artifacts. Resolving them
requires gene-discriminating read evidence or independent confirmation. The present
count records cannot supply that adjudication. No ambient-RNA correction has been
applied to these pinned counts. Donor overlap across collections remains a limitation.

Single molecules alone do not establish expression; reviewed gene-specific evidence
can provide independent confirmation. Low counts and absent findings are not proof
of safety. [Donor records, grades and reasons](audits/normal-cell-atlas-20261110/rare-cell/index.md).

## Protein staining and published normal expression

Reliable, identity-resolved normal protein findings can independently mark a protein
Unsafe under the tissue-tier policy. Antibody specificity, staining reliability,
cell type and source scope remain part of each finding. Conflicting or unresolved
assays cause review holds. Negative staining supplies only the existing limited
policy credit; it does not demonstrate absence or cancel a stronger finding.

## Published normal expression

Reviewed publications are adjudicated through the shared evidence records.
Accepted gene-specific normal expression in Vital or Protected tissues excludes;
other observations retain their tissue, scope and rationale. Findings can supersede
earlier findings without deleting history. MAGEC2's HPA062230/GTEx adjudication and
the DAZ gastric-RNA/paralog limitations use this same format.
[Literature rubric and sources](audits/target-selection-20261110/rebuild/literature.md).

## Home tissue the panels do not resolve

Olfactory epithelium, nasal mucosa, taste buds, hair follicle and dental pulp are not
resolved by the comparison panel. A recorded home in one of these tissues yields
Other:Tissue unless reproductive evidence qualifies. The home remains recorded
even when RA qualifies; a candidate label does not establish a measured margin there.
The optional normal-function filter can tighten that decision.

## Tumor quality and replication

Passing TGCT libraries contribute normally. The shared QC screen excludes libraries
with broad olfactory-receptor background, a coordinated late-spermatid signature,
or multiple unrelated tissue signatures. These are screening signatures, not proof
of gDNA contamination, normal-tissue fraction or a tumor activation mechanism.
[Every library measurement](audits/target-selection-20261110/rebuild/library-qc.csv.gz).

Catalog rows expose `samples_at_tenfold`, `independent_patients_at_tenfold`,
`anonymous_samples_at_tenfold` and `cancer_types_at_tenfold`. Unknown donors never
count as known independent patients. Counts refer to the catalog's 10× threshold.
`independent_patient_tumor_folds` stores each qualifying known patient's maximum
fold in descending order; filters use the k-th value to require k patients at
any selected fold ≥10. Legacy `*_above_normal`
columns count a weaker threshold and are descriptive only. Every tenfold witness
is downloadable, with source, sample, participant and expression value.

## API and peptide subtraction

```python
from oncoref import antigen_sets, target_filters, target_selection

ctas = antigen_sets.proteins("CTA")
ca = antigen_sets.proteins("CA")
ra = antigen_sets.proteins("RA")
union = antigen_sets.proteins("CA|RA")
background = antigen_sets.proteins("Unsafe|Other")
record = target_selection.results("MAGEC2")
replicated_ca = target_filters.select(target_filters.TargetFilter(
    membership="cancer_antigen", min_tumor_patients=2))
```

`Unsafe` and `Other` are also available separately. All annotated translations can
feed exact 8–15mer or optional bounded-distance subtraction, with occurrences and
sequence gaps retained. [Sequence API](antigen-sequences.md) · [Optional filters](target-filters.md).
The strict MAGE-prefix rule remains an explicit user-selected judgment; sequence
similarity screening complements it and does not establish TCR cross-reactivity.

## Reproduction and review status

The [rebuild](audits/target-selection-20261110/rebuild/index.md) publishes inputs,
receipts, every assignment and both sequential and sole-reason attribution.
`release_ready: false` describes outstanding scientific adjudications, including
read-level mapping review and the hair-follicle replay (#687). It is not the software
publication flag. Open holds stay excluded from positive sets.

Historical catalogs remain dated snapshots. The [1.38 threshold audit](audits/cta-threshold-review-20261116/index.md)
is not regenerated from today's membership rules; its generator rejects a mismatched policy.
