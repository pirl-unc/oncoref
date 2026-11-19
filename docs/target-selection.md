# Cancer and Reproductive Antigens

**CTA = CA ∩ RA. Candidate antigen proteins = CA ∪ RA.** Cancer-only proteins remain
Cancer Antigens. These are research screening memberships, not demonstrated antigen
presentation or clinical safety.

Catalog `oncoref.target-selection.2026-11-19.1`: **327 CTA** · **70 Cancer Antigen only** · **127 Reproductive Antigen only** · **1,116 Other** · **18,253 Unsafe**.

[Every protein and its evidence](audits/target-selection-20261110/rebuild/report.html) ·
[Funnel and Venn plots](audits/target-selection-20261110/rebuild/summary.md) ·
[Host choices](host-selection.md) · [Optional filters](target-filters.md)

[Pre-GTEx comparison for 78 CA-only proteins and 31 normal-cell holds](audits/antigen-evidence-review-20261118/index.md).

[Individual GTEx samples and symmetric reassessment of the ten testis leads](audits/gtex-reproductive-witnesses/index.md).

## Three separate questions

1. **Shared evidence:** Does the protein pass the normal-expression, coverage,
   identity, translation, sequence and expression-bounds checks?
2. **Biological membership:** Does cancer evidence qualify (CA), reproductive or
   developmental evidence qualify (RA), or both (CTA)?
3. **Host eligibility:** Which additional compartments must be protected for this
   research scenario, and is the evidence sufficient to assess them?

![Catalog and host selection](audits/target-selection-20261110/rebuild/summary/source-to-tiers-architecture.png)

## One threshold, two compartments

Both CA and RA independently require **≥10 × H₀**, where
**H₀ = max(non-reproductive normal N, Disposable D, 1 TPM)**, after shared eligibility.

| Set | Positive RNA evidence | Proteins |
|---|---|---:|
| CA | At least one QC-passing tumor sample ≥10 × H₀ | 397 |
| RA | Reproductive sample or supported embryo-lineage mean ≥10 × H₀ | 454 |
| CTA | Both routes | 327 |
| CA∪RA | Either route | 524 |

GTEx v11 contributes the **linearly interpolated donor P95 within each non-reproductive tissue site**. Complete member-locus TPMs are summed within a library first. Reproductive sites, including breast and prostate, retain individual samples as RA witnesses. HPA v23 contributes the within-tissue biopsy P95 and the within-population sorted immune donor P95. DICE contributes a P95 for single loci. Its ranks lack specimen IDs, so identical-protein groups with multiple loci use the sum of measured member maxima as an upper bound. GTEx medians, HPA tissue and brain means, and Monaco cell means retain their published statistics. The highest applicable value across these sources sets N or D. Single-cell and single-nucleus counts, nTPM and protein staining have separate normal-expression gates.

[All RNA sources, tissue coverage and set definitions](audits/target-selection-20261110/rebuild/normal-rna/index.md)

![Normal and reproductive TPM sources by tissue](audits/target-selection-20261110/rebuild/normal-rna/tpm-coverage.png)

A small cohort can have a P95 close to its maximum. The underlying measurements and sample QC remain available. [Download the source data and literature controls](antigen-data.md).

Missing required normal coverage does not become zero.
[Source receipt](audits/target-selection-20261110/rebuild/comparison.json).

Developmental means are calculated within embryo and lineage. A well-sampled
population requires ≥20 measured cells, ≥5 cells with ≥1 TPM and ≥10% positive cells.
Smaller populations require ≥5 measured cells, ≥50% positive cells and two independent
embryos of the same age and lineage reaching the tenfold threshold. Cells or lineages
from one embryo are not independent donors. These positive-support requirements do
**not** determine whether an observation can inform normal host exposure.

## Shared eligibility and Other

| Evidence or gap | Classification |
|---|---|
| Accepted normal RNA, normal protein or applicable clinical exclusion | Unsafe |
| Missing coverage, uncertain identity/translation or unresolved normal finding | Other:Review |
| Unsupported nested protein, ≥95% Unsafe 9-mer containment or read-through locus | Other:Sequence |
| Unmeasured home tissue without qualifying reproductive evidence | Other:Unscreened tissue |
| Neither enrichment route qualifies, or expression bounds fail | Other:Expression |

Insufficient coverage of a documented normal home blocks the CA route even when RA qualifies. Unsafe takes precedence, then Other:Review. Every failed gate remains on the row.
**Other describes unresolved or nonqualifying proteins, not absent normal expression.**
Disposable expression has separate 32 TPM and 32 nTPM bounds in the catalog; these
are not the limits for host-protected organs.

## Rare normal-cell RNA

Single-cell and single-nucleus measurements can identify populations diluted in bulk
RNA. Donor-resolved UMI evidence is graded independently of cluster nTPM and bulk TPM.

| Check | Rule |
|---|---|
| Signal | >= 1% of a donor's cells and >= 10 cells carry >= 2 molecules. |
| Replication | The signal in >= 2 independent studies of the same organ, which agree on a cell type or include >= 3 donors. Donors of one study share its libraries and reference; one donor in each of two studies, in different cell types, agrees on nothing but the organ. |
| Single study | Without replication, a signal counts only when >= 5% of a donor's cells carry multiple molecules. |
| Confirmation | Bulk RNA >= 1 native TPM in the same organ in both the HPA and GTEx sample sets, a different assay with its own read assignment. |
| Established | A replicated signal at a locus without mapping risk, or any signal with bulk confirmation. Unsafe. |
| Unresolved | A replicated signal with unresolved locus attribution, or a frequent single-study signal without independent confirmation. Other:Review pending adjudication. |
| Dismissed | A weak or inconsistent signal below replication and frequency requirements. Recorded with its reason; no category effect. |

A mapping-risk flag is not proof of read misassignment. Strong unresolved signals
remain Other:Review; low bulk expression alone cannot erase them.
[Donor records and rationale](audits/normal-cell-atlas-20261110/rare-cell/index.md).

## Protein staining

Identity-resolved, reliable normal protein evidence can independently mark a protein
Unsafe under the tissue-tier policy. Antibody reliability, specificity, cell type and
source scope remain part of each finding.

## Published normal expression

Reviewed gene-specific normal expression
retains sources, adjudication and supersession history.
[Literature rubric](audits/target-selection-20261110/rebuild/literature.md).

## Home tissue the panels do not resolve

Olfactory epithelium, nasal mucosa, taste buds, hair follicle and dental pulp are not
resolved by the native-TPM comparison panel. A protein with a documented home in one of these tissues can qualify as RA. CA and CTA require coverage of that home. A protein that reaches the cancer threshold alone is classified Other:Unscreened tissue.

## Tumor quality and replication

Library QC retains the olfactory-receptor background and unrelated-tissue screens.
Coordinated late-spermatid RNA is an **advisory note only**: normal reproductive
specificity does not establish absence from cancer or prove normal-cell admixture.
These notes never change sample QC status or remove tumor or reproductive evidence,
including when selecting only `sample_qc="pass"`. Independent QC exclusions still apply.
The [library audit](audits/target-selection-20261110/rebuild/library-qc.csv.gz) separates
blocking `reasons` from advisory `notes`; the public expression QC API exposes
`sample_qc_reasons` and `sample_qc_notes` separately.
**R** — Coordinated reproductive RNA: all six testis markers ≥0.1% of their reference; cell of origin unresolved. Advisory only.
Cohort availability and exclusions remain explicit
in the [source coverage table](audits/target-selection-20261110/rebuild/tumor-source-coverage.csv).
Catalog rows expose samples, independent patients, anonymous samples and cancer types
at tenfold. `independent_patient_tumor_folds` holds descending patient maxima; multiple
samples from one known patient count once. Unknown donors do not count as known patients.

## Biological origin and host protection

`reproductive_qualifying_compartments` names the observed origins. The legacy
`cta_tier` column mixes origin with an exposure convention and is retained only for
compatibility; it is **not a toxicity ranking**. `profile="core"` is rejected: choosing
an origin does not authorize gonadal injury.

Protect / Accept loss / Absent are separate host choices. Anatomy presets protect no
reproductive tissues by default; protection is explicit. Host rejection leaves biological membership
unchanged and adds the rejected candidate to the selection's subtraction background.
[Exact contract, evidence scope and migration](host-selection.md).

![Host choices and outcomes](audits/target-selection-20261110/rebuild/summary/cta-tiers-and-profiles.png)

## One selection API

```python
from oncoref import antigen_sets

selection = antigen_sets.select("CTA", mode="standard")
proteins = selection.proteins()
genes = selection.gene_ids()
sequences = selection.sequences()  # locally provisioned Ensembl 115
background = selection.background_index("Unsafe|Other")
```

[Sequence and peptide API](antigen-sequences.md) · [CLI and rebuild](curation.md).
The [rebuild](audits/target-selection-20261110/rebuild/index.md) records every input,
output and decision. Outstanding scientific adjudications remain explicit holds;
reproductive-content notes do not resolve malignant-cell localization.
