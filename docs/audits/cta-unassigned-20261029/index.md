# Common coverage, Reproductive identity and unassigned proteins

**Current local replay: 522 CTA proteins, 18,031 Unsafe,
20 Neither and 1,290 unassigned.** The maximum-versus-maximum
restriction rule has been removed. Complete GTEx core coverage is required. The
HPA diffuseness gate uses **10%** for adult/placental identity.
[Chosen method, changed proteins and cutoff sensitivity](diffuseness-impact.md).

[Current curation report](../cta-source-walkthrough-20261028/index.md) ·
[Paper titles and database versions](../cta-source-walkthrough-20261028/references.md) ·
[Every unassigned protein and its measurements](unassigned-protein-details.csv.gz) ·
[All protein scores](protein-diffuseness.csv.gz) · [Every ZNF730 calculation](znf730.md)

## One common GTEx core

The fixed [GTEx v11](https://gtexportal.org/home/downloads/adult-gtex) panel has **52 native
bulk contexts**: nine Reproductive and 43 other contexts. The latter represent 23
anatomical units after subdivisions share an organ's weight. Require a measured
value for **every member locus in every context**. A measured zero counts as coverage;
a missing or ambiguous identity does not. Read-through and retired identifiers need
identity repair, not substitution with a nearby gene.

**19,119 / 19,863 proteins have complete GTEx coverage.**
In this pinned matrix, all other proteins lack a complete member profile throughout
the panel; requiring 50%, 80% or 100% therefore gives the same answer. Coverage
failure excludes a protein from CTA and yields unassigned/incomplete evidence unless
independent safety evidence already establishes Unsafe. It does not establish Neither.
GTEx does not cover placenta, fetal stages, every rare cell or every Disposable tissue.

[All 52 contexts](gtex-v11-core-panel.csv) · [All GTEx profiles](gtex-v11-protein-profiles.csv.gz) ·
[23 previously retained proteins missing GTEx](gtex-coverage-losses.csv).
Removing the old maximum veto alone would add TGIF2-RAB5IF and ATF7-NPFF; both also
lack GTEx core coverage, so remain unassigned. The required HPA tissue/immune panels
remain in force. Optional DICE/cornea coverage does not shrink or expand the score's
denominator; measured counterevidence from those sources still contributes to safety.

## What the numerator and denominator mean

For each non-Reproductive context, calculate:

```text
tissue similarity = N / (N + R)
diffuseness = sum(fixed tissue weight × tissue similarity)
```

**N** is that tissue's RNA. **R** is the highest measured Reproductive RNA within
the **same source, version, unit and summary**. The numerator is normal RNA; the
denominator is that normal RNA plus the Reproductive reference. The outer average
has fixed weights summing to one. It is **not “10% of tissues”**, a count of positive
tissues, a percentage of total body RNA, or a probability of harm.

| Uniform normal profile | Score |
|---|---:|
| Every normal tissue is zero, with a positive reference | 0% |
| Every normal tissue is 1/10 of the reference | 9.09% |
| Every normal tissue equals the reference | 50% |
| Every normal tissue greatly exceeds the reference | Approaches 100% |

The **adult/placental identity boundary is >10%**. In the uniform example this means normals above
1/9 of the reference. One extreme GTEx organ contributes at most 1/23 = 4.35%; one
HPA organ at most 1/29 = 3.45%. Thus a single outlier does not veto identity, while
many comparable tissues raise the score. These weights correct panel sampling;
they are not organ sizes, cell counts or safety priorities.

Missing normals retain their fixed weight: report a lower/upper bound, never a
zero-filled or smaller-denominator score. A reference below **2 native units** gets
a descriptive value but no concentrated/diffuse verdict. This guard is a proposed
research parameter, not an assay detection limit. Zero reference gives no score.

## Consequences of diffuseness

**A supported diffuse adult/placental profile with no other qualifying route and
clear normal safety belongs in Neither, not Unsafe.** The HPA gate requires its
existing ≥5 nTPM positive reference and a complete same-source comparison.
Independent safety exclusions retain Unsafe precedence. Required coverage failures
remain unassigned. Missing optional concentration measurements remain unassessed,
without becoming a pass/fail value. [Full policy and impact](diffuseness-impact.md).

An adult comparison cannot reject supported embryonic evidence. Compatible
embryonic/adult measurements are needed to resolve its concentration. GTEx is a
separate descriptive comparison, not a second veto or a fallback denominator.
IHC affects safety but does not mathematically reduce RNA diffuseness.

## Measurements in familiar proteins

| Protein | category | source | reference_value | unit | Score (%) | review |
| --- | --- | --- | --- | --- | --- | --- |
| CTAG2 | CTA | hpa_rna_consensus | 10.1 | nTPM | 5.147 | concentrated |
| CT83 | CTA | hpa_rna_consensus | 62.2 | nTPM | 0.1494 | concentrated |
| FOXR2 | — | hpa_rna_consensus | 1.2 | nTPM | 0 | weak_reference |
| XAGE1A/B | CTA | hpa_rna_consensus | 98.8 | nTPM | 0.2247 | concentrated |
| MAGEA4 | CTA | hpa_rna_consensus | 21.5 | nTPM | 0.1515 | concentrated |
| ZNF730 | Unsafe | hpa_rna_consensus | 4.9 | nTPM | 9.532 | concentrated |
| CT55 | CTA | hpa_rna_consensus | 3.7 | nTPM | 0.334 | concentrated |
| PRAME | CTA | hpa_rna_consensus | 102.2 | nTPM | 0.304 | concentrated |
| TBXT | CTA | hpa_rna_consensus | 0.2 | nTPM | 15.15 | weak_reference |
| NY-ESO-1 | CTA | hpa_rna_consensus | 15.5 | nTPM | 0.0718 | concentrated |
| MTCP1 | Neither | hpa_rna_consensus | 5.1 | nTPM | 24.53 | diffuse |
| CTAG2 | CTA | gtex_tissue_rna | 14.33 | TPM | 0.0348 | concentrated |
| CT83 | CTA | gtex_tissue_rna | 113.4 | TPM | 0.04953 | concentrated |
| FOXR2 | — | gtex_tissue_rna | 2.671 | TPM | 0 | concentrated |
| XAGE1A/B | CTA | gtex_tissue_rna | 0.1847 | TPM | 0 | weak_reference |
| MAGEA4 | CTA | gtex_tissue_rna | 26.88 | TPM | 0 | concentrated |
| ZNF730 | Unsafe | gtex_tissue_rna | 4.914 | TPM | 1.646 | concentrated |
| CT55 | CTA | gtex_tissue_rna | 4.101 | TPM | 0.1538 | concentrated |
| PRAME | CTA | gtex_tissue_rna | 166.1 | TPM | 0.1416 | concentrated |
| TBXT | CTA | gtex_tissue_rna | 0.6035 | TPM | 5.597 | weak_reference |
| NY-ESO-1 | CTA | gtex_tissue_rna | 0.1067 | TPM | 2.062 | weak_reference |
| MTCP1 | Neither | gtex_tissue_rna | 9.746 | TPM | 16.82 | diffuse |

GTEx uses donor medians; HPA uses its consensus summary. These differ in sampling,
annotation and processing and are not independent replications (HPA incorporates
GTEx evidence). Do not average these scores or transfer their reference values.
NY-ESO-1 and XAGE1A/B illustrate why mandatory **coverage** must not silently become
a GTEx-only Reproductive expression floor. Both are well measured in GTEx, but their
GTEx Reproductive references are low. Their HPA evidence remains separately visible.

[HPA v23](https://v23.proteinatlas.org/about/download) provides a companion panel of
38 normal contexts / 29 anatomy units, including appendix, gallbladder and tonsil.
Thymus is separately **Permissive**, excluded from both sides of the identity
comparison. [All 50 context roles](hpa-v23-context-roles.csv). Disposable tissues
count in this identity assessment even though their expression
is permitted by the separate safety policy. [Panel weights](hpa-v23-normal-panel.csv) ·
[All HPA profiles](hpa-v23-protein-profiles.csv.gz).

## Reproductive is the umbrella

| Branch | Scope |
|---|---|
| Testicular | Testis and explicitly identified testicular cells |
| Other male | Other male reproductive organs and gametes |
| Ovarian | Ovary, follicles and oocytes |
| Female tract | Oviducts, uterus, cervix, vagina and other mapped tract tissues |
| Mammary | Mapped breast/mammary compartments |
| Placental | Mapped placenta and placental compartments |
| Other extraembryonic | Amnion, chorion, yolk sac, allantois and umbilical cord |
| Embryonic | Explicitly embryonic samples |
| Fetal | Explicitly fetal samples, including fetal organs |

These branches union to Reproductive. A protein can have evidence in multiple
branches. “Embryonic-only evidence” describes which measured route qualifies,
not proven absence from every other branch. Adult maternal organs do not change
category during pregnancy. The current developmental inputs cover embryonic stages;
they do not supply a comprehensive fetal panel. [Branch measurements](reproductive-branches.csv.gz).

## Open the original 617 / 649 buckets

The original release had **1,266 unassigned proteins**. Its statuses combined
observed failures with evidence gaps. They did not mean 617 had only missing data
or 649 were unsafe. Here are disjoint combinations of actual failed requirements:

### Incomplete evidence: 617 in release 1.22.0

| Failed requirements | Proteins |
| --- | --- |
| Reproductive expression + Adult restriction + Required coverage | 215 |
| Reproductive expression + Adult restriction | 193 |
| Reproductive expression + Cancer expression + Adult restriction + Required coverage | 87 |
| Required coverage | 48 |
| Reproductive expression + Adult restriction + Normal bulk | 19 |
| Reproductive expression + Adult restriction + HPA normal protein | 14 |
| HPA normal protein | 11 |
| Adult restriction + Required coverage | 11 |
| Reproductive expression + Cancer expression + Adult restriction | 5 |
| Cancer expression + Adult restriction + Required coverage | 4 |
| Reproductive expression + Adult restriction + Normal bulk + HPA tissue clusters | 3 |
| Reproductive expression + Adult restriction + HPA normal protein + Required coverage | 3 |
| Cancer expression + Required coverage | 2 |
| Reproductive expression + Adult restriction + Accepted normal-protein evidence | 1 |
| Reproductive expression + Cancer expression + Adult restriction + Required coverage + Normal bulk | 1 |

Missing developmental data or a multi-locus tumor profile can prevent the stricter
low-expression Neither determination without blocking an otherwise qualifying
adult-Reproductive CTA. The status tracked those gaps together. A failure of the
Reproductive expression floor means no qualifying RNA route in these inputs, not
proof that the protein is absent or that all reproductive IHC is negative.

### Outside definitions: 649 in release 1.22.0

| Failed requirements | Proteins |
| --- | --- |
| Reproductive expression + Adult restriction | 617 |
| Reproductive expression + Adult restriction + Normal bulk | 19 |
| Reproductive expression + Adult restriction + HPA normal protein + Normal bulk | 3 |
| Adult restriction | 2 |
| Reproductive expression + Adult restriction + HPA normal protein | 1 |
| Adult restriction + Normal bulk | 1 |
| Reproductive expression + Adult restriction + HPA tissue clusters | 1 |
| Adult restriction + HPA normal protein + Normal bulk | 1 |
| Reproductive expression + Adult restriction + HPA normal protein + Normal bulk + HPA tissue clusters | 1 |
| HPA tissue clusters | 1 |
| Reproductive expression + Cancer expression + Adult restriction | 1 |
| Normal bulk | 1 |

Most fail only the association floor. For example FOXR2 has 1.2 HPA nTPM despite
measured cancer expression; this failure is not normal-tissue danger. MAGEA3, TBATA,
SPRR4 and LY6G6D have qualifying association but separate normal-expression holds.
Those holds must be inspected independently of missingness or diffuseness.

Current unassigned status counts after removing the maximum veto and requiring
GTEx: **{'incomplete_evidence': 694, 'outside_category_definitions': 596}**. Download the
[current incomplete combinations](incomplete_evidence-current.csv) and
[current outside-definition combinations](outside_category_definitions-current.csv).
Every record in the detailed download retains its release-1.22 status and current
reason, so a coverage change cannot masquerade as newly measured toxicity.

## Method basis and reproducibility

[A benchmark of gene expression tissue-specificity metrics](https://doi.org/10.1093/bib/bbw008)
compares continuous tissue-specificity measures and discusses sensitivity to tissue
sampling. Our Reproductive-anchored bounded average is a **custom metric**, not
the paper's Tau statistic and not a validated 10% clinical or biological boundary.

Rebuild with `python scripts/audit_cta_profile_concentration.py` after the current
selection replay. `CANCERDATA_DATA_DIR` selects the shared source-cache root.
`oncoref.expression_diffuseness.assess()` exposes the score, weights, missingness
bounds and reference-quality status for downstream use. CSV source values remain
numeric and retain explicit units; labels in the report link to explanatory pages.
