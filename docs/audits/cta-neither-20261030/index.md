# Neither, missing categories, and targeting candidates

**Implemented update:** see [1.25 rules, category changes and review queues](update-1.25.0.md). The original audit below records the 1.24 baseline. Manageable excess remains a normal-evidence hold; it was not added to Unsafe.

This is an audit and proposal against **OncoRef 1.24.0**, dated 30 October 2026. No classification policy, released assignment, sequence mask or threshold was changed. The [numerical report](numbers.md), [reproducible script](https://github.com/pirl-unc/oncoref/blob/v1.25.0/scripts/audit_cta_neither_boundaries.py) and [input checksums](manifest.json) record the analysis.

**Neither is underspecified.** It currently recognizes three narrow profiles rather than all adequately assessed proteins that pass the normal screen but lack CTA identity. In particular, its low-expression route also requires low cancer expression. That makes a protein leave Neither when a cancer sample expresses it, even though nothing changed about normal-tissue evidence.

## Recommended categories

Use four disjoint operational categories, with biological identity and evidence strength recorded separately:

| Category | Meaning | Sequence-screening role |
|---|---|---|
| **Candidate** | Adequate normal coverage, no unresolved normal concern, passes all applicable normal ceilings/protein gates, and has a supported cancer-target rationale. | Supplies candidate sequences; still check each subsequence against exclusion evidence. |
| **Neither** | Adequate normal coverage and passes the normal screen, but does not meet a candidate definition. | Does not automatically blacklist every shared subsequence. Retain its measured normal-expression evidence. |
| **Unsafe** | Accepted normal-tissue evidence violates the chosen targeting policy. | Supplies the hard sequence-exclusion set, with the responsible tissue, measurement, tier and evidence attached. |
| **Unresolved** | Required normal coverage, sequence attribution or adjudication is insufficient to decide. | Supplies a separate uncertainty mask; a strict downstream screen can quarantine overlaps. Missing evidence is not a demonstrated injury signal. |

Candidate would have at least two biological subtypes: **CTA** (cancer plus reproductive/developmental association and the existing concentration requirement), and **non-reproductive tumor candidate** (a credible cancer-associated expression profile without the required reproductive association). Within the second subtype, reserve **tumor-restricted** for evidence that actually supports restriction. A single nonzero tumor RNA measurement does not establish tumor-only protein production. Cancer-specific altered sequences are a further distinct type and should not be confused with expression of the ordinary reference protein.

This does not propose a new tumor-versus-normal ratio for CTA admission. Keep cancer presence, reproductive identity, normal safety and diffuseness separate. The additional non-reproductive route needs its own transparent evidentiary definition before automatic admission; the present audit identifies a review queue, not a new threshold. RNA-supported and protein-confirmed candidates should remain distinguishable without imposing a new blanket protein-evidence requirement on established CTAs.

Neither sublabels can be **low**, **diffuse**, **Disposable**, or **other**. Low measured expression, broad low expression, and expression confined to permitted normal compartments should not require a special combination of tissue counts or low tumor RNA. A reproductive protein without qualifying cancer evidence can also be Neither: Neither does not mean that reproductive expression is necessarily absent.

There is one explicit policy decision in this proposal: if a Manageable ceiling is a targeting exclusion, a supported breach should join the operational Unsafe mask, with **Manageable** severity retained. The current implementation instead lets Manageable excess block CTA while sometimes assigning Neither. Calling all these proteins clinically dangerous would be unwarranted; calling them cleared for unrestricted targeting would also be misleading. An alternative is to name the operational bucket **Excluded** and retain Unsafe as an evidence label. Either convention must preserve the same exclusion behavior.

The decision order should be: accepted normal exclusion first; unresolved normal evidence second; supported candidate identity third; otherwise Neither. Missing cancer or reproductive evidence alone must not create normal-safety uncertainty. Optional missing sources remain visible and are not converted to zero.

## What the counts show

| Current result | Proteins | Audit finding |
|---|---:|---|
| CTA | 522 | Unchanged by this audit. |
| Unsafe | 18,031 | Predominantly supported normal expression; the size is not explained by the Neither gap. |
| Neither | 20 | Nine low-expression, ten “diffuse Manageable,” one diffuse reproductive. |
| Unassigned | 1,290 | 763 pass the recorded normal screen and required coverage; they lack a general fallback category. |

Of those 763, **757 fail only reproductive association** and six also lack qualifying cancer RNA. Mechanically adding a normal-pass/non-CTA fallback would produce **783 Neither and 527 Unassigned**, leaving CTA and Unsafe unchanged. This is a structural counterfactual, not 763 individual safety endorsements. Some belong in a non-reproductive candidate review queue, some may qualify as CTA through protein evidence, and some need rare-cell or sequence-attribution review.

There are **52** proteins in this queue with positive reproductive IHC for every member locus, including **38** with Medium/High staining and Supported/Enhanced reliability; 37 of those 38 also satisfy cancer presence. The current reproductive admission function is RNA-only, requiring >=5 nTPM or qualifying developmental TPM. This loses evidence that the earlier intended positive-protein route should consider. Examples include SAGE1, GAGE10, CDY1/B and CT45A2/8/9. An antibody/paralog review is still required; these are not 52 automatic rescues. See [the complete IHC review queue](reproductive-ihc-review.csv).

## USP17L3 and FAM90A9: interesting, but not demonstrated tumor-only proteins

| Pinned measurement | USP17L3 | FAM90A9 |
|---|---:|---:|
| Maximum supported non-permitted normal nTPM | 0.1 | 0.1 |
| Maximum supported non-permitted normal TPM | 0 | 0 |
| Reproductive bulk maximum, nTPM | 0 | 0 |
| Maximum measured embryonic lineage RNA, TPM | 0.000605 | 0.221 |
| OncoRef cancer maximum, TPM | 3.82 | 1.75 |
| Cancer associated with that OncoRef maximum | Head/neck squamous carcinoma | Stomach adenocarcinoma |
| Separate HPA tumor maximum/lower-bound summary, pTPM | 5.98 | 15.10 |
| Measured HPA normal tissue IHC | None | None |
| Optional DICE / corneal coverage | Complete / complete | Missing / missing |

These are measured source summaries, not the maximum over every normal donor or all possible cell states. TPM and HPA pTPM were not pooled or interpreted as replicate measurements. No measured IHC means missing protein evidence, not a panel of negative antibody stains. Maxima attest to individual samples, not prevalence. Tumor bulk RNA can also arise from nonmalignant cells.

**USP17L3.** The established USP17/DUB literature concerns a duplicated deubiquitinase family. [McFarlane et al., *The deubiquitinating enzyme USP17 is highly expressed in tumor biopsies, is cell cycle regulated, and is required for G1-S progression* (2010)](https://pubmed.ncbi.nlm.nih.gov/20388806/) reports tumor expression and cell-cycle biology. That is useful family context, but it does not establish that the specifically quantified USP17L3 product is a tumor-restricted antigen. [Burrows et al., *The DUB/USP17 deubiquitinating enzymes, a multigene family within a tandemly repeated sequence* (2005)](https://pubmed.ncbi.nlm.nih.gov/15780755/) establishes the duplicated architecture and explains why gene attribution matters.

A more specific counterexample is [Cross et al., *Characterisation of HBV and co-infection with HDV and HIV through spatial transcriptomics* (2024)](https://pmc.ncbi.nlm.nih.gov/articles/PMC11326438/): the authors report USP17L3 enrichment in the hepatitis B/delta coinfected liver specimen. This is diseased non-tumor tissue, not a healthy-donor baseline; it does not prove a protein-level hazard or give a comparable TPM ceiling. It does require an activation/inflammation review before asserting tumor restriction. My disposition is **normal-atlas low; non-reproductive candidate hypothesis; unresolved tissue-state and paralog attribution**, rather than a validated target.

**FAM90A9.** [Zheng et al., *Evaluation of PE5 Entries on the Human Proteome Project Target List* (2025)](https://pmc.ncbi.nlm.nih.gov/articles/PMC12459392/) directly discusses this locus. It reports no PeptideAtlas detection or HPA tissue detection, notes HPA subcellular staining, and recommends coding status based on homology rather than direct unique protein evidence. The authors emphasize the highly similar FAM90A copies and limited functional knowledge. [UniProt A6NNJ1](https://www.uniprot.org/uniprotkb/A6NNJ1) describes protein existence as inferred from homology. I found no target-specific clinical efficacy or safety evidence in the focused search. My disposition is **tumor-RNA candidate for investigation, with protein existence and paralog attribution unresolved**. It should not be called tumor-only or automatically admitted to a target sequence list.

## Individual review of all 20 Neither proteins

Numbers below are current supported **non-permitted normal maxima**, with nTPM and TPM kept separate. Complete source measurements and tumor values are in [the individual data table](neither-individual-data.csv). These are proposed review dispositions, not changes to released labels.

| Protein | Normal max nTPM / TPM | Interpretation and recommended handling |
|---|---:|---|
| **CRYGN** | 15 / 18.04 | Thyroid TPM exceeds its 16 ceiling. This is a Manageable hold, not a cleared low-expression profile. Review thyroid and ocular/neural coverage; animal crystallin biology must not be treated as proof of human toxicity. [Human gene annotation](https://www.ncbi.nlm.nih.gov/gene/155051); [gammaN-crystallin evolution study](https://pubmed.ncbi.nlm.nih.gov/15853812/). |
| **KRT34** | 64.5 / 0.932 | Skin RNA exceeds the ceiling. Hair-restricted protein evidence could support a permitted-compartment interpretation, but whole skin cannot simply be relabelled hair. [HPA hair staining](https://www.proteinatlas.org/ENSG00000131737-KRT34/tissue/hair). Hold pending compartment attribution. |
| **KRT85** | 167.2 / 0.219 | Same skin-versus-hair compartment issue, with much higher skin RNA than the missing-IHC ceiling. Review hair-specific localization before any rescue. [Human hair-keratin study](https://pubmed.ncbi.nlm.nih.gov/31580512/). |
| **KRT71** | 96.9 / 0.922 | A hair-associated keratin, but current whole-skin evidence violates the normal ceiling. Not equivalent to diffuse low expression. [Human hair-keratin study](https://pubmed.ncbi.nlm.nih.gov/31580512/). |
| **OXGR1** | 5.6 / 3.5 | Current normal screen passes. Neither is reasonable provisionally, but bulk maxima can miss rare renal cells: a [mouse study of renal acid-base regulation](https://pmc.ncbi.nlm.nih.gov/articles/PMC3696567/) locates functional OXGR1 in intercalated cells. This motivates human rare-cell review, not a human toxicity verdict. |
| **KRT9** | 11.8 / 2.519 | Current screen passes, but palm/sole expression is a major anatomical blind spot in ordinary skin panels. This is a specialized normal structural protein, not evidence of general silence. [Keratin 9 and palmoplantar epidermal integrity](https://pmc.ncbi.nlm.nih.gov/articles/PMC3923277/). Review volar skin before targeting. |
| **LEP** | 296.7 / 194.26 | Clearly substantial normal adipose expression. A poor unrestricted target candidate; it should not share a cleared Neither meaning with low-expression proteins. Human leptin deficiency establishes important normal physiology, not T-cell injury per se. [Montague et al., congenital leptin deficiency](https://pubmed.ncbi.nlm.nih.gov/9202122/). |
| **LCE2A** | 114.9 / 119.18 | Strong normal epidermal expression; not a Disposable hair-follicle protein. Exclude under the chosen Manageable ceiling unless that policy explicitly changes. [Human skin transcriptomics and antibody profiling](https://pmc.ncbi.nlm.nih.gov/articles/PMC4305515/). |
| **KRT25** | 108.8 / 0.263 | Another hair-associated keratin with a whole-skin ceiling failure. Needs exact compartment review, not the >=2-Manageable-tissue shortcut. [Human hair-keratin study](https://pubmed.ncbi.nlm.nih.gov/31580512/). |
| **MTCP1** | 4.1 / 7.319 | The clearest current diffuse-low Neither example: reproductive RNA 5.1, HPA diffuseness 24.53%, no recorded normal hold. Retain Neither at the current aggregate resolution. Isoform-specific cancer biology should be evaluated separately; [p13MTCP1 expression in translocation-bearing T-cell proliferations](https://www.sciencedirect.com/science/article/pii/S0006497120653022) cannot be inferred from all gene RNA. |
| **FAM90A5** | 0.1 / 0 | Low observed normal and cancer RNA; Neither is reasonable with a duplicated-family attribution flag. Lack of signal does not establish biological silence. [FAM90A family characterization](https://doi.org/10.1093/hmg/ddm209). |
| **USP17L5** | 0.1 / 0 | Low observed RNA; Neither is reasonable. Retain paralog/activation caveats rather than importing every USP17 family function into this locus. [USP17 family study](https://pubmed.ncbi.nlm.nih.gov/15780755/). |
| **USP17L12** | 0.5 / 0 | Same conclusion; its tumor maximum is 1.65 TPM, just below the low-expression definition's tumor cutoff. The cutoff should not determine normal safety. [USP17 family study](https://pubmed.ncbi.nlm.nih.gov/15780755/). |
| **FAM90A3** | 0.1 / 0 | Low observed normal and tumor RNA. Neither with family-attribution uncertainty; no locus-specific therapeutic inference. [FAM90A family characterization](https://doi.org/10.1093/hmg/ddm209). |
| **USP17L18** | 0.2 / 0.2 | Low observed normal RNA; tumor maximum 1.58 TPM. Keep Neither provisionally, with the same duplicated-family/state-coverage limitations. [USP17 family study](https://pubmed.ncbi.nlm.nih.gov/15780755/). |
| **TAF11L4** | 0 / 0 | Low in this panel, not demonstrated absent everywhere. External [Bgee expression calls](https://www.bgee.org/gene/ENSG00000284283) motivate source/identity review but are not quantitative TPM exclusions. Neither with explicit measurement scope. |
| **TMDD1** | 13.5 / 1.272 | Current screen passes; 13.5 nTPM in adipose fibroblasts is not “never expressed.” Neither is reasonable without reproductive association, but record that normal compartment and avoid a silence label. Functional evidence is limited; [HPA TMDD1](https://www.proteinatlas.org/ENSG00000284730-TMDD1). |
| **FAM90A22** | 0 / 0 | Low observed expression. Proteomics review reports only a non-unique peptide assignment: do not treat it as confirmed absence or unique protein detection. Neither with identity review. [Zheng et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC12459392/). |
| **FAM90A14** | 0.2 / 0 | Low observed RNA; retain Neither provisionally. The proteomics review does not establish unique tissue protein detection. [Zheng et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC12459392/). |
| **FAM90A7** | 0 / 0 | Same low-observed profile and highly similar family-member concern. Neither with explicit protein-evidence limitations. [Zheng et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC12459392/). |

Thus **seven** Neither proteins have normal-screen holds: CRYGN, KRT34, KRT85, KRT71, LEP, LCE2A and KRT25. Six exceed 30 TPM or nTPM in a non-permitted normal context. The current later policy allows Manageable failures in Neither; it does not enforce the earlier expectation that Neither would always remain below 30. This audit exposes that semantic conflict rather than silently reinstating a historical cap.

## Similar proteins in other categories

| Comparison | Evidence | Implication |
|---|---|---|
| USP17L3 vs USP17L5 | Both normal maxima 0.1 nTPM / 0 TPM; tumor maximum differs. | Cancer expression should not create an unexplained unassigned category. USP17L3 additionally needs the inflammation review above. |
| FAM90A9 vs FAM90A3/A5 | Similar very low normal summaries; the HPA cancer measurement moves FAM90A9 outside the low-expression template. | A category-definition artifact, with shared paralog uncertainty. |
| KRT35 vs KRT25/KRT34 | KRT35 has skin 56.3 nTPM and only one counted Manageable tissue; the others satisfy the two-tissue template. | Adding another normal tissue should not make a ceiling failure look more acceptable. |
| LCE3B, PLA2G2E, LCE3C | Tonsil 3.0, 2.4 and 23.4 nTPM; non-permitted normal maxima 0, 0.7 and 1.7 nTPM respectively. | Clear candidates for a permitted-normal/Disposable profile rather than unexplained Unassigned. This does not establish tumor-only production. |
| FAM90A15 | Brain 4.2 nTPM against a 4 ceiling. | RNA-only boundary review; not evidence of clinical injury. |
| KRT74 | DICE NK-cell RNA 8.23 TPM against an 8 ceiling. | Similar narrow margin, with assay and attribution review needed. |
| DAZL | DICE 8.01 TPM against an 8 ceiling. | A known reproductive protein excluded at a very small numerical margin; retain the raw value and uncertainty rather than present the label as established toxicity. |
| KRT72 / KRT73 | DICE naive-CD4 RNA 580.51 / 137.74 TPM against an 8 ceiling. | These are materially different from a 4.2-versus-4 boundary. An unexpected lineage warrants checking provenance; it is not grounds to ignore the measurement. |
| KRT75 | Pancreatic duct-cell RNA 210.4 nTPM against a 4 ceiling. | Stronger non-permitted normal concern, subject to source/identity validation. |

Across Unsafe there are **138 RNA-only cases** whose *every recorded Unsafe RNA exclusion* lies no more than 25% above its local ceiling. This is a sensitivity-review queue, not an argument to relax every threshold by 25%. Some also have Manageable excess or candidate-identity failures. [All boundary cases](unsafe-rna-boundaries.csv) and [source-specific normal peaks](review-normal-peaks.csv.gz) retain their reasons.

Conversely, passing the core panel is not an exhaustive safety assessment. POU4F3 is among the normal-pass Unassigned proteins, yet [human inherited-hearing-loss evidence](https://pubmed.ncbi.nlm.nih.gov/9506947/) highlights an inner-ear compartment that deserves explicit review. The same principle applies to palmoplantar skin, renal intercalated cells, inflammatory states and duplicated genes.

## What downstream sequence users need

The category table should not itself be the only sequence-exclusion mechanism. Export candidate sequences, accepted normal exclusion evidence, and unresolved evidence separately, with exact protein identity and tissue provenance.

An exact subsequence shared by a Candidate protein and an Unsafe protein remains excluded regardless of the candidate's gene label. A match only to Neither should not automatically be treated as an Unsafe match. Matches to Unresolved should be visible and optionally withheld by a strict screen. These policies concern exact sequence overlap; they do not establish absence of cross-reactivity to different sequences.

Keep full-length identical proteins aggregated, retain contributing genes and isoforms, and evaluate every substring against all possible normal origins. A low-expression locus cannot rescue an identical sequence produced by a different, normally expressed locus. Conversely, a family resemblance alone should not automatically condemn every unique region without evidence.

## Downloadable review tables

- [All 20 Neither proteins and measurements](neither-individual-data.csv)
- [All 1,290 Unassigned proteins](unassigned-all.csv)
- [763 normal-screen-pass Unassigned proteins](unassigned-normal-pass.csv)
- [52 reproductive-IHC review candidates](reproductive-ihc-review.csv)
- [138 RNA-only Unsafe boundary cases](unsafe-rna-boundaries.csv)
- [Related-family proteins across categories](family-neighbors.csv)
- [Disposable bulk examples](disposable-bulk-examples.csv)
- [USP17L3 and FAM90A9 measurements](tumor-positive-low-normal-examples.csv)
- [Numerical report and counterfactual](numbers.md)

Validation: the audit preserves all 19,863 protein identities and the released category partition; checks ensure that the 763-protein fallback queue contains no recorded Unsafe evidence. This is a local review artifact, not a merged policy change or a new release.
