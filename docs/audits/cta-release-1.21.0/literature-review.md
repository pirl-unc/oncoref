# Literature and implementation review

Reviewed 27 October 2026. The selection is an explicitly permissive **research prioritization policy**, calibrated using PRAME. Passing PRAME is not an independent validation result. The tier values, doubling factors and corneal scale factor are hypotheses; none was estimated against a prospective toxicity ground truth. Published expression, peptide presentation and product-specific clinical outcomes retain separate meanings.

## Decisions checked against primary evidence

| Topic | Evidence and interpretation | Release treatment |
|---|---|---|
| PRAME adrenal protein | Lezcano et al. report low normal adrenal-cortical expression, alongside reproductive/endometrial findings. The abstract supports a categorical finding, not a numerical intensity or stained-cell fraction. [Primary publication](https://pubmed.ncbi.nlm.nih.gov/33890816/) | Preserve `Low`, `author_description`, the source locator and abstract-level access. No invented H-score or HPA reliability grade. Low can pass if independent RNA gates pass. |
| PRAME clinical targeting | The IMA203 phase 1 publication reports 40 infused patients. Its safety and normal-peptide evaluation concern one receptor, PRAME004 and HLA-A*02:01, not every PRAME-targeting construct. [Wermke et al.](https://www.nature.com/articles/s41591-025-03650-6) | Clinical experience is contextual evidence, never an admission bypass or a gene-wide safety guarantee. |
| Cardiac fibroblasts | Controlled fibroblast reduction in adult mice was tolerated, with changes in extracellular matrix and cardiomyocyte physiology. Residual fibroblasts responded to injury. This does not establish the safety of eliminating all fibroblasts or targeting human heart stroma. [Kuwabara et al.](https://elifesciences.org/articles/69854) | A separate Manageable stromal tier is a research assumption. Heart bulk, cardiomyocytes and endothelial cells keep anatomical tiers. |
| Bowel epithelium | All three recipients in a CEA-directed TCR study developed severe transient colitis that was dose limiting. Direct epithelial targeting cannot be justified simply by analogy to accepted checkpoint-inhibitor toxicity. [Parkhurst et al.](https://pubmed.ncbi.nlm.nih.gov/21157437/) | The more permissive large-bowel tier remains explicit; measured undifferentiated-cell expression carries a stem/progenitor caution. |
| HPA IHC | HPA v23 distinguishes staining categories and reliability based on supporting evidence. Reliability is not a calibrated false-negative probability. [HPA assays and annotation](https://v23.proteinatlas.org/about/assays%2Bannotation) | Negative credit requires exact context, controls and resolved conflicts. Positive Low never grants negative credit. Author descriptions are not silently converted to HPA grades. |
| TBXT antibody disagreement | Clone 54-1 has specificity tests and a different normal-tissue pattern from the uncertain HPA reagent. Evidence is tissue-specific; a better antibody does not validate every untested tissue. [Hamilton et al.](https://www.oncotarget.com/article/3086/text/) | Retain categorical external observations and review flags; no blanket “TBXT is absent from all normal tissues” statement. TBXT passes the quantitative funnel. |
| Adult intervertebral disc | Adult nucleus-pulposus studies include brachyury RNA/protein observations, outside the routine bulk/IHC tissue panels. [Richardson et al.](https://www.nature.com/articles/s41598-017-01567-w) | Preserve the adult-disc evidence/coverage caveat. Selection is not a claim that the disc has no TBXT. |
| Developmental admission | The integrated Tyser atlas covers one CS7 human embryo, approximately days 16–19, with 1,195 cells. Cell support is not replication across donors. [Tyser et al.](https://www.nature.com/articles/s41586-021-04158-y) | Use the pinned E-MTAB-9388 gene-TPM matrix, supported lineages and explicit stage. Later fetal/DESCARTES coverage is not claimed. |

PRAME's retained normal measurements are adrenal bulk **3.6 nTPM** (ceiling 4), cardiac fibroblasts **12 nTPM** (16), and colon undifferentiated cells **20.3 nTPM** (32). These explain admission without hiding positive evidence. CTAG2/LAGE-1 remains **Unsafe**, including **24.9 nTPM in bronchial basal respiratory cells**, above the Vital stronger boundary of 8 with missing/unusable matched IHC. Its **5.1 nTPM heart bulk** also exceeds the candidate ceiling of 4, supplying a hold rather than that stronger boundary on its own. These measurements do not establish clinical injury or peptide presentation.

## Species source categories

Source availability and integration are separate fields in the [inventory](species-source-inventory.csv). Human selected groups, canine study nominations and an unavailable mouse snapshot are not interchangeable denominators.

| Resource | Biological categories | Unit / current integration |
|---|---|---|
| [BarkBase](https://pubmed.ncbi.nlm.nih.gov/31181663/) | Dog adult bulk, reproductive and embryonic expression | Native matrix, TMM-CPM; not used with human TPM cutoffs |
| [DoGA](https://www.nature.com/articles/s41467-024-52798-1) | Dog/wolf adult tissues and dog embryonic expression | Metadata registered; 5′ assay needs its own scale and dog/wolf/stage selection |
| [Pfizer canine atlas](https://doi.org/10.1038/s41597-026-07512-9) | Normal protein and normal bulk RNA | iBAQ protein matrix available; RNA metadata registered; pooled samples are not independent donors |
| [Tabula Muris/Senis](https://www.nature.com/articles/s41586-020-2496-1) | Mouse adult/aging normal cell and immune RNA | Identified, not integrated; assay-native counts/platforms |
| [MOCA](https://www.nature.com/articles/s41586-019-0969-x) | Mouse organogenesis | Identified, not integrated; embryonic UMIs |
| [Mouse gastrulation atlas](https://www.nature.com/articles/s41586-019-0933-9) | Mouse early developmental lineages | Identified, not integrated; embryonic UMIs |
| [GXD](https://www.informatics.jax.org/mgihome/GXD/aboutGXD.shtml) | Mouse developmental RNA, normal RNA and protein assays | Identified assay repository; retain assay and anatomical evidence, distinguish not detected from not analyzed |

Published CTA and placental lists annotate the human results; membership is not an escape hatch. Bradley et al.'s [cancer-placenta antigen study](https://doi.org/10.1038/s41467-020-19141-w) and Gong et al.'s [placental expression study](https://doi.org/10.1038/s41467-021-22695-y) have different evidentiary roles. A developmental or placental expression atlas is not automatically a validated oncofetal target list. Explicitly oncofetal resources identified in the earlier [developmental assessment](../../developmental-targets.md), including the Su liver program, remain prospective intake unless marked integrated in the source inventory.

## Corrections made during self-review

1. Grouped IHC previously flattened every positive member to Low. The new collapse preserves the strongest level, its own confidence, all observed levels and member completeness. It corrected 2,902 group/context labels (1,778 Medium; 1,124 High).
2. Group negative-IHC credit now requires positive profile controls for every expected member. A reliable weak stain cannot lend its confidence to an uncertain stronger stain.
3. Immune and corneal evidence is replayed from complete contexts. Rescaling only a maximum RNA ratio could lose an independent mapping or coverage hold. Unpaired group detection bounds remain explicit.
4. Zero expression no longer receives a stem/progenitor expression flag. Unknown anatomy and contradictory positive-IHC/negative-credit input fail validation.
5. Imported protein findings require an explicit reporting scale and provenance. Accepted positives with no extracted intensity remain exclusions, rather than being discarded for limited access.
6. Current results are packaged behind one species-scoped API. Historical decisions in locus tables are explicitly named as historical/baseline; they cannot masquerade as current classifications. Dog/mouse coverage gaps remain visible.

The corrected replay changes intermediate decisions but **no final group categories** relative to the prior 463-group proposal. The earlier proposal's exact implementation bytes are archived alongside its original hashes; the release has separate provenance. Full current decisions and retained known-CTA outcomes are in [the main report](index.md).
