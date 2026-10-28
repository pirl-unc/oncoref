# Policy simplification and developmental sources

This page separates the **active 1.21.0 policy** from proposed changes. The active output remains 463 CTA / 17,388 Unsafe / 2,012 Neither protein-equivalence groups. The [walkthrough](index.md) and [dog report](dog.md) describe what actually ran.

## Neither without a second RNA boundary

The released RNA rule has two boundaries: above the candidate ceiling gives a hold; at least twice that ceiling supplies Unsafe evidence. Older plots called the latter “Strong”. It is a policy multiplier, not measured toxicity, antibody confidence or a distinct biological state.

The requested simpler meaning does not need that intermediate RNA band. A possible replacement is:

| Category | Proposed definition |
|---|---|
| Unsafe | Credible normal expression exceeding the applicable single cutoff in vital/protected compartments, or independently accepted clinical injury counterevidence |
| CTA / oncofetal candidate | Tumor expression, reproductive/developmental association, adequate coverage and acceptable normal expression under the policy |
| Neither | All remaining genes: no qualifying tumor expression; no qualifying reproductive/developmental association; expression confined to other compartments; very low expression throughout; or unresolved coverage/evidence |

The category needs a separate reason/status field. **Neither because all measurements are low** differs from **Neither because measurements are missing**. A gene can lack cancer expression and still be Unsafe if a vital normal compartment supplies exclusion evidence. Every layer must evaluate every gene, including genes already outside candidate eligibility.

This proposal changes more than a plot label. The anatomy rule must apply consistently to protein and RNA. A positive result in a Manageable compartment should not automatically produce Unsafe under the proposed definition. Manageable expression can still prevent CTA admission under its ceiling and leave the gene Neither. Clinical counterevidence must retain its product, peptide and attribution scope. These changes have not been replayed or released here.

## Bowel and fibroblasts

Bowel epithelium can share the Manageable category. The numerical choice matters: the active Manageable base is 16, whereas the special bowel base is 32. Moving bowel to 16 alone gives **462 CTA / 17,391 Unsafe / 2,010 Neither**. PRAME moves to Neither; MAGE-A3, HOXA13 and FAM72B move from Neither to Unsafe. [Changed assignments](bowel-manageable16-comparison.csv).

Using 32 for all Manageable compartments is a different, broader change that requires replaying every context. The current comparison does not estimate that change. Calling a compartment Manageable expresses a tolerance assumption; it does not establish that injury is mild or reversible. Severe inflammatory colitis can occur during antigen-directed treatment ([Parkhurst et al., 2011](https://pubmed.ncbi.nlm.nih.gov/21157437/)).

Fibroblasts are resident cells of the organ even though they are not its cardiomyocytes or glandular epithelium. Partial depletion of PDGFRα-positive cardiac fibroblasts in mice was compatible with preserved cardiac function but altered tissue biology ([Kuwabara et al., 2022](https://elifesciences.org/articles/69854)). Targeting activated fibroblasts reduced fibrosis in a mouse injury model ([Aghajanian et al., 2019](https://doi.org/10.1038/s41586-019-1546-z)). Conversely, depletion of FAP-positive stromal cells in mouse skeletal muscle and bone marrow caused cachexia and anemia ([Roberts et al., 2013](https://doi.org/10.1084/jem.20122344)). These findings support distinguishing subsets and niches. They do not establish that all fibroblasts are disposable or that human cardiac fibroblast targeting is safe.

## PRAME has colon IHC

The frozen HPA v23 table reports PRAME **Not detected / Enhanced** in colon glandular cells, endothelial cells and peripheral nerve/ganglion cells. It also reports Not detected / Enhanced in adrenal glandular cells. [HPA PRAME colon](https://v23.proteinatlas.org/ENSG00000185686-PRAME/tissue/colon), [adrenal gland](https://v23.proteinatlas.org/ENSG00000185686-PRAME/tissue/adrenal+gland).

The RNA finding is **20.3 nTPM, 96 reads in colon undifferentiated cells**. The current matcher cannot establish that the broader IHC glandular-cell annotation covers that cluster. It therefore applies the unmatched-IHC ceiling of 32 while retaining the same-tissue negative evidence. At a Manageable base of 16, Approved negative credit would produce 32, and Enhanced credit would produce 128. Either would exceed 20.3, but neither credit should be assigned before verifying compartment coverage and the existing reagent/profile requirements. An independently curated Low adrenal protein finding also remains flagged and receives no negative-IHC credit.

CTAG2 has a similar matching limitation: HPA reports negative Enhanced respiratory epithelial staining, while the excluding RNA context is **bronchial basal respiratory cells, 24.9 nTPM, 58 reads**, against ceiling 4 and Unsafe boundary 8. The report retains the distinction between broad same-tissue negative staining and demonstrated coverage of the measured basal population.

## Developmental restriction and published nominations

The positive-context label should be **reproductive or developmental**, with distinct routes retained in the data. TBXT enters through embryonic expression, not testis expression: 245.088 TPM supported embryonic lineage mean versus 0.2 nTPM reproductive RNA. This makes it an oncofetal/developmental-route candidate within the combined set, rather than evidence that it is a classic cancer-testis gene.

The current gate measures **presence**, not developmental restriction. It admits 9,481 developmental-only groups; only 29 remain CTA after all filters. Another 1,427 groups pass both routes, of which nine remain CTA. The 12,139-group starting association set comes from whole-reference measurements, not a literature list.

Published oncofetal nominations can be curated symmetrically with published CTA lists: preserve every source row and its gene/exon/proteoform scope, map identities, report overlaps and apply independent normal/tumor filters. Relevant sources include:

| Source | Quantified or curated contribution | Scope and integration status |
|---|---|---|
| [Shaw et al. (2024)](https://doi.org/10.1038/s41467-024-47649-y), Supplementary Data 1 | 157 cancer-specific targets: 148 gene-level and nine alternatively spliced targets. Six have oncofetal annotations: TGFB2, WNT5A, GPC3, IGF2, and specific FN1/TNC splice targets | Explicit oncofetal annotations, derived through text mining and manual curation. Pediatric tumor comparison against adult GTEx. Newly identified source; not included in current counts. FN1/TNC annotations cannot establish whole-gene restriction |
| [Oncofetal TRIM71 study (2024)](https://pmc.ncbi.nlm.nih.gov/articles/PMC11388079/), Table S1 | Top-200 differential-expression oncofetal signature from fetal/adult/tumor liver comparisons | Liver-specific program; newly identified and not integrated. Differential expression does not establish absence in other normal organs |
| [Cardoso-Moreira et al. (2019)](https://pmc.ncbi.nlm.nih.gov/articles/PMC6658352/), supplementary developmental trajectories | Seven-organ developmental series, including prenatal, child and adult human samples; temporal specificity measurements | Provides a quantitative way to assess developmental enrichment and postnatal persistence. Organ/age coverage is incomplete; native CPM/RPKM cannot inherit TPM limits. Not integrated |
| [Tyser et al. (2021)](https://doi.org/10.1038/s41586-021-04158-y) | Supported cell-lineage means from one Carnegie-stage-7 human embryo | Used now for the expression route; not a curated list of developmentally restricted genes |

Later fetal atlases, including Cao/DESCARTES, can broaden stage and lineage coverage. Fetal expression alone cannot establish postnatal absence. A restricted-developmental nomination should record prenatal enrichment, observed postnatal measurements and missing organ/age coverage separately. No source above supplies a complete whole-body child/adult absence certificate. Adult HPA/GTEx filtering does not establish pediatric safety.

[Human source counts and funnel](index.md) · [Dog coverage](dog.md) · [Shared API](api.md).
