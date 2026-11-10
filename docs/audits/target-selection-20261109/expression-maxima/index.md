# Expression maxima and Disposable caps

**Requiring cancer or reproductive/developmental expression to exceed every measured nonreproductive value removes zero of the current 910 Target candidates.** This is an additional filter on the current catalog. The original cancer gate, reproductive association and shared eligibility rules still apply.

For each protein, compare the largest tumor-sample TPM, adult reproductive-tissue TPM or supported developmental-lineage mean TPM against the largest nonreproductive normal or Disposable value. The target value must be at least 1 TPM and strictly greater; ties fail. One tumor sample can establish the maximum. Protein equivalence records retain the rebuild's within-context member aggregation; maxima from different samples are never added.

## Absolute Disposable cap

The cap is inclusive: a value equal to it passes. Exceeding it moves a candidate to Other, without creating new Unsafe evidence. These are sensitivity analyses; no cap has been adopted.

| Disposable cap (native TPM) | Move to Other | CTA remaining | Cancer Antigen only remaining | Reproductive Antigen only remaining |
|---:|---:|---:|---:|---:|
| None | 0 | 632 | 250 | 28 |
| 1 | 10 | 627 | 246 | 27 |
| 4 | 2 | 631 | 249 | 28 |
| 8 | 1 | 632 | 249 | 28 |
| 16 | 1 | 632 | 249 | 28 |
| 32 | 0 | 632 | 250 | 28 |
| 64 | 0 | 632 | 250 | 28 |
| 128 | 0 | 632 | 250 | 28 |

LCE3C has the highest Disposable value, **21.4 TPM in tonsil**. NLRP4 is next at **4.6 TPM in tonsil**. The [complete protein table](protein-maxima.csv) retains each maximum, source, tissue and tumor sample. [Every changed protein at every cap](proteins-moved-to-other.csv) is listed separately.

## What this comparison measures

The TPM panel's Disposable tissues are **appendix, gallbladder and tonsil**. It contains no isolated hair-follicle or cutaneous-melanocyte TPM measurement. Thus zero losses means zero in this measured panel, not a complete all-tissue result. Normal values include tissue means and donor medians; this is not a maximum over every healthy individual. The underlying [tumor-source comparability concerns](../scrutiny/index.md) remain.

[HPA v23](https://v23.proteinatlas.org/about/download) has a **skin-melanocyte cluster** on a separate nTPM scale. 908 candidates have a mapped member observation there; 2 member loci are missing. No candidate has a positive signal meeting the existing 20-read requirement. The largest observed signal is **TSGA10IP: 24 nTPM, 14 reads**; it is retained as insufficiently supported, not converted to zero. HPA's source labels contain no hair-follicle cluster. [Cell observations and read support](hpa-disposable-cell-observations.csv), [missing members](hpa-disposable-cell-missing-members.csv) and [all source cell labels](hpa-cell-contexts.csv) are retained. TPM and nTPM are not pooled or treated as interchangeable. The cell result adds no supported cap exceedance on its own scale and does not fill the missing hair-follicle panel.

The proposed rule has two independent parts: **target expression must be highest**, and **Disposable expression must remain below a fixed ceiling**. A large tumor value cannot excuse unlimited Disposable expression. Required but missing normal measurements belong in Other:Review; an observed cap failure belongs in Other:Expression. The current required-coverage panel has not been expanded by this audit.

## Reproductive contexts and location labels

**Reproductive includes adult reproductive organs, placenta/trophoblast, extraembryonic membranes and structures, embryo and fetus.** The taxonomy already contains this union. The inventory now records selection scope separately from injury tier and uses the correct embryonic/fetal stage. Fourteen of the previously unresolved hazard entries belong to this Reproductive scope. Reproductive membership does not exempt the same immune, vascular or other shared cells at nonreproductive sites.

**Soft tissue, head and foot are location labels requiring localization.** They receive no new hazard tier and no Disposable allowance. Retain their observations and resolve the actual tissue/cell identity before using them for clearance. Among the currently registered specific human concepts, none remains without either a Reproductive selection scope or an adult hazard tier. This is a registry inventory, not a complete list of human tissues or source coverage.

Three concrete additions need separate curation:

| Compartment | Why the existing label is insufficient | Source to support mapping and expression review |
|---|---|---|
| Olfactory sensory neurons, sustentacular cells and basal progenitors | Nasopharynx and generic respiratory epithelium do not identify the olfactory neuroepithelium. | [Durante et al., 2020](https://www.nature.com/articles/s41593-020-0587-9), [GSE139522](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE139522) |
| Mineralized bone; osteoblasts, osteocytes and osteoclasts | Bone marrow and cartilage do not cover these cells. | [Human osteoblast PANX3 study](https://pubmed.ncbi.nlm.nih.gov/32982413/) provides a candidate-specific lead; its cell-line evidence does not establish an adult tissue hazard tier. |
| Dental pulp, odontoblasts and periodontal cells | Oral mucosa does not identify the tooth's living compartments. Adult cells belong in normal-tissue assessment; fetal specimens also retain developmental scope. | [Pagella et al., 2021](https://pubmed.ncbi.nlm.nih.gov/33997688/) |

These studies support distinct normal compartments; they do not establish tolerance of immune-mediated depletion. New hazard calls require the registry's functional-injury rationale and evidence.

[Current anatomy inventory](anatomy-inventory.csv) · [Scenario counts](cap-scenarios.csv) · [Summary](summary.json) · [Provenance manifest](manifest.json).

Reproduce: `python scripts/audit_target_expression_maxima.py --hpa-cell-file PATH_TO_V23_RNA_SINGLE_CELL_TYPE_TISSUE_TSV`.
