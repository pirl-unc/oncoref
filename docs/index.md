# oncoref

`oncoref` is the base-layer package for shared cancer reference data:
cancer-type ontology, cohorts, expression, clean-TPM normalization, TMB,
incidence/mortality, ICI response, HPA normal-tissue expression, and
HPA-derived cancer-testis antigen references.

Downstream packages such as [pirlygenes](https://github.com/pirl-unc/pirlygenes)
and [trufflepig](https://github.com/pirl-unc/trufflepig) should delegate
parity-clean shared primitives here, but they keep different kinds of work:
pirlygenes owns curated gene sets and panels, while trufflepig owns per-sample
interpretation and rule firing.

Use this page as the quick orientation. The [API guide](api.md) is the detailed
module map, with examples and notes about compatibility modules.

## Stack Boundary

- **oncoref**: empirical base facts and canonical identifiers. Use it for cancer
  codes, gene IDs, reference expression/normalization, epidemiology, TMB, ICI/aPD1
  response, HPA normal tissue, and source-anchored CTA facts. If a row has an
  `n`, confidence interval, source cohort, or PMID/DOI anchoring a measurement, it
  usually belongs here.
- **pirlygenes**: curated gene selections. Use it for lineage/family/compartment
  panels, discriminators, surfaceome, TME and stem-cell markers, response
  signature panels, target-to-therapy registries, and other purpose-specific
  gene sets keyed to oncoref IDs.
- **trufflepig**: per-sample application. Use it for sample QC narration,
  library-prep/source warnings, deconvolution, scoring, and tumor-sample rules.

## Start Here

- Need cancer-type codes, hierarchy, molecular subtypes, or matched normal
  tissues? Start with [Cancer Vocabulary](api.md#cancer-vocabulary).
- Need expression values or clean-TPM normalization? Start with
  [Expression And Normalization](api.md#expression-and-normalization).
- Need anti-PD-1 or broader checkpoint-inhibitor response estimates? Start with
  [ICI Response](api.md#ici-response).
- Need cancer-testis antigen definitions, patient coverage, or CTA-specific
  peptide load? Start with [CTA Antigens](api.md#cta-antigens).

```python
from oncoref import cancer_ontology, ici_response

crc_msi = cancer_ontology.cancer_type_records(subtype_group="MSI", under="CRC")
crc_msi[["code", "evidence_source_code", "normal_tissue_code"]]

# MMR/MSI classifier queries keep positive, negative, and confounder classes
# explicit, and can be restricted to direct expression-backed codes.
cancer_ontology.mmrd_cancer_codes(expression_only=True)
cancer_ontology.mmr_confounder_cancer_codes()

ici_response.best_available_ici_response("COAD_MSI")
```

## Data Domains

| Domain | Public modules | Use for |
| --- | --- | --- |
| Cancer vocabulary | [`oncoref.cancer_ontology`](api.md#cancer-vocabulary), [`oncoref.cohorts`](api.md#cancer-vocabulary) | Registry records, aliases, hierarchy, subtype axes, MMR/MSI classifier status, cohort IDs, matched normal tissues |
| ICI response | [`oncoref.ici_response`](api.md#ici-response) | Anti-PD-1 and broader ICI references, regimen-aware lookups, extracted endpoint estimates |
| CTA references | [`oncoref.cta`](api.md#cta-antigens), [`oncoref.cta_coverage`](api.md#cta-antigens), [`oncoref.cta_peptides`](api.md#cta-antigens) | HPA-derived CTA facts, patient coverage, CTA-specific 9-mer counts and load |
| Antigen panels | [`oncoref.antigen_coverage`](api.md#generic-antigen-panels) | Coverage calculations for caller-supplied non-CTA gene lists |
| Expression | [`oncoref.expression`](api.md#expression-and-normalization), [`oncoref.expression_builders`](api.md#expression-and-normalization) | Source-matrix ingestion, per-sample accessors, percentiles, representatives, within-sample summaries, and reference-expression accessors |
| Normalization | [`oncoref.normalization`](api.md#expression-and-normalization), [`oncoref.gene_families`](api.md#expression-and-normalization) | Clean TPM, housekeeping normalization, technical-RNA filtering, normalization/QC reference families |
| Genes and proteoforms | [`oncoref.gene_ids`](api.md#genes-and-proteoforms), [`oncoref.genome`](api.md#genes-and-proteoforms), [`oncoref.proteoforms`](api.md#genes-and-proteoforms) | Gene ID resolution, Ensembl lookup, proteoform grouping |
| Other references | [`oncoref.tmb`](api.md#burden-tmb-fusions-and-signatures), [`oncoref.incidence`](api.md#burden-tmb-fusions-and-signatures), [`oncoref.fusions`](api.md#burden-tmb-fusions-and-signatures) | TMB, incidence/mortality burden, defining fusions |
| Legacy compatibility | [`oncoref.response_signatures`](api.md#burden-tmb-fusions-and-signatures) | Transitional historical response-signature surface; new or extended therapy-signature panels belong in pirlygenes |
| Data management | [`oncoref.catalog`](api.md#data-management), [`oncoref.data_bundle`](api.md#data-management), [`oncoref.reference_data`](api.md#data-management), [`oncoref.hpa`](api.md#data-management) | Dataset inventory, download/cache status, HPA reference data |

## Install

```bash
pip install oncoref
```
