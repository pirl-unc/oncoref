# oncoref

`oncoref` is the base-layer package for shared cancer reference data:
cancer-type ontology, cohorts, expression, clean-TPM normalization, TMB,
incidence/mortality, ICI response, HPA normal-tissue expression, and
HPA-derived cancer-testis antigen references.

Downstream packages such as pirlygenes and trufflepig should delegate
parity-clean shared primitives here, but they may keep curated package-specific
tables, generated artifacts, and compatibility wrappers until a surface has a
clear oncoref contract.

Use this page as the quick orientation. The [API guide](api.md) is the detailed
module map, with examples and notes about compatibility modules.

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

ici_response.best_available_ici_response("COAD_MSI")
```

## Data Domains

| Domain | Public modules | Use for |
| --- | --- | --- |
| Cancer vocabulary | [`oncoref.cancer_ontology`](api.md#cancer-vocabulary), [`oncoref.cohorts`](api.md#cancer-vocabulary) | Registry records, aliases, hierarchy, subtype axes, cohort IDs, matched normal tissues |
| ICI response | [`oncoref.ici_response`](api.md#ici-response) | Anti-PD-1 and broader ICI references, regimen-aware lookups, extracted endpoint estimates |
| CTA references | [`oncoref.cta`](api.md#cta-antigens), [`oncoref.cta_coverage`](api.md#cta-antigens), [`oncoref.cta_peptides`](api.md#cta-antigens) | CTA gene sets, patient coverage, CTA-specific 9-mer counts and load |
| Antigen panels | [`oncoref.antigen_coverage`](api.md#generic-antigen-panels) | Coverage calculations for explicit non-CTA gene panels |
| Expression | [`oncoref.expression`](api.md#expression-and-normalization), [`oncoref.expression_builders`](api.md#expression-and-normalization) | Per-sample, percentile, representative, within-sample, and reference-expression accessors |
| Normalization | [`oncoref.normalization`](api.md#expression-and-normalization), [`oncoref.gene_families`](api.md#expression-and-normalization) | Clean TPM, housekeeping normalization, technical-RNA filtering, gene-family panels |
| Genes and proteoforms | [`oncoref.gene_ids`](api.md#genes-and-proteoforms), [`oncoref.genome`](api.md#genes-and-proteoforms), [`oncoref.proteoforms`](api.md#genes-and-proteoforms) | Gene ID resolution, Ensembl lookup, proteoform grouping |
| Other references | [`oncoref.tmb`](api.md#burden-tmb-fusions-and-signatures), [`oncoref.incidence`](api.md#burden-tmb-fusions-and-signatures), [`oncoref.fusions`](api.md#burden-tmb-fusions-and-signatures), [`oncoref.response_signatures`](api.md#burden-tmb-fusions-and-signatures) | TMB, incidence/mortality burden, defining fusions, response signatures |
| Data management | [`oncoref.catalog`](api.md#data-management), [`oncoref.data_bundle`](api.md#data-management), [`oncoref.reference_data`](api.md#data-management), [`oncoref.hpa`](api.md#data-management) | Dataset inventory, download/cache status, HPA reference data |

## Install

```bash
pip install oncoref
```
