# Download antigen evidence

Use `antigen_data.sources()` to find the exact inputs and observations behind the catalog. Each row gives the dataset selector, version, SHA-256 and source URL. An empty version selects the file shipped with the installed package. Use `antigen_data.download(dataset, version=...)` to verify any inventory file against the exact bytes used in curation. `catalog.ensure()` returns a local file and downloads external data on demand. `catalog.read_table()` reads delimited files, Excel sheets, GCT and Parquet in their native units.

```python
from oncoref import antigen_data, catalog, target_selection

sources = antigen_data.sources()
ct83 = target_selection.results("CT83").entity_id.item()
normal = antigen_data.expression("normal", entity_ids=[ct83])
tumor_hits = antigen_data.expression("tumor", entity_ids=[ct83])
staining = catalog.read_table("marcinkowski2019-ct83-expression")
protected_organs = antigen_data.tissues(tier="protected")
thymus_t_cells = antigen_data.tissue("thymus", "T cell")
lists = antigen_data.gene_lists()
positive_controls = antigen_data.gene_list("benchmark/CTA")
negative_controls = antigen_data.gene_list("benchmark/not-CTA")
coverage_controls = antigen_data.gene_list("benchmark/withhold")  # separate from negative controls
```

| Evidence | What the download contains |
|---|---|
| Normal and reproductive comparison | One protein group, source and site or sample per observation. Scope, native TPM, statistic, QC and locus coverage travel with each row. Use `normal` for the non-thymic, non-reproductive N scope, `thymic` for the separate thymus measurements, `disposable` for D, and `healthy` for the complete tissue comparison, including reproductive tissues. These views retain QC and missing observations. The `reproductive` panel includes embryo and lineage support. |
| Tumor RNA | The tumor witness table contains retained above-background observations. Use the inventory’s versioned `per-sample:<cohort>` entries for complete native matrices, including low and zero values. |
| Normal single-cell, single-nucleus RNA and IHC | `normal_cell_data.observations()` exposes supported assays with source identifiers, units, populations or staining levels. The inventory includes donor-level Census findings, HPA assays, supplementary tables and their upstream sources. Use the native file reader for RDS, Matrix Market, compound archives or study documents. |
| Literature, nomination tables and proteomics lists | `gene_lists()` enumerates published intake tables, reviewed claims and benchmark lists. `gene_list()` returns source-reported labels, resolved identities, evidence and scope. A proteomics nomination retains the paper’s identification scope. |
| CT83 tumor localization | Figure 4 specimen counts for gastric, triple-negative breast and other cancers, with histology, assay and denominator. Cohort positivity and native TPM remain separate measurements. [Marcinkowski et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC6712783/) |

The 452 contributing Census dataset versions and their original H5AD download addresses are available through `antigen_data.cell_atlases()`. OncoRef distributes the donor-level evidence used in grading alongside the version inventory.

## CT83 in the cohorts studied

The published RNA in situ measurements identify the cells expressing CT83. An independent lookup of QC-passing libraries for the exact CT83 locus (ENSG00000204019) gives the following native-TPM results. Each cohort is counted separately. The 10× threshold is 86.0307 TPM under this policy.

| Cohort | Measured libraries | At ≥1 TPM | At ≥10× | Peak TPM |
|---|---:|---:|---:|---:|
| Gastric adenocarcinoma | 413 | 225 | 2 | 88.92 |
| Triple-negative breast cancer | 157 | 96 | 1 | 90.35 |
| Cervical cancer | 306 | 61 | 10 | 193.17 |
| Lung adenocarcinoma | 515 | 257 | 6 | 142.12 |
| Lung squamous carcinoma | 498 | 157 | 1 | 90.57 |
| Pancreatic adenocarcinoma | 178 | 41 | 0 | 24.24 |

The RNAscope series counted tumors containing positive cells. The bulk 10× column counts whole libraries crossing a quantitative threshold. Individual samples from both series have their own denominators. [Study, Figure 4](https://pmc.ncbi.nlm.nih.gov/articles/PMC6712783/#Fig4)

## Literature controls

The benchmark is a test of the enrichment and tissue-coverage rules. CT83 and MAGEC1 have published positive evidence. OR6F1, OR52D1, OR2T4 and OR7E24 have direct adult olfactory observations. OR4M2B appears in the separate `benchmark/withhold` list. Its normal home lacks locus-resolved quantitative coverage. Published OR4M2 rows refer to a different human locus. Each control carries a citation and a qualification. The classifier reads expression, coverage and shared eligibility independently of these lists.

The healthy reference uses the linearly interpolated P95 within each measured tissue or sorted blood population where native TPM distributions exist. Replicate libraries from an identified donor receive one vote. Anonymous HPA biopsies receive one vote per biopsy. Small panels can have a P95 close to their largest observation. Summary-only inputs retain their published mean or median. DICE supplies unpaired specimen ranks, so a multi-locus group's conservative upper bound is recorded explicitly. Adult olfactory biopsy FPKM and donor qPCR copy estimates also have source-native P95 summaries through `normal_cell_data.olfactory_bulk_percentiles()`. Their units stay outside the TPM denominator. Reproductive and tumor samples can each establish positive enrichment at 10×. [Source coverage and measurements](audits/target-selection-20261110/rebuild/normal-rna/index.md)
