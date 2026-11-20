# Antigen lookup and curation builds

The end-user commands read the current catalog (`oncoref.target_selection`). The
maintainer command rebuilds the catalog, figures, and reports from the declared
evidence inputs in topological order.

```sh
# Default CTA membership across all reproductive compartments
oncoref antigens list --set CTA

# Host anatomy & tissue-tolerance profiles
oncoref antigens list --set CTA --profile male-patient
oncoref antigens list --set CTA --profile male-patient --absent-tissue prostate --absent-tissue "seminal vesicle"
oncoref antigens list --set CTA --profile female-patient
oncoref antigens list --set CTA --profile female-patient --lactating
oncoref antigens list --set CTA --profile female-patient --protect-tissue ovary

# Combine with filter modes, fold cutoffs, or patient replication
oncoref antigens list --set CTA --profile female-patient --mode standard
oncoref antigens list --set CA --min-tumor-fold 30 --min-patients 2 --format csv

# Inspect every gate, mode, and host profile for one gene
oncoref antigens explain LALBA --loss
oncoref antigens explain CTAG2 --format json

# Sequential catalog + filter attrition funnel
oncoref antigens funnel --profile male-patient --mode standard

# Exact or bounded-distance peptide subtraction against Unsafe|Other
oncoref antigens subtract protein.fa --background 'Unsafe|Other' --format csv
oncoref antigens subtract protein.fa --equiv ILV --max-distance 1
```

`--set` accepts `CTA`, `CA`, `RA`, `CA|RA`, `Unsafe`, `Other`, and `Unsafe|Other`.
`--profile` accepts `all` (no host screen), `female-patient` and `male-patient`
(assumed anatomy, with no reproductive tissue protected by default). Select
`--protect-tissue testis` or `--protect-tissue ovary` to request smaller sets.
`--accept-loss` (legacy alias
`--allow-tissue`) explicitly accepts tissue loss; `--absent-tissue` records absent
anatomy/state; `--protect-tissue` overrides an exception. `breast` expands to both
mammary states. `--pregnant` also protects the milk-secretory state, and pregnancy
remains unresolved because developmental coverage is incomplete. The old `core`
host shortcut is rejected. [Concepts and exact rules](host-selection.md).

`--ids` exports gene IDs. Background sets remain global when listed directly;
`subtract --profile ...` adds candidates excluded or unresolved under that host
context to the chosen background. List, explain and subtraction emit selection
provenance to stderr. `explain --profile ...` reports all host reasons. The legacy
`oncoref cta` command retains its historical snapshot and migration notice.

Subtraction requires locally installed Ensembl 115 translations and stops if any
background gene lacks a translation. Defaults are exact 8–15mers against
`Unsafe|Other`. Equivalence groups are explicit, disjoint assumptions; they are not a
universal recognition model. Weighted edits use `--substitution-cost I:L:0.5`,
`--insertion-cost`, and `--deletion-cost`. JSON provenance is written to stderr;
retained peptide occurrences go to stdout. Sequence retention does not demonstrate
peptide presentation or therapeutic safety.

## One maintainer entry point

From a complete source checkout:

```sh
python -m oncoref.curation fetch --include-tumors
python -m oncoref.curation build --reuse-tumor-cache
python -m oncoref.curation check
# Equivalent verification-only spelling:
python -m oncoref.curation build --check
python -m oncoref.curation diff --against prior-categories.csv.gz
```

`fetch` downloads and SHA-256-verifies every external reference dataset declared by the build (HPA v23, GTEx v11 individual tissue matrices, DICE, Monaco, and Petropoulos/Tyser embryo tables); pass `--include-tumors` to also fetch and verify all 132 pinned per-sample tumor matrices through `oncoref.catalog`. No out-of-band downloads or machine-local file paths are required.

The declared dependency graph (`oncoref/data/curation-build.json`) orders the main build as `rare_cell` → `catalog` → `host_exposure` → `plots` → `report` → `filter_audit` →
`filter_web_report`. The rebuilt catalog also feeds the GTEx witness comparison, literature audit and IHC audit. Failures stop dependent steps. Omitting `--reuse-tumor-cache`
re-extracts tumors; supplying it verifies the tumor cache receipt and input hashes first.
`check` verifies implementation and output hashes across all 210+ tracked files without
downloading or reprocessing raw data. `diff` uses stable protein identities and highlights
named regression controls (`PRAME`, `NY-ESO-1`, `MAGEA4`).

## Migration boundaries

The current build consumes the pinned normal-screen fact table and its replayed
adjudications. GTEx v11 individual gene TPM is registered under `gtex_individual_tissue_tpm` with a SHA-256 pin for each of the 52 primary tissue matrices. The build uses those primary matrices for two measurements. Non-reproductive tissue sites produce the derived native-TPM source `gtex_donor_p95:v11`. Each donor contributes one value, and its interpolated tissue percentile carries no individual sample ID. Reproductive tissues provide individual witnesses with donor and sample IDs. See the [source-by-tissue inventory and denominator figures](audits/target-selection-20261110/rebuild/normal-rna/index.md). Other inputs include checksum-verified HPA v23 individual tissue biopsies and
Petropoulos/Tyser human embryo single-cell datasets via `curation.reference_input(role)`.
Historical receipts and public compatibility modules (`cta_selection`, `target_sets`,
`cta`) remain frozen snapshots and are never overwritten by live policy updates.

Unresolved locus attribution remains an explicit research hold (`Other:Review`), not a
code failure. For reviewed tissue-loss consequences alongside observed reproductive RNA,
use `antigen_sets.loss_consequences(gene)` or `tissue_tiers.consequences(tissue,
cell_type, resolution="cell")`.

## Cell and isoform evidence

For receptor-level measurements from adult human nasal tissue, see [Adult olfactory expression](#adult-olfactory-expression).

Use `cta_literature.paper_gene_claims(gene)` for the observation, assay, tissue or developmental stage, and permitted use. For studies with downloadable artifacts, the claim’s `provenance["reference_sources"]` lists the keys for `reference_data.ensure(key)`. These studies supply cell localization and molecular identity. A scoped `unscreened_adult_home` in a reviewed human gene claim supplies anatomy to the coverage check. GFY has a documented adult olfactory home, so its tumor enrichment requires adequate native-TPM coverage there. The CA and RA fold calculations require native TPM.

| Gene | Primary source | What the data establish |
|---|---|---|
| GFY | [Adult olfactory biopsies, GSE139522](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE139522) | Detected in 94.1% of mature and 70.3% of immature olfactory neurons, pooled across the sampled cells. |
| AMELX | [Developing teeth, GSE184749](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE184749) | The primary dental epithelium count matrix covers 9–22 gestational weeks. At 20–22 weeks, 174 of 223 cells in the two ameloblast clusters expressed AMELX. |
| SP7 | [Embryonic skeleton, E-MTAB-14385](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-14385) · [Primary osteoblasts, GSE147390](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE147390) | Cranial osterix staining at 8.5 and 10 post-conception weeks. In a separate adult study, SP7 RNA was enriched in mature femoral-head osteoblasts from one donor with osteoarthritis and osteopenia. |
| CER1 | [Peri-implantation embryos, E-MTAB-8060](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-8060) | Days 9 and 11 yielded hypoblast CER1 RNA. Day-9 staining localized protein. The count matrix supplies UMI measurements, and the published lineage summary identifies the sampled embryos. |
| RPS4Y2 | [RPS4 tissue panel and primer supplement](https://doi.org/10.1186/1471-2199-11-33) | Transcript-directed assays distinguish splice products in testis and prostate. Attribution of the healthy naïve B-cell signal remains open. |
| GALP | [Alarin study](https://pubmed.ncbi.nlm.nih.gov/16954604/) · [HYPOMAP collection](https://cellxgene.cziscience.com/collections/d0941303-7ce3-4422-9249-cf31eb98c480) | Alarin was observed in neuroblastic tumors. HYPOMAP already provides adult hypothalamus coverage in the [atlas inventory](audits/normal-cell-atlas-20261110/source-registry.csv). The product expressed in the catalog tumor and the normal hypothalamus remains unresolved. |

### Adult olfactory expression

Use the public dataset catalog to retrieve the original files and the normal-cell API to compare exact source features. The two Olender tables use the same four biopsies. The two Verbeurgt tables describe the same 26 donors. Each pair counts as one study.

| Dataset ID | Data available downstream |
|---|---|
| `olfactory_olender_transcriptome` | Genome-wide biopsy expression from Olender Table S2, in FPKM. |
| `olfactory_olender_receptors` | Receptor-specific Table S7, with the four unique-read FPKM columns, hg19 locus and a flag for the authors' 0.003 zero replacement. |
| `olfactory_verbeurgt_cq`, `olfactory_verbeurgt_copies` | Receptor assay IDs and donor-level Cq or estimated copies per 20 ng of RNA. |
| `olfactory_durante_supplement` | Reported positive sensory neurons per listed receptor in supplementary Table 7. The original workbook includes additional cell populations and tables. |
| `olfactory_mssm_nuclei` | GSE317163. The full H5AD contains cell annotations and counts. The table reader summarizes the 448 author-labeled sensory-neuron nuclei by donor, retaining total nuclei, positive nuclei and UMI sum. Install `oncoref[olfactory]` for the H5AD reader. |

```python
from oncoref import catalog, normal_cell_data

genes = ["OR6F1", "OR52D1", "OR4M2B", "OR2T4", "OR7E24"]
bulk = normal_cell_data.observations("olfactory_olender_receptors", genes=genes)
nuclei = normal_cell_data.observations("olfactory_mssm_nuclei", genes=genes)
original = catalog.ensure("olfactory_mssm_nuclei")
```

| Receptor | Olender biopsies ≥0.01 FPKM / 4 | Peak FPKM | Verbeurgt donors ≥5 estimated copies / 26 | Durante positive neurons | MSSM positive sensory nuclei / 448 |
|---|---:|---:|---:|---:|---:|
| OR6F1 | 4 | 1.890 | 26 | 3 | — |
| OR52D1 | 3 | 0.115 | 26 | — | — |
| OR4M2B | — | — | — | — | — |
| OR2T4 | 3 | 0.321 | 15 | — | — |
| OR7E24 | 4 | 3.016 | 25 | 1 | 1 |

A dash means the source lacks an exact-name feature or assay in the table used. A measured zero retains its feature and denominator. OR4M2 has a distinct source label and cannot establish OR4M2B expression without locus-specific attribution. The four replicated mucosal RNA observations document normal tissue exposure. Their native units do not enter the 10× TPM denominator or establish protein-level toxicity.

[Olender bulk RNA](https://doi.org/10.1186/s12864-016-2960-3) · [Verbeurgt qPCR](https://doi.org/10.1371/journal.pone.0096333) · [Durante single-cell study](https://doi.org/10.1038/s41593-020-0587-9) · [MSSM single-nucleus atlas](https://doi.org/10.1038/s41467-026-75722-1)

### What downstream users can retrieve

`catalog.datasets()`, `catalog.versions(name)`, `catalog.ensure(name)` and `catalog.provenance(name)` cover registered assets. The data format and level of processing depend on the source. HPA normal-tissue IHC has scored observations, and the antibody XML adapter exposes antibody and reliability information. A paper’s curated gene list is also queryable evidence. For example, the da Silva tumor-proteomics list records nominations; its registry entry supplies no underlying raw mass spectra or peptide-level quantification. The HPA-distributed GTEx tissue TPM table is available as `catalog.read_table("rna_tissue_gtex", version="v23")`. The source inventory connects native assays, imported lists, exact cohort releases and the observations used in classification. [Download antigen evidence](antigen-data.md).

## Tissue and gene-list lookup

`antigen_data.tiers()` lists injury-tier definitions. `antigen_data.tissues(tier="protected")` lists organs, parent IDs and assessment IDs. `antigen_data.tissue("thymus", "T cell")` returns the matching cell rule, the expected loss and its sources. `antigen_data.gene_lists()` lists stable IDs for publication claims, benchmarks and marker sets. Pass one ID to `antigen_data.gene_list(list_id)` for its members and assay scope.

[Organ and cell injury hierarchy](tissue-hierarchy.md)
