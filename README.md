# oncoref

[![Tests](https://github.com/pirl-unc/oncoref/actions/workflows/tests.yml/badge.svg)](https://github.com/pirl-unc/oncoref/actions/workflows/tests.yml)
[![PyPI](https://img.shields.io/pypi/v/oncoref.svg)](https://pypi.org/project/oncoref/)

Curated cancer reference data — cancer-type ontology, tumor mutational burden
(TMB), incidence/mortality, checkpoint-inhibitor (ICI) response, per-cohort
RNA-seq expression, Human Protein Atlas (HPA) normal-tissue and cancer expression, and
HPA-derived cancer-testis antigen references, plus matched CPTAC RNA/protein
calibration sources — behind one small Python API, a data fetch/cache CLI, and a
set of reference plots.

The [holistic CTA/oncofetal membership set](docs/audits/cta-membership-20261026/index.md)
contains **1,416 candidate genes** after reprocessing all **20,097 canonical genes**.
Use `oncoref.cta_membership.gene_ids()` for this biological set, or
`oncoref.cta_membership.results()` for its evidence tiers and retained normal-expression
findings. PRAME, MAGE-A4 and both NY-ESO-1 loci (CTAG1A/B) are included, together
with additional established CTA, cancer-placental and oncofetal controls.
Membership describes biological association; it does not establish target safety,
HLA presentation, locus-specific protein evidence or clinical efficacy.

The set contains 21 specifically reviewed established controls, 1,047 candidates
supported by distinct cancer/reproductive publications, and 348 expression-supported
candidates. Classes overlap and are labelled separately. Oncofetal controls include
TBXT/brachyury, AFP, GPC3, CLDN6 and SALL4; developmental atlas measurements have
not yet been imported. Normal RNA, IHC, cell-specific expression, clinical concerns
and missing evidence remain visible, including 726 members with the earlier
expression-forbidden flag. That sequence reference has not been changed.

The earlier [215-gene expression/clinical shortlist](docs/audits/cta-research-20261026/index.md)
remains available through `oncoref.cta_research_policy`, and the
[213-gene expression policy](docs/audits/cta-policy-20261025/index.md) remains the
released `target_sets` classification. These are different views with different
purposes. [Normal-expression breadth](docs/audits/normal-breadth-20261026/index.md)
remains descriptive; the membership report adds its fixed anatomical score for
all canonical genes, including missing-coverage bounds. Donor variation is
individually assessed only where imported evidence exists, currently CTAG2.
[Reproductive tissue definitions](docs/normal-tissue-scopes.md) and
[developmental sources](docs/developmental-targets.md) document scope and gaps.
Comprehensive individual literature review remains incomplete.

The released `oncoref.target_sets` API partitions genes and exact proteins into
**CTA**, **Unsafe** and **Neither**. Only Unsafe supplies the forbidden protein
background. Admission requires reproductive RNA ≥5 nTPM, reliable reproductive
protein, or reviewed reproductive plus tumor/recognition evidence. Accepted
normal-protein findings exclude; material unresolved conflicts hold. The
[plain-language definitions](https://pirl-unc.github.io/oncoref/target-sets/)
give the full hierarchy and unchanged numerical thresholds. The current build
contains **213 CTA genes** and **359 distinct complete CTA protein sequences**.
Use the matching [v1.20.0 sequence assets](https://github.com/pirl-unc/oncoref/releases/tag/v1.20.0)
and their k-mer exclusion recipe. `oncoref.cta_policy.manifest()` identifies the active policy.

The [genome-wide expression groups](https://pirl-unc.github.io/oncoref/audits/normal-expression-categories-20261021-r2/)
classify Reproductive, Normal and Other independently of CTA nomination.
Other contains Disposable and Silent; incomplete evidence remains Unassessed.
Observed normal expression is preserved separately for subsequence comparisons.
[Disjoint gene and exact-protein sets](https://pirl-unc.github.io/oncoref/audits/normal-sequence-sets-20261021/)
include all complete Ensembl 115 isoforms, source records and FASTA downloads.
Shared sequences use Normal → Unassessed → Reproductive → Other precedence;
k-mers still require explicit sequence subtraction.

## Role in the stack

The openvax/PIRL tools are split by ownership boundary, not by file format.
`oncoref` is the base layer: downstream packages can depend on it, but it never
imports its consumers.

- **oncoref** is the empirical base: what is true, measured, or canonically named
  about cancer types, cohorts, genes, and reference datasets. It owns gene
  identity and canonicalization, the cancer-type ontology/registry, expression
  reference data and normalization, epidemiology, TMB, ICI/anti-PD-1 response, and
  source-anchored cancer-testis antigen (CTA) facts. If a row has an `n`, a
  confidence interval, a source cohort, or a PMID/DOI anchoring a measurement,
  it is usually an oncoref fact.
- [**pirlygenes**](https://github.com/pirl-unc/pirlygenes) owns curated gene
  sets and panels: which genes are useful for a purpose. That includes
  lineage/family/compartment/supertype panels, discriminators, surfaceome, tumor
  microenvironment (TME) and stem-cell markers, response-signature panels,
  target-to-therapy registries, and other opinionated selections keyed to
  oncoref cancer codes and gene IDs. An empty set can be a valid pirlygenes
  answer.
- [**trufflepig**](https://github.com/pirl-unc/trufflepig) owns per-sample
  interpretation: quality-control (QC) narration, library-prep/source warnings,
  deconvolution, scoring, and rule tables that fire against one tumor sample.

Adoption is staged: consumers can delegate parity-clean primitives while
retaining their own curated artifacts and compatibility APIs. Shared facts and
identifiers should be fixed here; purpose-specific panels and per-sample rules
stay downstream.

## Species

Use explicit species and dataset modules:

```python
import oncoref

oncoref.human.human_protein_atlas.rna_consensus()
oncoref.canine.barkbase.source()  # native TMM-CPM and provenance
oncoref.canine.references.registered_references()
oncoref.canine.antigens.coverage()

# Datasets also have direct top-level imports, with the same implementation.
from oncoref import barkbase, human_protein_atlas

dataset = barkbase.read(offline=True)  # verified cache or configured local overrides
```

`oncoref.human` and `oncoref.canine` are ordinary Python modules with named datasets
and references, expression, cancer, antigens and normalization groups. Dataset
modules expose source-specific operations; groups expose inventories and coverage.
Missing coverage is explicit; canine requests never fall back to human datasets.
Existing flat APIs retain their human defaults, and `for_species(...)` remains
available for dynamic selection. See the
[species and expression contract guide](docs/species.md).

## Core model

- **Canonical identities first.** Cancer facts key on the cancer-type registry;
  expression cohorts and evidence scopes are explicit rather than inferred from
  names. Gene APIs resolve to a canonical Ensembl space.
- **Small facts in the wheel, large matrices on demand.** Ontology, gene,
  burden, response, and provenance tables install with the package. The large
  expression bundle and HPA datasets use versioned download caches.
- **Provenance is part of the result.** Reference rows expose source scope,
  sample counts, versions, and review status. Missing or ineligible data stays
  distinguishable from a valid empty result.
- **Consumer policy stays downstream.** oncoref exposes HPA-derived CTA calls
  and empirical expression facts; it does not turn them into therapy panels or
  one-sample decisions.

## Install

Supported Python versions: 3.10–3.14.

```bash
pip install oncoref
```

Reference plotting is included in the standard installation. The pyensembl-backed
genome integration remains optional because it also requires a separately installed
Ensembl release:

```bash
pip install 'oncoref[genome]'  # pyensembl-backed transcript and gene lookup
```

The genome extra requires PyEnsembl 2.27.0 or newer. When upgrading from PyEnsembl
older than 2.10.17, rerun `pyensembl install --release N --species human` for each
installed GRCh38 release from 82 onward to build the complete patch/haplotype
annotation index. Existing FASTA files are reused.

## Quick start

The flat `oncoref` namespace remains available for compatibility and quick
interactive use. For new code, prefer the semantic submodules in
the [API guide](https://pirl-unc.github.io/oncoref/api/); they make it clearer
whether you are working with the cancer ontology, cohorts, ICI response, CTA
coverage, generic antigen-panel coverage, or CTA-specific peptides.

```python
import oncoref as od

od.resolve_cancer_type("prostate")        # -> "PRAD"
od.cancer_type_info("SARC_RMS_ARMS")      # full registry record + burden + tmb
od.cancer_tmb("FL")                       # 5.05 (source-checked panel median)
od.cancer_burden("pancreas", metric="us_mortality_pct")
od.burden_category("SARC_OS")             # -> "bone_and_joint" (incidence/mortality bucket)
od.cancer_ici_response("SKCM")            # 43.7% in the selected trial setting
od.cancer_ici_response("SKCM", regimen="PD-1+CTLA-4")   # 57.6  (pin a regimen)
od.therapy_benefit_toxicity_evidence(cancer_code="OV")  # sourced clinical facts
od.rna_protein_calibrations(gene="TP53")     # matched CPTAC cohort-specific models
od.rna_protein_hpa_priors(gene="TP53")       # HPA normal-tissue RNA/IHC weak prior
od.hpa_cancer_ihc_prevalence()               # scored-patient IHC counts, HPA 25.1
od.hpa_cancer_rna_prevalence(cohort="TCGA")   # genome-wide native-pTPM RNA summaries

# Cancer-testis antigens (HPA-derived tissue-restriction):
od.cta_gene_names()                       # expressed CTA symbols (MAGEA4, CT83, …)
od.cta_evidence()                         # full HPA restriction table
od.cta_clinical_target_evidence()         # explicit clinical/canonical tier + leak flags

# Per-cohort expression percentiles (downloads the data bundle on first use):
od.cohort_gene_percentiles("PRAD")        # per-gene p0…p100 vector (within-cohort)
od.within_sample_top_fraction("PRAD")     # per-gene frac of samples top-5% (within-sample)
```

## Domain map

| Question | Preferred modules |
| --- | --- |
| What does this cancer code mean? | `oncoref.cancer_ontology`, `oncoref.cohorts` |
| What expression reference is available? | `oncoref.expression`, `oncoref.source_matrices` |
| What matched RNA/protein calibration is available? | `oncoref.rna_protein` |
| What tumor-attributed or subtype summary is available? | `oncoref.tumor_references` |
| How is expression normalized or filtered? | `oncoref.normalization`, `oncoref.gene_families` |
| What is the canonical gene identity? | `oncoref.gene_ids`, `oncoref.genome`, `oncoref.proteoforms` |
| What is the TMB, burden, or ICI response evidence? | `oncoref.tmb`, `oncoref.incidence`, `oncoref.ici_response` |
| Which genes meet the HPA-derived CTA definition? | `oncoref.cta`, `oncoref.cta_coverage`, `oncoref.cta_peptides` |
| Where is a dataset cached or downloaded? | `oncoref.catalog`, `oncoref.data_bundle`, `oncoref.hpa` |

The [API guide](https://pirl-unc.github.io/oncoref/api/) explains each domain,
its provenance contract, and the distinction between reader, builder, and
compatibility APIs.

## Command line

See the [TMB/ICI/aPD1 curation audit](docs/curation-audit.md) for corrected values,
evidence gaps, and the complete row-by-row review ledger.

```bash
oncoref cancer-type prostate     # registry info as JSON
oncoref tmb FL                  # 5.05
oncoref ici SKCM                # 43.7 (--regimen to pin, --all-regimens to compare)
oncoref burden pancreas --metric us_mortality_pct
oncoref cta --count             # number of expressed CTAs
oncoref plot apd1-vs-tmb --out apd1_vs_tmb.png
oncoref plot patient-coverage --gene-set cta --out coverage_out
oncoref plot cta-curation --out cta_curation_out
```

### Data cache

```bash
oncoref data list               # every wheel/bundle/HPA/source dataset
oncoref data status bundle      # expression-bundle cache state (no download)
oncoref data metadata           # package/data/cache/release contract JSON
oncoref data fetch bundle       # download the large expression bundle
oncoref data fetch hpa          # HPA RNA / immunohistochemistry / single-cell data
oncoref data dir bundle         # where the data bundle is cached
oncoref data prune --yes        # delete stale bundle version caches
oncoref version
```

HPA files, expression bundles (including overlays and release metadata), and
ordinary expression build-source files use `datacache>=1.19.0` for streaming
and bounded HTTP retries. Per-cohort matrices retain their locked transfer and
strict HTTP length check. Existing cache paths, source versions,
archive/content checksums, and bundle validation remain in place.

Explicit fetches report the dataset/version, source URL, destination, cache
reuse, and completion size and average transfer rate. Datacache transfers in
terminals and notebooks also show byte progress and transfer rates; a
percentage/ETA appears only when the server supplies a usable total. Redirected
output gets status lines without animated bars. Per-cohort matrices report
status without a bar. Status and terminal bars go to stderr, keeping stdout
usable for paths and command results.

```bash
oncoref data fetch hpa --progress always  # force download bars
oncoref data fetch bundle --progress never  # status lines only
oncoref data fetch per-sample:LUAD --quiet  # results/errors, no status
ONCOREF_DOWNLOAD_PROGRESS=never python analysis.py  # applies to Python calls too
```

`--progress auto|always|never` overrides `ONCOREF_DOWNLOAD_PROGRESS` (default
`auto`). Existing `verbose=False` options, plus `reference_data.download(...,
verbose=False)`, suppress status and bars. Retry warnings retain datacache's
standard Python logging: transient HTTP failures get at most three attempts,
with attempt counts and delays reported. Failed transfers preserve previous
cache files. These controls do not change the specialized resumable SRA/FASTQ
builder or GDC/API query implementations.

## Development

```bash
./develop.sh   # create/refresh .venv and its editable install with dev extras
source .venv/bin/activate
oncoref version
./format.sh    # ruff format
./lint.sh      # ruff check + format --check
./test.sh      # lint + the full pytest suite with coverage and two workers
```

The repository virtualenv is the canonical development CLI. It imports the
checkout directly, so source changes are visible without reinstalling; rerun
`./develop.sh` when dependencies or packaging metadata change. To use the same
CLI outside the checkout, add this shell alias rather than maintaining a
separate global pip installation:

```bash
alias oncoref=/absolute/path/to/oncoref/.venv/bin/oncoref
```

## License

Apache 2.0.

### Annotation-scoped CTA identities

When using a particular Ensembl annotation, resolve CTA source IDs with
`oncoref.resolve_gene_identity(gene_id, genome=genome)` and obtain verified
negative-reference exclusion IDs with
`oncoref.cta_annotation_gene_ids(genome, unfiltered=True)`. Both use the same
versioned identity contract, retaining original source IDs and unresolved
mapping evidence. See [CTA gene identities](docs/cta-gene-identities.md).
