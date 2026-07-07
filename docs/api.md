# oncoref API Guide

`oncoref` keeps its historical flat top-level imports for compatibility, but new
code should prefer the semantic submodules below. They make the domain boundary
clear and avoid guessing whether a broad name such as `coverage` or `peptides` is
general or CTA-specific.

Package boundary: oncoref is the upstream home for empirical base facts and
canonical identifiers that are ready to be reused across the PIRL stack.
pirlygenes owns purpose-specific gene sets and panels; trufflepig owns
per-sample interpretation, QC narration, and rule firing. As a rule of thumb,
source-anchored measurements with denominators, confidence intervals, cohorts,
PMIDs/DOIs, or shared ontology implications belong in oncoref. Opinionated gene
selections and target-to-therapy registries belong in pirlygenes. One-sample
rules belong in trufflepig. When a missing data field, gene universe,
bundle-integrity rule, or source-QC decision affects shared reference artifacts,
the durable fix should live or be exposed here rather than only in a downstream
compatibility layer.

## Cancer Vocabulary

- `oncoref.cancer_ontology` — cancer-type registry, aliases, parent/child tree,
  lineage/family groupings, molecular subtype axes, matched normal tissues,
  source-scoped evidence resolution, and display helpers.
- `oncoref.cohorts` — expression/source cohort IDs, computed aggregate cohorts,
  source versions, and mixture-cohort flags.

Use these when asking "what cancer type or cohort does this code mean?"
Prefer the DataFrame-returning query helpers when code will be passed into
other oncoref domains; they keep the result type and columns stable.

```python
from oncoref import cancer_ontology, cohorts, expression

cancer_ontology.resolve_cancer_type("prostate")
cancer_ontology.cancer_type_tree("CRC")
cancer_ontology.cancer_type_path("COAD_MSI")

# CRC plus anatomical children and molecular leaves.
crc = cancer_ontology.cancer_type_records(under="CRC")
crc["code"].tolist()

# Cross-cutting molecular axes can be intersected with hierarchy or lineage.
msi_crc = cancer_ontology.cancer_type_records(subtype_group="MSI", under="CRC")
epithelial_msi = cancer_ontology.cancer_type_records(
    subtype_group="MSI", lineage_group="Epithelial"
)

# COAD_MSI / READ_MSI keep anatomical expression context but resolve evidence
# rows through CRC_MSI when published sources are colorectal-level.
msi_crc[["code", "evidence_source_code", "normal_tissue_code", "hpa_tissues"]]

# Join scalar references for the returned codes.
cancer_ontology.cancer_type_reference_data(msi_crc)

# Ask whether each ontology node has a direct expression reference or only
# molecular/no expression evidence.
cancer_ontology.expression_reference_coverage(subtype_group="MSI", under="CRC")
cancer_ontology.coverage_for_cancer_type("ASTB")

# Use codes directly with expression accessors.
codes = cancer_ontology.cancer_type_codes(subtype_group="MSI", under="CRC")
expression.cancer_reference_expression(codes)

# Matched normal RNA expression is an explicit HPA read.
cancer_ontology.matched_normal_tissue_expression("COAD", genes=["ENSG00000141510"])

cohorts.cohort_registry_df()
```

`expression_reference_coverage()` is the ontology-wide readiness table for
classifier consumers. It reports direct observed-bulk source-matrix coverage,
matched normal tissue availability, molecular/fusion-only definitions, canonical
gene/proteoform space, data/source-matrix versions, and a conservative
`consumer_recommendation`: `direct_reference`, `molecular_only`, or
`unsupported`. It intentionally does not synthesize marker-program or
discriminator fallbacks; those remain consumer-layer choices in packages such as
trufflepig.

## Gene Identity

- `oncoref.gene_ids` — canonical ENSG space, alt-haplotype / retired Ensembl ID
  migration, symbol/synonym resolution, and report-facing gene labels.
- `oncoref.genome` — pyensembl-backed transcript/gene lookup and transcript to
  gene aggregation for source matrices.

Use the gene-id helpers before building expression artifacts or joining
downstream gene sets to oncoref references. `canonical_gene_id()` is the primary
any-identifier entry point for the shipped ENSG + symbol/synonym space: it
normalizes versioned or case-varied Ensembl gene IDs, follows retired/alt ENSG
aliases into the canonical space, resolves symbols and synonyms, and returns
`None` for inputs that cannot be mapped to a canonical oncoref gene.
`canonical_gene_symbol()`, `display_gene_name()`, and `short_gene_name()` use
the same resolver so report code does not invent a separate symbol mapping.
`entrez_gene_mappings()` and `resolve_entrez_id()` expose the filtered NCBI
Entrez/GeneID table used by the resolver; it covers live IDs from NCBI dbXrefs
or current symbols plus discontinued IDs redirected through NCBI gene_history.
`gene_identifier_mapping_coverage()` and
`gene_identifier_mapping_summary()` make the shipped ENSG, symbol/synonym, and
Ensembl-alias coverage explicit for migration audits, including non-unique
symbols and missing-symbol rows. They do not claim that RefSeq or UniProt
coverage is complete.

```python
from oncoref import canonical_gene_id, canonical_gene_symbol, display_gene_name, gene_ids

canonical_gene_id("GNB2L1")        # previous symbol -> ENSG00000204628
canonical_gene_id("7157")          # Entrez/GeneID -> ENSG00000141510
canonical_gene_symbol("GNB2L1")    # previous symbol -> RACK1
display_gene_name("ENSG00000005955")  # retired Ensembl id -> GGNBP2
gene_ids.gene_identifier_mapping_summary()
```

## ICI Response

- `oncoref.ici_response` — checkpoint-inhibitor response anchors, anti-PD-1
  shortcuts, regimen-aware lookups, extracted endpoint estimates, and pooled
  response summaries.

`DEFAULT_ICI_REGIMEN_PRIORITY` is the unpinned regimen priority
(`PD-1`, then `PD-L1`, then `PD-1+CTLA-4`). The older
`REGIMEN_FALLBACK` name remains available in `oncoref.ici` for compatibility.

```python
from oncoref import ici_response

ici_response.apd1_response("SKCM")
ici_response.best_available_ici_response("SARC_ASPS")
ici_response.ici_response_by_regimen("SKCM")
ici_response.ici_response_estimates_df()
```

## CTA Antigens

- `oncoref.cta` — CTA definition, HPA restriction tiers, axes, aliases, and gene
  ID/name sets. Strict helpers such as `cta_gene_names()` and
  `cta_filtered_gene_names()` preserve the HPA reproductive-restriction default;
  `cta_clinical_target_evidence()` exposes a separate clinical/canonical tier for
  source-anchored CTA targets that may be strict-pass, HPA-excluded, or
  candidate-only. `cta_specificity_audit()` exposes machine-readable specificity
  demotion and candidate-only decisions for genes whose normal-tissue evidence
  makes strict-default inclusion unsafe or unresolved.
- `oncoref.cta_coverage` — CTA patient coverage over per-sample expression
  matrices.
- `oncoref.cta_peptides` — CTA-specific 9-mer counts and load.

`cta_specific_9mer_count_map()` returns a map from a join key to
`n_specific_9mers`; those counts are used as weights when computing
`cta_specific_9mer_load()`.

The CTA definition here is the HPA-derived tissue-restriction call. Broader
therapy-target curation, MS evidence, and downstream prioritization rules can
live in consumer packages while they remain package-specific.

```python
from oncoref import cta, cta_coverage, cta_peptides

cta.cta_gene_names()
cta.cta_clinical_target_evidence()
cta.cta_specificity_audit()
cta_coverage.cta_addressable_fraction("LUAD")
cta_peptides.cta_specific_9mer_count_map(by="proteoform_key")
```

## Generic Antigen Panels

- `oncoref.antigen_coverage` — coverage helpers for caller-supplied gene lists.

Use this when the panel is not necessarily CTA. The function names require
`gene_ids=` so a caller cannot accidentally rely on the CTA default. This module
computes coverage for a supplied list; it does not make oncoref the owner of
downstream panel curation.

```python
from oncoref import antigen_coverage

antigen_coverage.addressable_antigen_fraction("LUAD", gene_ids={"ENSG00000141510"})
antigen_coverage.greedy_antigen_coverage("LUAD", gene_ids={"ENSG00000141510"})
```

## Expression And Normalization

- `oncoref.expression` — read-time accessors for per-sample expression,
  percentile vectors, representative samples, within-sample top fractions, and
  pan-cancer reference tables. `sample_expression_qc` reports per-sample
  detected-gene counts, literal-zero fraction, top-gene/top-10 concentration,
  biological-housekeeping detection, source-scale class, and source-type caveats
  so sparse source-matrix artifacts can be audited before using absolute TPM
  floors or housekeeping normalization. `per_sample_expression(...,
  sample_qc="pass" | "pass_or_warn" | "all")` filters sample columns at read time; the raw per-sample accessor
  defaults to `"all"` for forensic access, while live summaries such as
  `cohort_stats` and `pooled_cohort_stats` default to QC-passing samples.
  `source_matrix_sample_qc_manifest`, `expression_artifact_build_metadata`, and
  `expression_artifact_build_summary` read the optional QC/build metadata emitted
  by regenerated expression bundles. Until a regenerated heavy bundle ships those
  files, they return schema-stable empty metadata by default; use
  `on_missing="raise"` when a downstream migration requires the manifests.
- `oncoref.expression_builders` — build-time ingestion and artifact cores used by
  data-bundle generation scripts. `GeoMatrixSource` /
  `build_source_matrices` own the generic supplementary-matrix path from raw
  source file to canonical per-code per-sample TPM parquet, mapping audit, parse
  diagnostics, and sample-QC sidecars. `geo_matrix_source_from_registry` and
  `scripts/build_geo_matrix.py` make `source_type: geo-matrix` entries in the
  packaged source registry directly buildable. `scripts/rebuild_expression_artifacts.py`
  then applies the same sample-QC policy to derived shards by default
  (`--sample-qc pass`) and emits `source-matrix-sample-qc.csv` plus
  `expression-artifact-build-metadata.*` in the staging directory so bundle
  releases record which source samples fed percentiles, representatives,
  proteoform summaries, and within-sample summaries. Representative
  sample selection uses `representative_sample_columns` / `cohort_medoids` on the
  biological clean-TPM view, then stores the selected samples' full
  clean_tpm_16_9_75 vectors. Release builds retain curated cohorts that have no
  strict QC-pass samples only through explicit source-aware fallbacks recorded in
  the build metadata, and clip invalid negative source expression values to zero
  with per-cohort counts.
- `oncoref.expression_engine` — reusable low-level builder primitives for
  expression tables: identity/value column detection, transcript-to-gene
  aggregation, source row ID-type detection, source gene-row mapping audits,
  missing-vs-non-parsing numeric diagnostics, and canonical ENSG aggregation in
  linear expression space. It is an explicit public module, so downstream
  builders can import `oncoref.expression_engine.map_source_gene_rows`,
  `canonicalize_source_gene_matrix`, and `coerce_source_expression_values`
  without reaching into scripts. Use these in builders before committing a
  source matrix so unresolved high-expression rows and duplicate canonical IDs
  are explicit artifacts rather than hidden cleanup. The source audit frames are
  intentionally unversioned public API objects: provenance belongs in build
  metadata, while the frames themselves use stable canonical columns such as
  `sample_qc_status`, `sample_qc_reasons`,
  `source_expression_nonzero_samples`, and
  `source_expression_sample_with_max`.
- `oncoref.source_matrices` — raw per-cohort source-matrix cache/fetch helpers.
  Use `source_matrices.sample_qc(code)` for the live source-matrix QC audit and
  `source_matrices.sample_qc_manifest(...)` for the optional generated-bundle QC
  manifest that records which samples fed derived artifacts.
- `oncoref.normalization` — TPM conversion, clean TPM, technical-RNA filtering,
  log transforms, percentile ranks, and housekeeping normalization.

The normalization helpers are intended to be reusable directly. Expression
accessors and bundles are also reusable, but downstream packages may keep their
own packaged expression artifacts until row-set, value, provenance, and QC
contracts are parity-clean for the specific accessor they want to replace.

`expression.pan_cancer_expression()` defaults to oncoref's entity-first schema:
HPA normal tissue columns are `<tissue>_nTPM_raw`, TCGA source/provenance
columns are `<CODE>_FPKM_raw`, deterministic TCGA TPM companions are
`<CODE>_TPM_raw`, and analysis columns append `_clean`, `_hk`, `_percentile`, or
`_log1p`. For migration code that needs pirlygenes' unsuffixed column names, use
`column_style="pirlygenes"`; the legacy `to_tpm=True` keyword is accepted as a
compatibility alias for that view and maps the default call to `normalize="tpm"`.

`expression.cancer_reference_expression()` returns cohort-level tumor reference
expression with stable long or wide output. It accepts canonical cancer codes,
aliases, and aggregate cohorts, resolves gene filters by ENSG or symbol, and can
return one or more normalization modes in one call:

- `normalize="tpm_clean"` / `"clean_tpm"` — shipped biological clean-TPM
  percentiles.
- `normalize="tpm_clean_biological"` — explicit name for that biological-only
  reference artifact.
- `normalize="tpm_clean_log1p"` — stored log1p biological clean-TPM percentiles.
- `normalize="tpm_raw"` / `"tpm"` — source-matrix raw TPM summaries recomputed
  through `cohort_stats`.

Long output includes source/provenance columns by default, including source
cohort, source type/unit, source scale class, reference method, `DATA_VERSION`,
and `SOURCE_MATRIX_VERSION`. This accessor is the compatibility surface for
reference-expression reads; expression artifact row-set/value parity is still
tracked separately in the upstream parity issues.
Raw source-matrix summaries default to `sample_qc="pass"` so sparse or
source-QC-failed samples do not shape newly computed reference rows. Use
`sample_qc="pass_or_warn"` or `"all"` for forensic audits or exact parity with
the current unfiltered source matrices. Long-form provenance and availability
rows report the live QC mode for raw source-matrix summaries and `artifact` for
existing shard-backed clean-TPM rows.

Use `expression.cancer_reference_expression_availability()` before delegating a
downstream reference-expression accessor that must distinguish unavailable
oncoref artifacts from empty gene filters. It returns one row per requested
code/mode with `requested_code`, expanded `cancer_code`, `request_kind`,
`available`, `missing_reason`, provenance fields, and the reference-expression
schema/data versions. `expression.cancer_reference_expression(...,
on_missing="empty")` returns a schema-stable empty frame and stores the same
missing rows in `df.attrs["missing_requests"]`; `on_missing="raise"` fails fast
for required cohorts. `include_request_metadata=True` adds request/availability
columns to long expression output, which is useful when a requested aggregate
expands to child expression cohorts.

Representative and percentile artifact readers have explicit downstream-facing
contracts:

- `expression.representative_cohort_samples(..., format="long",
  include_provenance=True)` includes the representative id, source cohort/project,
  source sample id, cohort sample count, deterministic selection rank/method/basis,
  artifact schema version, `DATA_VERSION`, and `SOURCE_MATRIX_VERSION`.
  Public representative ids default to pirlygenes-compatible `CODE_rep01`
  columns/values. Pass `representative_id_style="internal"` to expose the
  underlying bundle/provenance ids (`CODE__rep1`).
  Representatives are selected by central-medoid plus farthest-first traversal in
  log1p biological clean-TPM space, with stable sample-id tie-breaking; the
  persisted vectors remain full clean_tpm_16_9_75.
- `expression.cohort_gene_percentiles(..., include_provenance=True)` appends the
  cohort code, normalization, expression unit, percentile basis, artifact schema
  version, `DATA_VERSION`, and `SOURCE_MATRIX_VERSION`.
- Gene-level representative and percentile readers default to canonical oncoref
  ENSG IDs. For pirlygenes migration wrappers, pass
  `gene_id_style="pirlygenes"` to present known one-to-one
  `remapped_to_oncoref` rows with their legacy pirlygenes ENSG IDs. This is
  intentionally a presentation shim: it does not synthesize missing rows or alter
  expression values.
- Gene-level reference, representative, and percentile readers default to
  `gene_universe="artifact"`, which preserves the exact shipped row set. Pass
  `gene_universe="tumor_signal"` to drop only rows explicitly audited as
  oncoref-only technical extras for the requested artifact/cohort. Pass
  `include_gene_universe_flags=True` for long reference output or any
  representative/percentile output to append row-level `artifact_row_class`,
  `is_technical_extra`, `is_missing_biological`, and
  `recommended_consumer_action` columns. These options filter or label known
  artifact row classes; they never invent missing biological expression rows.
- Gene-level reference, representative, and percentile readers attach
  `df.attrs["gene_universe_delta_summary"]` and
  `df.attrs["gene_universe_delta_n"]` for the requested cohort/product. These
  attrs summarize the known pirlygenes/oncoref row-universe deltas that still
  apply to the returned artifact, so migration wrappers can separate remapped
  rows, missing upstream data, and intentional oncoref-only rows without
  reimplementing the audit-table matching logic.
- Missing percentile shards still raise by default. Use
  `on_missing="empty"` to return an empty but schema-stable frame with
  `df.attrs["missing_reason"]`, which is useful for compatibility adapters that
  need to distinguish unavailable upstream data from private downstream fallback
  data.

The QC-policy expression bundle ships representative, percentile,
within-sample, CTA-scope proteoform percentile, CTA-scope proteoform
within-sample, sample-QC, and build-metadata artifacts. Non-shipped proteoform
scopes can still recompute from cached source matrices. Row-set and value parity
with pirlygenes is still governed by the gene-universe and expression-artifact
parity issues.

`expression.expression_artifact_gene_universe_deltas()` exposes the known
pirlygenes/oncoref row-universe deltas from the current parity audit: canonical
remaps such as legacy `PAXX` to its oncoref ENSG, representative-sample rows
missing from oncoref source/artifact rows, and the full current set of oncoref-only representative
technical/noncoding, immune-receptor, Y-linked, and unresolved extra rows.
Use `expression.expression_artifact_gene_universe_delta_summary()` for counts by
product/cohort/status, or
`expression.expression_artifact_gene_universe_delta_report(product, cancer_types)`
for the compact request-scoped report used by accessor attrs. These tables include
`gene_biotype`, `artifact_row_class`, `is_technical_extra`,
`is_missing_biological`, and `recommended_consumer_action` so current-bundle row
classes do not have to be inferred from prose. Use
`expression.expression_artifact_technical_extra_gene_ids(...)` to get the
oncoref-only technical-extra ENSG IDs for a product/cohort filter. This surface is
intentionally provenance: it makes differences explicit for migration code, but
does not synthesize missing expression rows or alter artifact values.

Clean TPM has one public compartment contract:

- `clean-tpm-censored-genes.csv:category == "ribosomal_protein"` — 16%
  ribosomal compartment.
- `clean-tpm-censored-genes.csv:category == "technical"` — 9% other-technical
  compartment.
- genes absent from the censored table — 75% biological compartment.

The category-specific helper sets are available from `oncoref.gene_families`:

```python
from oncoref import gene_families

gene_families.clean_tpm_ribosomal_gene_ids()
gene_families.clean_tpm_other_technical_gene_ids()
gene_families.clean_tpm_censored_gene_ids()
```

For clean-TPM housekeeping denominators, use the biological HPA-stable panel:

```python
gene_families.clean_tpm_biological_housekeeping_gene_ids()
gene_families.clean_tpm_biological_housekeeping_genes()
gene_families.clean_tpm_biological_housekeeping_genes(primary_only=False)
```

Housekeeping normalization is defined as a median-of-ratios size factor against a
fixed, versioned per-gene reference profile:

```python
from oncoref import normalization

normalization.housekeeping_reference_profile()
normalization.tpm_to_housekeeping_normalized(matrix)
```

For each sample, oncoref computes:

```text
size_factor = median(housekeeping_clean_tpm[g] / reference_tpm[g])
normalized_expression[g] = clean_tpm[g] / size_factor
```

The default reference is the HPA v23-derived clean-TPM biological housekeeping
panel (`HOUSEKEEPING_REFERENCE_PROFILE_VERSION`). This is a sample-scale estimate
relative to a fixed biological HK profile, not the old "divide by the panel's
geometric mean" ratio. Prefer log1p clean TPM or percentile-rank clean TPM unless
the analysis specifically needs an HK-derived size factor.

The old geNorm-style denominator is deliberately buried behind
`method="legacy_geomean"` for explicit audits of historical outputs. The shorter
`method="geomean"` spelling is not accepted.

The legacy qPCR/reference-gene panel remains available as
`legacy_qpcr_housekeeping_*` and through the historical `housekeeping_*` helpers,
but it is not the clean-TPM biological denominator.

## Genes and Proteoforms

- `oncoref.gene_ids` — bundled canonical Ensembl gene space, cross-release alias
  resolution, symbol synonyms, and biotype checks.
- `oncoref.genome` — pyensembl-backed gene/transcript lookup against installed
  Ensembl releases.
- `oncoref.proteoforms` — identical-protein paralog grouping and expression
  collapse helpers.
- `oncoref.gene_qc` / `oncoref.gene_families` — technical-RNA and gene-family
  classification used by normalization. These are normalization/QC reference
  families, not the general home for pirlygenes marker panels.

## Burden, TMB, Fusions, and Signatures

- `oncoref.tmb` — tumor mutational burden reference values.
  `tmb.cancer_tmb_df()` includes evidence-schema columns (`estimate_type`,
  `source_scope`, `missing_reason`), and `tmb.cancer_tmb_record()` /
  `tmb.resolve_tmb_source()` preserve requested-code metadata for source-scoped
  lookups such as `COAD_MSI` or `READ_MSI` resolving through `CRC_MSI`. Direct
  audited gaps use `inheritance_kind="direct_missing"` so callers can distinguish
  “known no supported site-specific estimate” from an unmapped cancer code.
- `oncoref.incidence` — incidence/mortality burden and burden categories.
- `oncoref.fusions` — defining fusions and partner-family lookups.
- `oncoref.response_signatures` — legacy/compatibility response-signature
  surface used by oncoref plots. Treat it as transitional: new or extended
  therapy-response signature panels belong in pirlygenes unless they are recast
  as source-anchored empirical fact/provenance rows.

## Data Management

- `oncoref.catalog` — unified dataset inventory and fetch/status/path operations.
- `oncoref.data_bundle` — heavy expression bundle cache. Use
  `data_bundle.bundle_contract()` to inspect the downstream-stable package/data
  version linkage, release asset URLs, cache environment variables, completion
  marker policy, and expected artifact inventory for the active bundle. The
  inventory includes the generated sample-QC manifest, per-cohort build metadata,
  within-sample prevalence shards, and CTA-scope proteoform percentile/prevalence
  shards, not just the legacy pirlygenes expression tables. Use
  `data_bundle.bundle_release_manifest()` to fetch and validate only the small
  release manifest/checksum for the active `DATA_VERSION`, including tarball
  sha256 plus any artifact inventory, builder commit, source-matrix version, and
  sample-QC policy metadata published with the release. Use
  `data_bundle.bundle_metadata()` when a downstream package needs one
  no-heavy-download JSON object containing the static contract, local cache path
  and completeness state, local artifact inventory, and validated release
  manifest.
  CLI equivalents are available for CI/notebooks: `oncoref data contract` prints
  the static bundle contract, `oncoref data metadata [oncoref|pirlygenes]`
  prints the composed dependency state, and `oncoref data release-manifest
  [oncoref|pirlygenes]` prints only the validated release manifest/checksum
  metadata.
- `oncoref.reference_data` / `oncoref.hpa` — HPA reference-data cache and HPA
  tissue/cell-type accessors.

## Compatibility Modules

These modules remain importable but are less discoverable than the organized
facades above:

- `oncoref.apd1` — legacy anti-PD-1 response slice; prefer
  `oncoref.ici_response`.
- `oncoref.ici` — core ICI implementation; prefer `oncoref.ici_response` for the
  organized public surface.
- `oncoref.coverage` — original mixed CTA/generic antigen-panel coverage module;
  prefer `oncoref.cta_coverage` or `oncoref.antigen_coverage`.
- `oncoref.peptides` — original CTA-specific 9-mer module; prefer
  `oncoref.cta_peptides`.
