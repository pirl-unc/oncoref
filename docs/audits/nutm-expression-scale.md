# NUTM expression-scale incident

The published `source-v5.22.10/NUTM_per_sample_tpm.parquet` is not valid linear
TPM. Its 20,061 gene rows and three samples contain 2,163 negative measurements
(minimum -52.3163). The SHA-256 is
`ce4dfab53987f0b9492d74401b149bb240fca222bdec09297d188f5d06445ea7`.
A fresh release download matches the historical source cache exactly.

## Origin and consequences

The historical importer is
[`pirlygenes/scripts/build_nutm_tempus.py` at 6074ce3](https://github.com/pirl-unc/pirlygenes/blob/6074ce368078d22cbfbe06418f9dd206ad3840ec/scripts/build_nutm_tempus.py).
It reads `gene_tpm_cognizant_corrector`, then divides by the column sum and
multiplies by one million. Repeating those operations on the original input
reproduces the release matrix. The original field ranges from -3.6500 to
13.7014; the supplied dictionary calls it "Tempus external normalization
method v2" without defining the inverse transform. Rescaling these numbers
does not establish that they are linear TPM.

This corrupts absolute TPM thresholds and sums of paralog expression. In the
reported CTA analysis, HSFY1/2 meets p95 after its two transformed values are
summed. As a sensitivity check only, both common candidate inverses (`2**x`
and `max(2**x - 1, 0)`) remove that p95 hit. Neither inverse has been adopted
as corrected source data.

NUTM1 is present at the correct gene ID (`ENSG00000184507`). Its observed
ranks in the original supplied vectors are 69.11%, 86.41% and 80.44%. Those
numbers describe this three-sample, approximately 20,000-gene source; they
are not a validated biological conclusion. A monotone transform preserves
single-gene rank, but sample provenance, the measured gene universe, and
corrected versus uncorrected expression still require consideration. Missing
genes must not be padded with zeros to obtain a desired percentile.

The ranking populations are materially different across the plotted sources:

| Source matrix | Gene rows | Rows annotated protein-coding in current canonical gene space |
| --- | ---: | ---: |
| NUTM | 20,061 | 18,912 |
| BRCA_TNBC / GBM (Treehouse) | 34,571 | 19,101 |
| SARC_CHOR | 62,520 | 20,026 |

Consequently, a p90 label does not describe the same gene population across
these cohorts. This is a separate comparison limitation: declare and audit a
common ranking universe before interpreting cross-source percentile coverage.
The source records have three distinct RNA samples/analyses and no duplicated
sample/gene input rows. That rules out those mechanical duplication errors,
but does not establish complete independent diagnostic provenance for every
case or validate the expression scale.

## Related defects found during the audit

- Normalized source reads and release rebuilds silently clipped negative
  values to zero. This hid the scale defect without fixing it. Validation now
  occurs before canonicalization, where duplicate-row sums could also hide it.
- The generic declared-unit converter clipped negatives and converted NaN to
  zero before unit conversion. It now rejects incompatible values and
  preserves missing measurements, including completely unmeasured columns.
  The GDC STAR TPM reader also silently clipped negative and malformed values;
  it now rejects them before aggregation. Missing gene/sample values remain
  missing through GDC source assembly instead of becoming measured zeros.
- Count-to-TPM conversion silently assigned zero expression to rows lacking
  a valid gene length. It now requires finite positive lengths for every row.
- Source staging checked sample counts but not numerical validity, including
  when reusing a previous cache. It now validates all selections before
  copying any matrix into the cache or release directory.

Readers, builders and the proteoform-refresh workflow share the public
`oncoref.expression_engine.validate_nonnegative_expression(values, context=...)`
validator. It raises `ValueError` for negative, infinite or nonnumeric measurements
without modifying the input; missing values remain missing. Successful refresh
audits retain `n_negative_values_clipped=0` for compatibility with historical audit
files. No values are repaired by clipping. Raw and normalized reads both reject
the invalid NUTM source, including attempts to recompute stale protein summaries.

GDC duplicate-gene aggregation uses native `groupby.sum(min_count=1)` for TPM and
a separate first-symbol reduction. Entirely unmeasured groups remain missing,
partial observations sum the measured values, and measured zero remains zero.
This avoids a Python callback for every gene in every sample.

Minimum-value statistics from all 89 locally cached human source matrices
identified NUTM as the only negative-valued matrix. Every numerical column
had parquet minimum statistics. This does not validate uncached sources or
exclude all-positive transformed inputs, sample-label mistakes or other
provenance defects. It also does not retroactively repair derived artifacts.

## Repair status

The source registry and availability metadata flag NUTM as an invalid legacy
scale, with `linear_tpm_comparable=False`. The public per-sample reader rejects
the invalid artifact in every normalization mode.

### Derived artifacts are withheld (1.27.1)

Before 1.27.1 the per-sample reader rejected the source, but artifacts derived from
it were still served: `cancer_reference_expression("NUTM")` and
`cohort_gene_percentiles("NUTM")` returned the invalid summaries without warning.
Every reader now consults the availability manifest's `source_scale_class` for the
cohort's *selected* source, so the guard is driven by registry data and applies to any
source later marked `invalid_legacy_tpm`, not only NUTM:

| Reader | Named invalid cohort | Default or bulk listing |
|---|---|---|
| `cohort_gene_percentiles`, `gene_within_sample_top_fraction`, `proteoform_within_sample_top_fraction` | `ValueError` | cohort omitted from `available_*_cohorts()` |
| `gene_representative_samples`, `proteoform_representative_samples` | `ValueError` | cohort omitted from `available_representative_cohorts()` |
| `cancer_reference_expression` (artifact and `summary_rows`) | `ValueError`, whatever `on_missing` says | cohort omitted with a `UserWarning` |
| `cancer_reference_expression(reference_source="summary_rows_all")` | `ValueError` only if the invalid source is named in `source_cohort` | invalid source omitted with a `UserWarning`; a valid source such as Treehouse is returned |
| `cancer_reference_expression_availability` | describes the source (`available=False`) | still lists it, so the problem stays visible |

A named invalid cohort is never silently dropped as "unavailable", including in a
multi-cohort request. Valid cohorts are unaffected. Regression tests link the
registry's source scale to every derived-artifact reader
(`tests/test_invalid_source_derived_artifacts.py`).

### Still open

This does not repair the source. Issue #604 stays open until:

1. The exact source-field transform is established, or verified linear TPM and the
   quantification reference are obtained. Do not infer an inverse solely from the
   field name or from current publications using a similarly named output.
2. The selected source is rebuilt with explicit scale and sample provenance and
   published under a new per-cohort source-matrix version.
3. Every derived NUTM reference artifact is rebuilt, figures and fitted comparisons
   are regenerated, and NUTM1 ranks are reassessed on that verified input.
4. Cross-source p90/p95 comparisons declare a common ranking universe.

The controlled three-case input is not publicly downloadable and is not included in
this repository, so none of these steps can be done here. Until then NUTM stays
withheld. Consumers that read raw artifact files directly (for example the PirlyGenes
compatibility reader) bypass these guards and must apply the same rule.
