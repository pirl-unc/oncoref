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

## Related defects found during the audit

- Normalized source reads and release rebuilds silently clipped negative
  values to zero. This hid the scale defect without fixing it. Validation now
  occurs before canonicalization, where duplicate-row sums could also hide it.
- The generic declared-unit converter clipped negatives and converted NaN to
  zero before unit conversion. It now rejects incompatible values and
  preserves missing measurements, including completely unmeasured columns.
  The GDC STAR TPM reader also silently clipped negative and malformed values;
  it now rejects them before aggregation.
- Count-to-TPM conversion silently assigned zero expression to rows lacking
  a valid gene length. It now requires finite positive lengths for every row.
- Source staging checked sample counts but not numerical validity, including
  when reusing a previous cache. It now validates all selections before
  copying any matrix into the cache or release directory.

Minimum-value statistics from all 89 locally cached human source matrices
identified NUTM as the only negative-valued matrix. Every numerical column
had parquet minimum statistics. This does not validate uncached sources or
exclude all-positive transformed inputs, sample-label mistakes or other
provenance defects. It also does not retroactively repair derived artifacts.

## Repair status

The source registry and availability metadata flag NUTM as an invalid legacy
scale, with `linear_tpm_comparable=False`. The public per-sample reader rejects
the invalid artifact in every normalization mode. Historical derived NUTM
artifacts remain invalid and must not be used for quantitative analysis.

Remaining work before this incident is resolved:

1. Establish the exact source-field transform, or obtain verified linear TPM
   and the quantification reference. Do not infer an inverse solely from the
   field name or from current publications using a similarly named output.
2. Rebuild the selected source with explicit scale and sample provenance;
   publish it under a new per-cohort source-matrix version.
3. Rebuild every derived NUTM reference artifact and regenerate affected
   figures and fitted comparisons. Reassess NUTM1 ranks on that verified input.

Private case-level input and identifiers are not included in this audit.
