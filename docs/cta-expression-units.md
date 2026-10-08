# RNA expression units in CTA curation

The CTA normal-tissue evidence uses **HPA v23 consensus nTPM**. Keep that unit
when interpreting the tables, plots and thresholds; a tumor dataset reporting
TPM is not automatically on the same scale.

| Measure | What is normalized | Interpretation |
|---|---|---|
| Counts | No abundance normalization | Assigned reads/fragments; affected by sequencing depth and transcript length. |
| CPM | Counts divided by library size, multiplied by one million | Adjusts sequencing depth; does not adjust transcript length. |
| TPM | Length-adjusted abundance rescaled to sum to one million per sample | Relative RNA abundance over the included transcripts. |
| RPKM / FPKM | Reads / fragments divided by feature length and library size | Adjusts length and depth; unlike TPM, the column total is not fixed. |
| HPA nTPM | Protein-coding TPM followed by TMM normalization | HPA's normalized RNA expression scale, with a defined reference set and processing pipeline. |

See the [Salmon output definitions](https://salmon.readthedocs.io/en/latest/file_formats.html)
for TPM versus assigned reads, and the
[TPM methods paper](https://doi.org/10.1007/s12064-012-0162-3) for the distinction
from RPKM.

## Why HPA adds the “n”

HPA rescales TPM after removing noncoding transcripts, producing pTPM, then
applies TMM (trimmed mean of M values) normalization. Its normal-tissue consensus
uses the maximum of the HPA and GTEx estimates for each gene and tissue. The
result is an aggregated reference, not a measurement from one individual.
These definitions are pinned to the
[HPA v23 methods](https://v23.proteinatlas.org/about/assays%2Bannotation).

The extra scaling addresses RNA composition: a few highly abundant transcripts
can consume much of a sample's fixed total, making other transcripts look
relatively lower. TMM estimates a scaling factor from a trimmed set of gene
ratios to reduce this effect. It does not make RNA measurements absolute or
remove every biological or technical difference. See the
[original TMM paper](https://doi.org/10.1186/gb-2010-11-3-r25).

Consequently, do not substitute 10 TPM from an unrelated pipeline for 10 HPA
nTPM. There is no universal scalar conversion for an isolated gene value:
the included genes, normalization factors and aggregation matter. Neither
measure reports protein abundance or HLA peptide presentation.

## What the CTA fractions mean

`rna_*_ntpm` columns are expression measurements. The
`rna_deflated_reproductive_frac` column is a **dimensionless restriction score**:

```text
deflated signal[tissue] = max(nTPM[tissue] - 1, 0)
fraction = sum(deflated signal in selected tissues)
           / sum(deflated signal in all observed tissues)
```

`rna_reproductive_tissue_scope` identifies the numerator: `core` or `extended`.
Changing the scope does not change the expression unit. It changes which
tissues contribute to the numerator. For example, RNASE12 has a core fraction
of 0.0004 and an extended fraction of 0.9993 because its RNA is concentrated
in accessory reproductive tissues. Those two numbers are fractions, not nTPM.
See the [full panel audit](audits/cta-panel-curation.md).

When observed expression is entirely below the deflation threshold, the
deflated fraction follows the existing convention of 1.0; the separate
expression policy prevents that convention alone from admitting a gene to the
default panel. When a gene has **no RNA observations**, regeneration instead
clears its RNA measurements, assigns `NO_DATA` RNA axes and fails the RNA gate.
Missing values are not observed zeros, and old measurements are never reused
under a new scope label.

## Identical-protein expression and registry updates

The human CTA protein registry groups the longest Ensembl r112 protein sequences
from the current CTA source table and its explicitly seeded partners. This is
an identity reference, not a change to canonical CTA panel admission. The refresh
for [#567](https://github.com/pirl-unc/oncoref/issues/567) contains 37 groups and
93 member genes, including the identical 496-aa RBMY1F/RBMY1J pair. Both genes
now have the CTA-scope key `RBMY1F/J`.

Expression collapse sums independently allocated per-gene TPM contributions
within each sample **before** calculating percentiles or within-sample ranks.
Do not add previously computed prevalence fractions or percentiles, or duplicate
one ambiguous RNA observation across its possible source loci.

`scope="cta"` and `scope="genome"` select different collapse registries. Ranking
uses all remaining biological rows after the selected collapse; CTA scope does
not restrict the ranking denominator to CTA genes alone. Genome scope can merge
additional groups and therefore change ranks even for an unchanged CTA protein.

Generated proteoform summary parquets retain `proteoform_scope` and
`proteoform_registry_sha256` in their pandas attrs. Readers use a precomputed
summary only when both match the current registry. Legacy summaries without the
fingerprint, and summaries with an outdated fingerprint, are recomputed from the
cohort's per-sample matrix. If that matrix is unavailable, the accessor raises
an error explaining the rebuild requirement; the proteoform wrappers default to
`auto_fetch=True`. Local availability helpers exclude stale summaries unless
their source matrix is cached.

To regenerate the reference and summaries from biological clean-TPM matrices:

```bash
python scripts/generate_proteoform_groups.py --ensembl-release 112
python scripts/generate_cohort_percentiles.py --input <clean-matrix-dir> --proteoform --scope cta
python scripts/generate_within_sample_top5.py --input <clean-matrix-dir> --proteoform --scope cta
```

The summary inputs must use the intended sample-QC policy; drop technical genes
before building or supply the generators' `--drop-genes` file. The full
`rebuild_expression_artifacts.py` workflow also preserves the registry fingerprint.
Publishing regenerated bundles requires the usual separate data-bundle release;
runtime invalidation does not change the pinned bundle version.
