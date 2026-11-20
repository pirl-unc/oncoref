# Individual GTEx v11 tissue samples

The catalog uses 18,637 primary tissue samples from 946 GTEx donors. Each of the 52 tissue files has a versioned URL and SHA-256 checksum in the reference registry. Individual reproductive samples can qualify RA. GTEx non-reproductive tissue sites contribute interpolated donor P95 to the shared CA and RA denominator. Complete member-locus TPMs are summed within each sample first. The full matrices retain every measurement.

[GTEx and HPA contamination screen and donor sensitivity](../gtex-cross-tissue-signals/index.md) · [Source registry and rebuild instructions](../../curation.md) · [Sample inventory and QC](../target-selection-20261110/rebuild/gtex-library-qc.csv.gz) · [Comparison methods](../target-selection-20261110/rebuild/comparison.json)

## Symmetric comparison

The live comparison reports how many of the 414 testis donors pass each protein’s current threshold and how many tumor-fold units its best tumor reaches. The full table compares the pre-GTEx catalog with the current P95 policy.

[Ten pilot proteins and their donor counts](symmetric-comparison.csv) · [All 918 fold comparisons and category transitions](catalog-comparison.csv.gz) · [Baseline catalog checksum](baseline.json)

The earlier individual-maximum policy had 203 CTAs and moved 98 former CTAs to RA-only. Its prominent esophagus muscularis measurement and the alternative percentile calculations are documented in the [historical QC audit](../gtex-cross-tissue-signals/index.md). The [current source report](../target-selection-20261110/rebuild/normal-rna/index.md) gives the rebuilt totals.

## Original testis-only pilot

All ten CA-only leads identified in issue #716 have individual GTEx v11 testis samples at or above their pre-integration catalog threshold. Each of the 414 sample IDs in this matrix has a distinct GTEx participant prefix. Qualifying participant counts range from six for GOLGA6L1 to 376 for CTAGE1. The measured median for each of the eight leads originating from the v11 median file reproduces its reported value to the precision of that file. The two leads from the older HPA-hosted GTEx export also have qualifying v11 samples, although their v11 medians are below threshold.

[Protein results and highest samples](testis-summary.csv) · [All 4,140 individual measurements](testis-samples.csv.gz) · [Original 78-protein comparison](../antigen-evidence-review-20261118/ca-only-78.csv)

The original pilot applies the previous threshold of 10 × max(N, D, 1 TPM) to the provider's versioned Ensembl rows. Quantification of homologous loci retains the limitations of the original RNA-seq method.

The source is the [GTEx v11 testis TPM matrix](https://storage.googleapis.com/adult-gtex/bulk-gex/v11/rna-seq/tpms-by-tissue/gene_tpm_adult_gtex_v11_testis.gct.gz) from the [GTEx public downloads](https://gtexportal.org/home/downloads/adult-gtex/). The 84,558,723-byte gzip matches the provider MD5 of `937ce9a599ba69269a404cc7dea2351e`. The script verifies SHA-256 before reading, matches the exact source gene version in the catalog identity decisions and retains sample and participant IDs. Run from the repository root with the downloaded file:

```shell
PYTHONPATH=. python scripts/audit_gtex_testis_witnesses.py gene_tpm_adult_gtex_v11_testis.gct.gz
```
