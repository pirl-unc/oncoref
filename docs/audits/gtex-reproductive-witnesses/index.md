# Individual GTEx v11 tissue samples

The catalog uses 18,637 primary tissue samples from 946 GTEx donors. Each of the 52 tissue files has a versioned URL and SHA-256 checksum in the reference registry. Individual reproductive samples can qualify RA. Individual normal and Disposable samples contribute to the shared 10× denominator for CA and RA. The normal comparison retains the highest complete within-sample protein-group sum per GTEx tissue, with the sample and donor that attained it. The full matrices retain every measurement.

[GTEx and HPA contamination screen and donor sensitivity](../gtex-cross-tissue-signals/index.md) · [Source registry and rebuild instructions](../../curation.md) · [Sample inventory and QC](../target-selection-20261110/rebuild/gtex-library-qc.csv.gz) · [Comparison methods](../target-selection-20261110/rebuild/comparison.json)

## Symmetric comparison

All ten pilot proteins clear reproductive 10×. Eight qualify as CTAs. C1orf202 and TMEM95 qualify as RA-only because tumor expression peaks at 7.34× and 3.23×, respectively.

[Results for all ten, including donor counts and both fold values](symmetric-comparison.csv)

[Before and after for all 918 assessed groups, including every classification change](catalog-comparison.csv.gz) · [Previous catalog revision and checksum](baseline.json)

Across the inventory, 179 of the previous 303 CTAs retain membership and 24 former CA-only proteins gain it. Of the 124 leaving CTA, 98 become RA-only, 13 CA-only and 13 Other. One sample labeled esophagus muscularis (GTEX-1I1GR-1426-SM-CY8H8) sets the maximum normal TPM for 77 of the 124. It measures PRM2 at 3,047.26 TPM and carries the advisory R code. The denominator retains the sample under its provider tissue label. Its cellular origin is unresolved.

## Original testis-only pilot

All ten CA-only leads identified in issue #716 have individual GTEx v11 testis samples at or above their pre-integration catalog threshold. Each of the 414 sample IDs in this matrix has a distinct GTEx participant prefix. Qualifying participant counts range from six for GOLGA6L1 to 376 for CTAGE1. The measured median for each of the eight leads originating from the v11 median file reproduces its reported value to the precision of that file. The two leads from the older HPA-hosted GTEx export also have qualifying v11 samples, although their v11 medians are below threshold.

[Protein results and highest samples](testis-summary.csv) · [All 4,140 individual measurements](testis-samples.csv.gz) · [Original 78-protein comparison](../antigen-evidence-review-20261118/ca-only-78.csv)

The original pilot applies the previous threshold of 10 × max(N, D, 1 TPM) to the provider's versioned Ensembl rows. Quantification of homologous loci retains the limitations of the original RNA-seq method.

The source is the [GTEx v11 testis TPM matrix](https://storage.googleapis.com/adult-gtex/bulk-gex/v11/rna-seq/tpms-by-tissue/gene_tpm_adult_gtex_v11_testis.gct.gz) from the [GTEx public downloads](https://gtexportal.org/home/downloads/adult-gtex/). The 84,558,723-byte gzip matches the provider MD5 of `937ce9a599ba69269a404cc7dea2351e`. The script verifies SHA-256 before reading, matches the exact source gene version in the catalog identity decisions and retains sample and participant IDs. Run from the repository root with the downloaded file:

```shell
PYTHONPATH=. python scripts/audit_gtex_testis_witnesses.py gene_tpm_adult_gtex_v11_testis.gct.gz
```
