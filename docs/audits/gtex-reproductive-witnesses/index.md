# Individual GTEx testis witnesses

All ten CA-only leads identified in issue #716 have individual GTEx v11 testis samples at or above their existing catalog threshold. Each of the 414 sample IDs in this matrix has a distinct GTEx participant prefix. Qualifying participant counts range from six for GOLGA6L1 to 376 for CTAGE1. The measured median for each of the eight leads originating from the v11 median file reproduces its reported value to the precision of that file. The two leads from the older HPA-hosted GTEx export also have qualifying v11 samples, although their v11 medians are below threshold.

[Protein results and highest samples](testis-summary.csv) · [All 4,140 individual measurements](testis-samples.csv.gz) · [Original 78-protein comparison](../antigen-evidence-review-20261118/ca-only-78.csv)

This audit applies the existing threshold of 10 × max(N, D, 1 TPM) to the provider's versioned Ensembl rows. Quantification of homologous loci retains the limitations of the original RNA-seq method. Catalog admission requires the corresponding non-reproductive GTEx sample comparison, tissue and sample QC, and reassessment across the inventory. Adding individual normal measurements can raise the threshold.

The source is the [GTEx v11 testis TPM matrix](https://storage.googleapis.com/adult-gtex/bulk-gex/v11/rna-seq/tpms-by-tissue/gene_tpm_adult_gtex_v11_testis.gct.gz) from the [GTEx public downloads](https://gtexportal.org/home/downloads/adult-gtex/). The 84,558,723-byte gzip matches the provider MD5 of `937ce9a599ba69269a404cc7dea2351e`. The script verifies SHA-256 before reading, matches the exact source gene version in the catalog identity decisions and retains sample and participant IDs. Run from the repository root with the downloaded file:

```shell
PYTHONPATH=. python scripts/audit_gtex_testis_witnesses.py gene_tpm_adult_gtex_v11_testis.gct.gz
```
