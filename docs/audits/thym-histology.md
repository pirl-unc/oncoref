# TCGA-THYM WHO histology partition (`THYM` vs `THYMCA`) and `GSE125285` linear-TPM upgrade

Issue [#737](https://github.com/pirl-unc/oncoref/issues/737), reviewed 2026-11-22.

## 1. Disjoint `THYM` (Thymoma) and `THYMCA` (Thymic Carcinoma) cohorts

`TCGA-THYM` profiles 119 primary thymic epithelial tumors spanning WHO Type A/AB/B1/B2/B3
thymomas and WHO Type C (`ICD-O-3 8586/3`) thymic carcinomas. Previously, all 119
columns were bundled into `THYM` (`TREEHOUSE_POLYA_25_01_TCGA_SAMPLES`), leaving
`THYMCA` without a direct expression cohort (`parent_fallback_no_direct_cohort`)
and leaving `THYM_EPITHELIAL` with `partial_member_gap`.

`scripts/build_tcga_thym_histology.py` partitions the checksum-pinned released
`source-v5.22.10/THYM_per_sample_tpm.parquet` matrix (`34,571` genes × `119`
columns) without altering expression values or gene coverage, joining against
pinned `cBioPortal/datahub` clinical metadata at revision
`0cc9138746c08b304f8dac92c31983e0ef44af1d`:

- **`THYM` (`n = 108` primary tumors):** WHO Types A (`15`), AB (`35`), A|AB (`2`),
  B1 (`14`), B1|B2 (`1`), B2 (`27`), B2|B3 (`3`), and B3 (`11`), concordant with
  ICD-O-3 `8581`–`8585` and `8521/1`.
- **`THYMCA` (`n = 11` primary tumors):** WHO Type C (`Thymoma, Type C` /
  `Thymoma; Type C`) and ICD-O-3 `8586/3` (`TCGA-3S-A8YW-01`, `TCGA-5G-A9ZZ-01`,
  `TCGA-5U-AB0D-01`, `TCGA-XM-A8R8-01`, `TCGA-XU-A92Z-01`, `TCGA-XU-A933-01`,
  `TCGA-XU-A936-01`, `TCGA-ZB-A966-01`, `TCGA-ZB-A969-01`, `TCGA-ZB-A96V-01`,
  `TCGA-ZC-AAA7-01`), with zero cross-table disagreements.
- **`THYM_EPITHELIAL` (`n = 119` primary tumors):** Computed union of the two
  disjoint member cohorts (`THYM` + `THYMCA`), upgrading umbrella coverage to
  `complete_linear_tpm`.

See [thym-histology-membership.csv](thym-histology-membership.csv) and
[thym-histology-build.json](thym-histology-build.json) for full per-sample audit
records and input SHA-256 digests.

## 2. `GSE125285` (`BCC` & `cSCC`) upgrade from `log2(CPM+1)` proxy to linear RNA-seq TPM

`GSE125285` (`SRP180015`, 25 primary `BCC` and 10 primary `cSCC` tumors paired with
35 matched adjacent normal skin biopsies) previously used the author-supplied
`GSE125285_S35_log2GE.txt.gz` (`log2(CPM+1)` across `15,263` symbol-wrapped genes),
which lacked gene-length normalization (`bulk_rnaseq_cpm_proxy`,
`linear_tpm_comparable=False`).

`oncoref.expression_source_adapters.gse125285_skin_matrix` now ingests NCBI GEO's
checksum-pinned `GSE125285_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz`
(`sha256: 0a12577632d939426c6dd1d0d045b301c184c7351a834ac2db8c2a637b41d82f`,
exon-union gene-length-normalized TPM from raw counts), verifies 1-to-1 sample
metadata concordance against `GSE125285_family.soft.gz`, harmonizes NCBI `GeneID`s
to `29,936` canonical Ensembl 112 genes, and renormalizes per tumor to `1,000,000`
TPM (`linear_rnaseq_tpm`, `linear_tpm_comparable=True`, 35/35 `sample_expression_qc_v3`
`pass`). All `bulk_rnaseq_cpm_proxy` cohorts are eliminated from
`non_comparable_expression_cohorts()`.

## 3. Public RNA-seq audit of remaining non-comparable cohorts

`oncoref.non_comparable_expression_cohorts()` now exposes a machine-readable
`public_rnaseq_audit_status` column recording why the remaining 8 rank-only proxy
cohorts (`5` microarray TPM-proxies, `3` malignant-cell scRNA pseudobulks) and 1
withheld invalid cohort (`NUTM`) retain their explicit scale disclosures:

- `MTC`: primary human medullary thyroid RNA-seq (`EGAS00001001473`) is
  controlled-access EGA; public GEO RNA-seq (`GSE114069`/`GSE114070`) profiles
  `TT`/`MZ-CRC-1` cell lines rather than primary tumors.
- `SARC_MYXLPS`: public liposarcoma RNA-seq (`GSE352573`) uses Lexogen QuantSeq 3'
  tag sequencing without histological subtype labels (`DDLPS`/`WDLPS`/`MYXLPS`).
- `SARC_ESS_LG` / `SARC_ESS_HG`: public GEO RNA-seq has `0` HG-ESS cohorts and only
  `2` LG-ESS patients (`GSE196757`, technical replicates) versus `n = 9` LG-ESS and
  `n = 4` HG-ESS primary tumors in `GSE85383`.
- `CMN`: `0` public multi-patient bulk RNA-seq cohorts on GEO (`TREEHOUSE_RIBOD_25_01`
  has `n = 1` unselected sample).
- `CTCL`, `EPN`, `HCL`: purified malignant-cell scRNA pseudobulks prevent reactive
  T-cell / normal-parenchyma dilution and remain explicitly isolated from bulk TPM
  thresholds.
