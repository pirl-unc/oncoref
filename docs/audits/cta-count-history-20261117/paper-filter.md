# The 97 CTAs excluded at the standard paper-count step

This is the `f7847d51` snapshot. The [subsequent literature audit](../antigen-evidence-review-20261117/index.md)
corrects a missing SSU72L1 paper join: the same 97 exclusions become 48 with zero
counted tumor-support papers and 49 with one. Forty-three of the original 49
already have broader literature records when follow-ups are joined correctly.

These proteins pass the catalog and standard mode's earlier normal-RNA and curated-home gates. They leave only at the publication-count step. All remain CTAs in the broad catalog. The list is bound to catalog policy `oncoref.target-selection.2026-11-17.1`; default anatomy does not narrow it.

`min_papers` is an optional prioritization filter: all requires 0, loose 1, standard 2 and strict 3 counted sources. A source counts once per protein, across all eligible roles. A CTA catalog paper counts only proteins it actually nominates; oncofetal-program and tumor-expression papers count the proteins their records link. Reproductive-only biology papers do not count toward this tumor-support filter. Reviewed normal-expression literature is a separate curation layer and does not require two papers to establish a concern.

This is a count of distinct curated publications, not independent cohorts, patients, experimental validations, peptide presentation or clinical safety. Screens may reuse datasets. Zero means no eligible publication recorded in this registry, not no paper exists. Publication support is useful for prioritizing familiar targets but can penalize less-studied, expression-supported candidates.

Of 290 broad CTAs, 106 fail the two-paper threshold independently. Nine already fail an earlier standard gate, leaving the **97 sequential paper exclusions** below: **49 zero-paper and 48 one-paper proteins**. 59 of these 97 already have at least two known patients meeting the tumor 10× threshold.

[Download all 97 with fold, patient and publication counts](standard-paper-exclusions.csv) · [Counted publication links](standard-paper-links.csv)

## 49 proteins with 0 counted papers

`ANTXRL`, `CER1`, `CIROP`, `CST11`, `DEFB112`, `DEFB129`, `DEFB135`, `DPRX`, `ENSG00000283205`, `ERVV-1`, `ERVV-2`, `ETDA/B`, `FAM90A27P`, `FOXL3`, `HMGB1P1`, `HNRNPCL3/4`, `HOXB1`, `IZUMO1R`, `KLF18`, `KRTAP20-4`, `LGALS16`, `MBD3L4`, `MPC1L`, `MSGN1`, `NANOGNB`, `OR10J1`, `PABPN1L`, `PNMA6E`, `PRR23D1/2`, `PRR23E`, `RNASE9`, `SCGB1D4`, `SMIM42`, `SPADH`, `SPANXN4`, `SPATA31D3`, `SPATA31F3`, `SSU72L1`, `SSU72L2`, `SSU72L3`, `SSU72L4`, `SSU72L5`, `SSU72L6`, `TCP11X2`, `TEX28`, `TMEM239`, `TPRX2`, `TSPY8`, `UBTFL1`.

## 48 proteins with 1 counted paper

| Protein group | Counted publication |
|---|---|
| ACCSL | [Xiang C et al. (2021)](https://pubmed.ncbi.nlm.nih.gov/34721544/) |
| ANHX | [Chang et al. (2019)](https://doi.org/10.1002/cam4.2223) |
| C16orf82 | [da Silva et al. (2017)](https://doi.org/10.18632/oncotarget.21715) |
| C19orf85 | [Smith et al., 2023](https://pmc.ncbi.nlm.nih.gov/articles/PMC10578318/) |
| C20orf144 | [da Silva et al. (2017)](https://doi.org/10.18632/oncotarget.21715) |
| CALR3 | [da Silva et al. (2017)](https://doi.org/10.18632/oncotarget.21715) |
| CCDC179 | [Chang et al. (2019)](https://doi.org/10.1002/cam4.2223) |
| CIMIP4 | [Chang et al. (2019)](https://doi.org/10.1002/cam4.2223) |
| CSAG1 | [Wang et al. (2016)](https://doi.org/10.1038/ncomms10499) |
| CT45A10 | [Loriot et al. (2025)](https://doi.org/10.1371/journal.pgen.1011734) |
| CXorf51A/B | [Loriot et al. (2025)](https://doi.org/10.1371/journal.pgen.1011734) |
| DDI1 | [da Silva et al. (2017)](https://doi.org/10.18632/oncotarget.21715) |
| EBLN1 | [Chang et al. (2019)](https://doi.org/10.1002/cam4.2223) |
| EDDM3B | [da Silva et al. (2017)](https://doi.org/10.18632/oncotarget.21715) |
| GOLGA6D | [Wang et al. (2016)](https://doi.org/10.1038/ncomms10499) |
| H2AB1 | [da Silva et al. (2017)](https://doi.org/10.18632/oncotarget.21715) |
| H2BC1 | [Chang et al. (2019)](https://doi.org/10.1002/cam4.2223) |
| HNRNPCL1 | [Chew et al., 2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6736738/) |
| HNRNPCL2 | [da Silva et al. (2017)](https://doi.org/10.18632/oncotarget.21715) |
| IQCF2 | [Chang et al. (2019)](https://doi.org/10.1002/cam4.2223) |
| IQCF6 | [Chang et al. (2019)](https://doi.org/10.1002/cam4.2223) |
| LALBA | [Wang et al. (2016)](https://doi.org/10.1038/ncomms10499) |
| LYPD4 | [da Silva et al. (2017)](https://doi.org/10.18632/oncotarget.21715) |
| MAGEB4 | [da Silva et al. (2017)](https://doi.org/10.18632/oncotarget.21715) |
| MBD3L2B | [Smith et al., 2023](https://pmc.ncbi.nlm.nih.gov/articles/PMC10578318/) |
| MEIOSIN | [Chang et al. (2019)](https://doi.org/10.1002/cam4.2223) |
| OR1N1 | [da Silva et al. (2017)](https://doi.org/10.18632/oncotarget.21715) |
| OR4M2 | [da Silva et al. (2017)](https://doi.org/10.18632/oncotarget.21715) |
| PAGE3 | [Chang et al. (2019)](https://doi.org/10.1002/cam4.2223) |
| PATE1 | [PATE, a gene expressed in prostate cancer, normal prostate, and testis, identified by a functional genomic approach](https://pubmed.ncbi.nlm.nih.gov/11880645/) |
| PPP1R2B | [Bruggeman et al. (2018)](https://doi.org/10.1038/s41388-018-0357-2) |
| PRAMEF15 | [Chew et al., 2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6736738/) |
| PRAMEF25/6 | [Chew et al., 2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6736738/) |
| RBMY1F/J | [Wang et al. (2016)](https://doi.org/10.1038/ncomms10499) |
| RNASE12 | [da Silva et al. (2017)](https://doi.org/10.18632/oncotarget.21715) |
| SEMG2 | [SEMG1/2 augment energy metabolism of tumor cells](https://pmc.ncbi.nlm.nih.gov/articles/PMC7733513/) |
| SPAG11B | [da Silva et al. (2017)](https://doi.org/10.18632/oncotarget.21715) |
| SPO11 | [Chang et al. (2019)](https://doi.org/10.1002/cam4.2223) |
| TBC1D21 | [Chang et al. (2019)](https://doi.org/10.1002/cam4.2223) |
| TMEM225 | [Chang et al. (2019)](https://doi.org/10.1002/cam4.2223) |
| TRIM43B | [Smith et al., 2023](https://pmc.ncbi.nlm.nih.gov/articles/PMC10578318/) |
| TXNDC8 | [Chang et al. (2019)](https://doi.org/10.1002/cam4.2223) |
| WFDC10A | [da Silva et al. (2017)](https://doi.org/10.18632/oncotarget.21715) |
| ZCCHC13 | [Chang et al. (2019)](https://doi.org/10.1002/cam4.2223) |
| ZNF722 | [Loriot et al. (2025)](https://doi.org/10.1371/journal.pgen.1011734) |
| ZNF723 | [Loriot et al. (2025)](https://doi.org/10.1371/journal.pgen.1011734) |
| ZNF729 | [Loriot et al. (2025)](https://doi.org/10.1371/journal.pgen.1011734) |
| ZP4 | [Somes LK et al. (2025)](https://pubmed.ncbi.nlm.nih.gov/39980195/) |
