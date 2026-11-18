# GTEx esophagus muscularis and the 98 former CTAs

Historical sensitivity analysis against the individual-maximum catalog. Reproduce these tables from commit `819111a9`. The [current source report](../target-selection-20261110/rebuild/normal-rna/index.md) uses interpolated donor P99 for GTEx normal tissue sites.

**GTEX-1I1GR-1426-SM-CY8H8 is a pathology-confirmed muscularis specimen with convincing evidence of testis RNA admixture in its sequencing library.** The GTEx pathologist recorded “6 pieces, confirmed target muscularis.” The slide record explicitly links the specimen to this RNA library from a male donor aged 60–69. The note records no tumor finding. It gives no microscopic disease-free certification of every piece or independent identity check on the sequenced RNA. [GTEx slide record](https://gtexportal.org/api/v2/histology/image?tissueSampleId=GTEX-1I1GR-1426) · [GTEx v11 specimen annotations](https://storage.googleapis.com/adult-gtex/annotations/v11/metadata-files/GTEx_Analysis_v11_Annotations_SampleAttributesDS.txt)

## Tissue and batch evidence

* Across 36,864 genes reaching 1 TPM in at least one primary tissue median, the sample most closely resembles esophagus muscularis (Spearman 0.872), followed by gastroesophageal junction (0.870). Its MYH11 expression is 5,639 TPM, versus a median of 8,199 in the other 560 muscularis samples and 140 in testis. [All 52 tissue correlations](tissue-correlations.csv) · [Marker TPM](marker-comparison.csv)
* PRM2 reaches 3,047 TPM. The next highest muscularis sample reaches 49.9 TPM and the peer median is 0.74 TPM. Excluding all six original warning markers, 702 of 717 additional testis-enriched loci exceed **every** other muscularis library. The median excess above the muscularis baseline is 10.9% of the corresponding testis median across these genes. This describes RNA levels. A tissue or cell fraction would require an independently calibrated mixture model. The donor’s esophageal mucosa and gastroesophageal junction samples both fall below the coordinated reproductive signature threshold. [Same-donor tissue comparison](same-donor-tissues.csv)
* Comparing the **pattern of deviations** across those 717 genes to 414 individual testis libraries identifies GTEX-1ICLY-2926-SM-C1YS5 as the closest match (Spearman 0.826). The runner-up reaches 0.449. The donor’s own testis reaches −0.056. The closest match comes from a different donor and shares the suspect library’s RNA extraction batch, BP-72061. The same sample ranks first in each disjoint half of the panel (0.804 and 0.846). This supports technical carryover. The exact source and processing step remain unconfirmed without read-level genotype concordance or an independent extraction. Cross-tissue carryover has also been documented in other GTEx libraries. [All testis comparisons and batch IDs](testis-profile-comparison.csv) · [Nieuwenhuis et al.](https://www.nature.com/articles/s41467-020-15821-9)

## Tumor evidence

All 98 proteins changed from CTA to RA-only when individual GTEx normal measurements entered the comparison. All have a measured tumor peak of at least 10 TPM, and 13 reach 100 TPM. At their earlier protein-specific tenfold thresholds, 35 had one qualifying tumor library, 40 had two to four, 11 had five to nine and 12 had at least ten. Forty-two have at least two identified patients at the earlier threshold. The same 98 now have zero libraries at the updated threshold. Bulk tumor expression leaves the producing cell type unresolved.

[All 98 proteins with peak tumor type, sample and TPM, qualifying libraries and known patients, old and new thresholds, and sensitivity to this specimen](tumor-evidence-98.csv)

This esophageal library sets the normal maximum for 77 of the 98. Omitting only that sample restores the tumor tenfold margin for **nine**: CCER1, DUSP21, FSCN3, H2AP, MPC1L, PRSS37, SPANXN4, SPATA31D1 and SPATA31D3. Among the 77, the next maximum comes from the same three libraries for 73 proteins: thyroid GTEX-1HSEH-0826-SM-B2LXE (47), tibial nerve GTEX-13NYS-0426-SM-5MR4X (16), and sigmoid colon GTEX-16GPK-2526-SM-GQ1DH (10). All three carry the advisory reproductive-program flag. Their provenance requires separate review. This sensitivity analysis leaves the catalog unchanged.

## Reproduction

The audit uses the catalog and protein-group identities at commit `47d59b38`, the verified tumor matrix and sample identity receipt, the pinned GTEx v11 tissue matrices, and the linked v11 specimen annotations. Run the following from a source checkout with the tumor cache and the GTEx v11 median matrix staged. The samples directory contains the registered `esophagus_muscularis.gct.gz` and `testis.gct.gz` files.

```sh
PYTHONPATH=. python scripts/audit_gtex_esophagus_library.py GTEx_Analysis_v11_Annotations_SampleAttributesDS.txt /path/to/samples-directory
```

The independent testis panel requires a testis median of at least 10 TPM and at least 20 times the greater of 1 TPM and the highest other primary tissue median. For the sample-to-sample comparison, each gene’s testis measurement is divided by its testis median. The suspect library is compared on that scale after subtracting its muscularis peer median and flooring at zero. Donor-specific expression similarity is exploratory and depends on the panel and normalization. The tumor sensitivity calculation recomputes each protein-group maximum among the other 560 muscularis samples after summing member loci within the same sample. It retains all remaining tissue summaries, individual observations, Disposable bounds, and the 1 TPM floor.
