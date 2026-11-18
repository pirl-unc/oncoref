# Cross-tissue expression in individual GTEx v11 libraries

The survey ranked all 18,637 primary tissue libraries. At a relative expression threshold of 1%, 384 libraries showed an unusual foreign tissue program. Twenty-one had a reproducible best match to a different donor's source tissue in the same technical batch, supporting probable carryover. Six of the 21 included a testis signal. These are research flags. The measurements cannot certify contamination or estimate prevalence across all types of contamination.

[Ranked batch matches with sample and batch IDs](batch-matched.csv) · [Counts by threshold](thresholds.csv) · [All sample–panel flags](signals.csv.gz) · [Marker panel sizes](panels.csv)

## Scoring and support

The screen aggregates GTEx sites within their named organ families and compares each library with other organ families. A marker must reach at least 20 times the greater of 1 TPM and its highest median in any other organ family. A panel needs at least 10 such loci. Fifteen source panels met that requirement. For each marker, subtract the recipient site's median from sample TPM, floor at zero, and divide by the marker's source-family reference abundance (the highest native site median). The score is the median across the panel. At least 75% of its markers must also exceed their 95th percentile within the recipient site. A score of 0.10 denotes excess RNA at 10% of source-reference expression, with no calibration to mass, volume or percentage of cells.

For panels of at least 20 loci and scores of at least 1%, rank every native library by Spearman correlation across normalized marker levels with the recipient's excess. The more stringent classification requires a unique best library in the full panel and in both halves (alternating sorted Ensembl IDs), overall correlation ≥0.5, different donor IDs and overlap of a recorded RNA extraction or expression batch. The two halves can include co-regulated genes. Batch effects and shared biology can affect correlations. Read-level donor variants or an independent extraction would provide stronger attribution.

True tissue composition can yield large screening scores. For example, a minor salivary gland sample with a large muscle signal has a pathology note describing mostly connective tissue and muscle. Conversely, same-donor carryover, mixtures of several sources, low-level contamination, and donors or tissues lacking adequate reference data can evade the stringent screen. The six-marker testis warning has its own denominator of 18,223 non-testis libraries and uses the least expressed of six markers relative to its reference. It answers a different question from the broader marker panel.

## Cases relevant to the antigen catalog

* Esophagus muscularis **GTEX-1I1GR-1426-SM-CY8H8** matches both testis and liver from donor **GTEX-1ICLY** in the same extraction batch. The testis score is 10.9% and liver is 2.7%.
* Sigmoid colon **GTEX-16GPK-2526-SM-GQ1DH** matches testis from donor **GTEX-15UF7** in its extraction batch. A cerebellum library from donor GTEX-16GPK independently matches liver from GTEX-15UF7 in a shared extraction batch.
* Thyroid **GTEX-1HSEH-0826-SM-B2LXE** and tibial nerve **GTEX-13NYS-0426-SM-5MR4X** show strong testis signals. Their best reference-library matches lack the required batch support. These are separate from the thyroid/muscle and tibial nerve/pancreas pairs in the ranked table.

The inputs are the checksum-pinned GTEx files used by the catalog and the [v11 sample annotations](https://storage.googleapis.com/adult-gtex/annotations/v11/metadata-files/GTEx_Analysis_v11_Annotations_SampleAttributesDS.txt). To rebuild, run `PYTHONPATH=. python scripts/audit_gtex_cross_tissue_signals.py /path/to/tissue-shards /path/to/GTEx_Analysis_v11_Annotations_SampleAttributesDS.txt` from a source checkout with the GTEx median matrix staged.

## Independent donors and HPA biopsies

Each of the 52 GTEx primary sites has one library per donor within that site. Site counts range from 11 (kidney medulla) to 818 (skeletal muscle). A person may contribute several different tissues. The 186 HPA v23 biopsies span 40 tissue types, with one to 13 samples per tissue. This HPA matrix omits donor and batch identifiers. [Counts by source and tissue](sample-coverage.csv)

Of the HPA biopsies, 176 from 34 tissues had at least two same-tissue peers. For each marker, the screen subtracts the median of those other biopsies and scales positive excess by its GTEx source-organ reference. An advisory flag requires a panel median of at least 1% and at least 75% of the markers above every peer. Panels require ten measured, identity-mapped loci and at least half their GTEx markers. Five biopsies flagged. The strongest signal, whole blood in bone marrow, illustrates how tissue composition can explain a positive result. HPA's small peer groups and cross-atlas scaling limit comparisons with GTEx flag rates. [HPA hits](hpa-signals.csv) · [Panel coverage](hpa-panels.csv)

**O** denotes an expression outlier. **M** denotes a GTEx outlier with a different-donor match supported by batch metadata. Four M libraries attain the combined normal/Disposable maximum for 162 of the 918 assessed proteins. The esophagus muscularis specimen accounts for 153 of them. A sample flag describes the whole library. The affected gene's cellular origin requires separate evidence. [Flagged maxima, second donor in the same GTEx site, and sensitivity to M libraries](flagged-maxima.csv)

## Sensitivity of the tumor 10× threshold

The analysis reconstructs within-sample protein-group sums and reproduces every assessed group's catalog baseline. Each alternative is computed within each healthy tissue site before taking the maximum across sites and the remaining normal and Disposable references. Reproductive biopsies remain in their assigned compartments. Individual GTEx samples contribute once per donor within a site. P99 is reported both as an observed donor rank and with linear interpolation. P95 and P75 use linear interpolation. When both sample sources vary, HPA retains its maximum for tissues with fewer than three measured biopsies. Every scenario retains the 1 TPM floor. The outcome below is the tumor expression gate alone.

| Healthy tissue statistic | Of 918 over tumor 10× | Of the 98 former CTAs now RA-only over tumor 10× |
|---|---:|---:|
| Current maximum | 356 | 0 |
| Maximum after omitting GTEx M libraries | 366 | 10 |
| Maximum after omitting only source-panel observations in matched GTEx libraries | 365 | 9 |
| Second-highest GTEx donor in each site | 410 | 33 |
| GTEx P99, nearest rank | 473 | 84 |
| GTEx P99, interpolated | 481 | 89 |
| GTEx 95th percentile | 515 | 95 |
| GTEx 75th percentile | 536 | 98 |
| GTEx and HPA P99, interpolated | 482 | 89 |
| GTEx and HPA 95th percentile | 518 | 95 |
| GTEx and HPA 75th percentile | 552 | 98 |

Nearest-rank P99 selects the observed value at rank `ceil(0.99 × donors)` in ascending TPM order. This retains the maximum in sites with fewer than 100 donors. Across the 52 GTEx sites, five have fewer than 100. In the GTEx-only scenario, HPA biopsy maxima and other reference sources still contribute. For PRM2, nearest-rank P99 is 18.62 TPM in esophagus muscularis and 25.05 TPM in tibial nerve.

For an automated enrichment rule, use the highest tissue-specific nearest-rank P99 from GTEx alongside the existing other normal references, Disposable reference and 1 TPM floor. The existing orthogonal normal-expression gates act independently. This rule always resolves the influence of a single outlier when donor measurements are available. The catalog's existing required-coverage rule still applies to missing tissues. The maximum and the sample flags remain available to explain the decision. This is an evaluated policy option.

The source-panel maximum scenario tests a separate carryover hypothesis. It omits a GTEx group observation only when its sample has a batch-supported different-donor match and every member locus is a marker of that matched source organ. In a site with no remaining observation for the group, it uses the original maximum. It leaves outlier-only samples in the calculation.

For PRM2, the GTEx esophagus muscularis maximum is 3,047.26 TPM, the second-highest donor reaches 49.93, P95 is 5.65, and P75 is 2.12. Ten of 561 donors reach 10 TPM. In tibial nerve, the corresponding values are 1,495.04, 214.78, 6.04 and 1.72 TPM, and 19 of 670 donors reach 10 TPM. These are bulk RNA measurements.

[Scenario totals](normal-sensitivity.csv) · [All 918 per-protein ceilings and tumor peaks](protein-sensitivity.csv.gz) · [Original 98 transitions](../gtex-esophagus-muscularis/tumor-evidence-98.csv)

After running the GTEx screen above, run `PYTHONPATH=. python scripts/audit_healthy_sample_qc.py /path/to/tissue-shards` using the pinned HPA biopsy archive and the same catalog checkout.
