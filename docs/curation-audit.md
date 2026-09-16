# TMB, ICI and aPD1 curation audit

Reviewed 2026-09-15, against repository baseline `52f423e`.

The former `UCEC_POLE` 100% response anchor was a selected exceptional-responder
case, not a population response rate. It has been withdrawn. The broader audit
found unrelated citations, means entered as medians, mismatched source populations,
combined treatment arms, control arms labeled as ICI, and incorrect response
denominators. This change corrects confirmed failures and makes unresolved source
verification explicit; it does **not** certify every remaining number.

Issues were filed before their corresponding fixes:
[#523](https://github.com/pirl-unc/oncoref/issues/523) (POLE),
[#524](https://github.com/pirl-unc/oncoref/issues/524) (complete audit tracker),
[#525](https://github.com/pirl-unc/oncoref/issues/525) (regimens/populations/pooling),
[#526](https://github.com/pirl-unc/oncoref/issues/526) (TMB statistics/provenance),
[#527](https://github.com/pirl-unc/oncoref/issues/527) (incorrect citations), and
[#528](https://github.com/pirl-unc/oncoref/issues/528) (anchors/denominators), and
[#530](https://github.com/pirl-unc/oncoref/issues/530) (PR review and source recovery).

## Scope and reproducible evidence

Every row in the four numerical reference tables is inventoried, including blanks
and contextual estimates: 130 TMB, 110 ICI anchors, 84 aPD1 anchors, and 786 endpoint
estimates (1,110 total). No estimate IDs were deleted. The
[row ledger](audits/tmb-ici-apd1.csv) records each value, citation, locator status and
unresolved checks. Run:

```bash
python scripts/audit_tmb_ici_apd1.py --output docs/audits/tmb-ici-apd1.csv
```

The audit checks numeric validity, response/count arithmetic, confidence-interval
ordering, drug/regimen consistency, compact-anchor agreement, missing denominators,
and provenance. Its output is tested for reproducibility. These checks cannot prove
that a source supports the biological population assigned to a row.

The [citation inventory](audits/source-citations.csv) includes 183 distinct original
and current references. PMID/DOI metadata identifies the actual article behind a
citation; resolving an identifier is **not** numerical source validation. The
[TMB review table](../oncoref/data/cancer-tmb-source-audit.csv) records a disposition
for all 130 rows, including rejected legacy values, assay, locator, and review notes.

After source recovery there are 38 source-checked numeric TMB rows, 70 numeric entries
still requiring exact source verification, 10 withdrawn population-median claims,
and 12 other explicit gaps. Ten of the initial 20 withdrawals now have a verified
replacement or repaired citation; an existing MPN gap was also filled. The endpoint
inventory still flags 338 rows/anchors
whose numeric source locator is not verified, 33 without response denominators,
and 46 with denominators below 10. These counts overlap; they are review flags, not
338 independently demonstrated errors. Earlier `source_verified` and locator
`verified` fields remain historical curation metadata, not a fresh scientific
certification. #524 remains open for this source work.

## POLE: strong biology does not establish a response-rate ranking

[Mehnert 2016](https://www.jci.org/articles/view/84940) describes one patient's
partial response. The retained estimate now records a patient response count as
context, with no population ORR or binomial CI. Both ICI and aPD1 lookups for
`UCEC_POLE` return an explicit gap. Its unsourced TMB value of 100 mut/Mb is replaced
by a separately supported genomic median of **150.8 mut/Mb in 61 patients** from
[Nero 2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC11771542/), Results and Figure 1B.
That study uses a tumor-only TSO500 panel and the 11 pathogenic hotspot definition,
including 16 multiple-classifier cases. It is not an ICI response trial.

The supplied critique correctly rejects the 100% population claim, but several of
its supporting statements also need qualification:

* Prospective POLE-selected evidence exists. [Rousseau 2022, AcSé nivolumab](https://pmc.ncbi.nlm.nih.gov/articles/PMC9167784/)
  studied MMR-proficient advanced solid tumors, reporting responses in 5/11
  assessable patients with proofreading-deficient POLE variants. This is a
  pan-cancer result; it does not establish UCEC_POLE as the most responsive
  indication, and progression occurred in endometrial cancer despite pathogenic
  variants.
* [GARNET](https://pmc.ncbi.nlm.nih.gov/articles/PMC10643997/), Table 2, reports 2/5
  responses in its exploratory POLE-mutant category. Both responders were dMMR.
  That small, overlapping biomarker category is not a clean pathogenic-EDM UCEC
  population estimate.
* The final [avelumab endometrial trial](https://pmc.ncbi.nlm.nih.gov/articles/PMC9798913/)
  enrolled **no POLE-mutated tumors**. Its dMMR result cannot be treated as a
  measured POLE result. The [toripalimab basket trial](https://pmc.ncbi.nlm.nih.gov/articles/PMC11366758/)
  had only three assessable EDM cases; its two EDM responders had colon cancer.
* The original [2022 rectal dostarlimab report](https://doi.org/10.1056/NEJMoa2201445)
  reported 12 patients, not 42. The later [2025 expansion](https://doi.org/10.1056/NEJMoa2404512)
  reported 49/49 rectal and 35/54 nonrectal clinical complete responses. Neoadjuvant
  clinical CR and advanced-disease RECIST ORR are different endpoints and settings.
  These results cannot establish a universal numerical response ceiling.
* The [POLE pathogenicity study](https://pubmed.ncbi.nlm.nih.gov/31829442/)
  supports variant-level interpretation. The quoted hotspot list is not a complete
  or accurate whitelist: the recognized endometrial set includes L424I, M444K and
  D368Y. Exonuclease-domain location alone does not establish pathogenicity, and a
  non-EDM variant is not automatically proven biologically inert. Multiple-classifier
  assignment and treatment-response evidence are separate questions.

The cHL response anchor remains 71.2% (173/243), now citing the corresponding
[CheckMate 205 five-year report](https://pubmed.ncbi.nlm.nih.gov/37530622/), rather
than the older 69% report. We do not adopt the accompanying assertion that cHL is
uniformly low-TMB: the cited [Hodgkin genomic cohort](https://pmc.ncbi.nlm.nih.gov/articles/PMC6369943/)
uses clinical panels and does not support the former approximate WES median of 2.
No substitution-spectrum or netMHCpan binding claim was used to assign clinical
response rates.

## Confirmed numerical and source corrections

| Area | Previous claim | Corrected treatment |
| --- | --- | --- |
| POLE ICI/aPD1 | 100% ORR from a selected case | explicit evidence gap; case response retained as context |
| COAD, READ, UCEC | prevalence-weighted modeled ORRs | models retained as context; canonical anchors blank |
| UCEC_CNH/CNL | pMMR/MSS response assigned to molecular subtypes | gap: source does not isolate those subtypes |
| LUAD_STK11 | KRAS/STK11 evidence assigned to STK11/KEAP1 union; genomic n=637 used as response n | gap; invalid denominator removed |
| UVM ICI | mixed anti-PD-1/PD-L1 cohort labeled PD-1 | mixed-agent context; no PD-1 anchor |
| Ependymoma | 1/22 across mono and combo arms | nivolumab 1/12 = 8.3%; combination n=10 kept separate |
| DIPG, medulloblastoma | inferred zero ORRs with ordinary anchor confidence | formal-ORR gaps; outcome inferences remain context |
| GBM | enrolled n=184 used for ORR | response-evaluable n=153; bevacizumab comparator n=156 |
| MDS KEYNOTE-013 | enrolled n=28 used for response categories and calculated CIs | response-evaluable n=27; enrollment preserved separately |
| ESCA CheckMate 648 | 89/325 paired with 28% | 90/325, reported rounded 28%; correct source CI |
| ASPS disease control | 91% paired with 47/52 | descriptive count-derived 90.4%, explicitly distinguished from a reported DCR endpoint |
| FL TMB | 1.35, described as WES/mean | panel median 5.05 mut/Mb, n=119 |
| Pancreatic / midgut NET TMB | 1.9 / gap | advanced-site genome-wide WGS medians 1.35 / 1.05 |
| BCC / cSCC TMB | 65 / 45 with unrelated citations | Chalmers panel medians 47.3 (n=92) / 45.2 (n=266) |
| LCNEC TMB | George mean 8.6 entered as median | Chalmers panel median 9.9 (n=288), explicitly a different assay/cohort |
| MDS TMB | inferred conversion to 0.3/Mb | explicitly reported Chalmers median 0.8/Mb |
| Penile SCC TMB | 4.5 attributed to an erdafitinib trial | same value supported by Chalmers Table 1, n=60 |

Primary sources for these corrections:
[CheckMate 908](https://pmc.ncbi.nlm.nih.gov/articles/PMC10398811/),
[CheckMate 143](https://pubmed.ncbi.nlm.nih.gov/32437507/),
[KEYNOTE-013 Table 3](https://doi.org/10.1080/10428194.2022.2034155),
[CheckMate 648 regulatory report, Figure 2.3.2.2-6](https://www.fda.gov/media/182144/download),
[ASPS trial](https://pmc.ncbi.nlm.nih.gov/articles/PMC10729808/),
[FL Results 3.2](https://pmc.ncbi.nlm.nih.gov/articles/PMC12985221/),
[advanced NEN WGS](https://pmc.ncbi.nlm.nih.gov/articles/PMC8322054/), and
[Chalmers Table 1 and Results](https://pmc.ncbi.nlm.nih.gov/articles/PMC5395719/).

The rejected ADCC and RB values were means, not medians. ADCC now has a replacement
panel median; RB remains a gap. Other directly checked medians include rectal NET, GIST, DSRCT,
BTC, NSCLC and thoracic SMARCA4-deficient tumors; their source populations, assay
methods and specimen-versus-patient distinctions remain explicit.

## PR review: original sources and recovered medians

An unsupported citation is not evidence that no median exists. The
[20-row recovery ledger](audits/tmb-source-recovery.csv) preserves the original
value/citation, what the claimed source actually reports, replacement evidence,
and the reason for retaining any gap. Review findings were filed in #530 before
the follow-up correction.

**P1: retained Chalmers claims disagreed with its actual supplementary table.**
Additional file 3, Table S1 supports ACC 2.7 (n=204), CHOL 2.5 (liver cohort,
n=1456), KIRP 2.7 (n=152), LIHC 3.6 (n=602), PAAD 1.8 (n=2483), THCA 1.8
(n=350), THYM 1.3 (n=108), THYMCA 2.5 (n=168), and UCS 3.6 (n=245).
The old numbers did not match. Those rows are corrected, and sample counts are
filled for 25 checked Chalmers rows, including an MPN median of 0.8 (n=138).
The [extracted source rows](audits/tmb-chalmers-source-rows.csv) record exact workbook
cell ranges, the download URL and SHA-256. No published medians were averaged.

**P2: withdrawal was unnecessarily final for recoverable claims.**

| Code | Recovered median, mut/Mb | n | Source and limit |
| --- | ---: | ---: | --- |
| UVM | 0.34 | 80 | [Wu 2019 Table 1](https://pmc.ncbi.nlm.nih.gov/articles/PMC6944566/); old number was correct, DOI was wrong |
| NUTM | 1.0 | 71 | [Kim 2025](https://pubmed.ncbi.nlm.nih.gov/40704901/); known-TMB subset of 116 registry patients; heterogeneous panels |
| UCEC_POLE | 150.8 | 61 | Nero 2025; pathogenic hotspot cohort including multiple classifiers |
| MCL | 3.3 | 75 | Chalmers Table S1, lymph-node mantle-cell cohort |
| MTC | 1.8 | 96 | Chalmers Table S1, medullary thyroid |
| ACINIC | 1.8 | 81 | Chalmers Table S1, salivary acinic-cell cohort |
| ADCC | 1.8 | 184 | Chalmers Table S1, salivary adenoid cystic cohort |
| VSCC | 5.2 | 72 | Chalmers Table S1, vulvar SCC |
| ANSC | 5.4 | 232 | Chalmers Table S1, anal SCC |
| SARC_CHON | 1.7 | 93 | Chalmers Table S1, bone cohort; not pooled with soft-tissue chondrosarcoma |

All Chalmers replacements above come from the actual
[published supplement](https://static-content.springer.com/esm/art%3A10.1186%2Fs13073-017-0424-2/MediaObjects/13073_2017_424_MOESM3_ESM.xlsx).
Its coding-panel definition includes synonymous substitutions and indels, then
filters drivers and germline variants. Wu uses nonsynonymous coding mutations
divided by 38 Mb. These medians cannot be interpreted as assay-harmonized estimates.

For UCEC_MSI, [Leon-Castillo 2020 Table 2](https://pmc.ncbi.nlm.nih.gov/articles/PMC7065171/)
provides a median of 21.5 in 127 MSI-H, POLE-wild-type tumors. This replaces 18,
which the old TCGA citation did not establish as a median. The reanalysis includes
synonymous coding variants and uses a 38 Mb denominator.

Ten original gaps remain: RB, HL, CTCL, HCL, BRCA_Normal, LUAD_EGFR, UCEC_CNL,
UCEC_CNH, NBL_MYCNamp and NBL_MYCNnonamp. These are specific unresolved claims,
not declarations that the diseases have no TMB literature. For example:

* The originally cited salivary paper reports alterations per tumor and fractions
  above a TMB threshold, not the claimed ACINIC median of 2.6.
* The [neuroblastoma paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC6624164/) does
  mention an overall approximately 0.6 exonic mut/Mb median in its introduction,
  citing earlier work. It does not establish separate MYCN subgroup medians.
* EGFR evidence is available: [Hastings 2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6683857/)
  reports 3.8 in 383 EGFR-mutant lung cancers, and
  [Offin 2019](https://pubmed.ncbi.nlm.nih.gov/30045933/) reports 3.77 in 153
  metastatic exon19del/L858R cases. These are recorded as candidates; assigning
  them to the complete LUAD_EGFR entity requires resolving histology/allele scope.
* The primary Lawrence 2013 text establishes the median statistic and reports AML
  0.37/Mb, but the precise remaining per-type claims require its sample supplement.
  Download attempts returned access challenges, not the data. These remain
  explicitly unverified; reading the abstract does not validate the numbers.

Runtime review found no additional failures in the current regimen-fallback,
overlap-exclusion and gap-resolution cases. Passing software tests does not close
the remaining scientific verification work in #524.

Canonical response anchors now agree numerically with their selected endpoint
rows. This includes 5.3% for pretreated TNBC, 14.6% for PD-L1-positive cervical
cancer, 19.3% for response-evaluable ESCC, 7.8% for GBM, 41.6% for intermediate/poor
risk RCC on nivolumab/ipilimumab, 8.5% for mesothelioma, 25.9% for PD-L1-positive
NPC, 5% for the selected PD-L1-positive prostate cohort, 43.7% for melanoma, and
11.6% for gastric/GEJ cancer. These are **selected source settings**, not universal
unselected rates for each ontology label. Prior anchors that rounded or blended
these values are preserved in git history and described in row notes.

## Interpretation and safeguards

* Default pooling selects one regimen and primary evidence only. Opting into
  alternates retains source rows, but shared citation/NCT/trial identities block
  a numerical pool because overlap has not been excluded. Invalid count/rate
  pairs also block pooling. Independent trials can still differ clinically;
  this guard is not a meta-analysis or a proof of cohort independence.
* Drug labels now separate monotherapy, PD-1/CTLA-4, PD-L1 combinations,
  chemotherapy-containing arms, CTLA-4 alone, non-ICI controls, and mixed-agent
  cohorts. Controls, inferred outcomes, modeled blends and out-of-population
  context do not contribute to pools. The forest plot explicitly requests
  alternate trial points while retaining a primary-only numerical summary.
* DART high-grade pan-NEN evidence is not LCNEC-specific; differentiated thyroid
  context is not ATC evidence; overall soft-tissue sarcoma is not synovial sarcoma.
  These rows remain inspectable, but cannot supply those subtype pools.
* The legacy `apd1` table includes explicitly tagged fallback regimens. Filter
  `drug_target == 'PD-1'` for monotherapy analyses. Inspect endpoint population,
  denominator, biomarker selection and evidence inheritance alongside every value.
* TMB `published_median` provenance requires a checked source. Other retained
  numeric values are estimates awaiting verification, even if their old
  confidence label was high. For the reviewed subset:

```python
from oncoref import tmb

reviewed = tmb.cancer_tmb_df().query("source_review_status == 'source_checked'")
# Keep tmb_assay, source_scope, source_review_notes and n_samples in comparisons.
```

Open source work under #524 includes reproducing exact Lawrence supplementary
medians, resolving unmatched Chalmers populations, validating pediatric and molecular-subtype proxies,
resolving review-derived citations, and checking all remaining ICI endpoint
populations and denominators against their original tables. No inference from TMB,
mutation mechanism, binding predictions, or a response-rate ranking fills those gaps.
