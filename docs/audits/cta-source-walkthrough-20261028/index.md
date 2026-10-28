# CTA and oncofetal selection

**466 CTA groups (484 loci), 17,718 Unsafe, 17 Neither and 1,662 unassigned**, across 19,863 protein-equivalence groups. Policy: `oncoref.cta-selection.2026-10-28`. This is the current PR policy; the [1.21.0 release](../cta-release-1.21.0/index.md) remains a frozen comparison.

[Slides](slides/cta-evidence-walkthrough-expanded.pptx) · [Dog evidence](dog.md) · [Shared API](api.md) · [Policy and IHC matching](policy-review.md) · [All decisions](../cta-selection-20261028/protein-group-classifications.csv.gz) · [Changed decisions](../cta-selection-20261028/classification-changes.csv.gz)

## Categories

- **CTA:** reproductive or embryonic association, measured tumor expression, required coverage and acceptable normal expression all pass.
- **Unsafe:** above the single RNA ceiling in Vital/Protected, accepted non-Low protein in those compartments, or retained clinical exclusion evidence. This is a research-policy exclusion, not demonstrated toxicity.
- **Neither:** an affirmative low-expression profile or diffuse-Manageable profile, with complete required RNA coverage and no unresolved evidence. Failed candidate gates alone do not establish Neither.
- **Unassigned:** missing/uncertain evidence or an observed profile that fits none of those definitions. In the API, category is null and assessment_status explains why.

The separate **2× Strong boundary is removed**. Every layer evaluates every group; Neither is assigned only after its complete-profile check. DICE and corneal coverage remain optional, but observed counterevidence contributes.

## Neither profiles

The shared appreciable-expression floor is **2** in native TPM/nTPM, matching the tumor gate. Corneal CP10K uses the existing nominal /100 policy. This is an operational threshold, not a source-validated detection limit.

| Profile | Required evidence | Groups |
|---|---|---:|
| Low expression | All measured normal, reproductive, embryonic and tumor RNA below 2; required normal panels and mapped tumor/developmental measurements present; no detected protein | 0 |
| Diffuse Manageable | Fails both reproductive/developmental association routes; RNA ≥2 in at least two distinct Manageable tissue concepts; required normal coverage; no above-cutoff or accepted non-Low Vital/Protected protein evidence | 17 |

Colon and rectum count as one tissue for breadth. A high reproductive measurement cannot conceal broad normal expression. Missing optional measurements are disclosed, not set to zero. No group currently meets the complete low-expression profile; this does not establish that such genes do not exist. Unassigned statuses: {'incomplete_evidence': 1329, 'outside_category_definitions': 333}.

## What the assessment statuses mean

| API status | Meaning | Human groups |
|---|---|---:|
| `screened` | A CTA, Unsafe or Neither category has been assigned by this policy. It does not mean safe, clinically validated or complete coverage of every source. | 18,201 |
| `incomplete_evidence` | No category was assigned, and required normal coverage, mapped tumor/developmental measurements or evidence resolution is missing. | 1,329 |
| `outside_category_definitions` | No category was assigned; the tracked coverage/resolution checks pass but the observed profile meets none of the three definitions. | 333 |
| `unassessed` | No species-specific quantitative policy has been applied; used for native dog nominations. | 0 |

An available exclusion can establish Unsafe despite another missing assay. A reproductive-route CTA need not have embryonic measurements. The low-expression Neither profile needs both tumor and embryonic coverage; the diffuse-Manageable profile has different requirements.

Among the 1,329 incomplete assessments, gaps overlap:

| Evidence gap | Groups |
| --- | --- |
| normal panel missing | 574 |
| tumor measurement missing | 71 |
| embryonic measurement missing | 1011 |
| unresolved evidence | 33 |

[Per-group gaps, including all member gene IDs](unassigned-evidence-gaps.csv). MAGE-A1 fails the reproductive association gate (4.4 nTPM) and lacks embryonic data. MAGE-A3 has colon undifferentiated-cell RNA of 53.1 nTPM above its Manageable ceiling of 16, and also lacks embryonic data. That missing measurement is not the cause of its colon hold; it limits further profile assessment.

## Sources of candidate annotations

These published lists annotate the selection; they do not provide a gene-specific admission bypass. The union is **2,482 mapped groups / 2,538 loci**. **421** selected groups have a listed-paper annotation; **45** do not. Source contribution counts overlap. **Mapped groups means curated nominations linked to the current reference, not newly discovered genes or selected CTAs.** The CTA column reports the subset passing the current funnel.

| citation | mapped_entities | CTA | Unsafe | Neither | Unassigned | source_url |
| --- | --- | --- | --- | --- | --- | --- |
| da Silva et al. (2017) | 1017 | 358 | 327 | 1 | 331 | https://doi.org/10.18632/oncotarget.21715 |
| Wang et al. (2016) | 988 | 266 | 587 | 1 | 134 | https://doi.org/10.1038/ncomms10499 |
| Chang et al. (2019) | 985 | 313 | 588 | 0 | 84 | https://doi.org/10.1002/cam4.2223 |
| Gong et al. (2021) | 729 | 361 | 250 | 0 | 118 | https://doi.org/10.1038/s41467-021-22695-y |
| Bruggeman et al. (2018) | 692 | 27 | 643 | 2 | 20 | https://doi.org/10.1038/s41388-018-0357-2 |
| Jamin et al. (2021) | 469 | 33 | 422 | 0 | 14 | https://doi.org/10.1002/1878-0261.12900 |
| Loriot et al. (2025) | 202 | 82 | 87 | 0 | 33 | https://doi.org/10.1371/journal.pgen.1011734 |
| Carter et al. (2023) | 103 | 59 | 27 | 0 | 17 | https://doi.org/10.1136/jitc-2023-007935 |
| Seager et al. (2024) | 15 | 9 | 2 | 0 | 4 | https://doi.org/10.1186/s12967-024-04918-0 |
| Bradley et al. (2020) | 9 | 1 | 8 | 0 | 0 | https://doi.org/10.1038/s41467-020-19141-w |
| Garg et al. (2007) | 5 | 2 | 2 | 0 | 1 | https://doi.org/10.1158/1078-0432.CCR-06-2340 |
| Boël et al. (1995) | 2 | 0 | 0 | 0 | 2 | https://doi.org/10.1016/s1074-7613(95)80053-0 |
| Bai et al. (2016) | 1 | 0 | 1 | 0 | 0 | https://doi.org/10.1158/0008-5472.CAN-16-0225 |

### Expanded author-recognized CTA sets

The original targeted intake recorded only each paper's principal target. The current paper counts also include explicitly recognized prior CT antigens, with a separate evidence role and source location. Garg therefore contributes five mapped groups and Boël two, plus its still-unresolved BAGE target. These additions are already in other lists: the overall union and selection counts do not change.

| Paper | Principal target | Previously recognized CTAs | Mapped groups | Selected CTA | Source |
| --- | --- | --- | --- | --- | --- |
| Bai et al. (2016) | EGFL6 | None curated | 1 | 0 | https://doi.org/10.1158/0008-5472.CAN-16-0225 |
| Garg et al. (2007) | SPAG9 | SPA17; SYCP1; NY-ESO-1 (CTAG1A/B); SSX4/B | 5 | 2 | https://doi.org/10.1158/1078-0432.CCR-06-2340 |
| Boël et al. (1995) | BAGE | MAGE-A1; MAGE-A3 | 2 | 0 | https://doi.org/10.1016/s1074-7613(95)80053-0 |

Bai studies EGFL6 as an ovarian-cancer therapeutic target, without establishing a CTA-restricted set. Its larger analysis contains **538 EGFL6-correlated genes**, not 538 proposed CTAs; signaling partners such as SHP2/ERK are not added as CT antigens. Garg explicitly compares SPAG9 with SPA17, SYCP1, NY-ESO-1 and SSX4. Boël recognizes MAGE-1/MAGE-3 as prior testis-restricted tumor antigens alongside its BAGE discovery. Family mentions without resolved members are not expanded into all paralogs. [Bai et al.](https://pubmed.ncbi.nlm.nih.gov/27803106/), [Garg et al.](https://pubmed.ncbi.nlm.nih.gov/17332284/), [Boël et al.](https://pubmed.ncbi.nlm.nih.gov/7895173/).

EGFL6 and SPAG9 occur in no other currently registered paper membership, though EGFL6 also has a separate CTpedia annotation. BAGE occurs in Seager et al. (2024) as an unresolved target label: that textual overlap does not enter the mapped-gene overlap matrix. These are coverage statements about our curation, not claims that no other literature studies them. Bradley et al. (2020) is different: its ten curated genes collapse to nine protein groups, of which **PLAC1 is the one selected CTA**; the other eight are Unsafe.

| Paper | Target | Category | Other registered papers |
| --- | --- | --- | --- |
| Bai et al. (2016) | EGFL6 | Unsafe | None |
| Garg et al. (2007) | SSX4/B | Unassigned | Jamin et al. (2021); da Silva et al. (2017) |
| Garg et al. (2007) | SPA17 | Unsafe | Chang et al. (2019); Wang et al. (2016); da Silva et al. (2017) |
| Garg et al. (2007) | SYCP1 | CTA | Chang et al. (2019); Wang et al. (2016); Gong et al. (2021); Loriot et al. (2025); da Silva et al. (2017) |
| Garg et al. (2007) | NY-ESO-1 | CTA | Jamin et al. (2021); Wang et al. (2016); Seager et al. (2024); Loriot et al. (2025); da Silva et al. (2017) |
| Garg et al. (2007) | SPAG9 | Unsafe | None |
| Boël et al. (1995) | MAGEA3 | Unassigned | Jamin et al. (2021); Chang et al. (2019); Wang et al. (2016); Bruggeman et al. (2018); Gong et al. (2021); Carter et al. (2023); Seager et al. (2024); Loriot et al. (2025); da Silva et al. (2017) |
| Boël et al. (1995) | MAGEA1 | Unassigned | Jamin et al. (2021); Chang et al. (2019); Wang et al. (2016); Gong et al. (2021); Carter et al. (2023); Seager et al. (2024); Loriot et al. (2025); da Silva et al. (2017) |

[Expanded memberships and group IDs](targeted-paper-overlaps.csv) · [Evidence roles and exact source locations](publication-context.csv) · [Every registered paper table and its extraction scope](paper-source-scope.csv).

The larger papers already contribute their registered gene sets, including Wang's 1,019 CT-pattern genes, Bruggeman's 756 germ-cell genes, da Silva's 1,103 testis-biased and nested 745 CT genes, Chang's 1,036 TGCT genes, Carter's 103 final CT genes, and Loriot's 146 strict plus 134 preferential annotations. Input rows, gene loci and collapsed protein groups differ; these numbers must not be summed. Gong contributes reproductive/placental expression nominations rather than a claimed CTA list. This audit extends the three targeted entries; it does not claim a fresh exhaustive full-text review of every large source.

EGFL6 is excluded by lung bulk RNA (24.8 nTPM), lung smooth-muscle RNA (363.3 nTPM) and corneal evidence. SPAG9 has accepted bronchial protein staining and broad RNA counterevidence, including spinal cord (66.3 nTPM) and brain oligodendrocytes (287.1 nTPM). These are current quantitative-policy exclusions, irrespective of the original paper's terminology.

![Published source contributions](figures/published-list-contributions.png)
![Published source overlaps](figures/published-list-overlaps.png)

[All pairwise overlaps](human-overlaps.csv.gz) · [Memberships and gene IDs](human-memberships.csv.gz)

## Why roughly 12,000 groups have reproductive or developmental association

The starting association set comes from whole-reference **RNA measurements**, not a 12,000-gene literature list. Reproductive association requires HPA ≥5 nTPM and greater than adult nonreproductive bulk maximum. The developmental route requires a Tyser lineage mean ≥5 TPM, ≥20 cells, ≥5 positive cells and ≥10% positive at ≥1 TPM. It establishes embryonic presence, not embryonic restriction. Later fetal and child/adult developmental atlases are not yet integrated.

| route | entrants | CTA |
| --- | --- | --- |
| Reproductive only | 1231 | 428 |
| Developmental only | 9481 | 29 |
| Both | 1427 | 9 |

[Tyser et al. (2021)](https://doi.org/10.1038/s41586-021-04158-y), [HPA v23 downloads](https://v23.proteinatlas.org/about/download).

## Measurement sources and filtering rules

| name | evidence_category | native_unit | status | source_url |
| --- | --- | --- | --- | --- |
| HPA v23 / GTEx consensus | normal_bulk_rna | nTPM | used_in_selection | https://v23.proteinatlas.org/about/download |
| HPA v23 / GTEx consensus | reproductive_expression | nTPM | used_in_selection | https://v23.proteinatlas.org/about/download |
| Tyser human gastrulation atlas | developmental_expression | TPM | used_in_selection | https://doi.org/10.1038/s41586-021-04158-y |
| OncoRef source-mapped tumor union | tumor_expression | TPM | used_in_selection | https://github.com/pirl-unc/oncoref/tree/v1.21.0/docs/audits/cta-expanded-funnel-20261027/tumor-union |
| HPA cancer sample measurements | tumor_expression | pTPM | used_in_selection | https://v23.proteinatlas.org/about/download |
| HPA tissue/cell clusters | normal_cell_rna | nTPM | used_in_selection | https://v23.proteinatlas.org/about/download |
| HPA immune summaries | immune_rna | nTPM | used_in_selection | https://v23.proteinatlas.org/about/download |
| Monaco immune summaries | immune_rna | TPM | used_in_selection | https://v23.proteinatlas.org/about/download |
| DICE immune distributions | immune_rna | TPM | used_in_selection | https://dice-database.org/ |
| Collin adult corneal atlas | normal_cell_rna | CP10K | used_in_selection | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE155683 |
| HPA normal-tissue IHC | normal_protein | HPA categorical level and reliability | used_in_selection | https://v23.proteinatlas.org/about/assays%2Bannotation |
| Curated external staining observations | normal_protein | native categorical findings | used_in_selection | https://github.com/pirl-unc/oncoref/tree/v1.21.0/docs/audits/external-ihc-20261027 |
| Reviewed normal protein and clinical exclusions | normal_protein | categorical reviewed findings | used_in_selection | https://github.com/pirl-unc/oncoref/blob/v1.21.0/oncoref/data/cta-policy-findings.csv |
| Reviewed normal protein and clinical exclusions | clinical_safety | categorical reviewed findings | used_in_selection | https://github.com/pirl-unc/oncoref/blob/v1.21.0/oncoref/data/cta-policy-findings.csv |
| Gong et al., Nature Communications (2021) | nomination | gene-list membership | annotation | https://doi.org/10.1038/s41467-021-22695-y |
| Gong et al., Nature Communications (2021) | nomination | gene-list membership | annotation | https://doi.org/10.1038/s41467-021-22695-y |
| Bradley et al., Nature Communications (2020) | nomination | gene-list membership | annotation | https://doi.org/10.1038/s41467-020-19141-w |
| Wang 2016 | nomination | gene-list membership | annotation | https://doi.org/10.1038/ncomms10499 |
| Bruggeman 2018 | nomination | gene-list membership | annotation | https://doi.org/10.1038/s41388-018-0357-2 |
| da Silva 2017 · testis-biased | nomination | gene-list membership | annotation | https://doi.org/10.18632/oncotarget.21715 |
| da Silva 2017 · CT | nomination | gene-list membership | annotation | https://doi.org/10.18632/oncotarget.21715 |
| Jamin 2021 | nomination | gene-list membership | annotation | https://doi.org/10.1002/1878-0261.12900 |
| Jamin 2021 · core | nomination | gene-list membership | annotation | https://doi.org/10.1002/1878-0261.12900 |
| Chang 2019 · TGCT | nomination | gene-list membership | annotation | https://doi.org/10.1002/cam4.2223 |
| Carter 2023 | nomination | gene-list membership | annotation | https://doi.org/10.1136/jitc-2023-007935 |
| Seager 2024 · panel | nomination | gene-list membership | annotation | https://doi.org/10.1186/s12967-024-04918-0 |
| Loriot 2025 · strict | nomination | gene-list membership | annotation | https://doi.org/10.1371/journal.pgen.1011734 |
| Loriot 2025 · preferential | nomination | gene-list membership | annotation | https://doi.org/10.1371/journal.pgen.1011734 |
| Gong 2021 · reproductive coding | nomination | gene-list membership | annotation | https://doi.org/10.1038/s41467-021-22695-y |
| Gong 2021 · reproductive noncoding | nomination | gene-list membership | annotation | https://doi.org/10.1038/s41467-021-22695-y |
| da Silva 2017 · tumor proteomics | nomination | gene-list membership | annotation | https://doi.org/10.18632/oncotarget.21715 |
| Bai 2016 | nomination | gene-list membership | annotation | https://doi.org/10.1158/0008-5472.CAN-16-0225 |
| Boël et al., Immunity (1995) | nomination | gene-list membership | annotation | https://doi.org/10.1016/s1074-7613(95)80053-0 |
| Garg et al., Clinical Cancer Research (2007) | nomination | gene-list membership | annotation | https://doi.org/10.1158/1078-0432.CCR-06-2340 |

| stage | rule |
| --- | --- |
| Reproductive or developmental expression | Reproductive ≥5 nTPM and > adult bulk maximum OR Tyser embryonic lineage mean ≥5 TPM; ≥20 cells, ≥5 positive and ≥10% positive at ≥1 TPM. |
| Tumor sample ≥2 TPM/pTPM | At least one measured compatible tumor sample ≥2 in native TPM or HPA pTPM; do not add maxima from different samples. |
| HPA normal protein | Reliable non-Low positives in Vital/Protected supply Unsafe evidence; Manageable positives block CTA only; uncertain evidence remains unresolved. Low is flagged. |
| Required RNA coverage | Required bulk and tissue-cell/immune panels; missing members are not zeros. DICE and cornea are optional coverage. |
| Normal bulk | One inclusive tissue/IHC ceiling: above it in Vital/Protected supplies Unsafe evidence; above it in Manageable blocks CTA. No second RNA boundary. |
| HPA tissue clusters | RNA requires ≥20 reads and uses the single ceiling. IHC credit requires an exact or declared compartment match; unrelated cells in the same tissue cannot inherit it. |
| HPA immune | Protected-tier native nTPM limits without matched negative-IHC credit; reproductive/placental location remains neutral. |
| Monaco immune | Protected-tier native TPM limits without matched negative-IHC credit. |
| DICE immune | Third-highest reported TPM value; native protected-tier limits. Unpaired group distributions retain bounds. Missing coverage alone does not exclude. |
| Adult cornea | Recurrent native CP10K: ≥5 positive cells and ≥20 counts in ≥2 sample groups. Nominal RNA limits /100; not a TPM conversion. Missing coverage alone does not exclude. |
| Accepted normal-protein evidence | Accepted positive protein findings exclude unless extracted Low intensity retains a flag; unquantified does not mean negative. |
| Clinical safety exclusion | Retain reviewed product-specific clinical injury counterevidence and causal-attribution limits; no admission shortcut from clinical use. |
| External validated normal protein | Use categorical/quantitative staining only within verified reagent, tissue and cell scope; no inherited HPA reliability label. |
| Neither expression profiles | Neither requires complete required RNA coverage and either all measured RNA below the shared 2-unit floor with no detected protein, or failed reproductive/developmental association with appreciable RNA in >=2 Manageable tissue concepts and no Vital/Protected exclusion. Missing/unresolved evidence and other profiles are unassigned. |

![RNA cutoffs](figures/active-rna-cutoffs.png)

The negative-IHC columns require a usable profile control, complete group-member coverage and a supported tissue/cell match. Same-tissue IHC for an unrelated cell cannot increase a ceiling. Values retain their native units; pTPM denotes protein-coding **RNA** TPM. Tumor eligibility uses ≥2 in at least one compatible OncoRef TPM or HPA pTPM sample, without summing maxima from different samples.

## Filtering funnel

| layer | new_candidate_removals | new_Unsafe | CTA_eligible_so_far | Unsafe | Neither | Unassigned |
| --- | --- | --- | --- | --- | --- | --- |
| Reproductive / embryonic association | 7724 | 0 | 12139 | 0 | 0 | 7724 |
| Tumor sample ≥2 TPM/pTPM | 7 | 0 | 12132 | 0 | 0 | 7731 |
| HPA normal protein | 8850 | 10094 | 3282 | 10094 | 0 | 6487 |
| Required RNA coverage | 90 | 0 | 3192 | 10094 | 0 | 6577 |
| Normal bulk | 2288 | 6441 | 904 | 16535 | 0 | 2424 |
| HPA tissue clusters | 393 | 988 | 511 | 17523 | 0 | 1829 |
| HPA immune | 8 | 48 | 503 | 17571 | 0 | 1789 |
| Monaco immune | 11 | 38 | 492 | 17609 | 0 | 1762 |
| DICE immune | 21 | 88 | 471 | 17697 | 0 | 1695 |
| Adult cornea | 4 | 18 | 467 | 17715 | 0 | 1681 |
| Accepted normal-protein evidence | 1 | 2 | 466 | 17717 | 0 | 1680 |
| Clinical safety exclusion | 0 | 1 | 466 | 17718 | 0 | 1679 |
| External validated normal protein | 0 | 0 | 466 | 17718 | 0 | 1679 |
| Neither expression profiles | 0 | 0 | 466 | 17718 | 17 | 1662 |

![Human funnel](figures/human-funnel.png)

[Independent layer effects](../cta-selection-20261028/independent-layer-effects.csv) distinguish overlapping evidence from order-dependent first removals. All counts use protein-equivalence groups; gene-locus counts are explicitly labelled.

## Selected CTA expression

These plots compare reproductive nTPM against the highest measured embryonic lineage mean TPM. **204 of 466 selected groups (44%) have both coordinates; 262 (56%) lack an embryonic coordinate.** The plotted source is the EBI filtered gene matrix from one Carnegie-stage-7 embryo, approximately 16–19 days after fertilization—not a broad fetal atlas. PRAME, MAGE-A4, CT83, CTAG1A/B and XAGE1A/B have no gene rows in that downloaded matrix. Missing rows alone cannot establish zero or absent developmental expression.

The follow-up checked the study authors' **57,490-column raw-count matrix**, separately from the atlas's 25,330-row filtered matrix. The atlas download labelled “Raw counts” is also filtered and retains the same gene rows. Direct author observations are:

| source_symbol | gene_column_present | cells | total_reported_reads | positive_cells |
| --- | --- | --- | --- | --- |
| PRAME | True | 1195 | 0 | 0 |
| MAGEA4 | True | 1195 | 14 | 1 |
| CTAG1A | True | 1195 | 0 | 0 |
| CTAG1B | True | 1195 | 0 | 0 |
| XAGE1A | False | 1195 | — | — |
| XAGE1B | True | 1195 | 0 | 0 |
| CT83 | True | 1195 | 0 | 0 |
| TBXT | True | 1195 | 3.115e+04 | 226 |

PRAME, CT83, CTAG1A/B and XAGE1B have zero reported reads across these 1,195 cells; MAGE-A4 has 14 reads in one cell. TBXT has 31,149 reported reads in 226 cells. XAGE1A has no separate column, so complete XAGE1A/B group coverage is not established. These are observed nondetections or sparse signals in one early embryo, **not absence throughout fetal development**. The exact per-gene reason for exclusion from the atlas remains unavailable; the authors' counts support low detection but come from a separate quantification. They are not converted to TPM or inserted into the selection funnel. [Author data](https://github.com/ScialdoneLab/human-gastrula-shiny), [count evidence and gene IDs](embryonic-comparator-coverage.csv), [source hashes](embryonic-coverage-receipt.json).

Each figure now shows the 204 measured groups above and the 262 missing-embryonic groups in a separate categorical band below. All selected groups are represented; the band has no embryonic-expression value. Gray points mean the healthy-expression color is unavailable. Labels highlight the widely studied antigens, including their gene-group names, and every measured group ≥5 on both axes. The latter is a descriptive expression criterion; developmental admission also requires cell support.

![Selected CTA reproductive versus embryonic expression](figures/selected-reproductive-embryonic.png)
![Selected CTA colored by normal immune TPM](figures/selected-reproductive-embryonic-immune.png)
![Selected CTA colored by adult corneal CP10K](figures/selected-reproductive-embryonic-cornea.png)

Color uses the maximum supported value across non-disposable healthy compartments within each unit: HPA bulk/cell/immune nTPM, Monaco/DICE TPM, or corneal CP10K. No maximum mixes units. Reproductive/developmental locations are excluded from the color scale. Missing embryonic measurements are not plotted at zero. [Every selected group's plot data](selected-expression-scatter.csv.gz).

## Widely studied antigens

| protein_group_label | category | assessment_status | reproductive_ntpm | developmental_all_lineage_max_tpm | adult_bulk_max_ntpm | exclusion_reasons |
| --- | --- | --- | --- | --- | --- | --- |
| PAGE1 | CTA | screened | 16.7 | — | 0 |  |
| SMC1B | CTA | screened | 30.2 | 0.4186 | 0.3 |  |
| CTAG2 | CTA | screened | 10.1 | — | 5.1 |  |
| MAGEA4 | CTA | screened | 21.5 | — | 0.6 |  |
| TBXT | CTA | screened | 0.2 | 245.1 | 2.5 |  |
| PRAME | CTA | screened | 102.2 | — | 3.6 |  |
| MAGEA1 | Unassigned | incomplete_evidence | 4.4 | — | 0 | Reproductive 4.4 nTPM versus adult bulk 0 nTPM; supported embryonic mean unavailable TPM. Neither association route passes. |
| CT83 | CTA | screened | 62.2 | — | 2.6 |  |
| MAGEA12 | Unsafe | screened | 2.7 | — | 0.5 | Reproductive 2.7 nTPM versus adult bulk 0.5 nTPM; supported embryonic mean unavailable TPM. Neither association route passes.  /  Retain reviewed product-specific clinical injury counterevidence and causal-attribution limits; no admission shortcut from clinical use. |
| MAGEA3 | Unassigned | incomplete_evidence | 7.8 | — | 0.2 | colon: undifferentiated cells; 53.1 nTPM (ceiling 16) |
| NY-ESO-1 | CTA | screened | 15.5 | — | 0.2 |  |
| SSX2/B | CTA | screened | 32.7 | 0.4237 | 0.5 |  |
| XAGE1A/B | CTA | screened | 98.8 | — | 1.5 |  |

Each comparison has its own full-width figure. The embryonic comparison retains missing-measurement antigens in a separate band.

![Widely studied antigens: reproductive expression](figures/widely-studied-antigens-reproductive.png)
![Widely studied antigens: embryonic expression](figures/widely-studied-antigens-embryonic.png)

[Detailed tissue/cell/value explanations](human-explanations.csv.gz), [all comparator contexts](../cta-selection-20261028/widely-studied-contexts.csv.gz), and [every previously annotated CTA and its current category](../cta-selection-20261028/known-cta-decisions.csv.gz). Previous CTA labels are comparison annotations only. Historical-cohort counts, where available: {'CTA': 162, 'Unsafe': 30, 'Unassigned': 21}.

## Evidence limits

Adult normal panels do not establish pediatric safety. Tyser contains one early embryo. RNA expression does not measure peptide-HLA presentation; HPA negative staining has finite sensitivity. PRAME informed earlier parameter choices and is not independent validation. All gene IDs, equivalence-group members, source versions and native units remain in the exports.
