# CTA release 1.19.1: the list, filtering rules and caveats

[Download this review PDF](cta-release-review.pdf) · [All supporting tables](cta-release-tables.zip) · [Search every gene](../cta-atlas-20261024/explore.html) · [Complete gene catalogs](../cta-atlas-20261024/index.md)

[Retained-gene review flags](retained-review.csv) · [Immune audit](immune-summary.csv) · [Activated-only review](activated-only-review.csv) · [Input provenance](provenance.json) · [Output checksums](artifacts.json) · [Review record](review.md)

## Release snapshot

oncoref 1.19.1 documents the current CTA list and its caveats. It adds this detailed review and an immune-state audit; it changes no gene, protein, tissue category or threshold. The 399 retained genes remain the frozen screening result under policy 2026-10-24.1, not a list with established clinical safety.

| Scope | Retained CTA | Neither | Excluded |
|---|---|---|---|
| Mapped CTA nominations | 399 | 561 | 1,580 |
| Canonical human genes | 399 | 1,999 | 17,699 |
| Complete unique proteins | 723 | 2,878 | 126,384 |

CTA means selected by this research screen. Unsafe means membership in the forbidden background, not necessarily observed clinical toxicity. Neither supplies neither targets nor default forbidden proteins. It includes uncertainty, weak target evidence and reconsiderable findings; it is not a safe reserve list. [Current high-normal-RNA policy and candidate follow-up](https://pirl-unc.github.io/oncoref/audits/cta-neither-review-20261024/)

The most consequential new literature finding is that 37 retained genes have normal-tissue or isoform conflicts. ABCC11 and MAGEC2 need priority neural protein review; YBX2 has vascular/adrenal staining to adjudicate. Their flags are visible beside the frozen decisions. This report does not silently clear those conflicts. [Complete CTA evidence atlas and gene catalogs](https://pirl-unc.github.io/oncoref/audits/cta-atlas-20261024/)

Read this review for rules and interpretation; use the companion atlas for every gene's full evidence. Its 138-page retained, 189-page Neither and 588-page excluded catalogs cover all 2,540 mapped nominations, with bookmarks and PDF page links. [Complete CTA evidence atlas and gene catalogs](https://pirl-unc.github.io/oncoref/audits/cta-atlas-20261024/)

Evidence labels matter: a trial registration is a nomination, human tumor recognition is evidence of recognition, and a clinical response applies to the studied product and patient selection. None automatically proves reproductive restriction or general safety.

## How a gene enters the retained set

Start from a mapped curated CTA nomination with usable normal-tissue coverage. Require measured reproductive expression and either reproductive bulk RNA >=5 nTPM, reliable positive reproductive protein, or cited reproductive expression plus tumor/antigen-recognition support. Good reproductive protein evidence can satisfy admission below 5 nTPM. CT55, with reproductive bulk RNA 3.7 nTPM, qualifies through published support. Tumor expression alone cannot override an exclusion. [Normal-cell screen: source-specific rules and measurements](https://pirl-unc.github.io/oncoref/audits/cta-normal-cells-20261023/)

All normal-tissue, donor, antibody, immune/cellular and safety gates still apply. Insufficient reproductive evidence or selectivity alone leads to Neither, not automatically Unsafe. Unresolved identity and missing observations remain explicit gaps. BAGE and GAGE2B remain recorded in intake but outside the 2,540 mapped-gene denominator.

| Tissue category | Examples | Treatment in this screen |
|---|---|---|
| Reproductive | Testis and male tract; ovarian/female tract; mammary; placenta, extraembryonic and explicitly fetal/embryonic compartments | Required positive expression; permitted normal context |
| Disposable | Mapped hair follicles, cutaneous melanocytes, tonsil, appendix, gallbladder | Permitted under conditional project assumptions |
| Manageable | Other skin, spleen, thymus, thyroid, adipose | Still screened; positive normal IHC excludes |
| Protected | Salivary glands, lymph nodes, cartilage, several mucosal and urinary compartments | Stricter screened normal contexts |
| Vital | Heart, lung, CNS, liver, kidney, pancreas, eye, muscle, blood/marrow and other critical compartments | Strictest normal RNA tier |

These examples are not the full anatomy registry. A fetal liver sample is reproductive by developmental scope; adult maternal liver is not. Immune cells inside a permitted organ retain their immune identity. The broad set is not a pregnancy- or fertility-sparing selection. Surgical removability does not establish tolerability of systemic immune injury.

## RNA ceilings, IHC confidence and selectivity

Accepted positive normal IHC, including Low, vetoes selection in Vital, Protected or Manageable compartments. Uncertain/cross-reactive primary staining is handled as a separate antibody review, not accepted as a definitive protein finding.

RNA relaxation requires complete negative IHC across the tissue's routine measured compartments plus a reliable positive reproductive control in the curated gene profile. Missing/uncertain or incomplete local staining, or no reproductive control, uses baseline limits. The export does not prove matched antibody/batch/donor across tissues or calibrate RNA-to-protein conversion.

| Category | Baseline | Approved negative | Supported negative | Enhanced negative |
|---|---|---|---|---|
| Vital | 2 / 80% | 4 / 75% | 6 / 70% | 12 / 65% |
| Protected | 6 / 75% | 8 / 70% | 12 / 65% | 24 / 60% |
| Manageable | 10 / 70% | 20 / 65% | 30 / 60% | 30* / 55% |

Cells show effective maximum nTPM / selectivity setting. *The underlying Enhanced Manageable ceiling is 60; the independent >30 rule makes the effective ceiling 30. The thresholds are project research settings, not clinically validated safety limits.

For each screened bulk tissue: N <= C; N <= R(1-s)/s + 0.1; and N < R. R is the highest reproductive nTPM, N the normal tissue's nTPM, C its absolute ceiling, and s its selectivity setting. All three comparisons must pass. For example, with R=5 and Vital baseline s=0.80, the selectivity ceiling is 1.35 nTPM. A normal value of 1.5 passes the absolute ceiling of 2 but fails selectivity.

This is a maximum-based comparison, not a global expression-breadth statistic. One high reproductive compartment can dominate R. The fixed 0.1 nTPM allowance softens near-zero ratios; it is not a fitted uncertainty interval. Broad low-level normal expression still merits scrutiny even when each individual comparison passes.

## Cellular evidence and the independent hard ceiling

Any observed >30 TPM or >30 nTPM outside explicitly mapped Reproductive/Disposable contexts yields Unsafe / High RNA. Exactly 30 does not trigger this gate. A single high observation suffices; no average, strong reproductive signal or negative IHC can waive it. Unresolved anatomy does not establish permission. This also overrides a prior Variable hold. [Current high-normal-RNA policy and candidate follow-up](https://pirl-unc.github.io/oncoref/audits/cta-neither-review-20261024/)

| Source | Lower review trigger | Usual result without a stronger exclusion |
|---|---|---|
| HPA tissue/cell clusters | >=5 nTPM and >=20 reads; normal cluster at least as high as the reproductive-cell maximum within that source. Substantial unresolved context or absent reproductive-cell coverage also prompts review. | Neither / Cellular |
| Original DICE | >=5 TPM in at least three reported values in any immune cell/state distribution | Neither / Cellular |
| HPA/Monaco sorted immune summaries | >=5 in the source's native nTPM or TPM; HPA's DICE summary is descriptive, not independent replication | Neither / Cellular |
| Adult corneal atlas | Mean per-cell CP10K >=1, >=5 positive cells and >=20 counts within a group; support in at least two groups | Neither / Cellular |
| Primary antibodies | Unresolved positive normal staining / curated-profile conflicts | Neither / Antibody |
| Curated normal protein literature | Accepted production by nonmalignant cells outside permitted contexts, including activation-dependent production | Unsafe / Normal |

Bulk RNA, sorted-cell RNA and cluster RNA are separate evidence sources. Even when labelled nTPM, their aggregation and normalization differ. TPM, nTPM and corneal CP10K are never pooled or numerically converted. The hard ceiling deliberately applies the same number separately to native TPM/nTPM, not to CP10K. Low read support, contamination and multimapping require explicit adjudication rather than silent exemptions. [Normal-cell screen: source-specific rules and measurements](https://pirl-unc.github.io/oncoref/audits/cta-normal-cells-20261023/)

## What remains: named targets and cancer occurrence

The following are examples, not a ranking of clinical benefit. Cancer RNA shows the highest reported cohort fraction at >=5 native pTPM among cohorts with at least 20 measured samples. TCGA and validation cohorts remain separate. Counts are samples, not the fraction of tumor cells expressing a target.

| Gene | Screen result | Example cancer RNA occurrence |
|---|---|---|
| MAGEA4 | CTA / Pass | Lung Squamous Cell Carcinoma (TCGA): 270/489 (55.2%) |
| PRAME | CTA / Pass | Skin Cuteneous Melanoma (TCGA): 94/99 (94.9%) |
| CTAG1A | CTA / Pass | Bladder Urothelial Carcinoma (TCGA): 24/169 (14.2%) |
| CTAG1B | CTA / Pass | Bladder Urothelial Carcinoma (TCGA): 16/169 (9.5%) |
| XAGE1A | CTA / Pass | Lung Adenocarcinoma (validation): 81/105 (77.1%) |
| XAGE1B | CTA / Pass | Lung Adenocarcinoma (validation): 71/105 (67.6%) |
| CT83 | CTA / Pass | Stomach Adenocarcinoma (TCGA): 174/346 (50.3%) |
| CT55 | CTA / Supported | Skin Cuteneous Melanoma (TCGA): 5/99 (5.1%) |
| SAGE1 | CTA / Pass | Lung Squamous Cell Carcinoma (validation): 8/68 (11.8%) |
| CTCFL | CTA / Pass | Ovary Serous Cystadenocarcinoma (TCGA): 187/349 (53.6%) |
| GAGE1 | CTA / Pass | Lung Squamous Cell Carcinoma (validation): 6/68 (8.8%) |
| SSX4 | CTA / Pass | Ovary Serous Cystadenocarcinoma (TCGA): 39/349 (11.2%) |
| SSX5 | CTA / Pass | Lung Adenocarcinoma (validation): 6/105 (5.7%) |

RNA occurrence is not protein detection, peptide presentation, HLA eligibility or a response rate. Bulk tumor signals can originate from tumor-associated normal cells. The full tables include smaller cohorts, thresholds 0.1/1/2/5 pTPM, and separate Medium/High cancer-IHC counts. RNA and IHC subjects are unpaired; the aggregate cancer-IHC source does not supply antibody reliability. [Complete CTA evidence atlas and gene catalogs](https://pirl-unc.github.io/oncoref/audits/cta-atlas-20261024/)

## Why candidates are excluded or left in Neither

| Result | Principal recorded reason | Genes |
|---|---|---|
| Neither | Unresolved antibody evidence | 199 |
| Neither | Cell / immune RNA concerns | 70 |
| Neither | Weak affirmative reproductive evidence | 191 |
| Neither | Other admission / curation requirements | 59 |
| Neither | Insufficient reproductive enrichment | 10 |
| Neither | Coverage / identity gap | 31 |
| Neither | Vital donor zero/high pattern | 1 |
| Unsafe | Additional >30 normal RNA rule | 125 |
| Unsafe | Published normal-cell protein | 3 |
| Unsafe | Accepted positive normal IHC | 845 |
| Unsafe | Bulk RNA above applicable ceiling | 606 |
| Unsafe | Independent safety hold | 1 |

These are mutually exclusive final labels, not counts of every failing observation. A gene may fail several gates. Existing exclusions were preserved: a gene already excluded for normal protein need not be relabelled High RNA when another source is high. Earlier Unsafe assignments were not all re-audited across every expanded source.

| Examples | Why not retained |
|---|---|
| MAGEA3 | Colon undifferentiated cluster 53.1 nTPM; >30 rule |
| ROPN1 | Salivary serous cells 32.9 nTPM; near the hard ceiling but still excluded |
| ACRBP / TEX101 | High normal-cell/immune RNA; source-native maxima are in the follow-up table |
| IL5 / IL13 / PAPPA | Published normal activated immune or vascular protein production |
| MAGEA12 | Independent reviewed safety hold |
| SPAG9 | Mapped candidate with positive normal protein evidence |
| FOXR2 | Strong cancer biology supports interest, but uncertain primary antibodies and normal staining remain unresolved |
| CTAG2 | Variable: donor-level vital RNA fails despite negative curated local IHC |

## Retained genes with higher normal RNA

| Gene | Highest non-permitted nTPM observation | Interpretation |
|---|---|---|
| TEX14 | 29.7 nTPM; breast/macrophages; hpa_single_cell_tissue; vital | Below the hard ceiling; not evidence of absent normal protein |
| MORC1 | 28.5 nTPM; pbmc/platelets; hpa_single_cell_tissue; vital | Below the hard ceiling; not evidence of absent normal protein |
| KIAA1210 | 27.7 nTPM; prostate/t-cells; hpa_single_cell_tissue; vital | Below the hard ceiling; not evidence of absent normal protein |
| KLF17 | 27.4 nTPM; unsited/Oligodendrocyte precursor cells; hpa_single_cell; unresolved | Below the hard ceiling; not evidence of absent normal protein |
| YBX2 | 27 nTPM; colon/distal enterocytes; hpa_single_cell_tissue; protected | Below the hard ceiling; not evidence of absent normal protein |
| CGB3 | 25.9 nTPM; pancreas/mixed cell types; hpa_single_cell_tissue; vital | Below the hard ceiling; not evidence of absent normal protein |
| C10orf62 | 24.8 nTPM; testis/macrophages; hpa_single_cell_tissue; vital | Below the hard ceiling; not evidence of absent normal protein |
| MEIOB | 24.2 nTPM; placenta/mixed immune cells; hpa_single_cell_tissue; vital | Below the hard ceiling; not evidence of absent normal protein |
| TEX15 | 24 nTPM; brain/inhibitory neurons; hpa_single_cell_tissue; vital | Below the hard ceiling; not evidence of absent normal protein |
| TMEM210 | 21.6 nTPM; pbmc/dendritic cells; hpa_single_cell_tissue; vital | Below the hard ceiling; not evidence of absent normal protein |

These are cell-context observations, not comparable bulk-tissue concentrations. The lower tissue-cluster review also requires a comparison with reproductive-cell RNA, so a substantial normal signal can remain when the reproductive-cell maximum is larger. This is an explicit limitation of the current rule, not a conclusion that those cells could safely be targeted.

231 retained genes have a reliable curated reproductive IHC positive control; 168 have neither that control nor additional reliable reproductive protein recorded. The latter can qualify through RNA or cited support, but have a weaker direct protein foundation. The retained-review table exposes these flags per gene.

## Retained literature conflicts: priority cases

All 37 flagged retained genes are listed in the companion table. These flags are independent of the frozen pass/fail result. A conflicting annotation requires checking the original tissue, antibody/probe, isoform, developmental stage and patient context. A shared family record is not independent evidence for every locus.

ABCC11: Published human brain immunoblotting and axonal immunofluorescence conflict with a reproductive-restricted target label. Action: Reassess neural protein specificity before treating this screen pass as a usable target; review should enter the existing normal-protein adjudication process.

Source evidence. [ABCC11 normal-tissue evidence](https://pubmed.ncbi.nlm.nih.gov/16359813/)

MAGEC2: A 2003 study reported cerebellar Purkinje-cell staining with the TC-1 polyclonal anti-HCA587 reagent. Action: Resolve antibody/paralog specificity and modern neural protein evidence before relying on the retained label; the report is not a demonstrated clinical injury.

Source evidence. [MAGEC2 normal-tissue evidence](https://pubmed.ncbi.nlm.nih.gov/12920247/)

YBX2: Published staining includes vascular smooth muscle in pulmonary artery, myocardium and skeletal muscle, plus adrenal cortical cells. Action: Review the original staining and antibody specificity; distinguish vascular cells from cardiomyocytes and reconcile with the current negative profile.

Source evidence. [YBX2 normal-tissue evidence](https://pmc.ncbi.nlm.nih.gov/articles/PMC2361212/)

PATE2: UniProt assigns normal male/female brain expression to isoform 2, while reproductive expression includes both forms. Action: Resolve transcript and protein identities; a gene-wide CTA decision cannot establish an isoform-safe target.

Source evidence. [PATE2 normal-tissue evidence](https://www.uniprot.org/uniprotkb/Q6UY27/entry)

TPTE2: Primary work describes TPIP splice variants in brain and stomach as well as testis. Action: Resolve historical TPIP nomenclature and the normal-tissue isoforms against the exact protein inventory.

Source evidence. [TPTE2 normal-tissue evidence](https://pubmed.ncbi.nlm.nih.gov/11716755/)

SUN5: UniProt explicitly records conflicting reports of widespread versus testis-specific expression. Action: Adjudicate original probes, transcript identities and tissue panels rather than selecting the preferred annotation.

Source evidence. [SUN5 normal-tissue evidence](https://www.uniprot.org/uniprotkb/Q8TC36/entry)

Other flags include ACTL8/pancreas, FAM170A/brain, FKBP6/broad normal expression, PRAME/low normal annotations, MAGEA11/adrenal, CT47 shared-family brain attribution, SSX1/2 thyroid annotations and developmental-versus-adult compartment ambiguity. These findings are not all the same evidentiary strength and should not be collapsed into a single toxicity claim.

## Activated immune expression: what the audit actually shows

The current mapping sends sorted immune populations to blood/Vital. The additional immune RNA triggers are also shared across lineages. They do not estimate cell rarity, depletion tolerance or reversibility. This is conservative and biologically coarse. DICE's activated CD4/CD8 measurements are stimulated normal T cells; activated is a state, not a fixed disposable lineage. [Schmiedel et al., DICE, Cell 2018](https://pubmed.ncbi.nlm.nih.gov/30449622/)

| Lower-threshold immune finding | Nominated genes | Meaning |
|---|---|---|
| Activated DICE trigger | 33 | Includes genes with other immune flags |
| Unstimulated DICE trigger | 99 | Can overlap activated and summary flags |
| Other HPA/Monaco immune trigger | 70 | Source-native summaries; not independent replication |
| Activated DICE triggers only | 7 | No other lower-threshold immune trigger in assessed tables; not absence elsewhere |

| Gene | Current result | Observed activated distribution(s) |
|---|---|---|
| CFAP141 | Unsafe / High RNA | DICE CD4_STIM (activated): 32/88 values ≥5 TPM, max 11.3571; DICE CD8_STIM (activated): 35/88 values ≥5 TPM, max 10.1025 |
| ENTHD1 | Neither / Antibody | DICE CD4_STIM (activated): 11/88 values ≥5 TPM, max 9.93135 |
| IL21 | Unsafe / High RNA | DICE CD4_STIM (activated): 88/88 values ≥5 TPM, max 247.982 |
| IL3 | Unsafe / High RNA | DICE CD4_STIM (activated): 48/88 values ≥5 TPM, max 78.673; DICE CD8_STIM (activated): 24/88 values ≥5 TPM, max 55.1332 |
| IL31 | Neither / Antibody | DICE CD4_STIM (activated): 18/88 values ≥5 TPM, max 16.2381 |
| SLFNL1 | Neither / Cellular | DICE CD4_STIM (activated): 54/88 values ≥5 TPM, max 9.18486; DICE CD8_STIM (activated): 19/88 values ≥5 TPM, max 7.51181 |
| SMIM47 | Neither / Cellular | DICE CD4_STIM (activated): 3/88 values ≥5 TPM, max 6.76639 |

This seven-gene flag is not a rescue list. Three remain Unsafe / High RNA, two have antibody holds, and two are Neither / Cellular. Removing activation-related review would not automatically admit them. All other gates need a complete recomputation.

A fraction such as 3/88 is a fraction of reported expression values meeting a threshold, not a fraction of cells. DICE columns are anonymous per-gene distributions here and cannot be paired across genes as stable donor IDs. Rarity across people, rarity of a cell lineage and rarity of expressing cells within that lineage are different questions.

## Immune depletion: a caveat, not a new permission rule

Clinical experience supports distinguishing immune lineages rather than equating all blood expression with the same consequence. It does not establish that a gene expressed in a depleted lineage is reproductive-restricted. The following interpretations are recommendations for future review, not implemented exemptions.

Eosinophils: benralizumab follow-up in selected asthma patients found no new safety signals during treatment extending up to five years. This supports considering conditional manageability of eosinophil depletion; it does not validate permanent elimination in every tissue, population or treatment context. [Korn et al., benralizumab treatment for up to five years, 2021](https://pubmed.ncbi.nlm.nih.gov/34487870/)

Basophils: a small population with observed reduction in benralizumab studies, including blood and marrow findings. These observations support separate review rather than an unconditional Disposable label. Depletion extent, tissue compartment and duration matter. [Benralizumab randomized allergic-asthma study, 2024](https://publications.ersnet.org/content/erj/64/3/2400512)

B cells: depletion is an accepted therapeutic trade-off in some settings, but hypogammaglobulinemia and infection can require monitoring and immunoglobulin replacement. A clinical cohort links post-CD19-CAR-T antibody deficiency to infection; associated conditioning, prior treatment and other cytopenias complicate attribution. Conditionally manageable is more defensible than harmless or disposable. [CD19 CAR T therapy, hypogammaglobulinemia and infection, clinical cohort](https://pmc.ncbi.nlm.nih.gov/articles/PMC11805655/)

Activated CD4/CD8 cells: ordinary T cells can enter these states during immune responses. A low resting signal cannot establish that depletion during activation is tolerable. The relevant lineage and state distribution need review; the RNA assay alone does not show presentation of a target peptide or actual cell killing. [Schmiedel et al., DICE, Cell 2018](https://pubmed.ncbi.nlm.nih.gov/30449622/)

Regulatory T cells: small numbers do not imply dispensability. Experimental depletion in neonatal and adult mice showed a critical role in preventing autoimmune pathology. That is mechanistic mouse evidence, not a numerical human safety limit. [Kim et al., regulatory T-cell depletion and autoimmunity in mice, 2007](https://www.nature.com/articles/ni1428)

IL5 illustrates a modality distinction: neutralizing a secreted cytokine, depleting IL-5-receptor-bearing eosinophils, and killing IL5-producing cells are different actions. The latter includes normal activated T-helper biology; tolerating eosinophil depletion does not by itself clear an IL5-producing-cell target. [Normal human IL-5-producing T-helper-cell biology](https://pmc.ncbi.nlm.nih.gov/articles/PMC3445433/)

A future revision should record lineage, activation state, evidence of normal protein or peptide, depletion extent/duration, clinical manageability and nonimmune expression separately. This release leaves blood/Vital mapping and all immune thresholds unchanged.

## Neither follow-up and numerical boundaries

CTAG2 - Neither / Variable. When and why do some normal heart specimens contain CTAG2 RNA? Resolve donor variability and obtain paired RNA/protein information from the implicated cardiac context.

FOXR2 - Neither / Antibody. Which primary antibody staining is specific to FOXR2? Resolve discordance between HPA034487 and HPA057358 using independently validated protein evidence.

HORMAD1 - Neither / Cellular. Why is RNA repeatedly detected in monocytes and naive B cells despite reproductive-positive curated IHC? Review gene-specific read mapping and obtain reliable protein evidence in the implicated immune populations.

TDRD1 - Neither / Cellular. Are the few activated-CD4 values and Monaco plasmacytoid-DC summary reproducible? Resolve the borderline immune RNA with source-specific read review and reliable local protein evidence.

BRDT - Neither / Cellular. Does the Monaco-only borderline signal reflect genuine immune expression? Reconcile Monaco CD8/basophil RNA with the much lower DICE values and obtain relevant protein evidence.

DAZL - Neither / Cellular. Is the naive-B-cell RNA translated in normal B cells? Review gene-specific mapping and reliable protein evidence in naive B cells.

CTAG2's pinned heart RNA values are 13.5, 0, 0 and 7.1 nTPM. The applicable absolute ceiling is 12 and the selectivity-derived ceiling is 5.538; both nonzero donors fail at least one rule. Negative IHC specimens are different people. Variable requires reported zero and failing RNA in the same tissue/cohort/version/units and appropriate negative local IHC, not a claim that RNA-zero patients are treatment-eligible. [Current high-normal-RNA policy and candidate follow-up](https://pirl-unc.github.io/oncoref/audits/cta-neither-review-20261024/)

The atlas flags 55 genes near a numerical boundary: retained maxima 24-30, high-RNA exclusions with all audited native-unit maxima <=33, or failed bulk numeric comparisons within 20% of their limits. These descriptive bands are arbitrary review aids, not new selection thresholds. ROPN1 at 32.9 is an example; it remains excluded.

Reconsideration should resolve the particular uncertainty: assay/paralog mapping, ambient contamination, cell identity, isoform, activation state or donor dependence. Familiarity, a trial registration or tumor abundance alone cannot waive a normal finding.

## Clinical efficacy, safety and biological evidence

Nineteen mapped genes received targeted primary-source reviews in the atlas; 30 have curated trial/therapy links. This is not a systematic clinical review of all 2,540 genes. The source-linked biological annotations cover 2,505 genes, with function annotations for 1,873. Missing review is not evidence of no publications, no benefit or no risk.

FDA TECELRA label, revised June 2026: Selected synovial sarcoma patients aged 12+: 60/137 responses (43.8%), median response duration 5.3 months. CRS 72% (2% grade >=3), ICANS 4%, grade >=3 infections 12%. Product-specific evidence; not blanket MAGEA4 safety. [FDA TECELRA label, revised June 2026](https://www.fda.gov/media/180565/download)

IMA203 phase 1, Nature Medicine 2025: Confirmed responses: 11/38 (28.9%); unconfirmed-plus-confirmed: 21/40 (52.5%). Severe CRS: 2/41. No severe neurotoxicity or treatment-related deaths reported. This does not establish safety of other PRAME-directed products or rule out rare adrenal injury. [IMA203 phase 1, Nature Medicine 2025](https://www.nature.com/articles/s41591-025-03650-6)

Robbins et al., NY-ESO-1 TCR follow-up, 2015: Objective responses in 11/18 synovial sarcomas and 11/20 melanomas. HLA- and antigen-selected, uncontrolled cohort; shared NY-ESO-1 protein attribution does not distinguish CTAG1A from CTAG1B. [Robbins et al., NY-ESO-1 TCR follow-up, 2015](https://pubmed.ncbi.nlm.nih.gov/25538264/)

Lete-cel myxoid/round-cell liposarcoma pilot, 2025: Responses: 2/10 and 4/10 across two cohorts. A potentially related fatal cardiac arrest occurred on day 161 amid substantial competing illness/treatment. It does not establish CTAG2-mediated myocardial injury; the independent donor-expression hold remains. [Lete-cel myxoid/round-cell liposarcoma pilot, 2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC12084024/)

## Clinical safety and biological interpretation

Linette et al., cardiovascular toxicity and titin cross-reactivity, 2013: Two fatal cardiac events followed affinity-enhanced MAGE-A3 TCR treatment. Cross-recognition of unrelated cardiac titin explains the hazard; it does not demonstrate native MAGEA3 expression in heart. [Linette et al., cardiovascular toxicity and titin cross-reactivity, 2013](https://pubmed.ncbi.nlm.nih.gov/23770775/)

Morgan et al., anti-MAGE-A3 neurological toxicity, 2013: Five of nine patients had tumor regression; three developed neurological toxicity and two died. The receptor cross-recognized MAGE-A12, with brain RNA evidence; pan-MAGE neuronal staining did not uniquely establish the causal paralog. [Morgan et al., anti-MAGE-A3 neurological toxicity, 2013](https://pubmed.ncbi.nlm.nih.gov/23377668/)

MAGRIT randomized phase 3 trial, 2016: MAGE-A3 vaccination did not improve disease-free survival after NSCLC resection (hazard ratio 1.02; 95% CI 0.89-1.18). Trial nomination is therefore not evidence of benefit. [MAGRIT randomized phase 3 trial, 2016](https://pubmed.ncbi.nlm.nih.gov/27132212/)

CXorf48 as a CML stem-cell immune target, 2017: Human antigen-specific CTLs recognized patient-derived CML stem/progenitor cells. Endogenous immune associations after imatinib cessation support interest, but are not an interventional efficacy or safety trial. [CXorf48 as a CML stem-cell immune target, 2017](https://pubmed.ncbi.nlm.nih.gov/28862699/)

FOXR2 is an epigenetically regulated pan-cancer oncogene, 2022: Activation occurred in 8% of individual tumors across 70% of cancer types; experimental work supports ETS-circuit-driven growth. Testis restriction supports reconsideration, but the antibody conflict still needs resolution. No therapeutic outcome is established here. [FOXR2 is an epigenetically regulated pan-cancer oncogene, 2022](https://pubmed.ncbi.nlm.nih.gov/35802025/)

Mechanism matters: an engineered receptor's recognition of an unrelated normal protein does not establish native CTA expression in that organ. Conversely, a tolerated product does not establish safety for every receptor, epitope, modality, dose or population. Cytokine-release syndrome and conditioning toxicity must remain separate from established target-mediated destruction of a normal compartment.

## Limits of this release and downstream use

Frozen selection versus new concern: the 37 retained literature conflicts do not automatically modify the target API or protein FASTAs. Users must inspect the separate review flags; the raw 399-gene label alone does not resolve them.

Expression breadth: the current selectivity rule uses a reproductive maximum. It has no continuous whole-profile breadth penalty, no fitted RNA uncertainty model and no empirical false-negative IHC calibration.

Measurement coverage: no missing observation is zero. Adult, fetal, stimulated, healthy-volunteer and tumor-adjacent contexts are not interchangeable. Earlier Unsafe genes have incomplete expanded-source audit coverage because the review focused on potential targets and Neither.

Donor scope: original DICE/immune distributions are available, but the explicit reviewed vital-organ zero/high donor rule is currently populated only for CTAG2 heart. Other organs and genes are not thereby donor-negative.

Gene versus protein: gene-level RNA does not allocate expression to isoforms. Exact complete sequences shared by multiple genes use Unsafe > Neither > CTA. Full protein sets are disjoint but can share k-mers; the published exclusion recipe still requires sequence-level subtraction.

Species: anatomy definitions include human, dog and mouse mappings, but this release's HPA/DICE/corneal measurements and numerical limits are human-specific. Human absence cannot establish absence in another species.

Clinical scope: reproductive and Disposable permissions are research assumptions. The broad list does not protect pregnancy or fertility; separate sex/developmental selections are needed. Manageable still means screened, not permitted.

Release contents: this PDF, retained-review.csv, immune-summary.csv, immune-trigger-details.csv, activated-only-review.csv, summary.json and a provenance manifest. The all-gene atlas and v1.18.0 exact-protein assets remain the authoritative frozen evidence/sequence exports. The release manifest hashes the inputs and outputs; no raw source download or threshold change is needed to reproduce this supplement.

Reproduce: python scripts/cta_release_report.py. Use oncoref.target_sets.results() for current gene decisions; earlier cta_screen exports preserve their dated policy.

## Sources and further reading

- [Complete CTA evidence atlas and gene catalogs](https://pirl-unc.github.io/oncoref/audits/cta-atlas-20261024/)
- [Current high-normal-RNA policy and candidate follow-up](https://pirl-unc.github.io/oncoref/audits/cta-neither-review-20261024/)
- [Normal-cell screen: source-specific rules and measurements](https://pirl-unc.github.io/oncoref/audits/cta-normal-cells-20261023/)
- [Schmiedel et al., DICE, Cell 2018](https://pubmed.ncbi.nlm.nih.gov/30449622/)
- [Korn et al., benralizumab treatment for up to five years, 2021](https://pubmed.ncbi.nlm.nih.gov/34487870/)
- [Benralizumab randomized allergic-asthma study, 2024](https://publications.ersnet.org/content/erj/64/3/2400512)
- [CD19 CAR T therapy, hypogammaglobulinemia and infection, clinical cohort](https://pmc.ncbi.nlm.nih.gov/articles/PMC11805655/)
- [Kim et al., regulatory T-cell depletion and autoimmunity in mice, 2007](https://www.nature.com/articles/ni1428)
- [Normal human IL-5-producing T-helper-cell biology](https://pmc.ncbi.nlm.nih.gov/articles/PMC3445433/)
- [ABCC11 normal-tissue evidence](https://pubmed.ncbi.nlm.nih.gov/16359813/)
- [MAGEC2 normal-tissue evidence](https://pubmed.ncbi.nlm.nih.gov/12920247/)
- [YBX2 normal-tissue evidence](https://pmc.ncbi.nlm.nih.gov/articles/PMC2361212/)
- [PATE2 normal-tissue evidence](https://www.uniprot.org/uniprotkb/Q6UY27/entry)
- [TPTE2 normal-tissue evidence](https://pubmed.ncbi.nlm.nih.gov/11716755/)
- [SUN5 normal-tissue evidence](https://www.uniprot.org/uniprotkb/Q8TC36/entry)
- [FDA TECELRA label, revised June 2026](https://www.fda.gov/media/180565/download)
- [IMA203 phase 1, Nature Medicine 2025](https://www.nature.com/articles/s41591-025-03650-6)
- [Robbins et al., NY-ESO-1 TCR follow-up, 2015](https://pubmed.ncbi.nlm.nih.gov/25538264/)
- [Lete-cel myxoid/round-cell liposarcoma pilot, 2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC12084024/)
- [Linette et al., cardiovascular toxicity and titin cross-reactivity, 2013](https://pubmed.ncbi.nlm.nih.gov/23770775/)
- [Morgan et al., anti-MAGE-A3 neurological toxicity, 2013](https://pubmed.ncbi.nlm.nih.gov/23377668/)
- [MAGRIT randomized phase 3 trial, 2016](https://pubmed.ncbi.nlm.nih.gov/27132212/)
- [CXorf48 as a CML stem-cell immune target, 2017](https://pubmed.ncbi.nlm.nih.gov/28862699/)
- [FOXR2 is an epigenetically regulated pan-cancer oncogene, 2022](https://pubmed.ncbi.nlm.nih.gov/35802025/)
