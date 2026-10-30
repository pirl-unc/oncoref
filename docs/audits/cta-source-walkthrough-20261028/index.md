# CTA and oncofetal selection

**698 CTA proteins (741 loci), 18,031 Unsafe, 600 Other and 534 unassigned**, across 19,863 proteins (identical selected reference sequences collapsed). Policy: `oncoref.cta-selection.2026-10-30.2`. This is the current policy; the [1.21.0 release](../cta-release-1.21.0/index.md) remains a frozen comparison.

[PDF slides](cta-selection-1.25.0.pdf) · [Editable slides](slides/cta-evidence-walkthrough-1.25.0-final.pptx) · [TCR-T serious-event review and PDF](../cta-tcr-trials-20261029/index.md) · [Unassigned proteins and diffuseness](../cta-unassigned-20261029/index.md) · [Papers and database versions](references.md) · [Dog evidence](dog.md) · [Shared API](api.md) · [Policy and IHC matching](policy-review.md) · [All decisions](../cta-selection-20261028/protein-group-classifications.csv.gz) · [Changed decisions](../cta-selection-20261028/classification-changes.csv.gz)

[1.25 changes: Reproductive RNA/IHC/cell routes and categorical Other assignments](../cta-neither-20261030/update-1.25.0.md).

## Categories

**Safety, cancer expression and Reproductive
evidence are separate requirements. Only safety evidence can establish Unsafe.**
Cancer needs one qualifying measured sample, without a tumor-versus-normal ratio.
The reproductive-versus-normal maximum veto has been removed. Whole-profile
HPA diffuseness above 10% blocks a completely measured adult/placental identity
route. Thymus is Permissive and excluded; Disposable tissues count. Supported
embryonic evidence remains independent and its concentration is unassessed.
[Rule, changed proteins and sensitivity](../cta-unassigned-20261029/diffuseness-impact.md). [Definitions and logic](../../cta-evidence-axes.md).

![Three separate CTA requirements](../../images/cta-evidence-axes.png)

- **CTA:** qualifying association after the concentration check, one measured cancer sample, required coverage and acceptable normal expression all pass.
- **Unsafe:** above the single RNA ceiling in Vital/Protected, accepted non-Low protein in those compartments, or retained clinical exclusion evidence. This is a research-policy exclusion, not demonstrated toxicity.
- **Other:** complete required normal coverage, resolved evidence and every normal gate passing, but no candidate subtype qualifies. Cancer or reproductive evidence may be absent or missing.
- **Unassigned:** required normal coverage or evidence resolution is missing, or a normal-screen hold remains. In the API, category is null and assessment_status explains why.

The separate **2× Strong boundary is removed**. Every layer evaluates every protein; Other is assigned only after its complete-profile check. The fixed [GTEx v11](https://gtexportal.org/home/downloads/adult-gtex) core requires all 52 bulk contexts for every protein member. DICE and corneal coverage remain optional, but observed counterevidence contributes.

## Other categories

Every normal-pass protein has one categorical assignment: CTA or a named Other category. Safety exclusions remain Unsafe; incomplete normal evidence remains Unassessed. The API exposes `classification` (for example `Other:Low`) as well as `category` and `subtype`.

| Category | Meaning | Proteins |
| --- | --- | --- |
| Other:Diffuse | Complete HPA profile is dispersed (score >10% with a reference ≥2 nTPM), or the active concentration gate fails; no CTA route remains. | 8 |
| Other:Reproductive | Qualifying Reproductive evidence, but no qualifying cancer expression; the normal screen passes. | 2 |
| Other:Inducible | Cited expression in activated nonreproductive normal cells; a localized indirect atlas signal cannot establish Reproductive identity. | 5 |
| Other:Tumor-supported | Two distinct tumors ≥10 native TPM or reviewed reliable tumor protein, with low measured normal/developmental expression; not established tumor restriction. | 66 |
| Other:Disposable | Positive Disposable RNA/protein; measured RNA elsewhere below the descriptive 2-unit floor and no detected protein elsewhere, including low measured developmental RNA. | 2 |
| Other:Low | All measured normal RNA below the descriptive 2-unit floor and no detected normal protein. Cancer expression need not be low. | 254 |
| Other:Normal | Supported nonreproductive normal RNA reaches the descriptive 2-unit floor or reliable protein is detected; no preceding subtype or CTA definition qualifies. | 171 |
| Other:Weak | Reproductive signals are below admission requirements or lack support; no preceding subtype qualifies and no supported appreciable nonreproductive expression is established. | 92 |

The categories are disjoint and applied in documented priority order. Low and Disposable use the descriptive <2 native TPM/nTPM convention (CP10K retains its separate nominal /100 floor), not an assay detection limit. Tumor abundance need not be low. Every Other category still requires a complete normal-screen pass. No unnamed residual category remains.

[Full categorical table and individual changes](../cta-neither-20261030/update-1.25.0.md). The current admission fixes use the full mapped Reproductive bulk panel plus supported reproductive-lineage cells. Generic/immune cells can provide indirect sample-site evidence only with a complete same-source outside panel, every outside cluster <5 nTPM, and no curated normal-induction conflict. This does not identify the producing cell; macrophages may carry passenger RNA. All safety and concentration gates remain independent.

## What the assessment statuses mean

| API status | Meaning | Human proteins |
|---|---|---:|
| `screened` | A CTA, Unsafe or named Other category has been assigned by this policy. It does not mean safe, clinically validated or complete coverage of every source. | 19,329 |
| `incomplete_evidence` | Required normal coverage or evidence resolution is missing. Additional missing cancer/developmental measurements are reported separately and do not themselves prevent Other. | 475 |
| `normal_evidence_hold` | Required coverage is present, but a normal-screen hold remains, including Manageable ceiling failures. | 59 |
| `unassessed` | No species-specific quantitative policy has been applied; used for native dog nominations. | 0 |

An available exclusion can establish Unsafe despite another missing assay. A CTA supported by adult/placental Reproductive evidence need not have embryonic measurements. Other does not require tumor or embryonic coverage when the required normal screen is complete and resolved.

Among the 475 incomplete assessments, gaps overlap:

| Evidence gap | Proteins |
| --- | --- |
| normal panel missing | 453 |
| tumor measurement missing | 34 |
| tumor group profile incomplete | 3 |
| embryonic measurement missing | 86 |
| unresolved evidence | 28 |

[Exact data sources, missing tissue/cell contexts, member identities and cohort coverage](coverage.md) · [Per-protein gap flags](unassigned-evidence-gaps.csv). The previous 71 tumor gaps combined 41 missing measurements and 30 measured groups with incomplete low-expression-profile coverage. The table above reflects the subsequent identity repair. The comparator table below reflects the new developmental and normal-bulk sources; MAGE family members are not admission controls.

## Sources of candidate annotations

These published lists annotate the selection; they do not provide a gene-specific admission bypass. The union is **2,514 mapped proteins / 2,570 loci**. **563** selected proteins have a listed-paper annotation; **135** do not. Source contribution counts overlap. **Mapped proteins means curated nominations linked to the current reference, not newly discovered genes or selected CTAs.** The CTA column reports the subset passing the current funnel.

| citation | mapped_entities | CTA | Unsafe | Other | Unassigned | source_url |
| --- | --- | --- | --- | --- | --- | --- |
| [da Silva et al. (2017)](https://doi.org/10.18632/oncotarget.21715) | 1017 | 473 | 401 | 119 | 24 | [Source](https://doi.org/10.18632/oncotarget.21715) |
| [Wang et al. (2016)](https://doi.org/10.1038/ncomms10499) | 988 | 321 | 640 | 18 | 9 | [Source](https://doi.org/10.1038/ncomms10499) |
| [Chang et al. (2019)](https://doi.org/10.1002/cam4.2223) | 985 | 332 | 639 | 4 | 10 | [Source](https://doi.org/10.1002/cam4.2223) |
| [Gong et al. (2021)](https://doi.org/10.1038/s41467-021-22695-y) | 729 | 399 | 306 | 17 | 7 | [Source](https://doi.org/10.1038/s41467-021-22695-y) |
| [Bruggeman et al. (2018)](https://doi.org/10.1038/s41388-018-0357-2) | 692 | 24 | 663 | 0 | 5 | [Source](https://doi.org/10.1038/s41388-018-0357-2) |
| [Jamin et al. (2021)](https://doi.org/10.1002/1878-0261.12900) | 469 | 39 | 427 | 1 | 2 | [Source](https://doi.org/10.1002/1878-0261.12900) |
| [Loriot et al. (2025)](https://doi.org/10.1371/journal.pgen.1011734) | 202 | 99 | 101 | 0 | 2 | [Source](https://doi.org/10.1371/journal.pgen.1011734) |
| [Carter et al. (2023)](https://doi.org/10.1136/jitc-2023-007935) | 103 | 65 | 37 | 0 | 1 | [Source](https://doi.org/10.1136/jitc-2023-007935) |
| [Hirbo et al. (2015)](https://doi.org/10.1371/journal.pone.0144155) | 78 | 13 | 64 | 0 | 1 | [Source](https://doi.org/10.1371/journal.pone.0144155) |
| [Seager et al. (2024)](https://doi.org/10.1186/s12967-024-04918-0) | 15 | 12 | 2 | 0 | 1 | [Source](https://doi.org/10.1186/s12967-024-04918-0) |
| [Sasaroli et al. (2011)](https://doi.org/10.4161/cbt.12.3.16260) | 13 | 0 | 13 | 0 | 0 | [Source](https://doi.org/10.4161/cbt.12.3.16260) |
| [Garg et al. (2007)](https://doi.org/10.1158/1078-0432.CCR-06-2340) | 5 | 3 | 2 | 0 | 0 | [Source](https://doi.org/10.1158/1078-0432.CCR-06-2340) |

### Broader nominations and a compact primary-source list

[Sasaroli et al. (2011)](https://pmc.ncbi.nlm.nih.gov/articles/PMC3230481/) contributes its 13 detailed tumor-vascular candidates, all mapped. [Hirbo et al. (2015)](https://doi.org/10.1371/journal.pone.0144155) contributes 86 placenta-enriched nominations: 82 mapped identities and four unresolved names. Its 30 fast-evolving genes are a nested subset, not 30 additional genes. Noncoding mapped rows remain in the intake audit but do not enter protein-group counts. These are vascular and placental nomination scopes, not proof of CTA restriction.

EGFL6 now occurs in Sasaroli, Hirbo and Bai. The main source table and plots omit fully redundant single-digit sets: Bai, Boël and Bradley. Their memberships and citations remain in the downloadable provenance. Garg stays because it uniquely contributes SPAG9 in the current registry. A source is retired from the main view only if the retained sources cover its entire mapped set; the union is unchanged by this presentation choice.

| citation | mapped_entities | covered_by | source_url |
| --- | --- | --- | --- |
| [Bradley et al. (2020)](https://doi.org/10.1038/s41467-020-19141-w) | 9 | ('10.1002/1878-0261.12900', '10.1002/cam4.2223', '10.1038/ncomms10499', '10.1038/s41388-018-0357-2', '10.1038/s41467-021-22695-y', '10.1371/journal.pone.0144155', '10.18632/oncotarget.21715', '10.4161/cbt.12.3.16260') | [Source](https://doi.org/10.1038/s41467-020-19141-w) |
| [Boël et al. (1995)](https://doi.org/10.1016/s1074-7613(95)80053-0) | 2 | ('10.1002/1878-0261.12900', '10.1002/cam4.2223', '10.1038/ncomms10499', '10.1038/s41388-018-0357-2', '10.1038/s41467-021-22695-y', '10.1136/jitc-2023-007935', '10.1186/s12967-024-04918-0', '10.1371/journal.pgen.1011734', '10.18632/oncotarget.21715') | [Source](https://doi.org/10.1016/s1074-7613(95)80053-0) |
| [Bai et al. (2016)](https://doi.org/10.1158/0008-5472.CAN-16-0225) | 1 | ('10.1371/journal.pone.0144155', '10.4161/cbt.12.3.16260') | [Source](https://doi.org/10.1158/0008-5472.CAN-16-0225) |

Bai's 538 EGFL6-correlated genes are not a CTA list. Garg recognizes SPA17, SYCP1, NY-ESO-1 and SSX4 alongside SPAG9; Boël recognizes MAGE-A1/A3 alongside the unresolved BAGE label. BAGE's textual overlap with Seager remains outside mapped overlap counts. [Original targeted sets and overlaps](targeted-paper-overlaps.csv) · [Evidence roles](publication-context.csv) · [Every paper's extraction scope](paper-source-scope.csv) · [New source rows, aliases and unresolved names](publication-intake/mapped-source-rows.csv).

The larger papers already contribute their registered gene sets, including Wang's 1,019 CT-pattern genes, Bruggeman's 756 germ-cell genes, da Silva's 1,103 testis-biased and nested 745 CT genes, Chang's 1,036 TGCT genes, Carter's 103 final CT genes, and Loriot's 146 strict plus 134 preferential annotations. Input rows, gene loci and collapsed protein groups differ; these numbers must not be summed. Gong contributes reproductive/placental expression nominations rather than a claimed CTA list. This audit adds the declared vascular and placental sets; it does not claim a fresh exhaustive full-text review of every large source.

EGFL6 is excluded by lung bulk RNA (24.8 nTPM), lung smooth-muscle RNA (363.3 nTPM) and corneal evidence. SPAG9 has accepted bronchial protein staining and broad RNA counterevidence, including spinal cord (66.3 nTPM) and brain oligodendrocytes (287.1 nTPM). These are current quantitative-policy exclusions, irrespective of the original paper's terminology.

![Published source contributions](figures/published-list-contributions.png)
![Published source overlaps](figures/published-list-overlaps.png)

[All pairwise overlaps](human-overlaps.csv.gz) · [Memberships and gene IDs](human-memberships.csv.gz)

## Reproductive evidence and its branches

**Reproductive** includes adult reproductive, placental, extraembryonic, embryonic and fetal tissues. Branch labels describe the measured anatomical/stage evidence; “embryonic-only evidence” does not imply absence from unmeasured branches. [Branch definitions and measurements](../cta-unassigned-20261029/index.md#reproductive-is-the-umbrella). The association set uses whole-reference RNA and reliable reproductive IHC. Reproductive association requires HPA ≥5 nTPM **or mapped positive reproductive IHC with Approved, Supported or Enhanced reliability for every member, without an antibody conflict**. Protein evidence does not bypass cancer, normal-safety or concentration checks. The developmental route requires a supported source-specific stage/lineage mean ≥5 TPM in Tyser or Petropoulos, ≥20 cells, ≥5 positive cells and ≥10% positive at ≥1 TPM. These establish presence in a qualifying context. This presence check does not compare with normal expression; normal safety and profile diffuseness answer separate questions. Later fetal and child/adult developmental atlases are not yet integrated.

| route | entrants | CTA |
| --- | --- | --- |
| Adult/placental evidence only | 4685 | 467 |
| Embryonic-only evidence | 459 | 92 |
| Both measured routes | 12546 | 139 |

[Tyser et al. (2021)](https://doi.org/10.1038/s41586-021-04158-y), [HPA v23 downloads](https://v23.proteinatlas.org/about/download).

## Measurement sources and filtering rules

| name | evidence_category | native_unit | status | source_url |
| --- | --- | --- | --- | --- |
| [Published human normal immune induction contexts](https://pmc.ncbi.nlm.nih.gov/articles/PMC6457222/) | reproductive_expression | qualitative reported RNA induction | used_in_selection | [Source](https://pmc.ncbi.nlm.nih.gov/articles/PMC6457222/) |
| [Published human normal immune induction contexts](https://pmc.ncbi.nlm.nih.gov/articles/PMC6457222/) | normal_cell_rna | qualitative reported RNA induction | used_in_selection | [Source](https://pmc.ncbi.nlm.nih.gov/articles/PMC6457222/) |
| [HPA v23 / GTEx consensus](https://v23.proteinatlas.org/about/download) | normal_bulk_rna | nTPM | used_in_selection | [Source](https://v23.proteinatlas.org/about/download) |
| [HPA v23 / GTEx consensus](https://v23.proteinatlas.org/about/download) | reproductive_expression | nTPM | used_in_selection | [Source](https://v23.proteinatlas.org/about/download) |
| [Tyser human gastrulation atlas](https://doi.org/10.1038/s41586-021-04158-y) | developmental_expression | TPM | used_in_selection | [Source](https://doi.org/10.1038/s41586-021-04158-y) |
| [Petropoulos preimplantation atlas](https://www.ebi.ac.uk/gxa/sc/experiments/E-MTAB-3929/Results) | developmental_expression | TPM | used_in_selection | [Source](https://www.ebi.ac.uk/gxa/sc/experiments/E-MTAB-3929/Results) |
| [HPA v23 human brain, 193 subregions](https://v23.proteinatlas.org/about/download) | normal_bulk_rna | nTPM | used_in_selection | [Source](https://v23.proteinatlas.org/download/rna_brain_hpa.tsv.zip) |
| [Direct GTEx v11 tissue medians](https://gtexportal.org/home/downloads/adult-gtex) | normal_bulk_rna | TPM | used_in_selection | [Source](https://gtexportal.org/home/downloads/adult-gtex) |
| [OncoRef source-mapped tumor union](https://github.com/pirl-unc/oncoref/tree/v1.21.0/docs/audits/cta-expanded-funnel-20261027/tumor-union) | tumor_expression | TPM | used_in_selection | [Source](https://github.com/pirl-unc/oncoref/tree/v1.21.0/docs/audits/cta-expanded-funnel-20261027/tumor-union) |
| [HPA v23 cancer sample measurements](https://v23.proteinatlas.org/about/download) | tumor_expression | pTPM | used_in_selection | [Source](https://v23.proteinatlas.org/about/download) |
| [HPA v23 tissue/cell clusters](https://v23.proteinatlas.org/about/download) | normal_cell_rna | nTPM | used_in_selection | [Source](https://v23.proteinatlas.org/about/download) |
| [HPA v23 tissue/cell clusters](https://v23.proteinatlas.org/about/download) | reproductive_expression | nTPM | used_in_selection | [Source](https://v23.proteinatlas.org/about/download) |
| [HPA v23 immune summaries](https://v23.proteinatlas.org/about/download) | immune_rna | nTPM | used_in_selection | [Source](https://v23.proteinatlas.org/about/download) |
| [Monaco immune summaries via HPA v23](https://v23.proteinatlas.org/about/download) | immune_rna | TPM | used_in_selection | [Source](https://v23.proteinatlas.org/about/download) |
| [DICE immune distributions](https://dice-database.org/) | immune_rna | TPM | used_in_selection | [Source](https://dice-database.org/) |
| [Collin adult corneal atlas](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE155683) | normal_cell_rna | CP10K | used_in_selection | [Source](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE155683) |
| [HPA v23 normal-tissue IHC](https://v23.proteinatlas.org/about/assays%2Bannotation) | normal_protein | HPA categorical level and reliability | used_in_selection | [Source](https://v23.proteinatlas.org/about/assays%2Bannotation) |
| [HPA v23 normal-tissue IHC](https://v23.proteinatlas.org/about/assays%2Bannotation) | reproductive_expression | HPA categorical level and reliability | used_in_selection | [Source](https://v23.proteinatlas.org/about/assays%2Bannotation) |
| [Curated external staining observations](https://github.com/pirl-unc/oncoref/tree/v1.21.0/docs/audits/external-ihc-20261027) | normal_protein | native categorical findings | used_in_selection | [Source](https://github.com/pirl-unc/oncoref/tree/v1.21.0/docs/audits/external-ihc-20261027) |
| [Reviewed normal protein and clinical exclusions](https://github.com/pirl-unc/oncoref/blob/v1.21.0/oncoref/data/cta-policy-findings.csv) | normal_protein | categorical reviewed findings | used_in_selection | [Source](https://github.com/pirl-unc/oncoref/blob/v1.21.0/oncoref/data/cta-policy-findings.csv) |
| [Reviewed normal protein and clinical exclusions](https://github.com/pirl-unc/oncoref/blob/v1.21.0/oncoref/data/cta-policy-findings.csv) | clinical_safety | categorical reviewed findings | used_in_selection | [Source](https://github.com/pirl-unc/oncoref/blob/v1.21.0/oncoref/data/cta-policy-findings.csv) |
| [Letetresgene Autoleucel in Advanced/Metastatic Myxoid/Round Cell Liposarcoma](https://doi.org/10.1200/JCO-24-01466) | clinical_safety | reviewed clinical attribution | used_in_selection | [Source](https://doi.org/10.1200/JCO-24-01466) |
| [Autologous T cell therapy for PRAME+ advanced solid tumors in HLA-A*02+ patients: a phase 1 trial](https://doi.org/10.1038/s41591-025-03650-6) | clinical_safety | reviewed clinical attribution | used_in_selection | [Source](https://doi.org/10.1038/s41591-025-03650-6) |
| [Autologous T cell therapy for MAGE-A4+ solid cancers in HLA-A*02+ patients: a phase 1 trial](https://doi.org/10.1038/s41591-022-02128-z) | clinical_safety | reviewed clinical attribution | used_in_selection | [Source](https://doi.org/10.1038/s41591-022-02128-z) |
| [SPEARHEAD-1: a single-arm phase 2 trial of afamitresgene autoleucel in advanced synovial sarcoma and myxoid/round cell liposarcoma](https://doi.org/10.1016/S0140-6736(24)00319-2) | clinical_safety | reviewed clinical attribution | used_in_selection | [Source](https://doi.org/10.1016/S0140-6736(24)00319-2) |
| [Identification of a Titin-Derived HLA-A1-Presented Peptide as a Cross-Reactive Target for Engineered MAGE A3-Directed T Cells](https://doi.org/10.1126/scitranslmed.3006034) | clinical_safety | reviewed clinical attribution | used_in_selection | [Source](https://doi.org/10.1126/scitranslmed.3006034) |
| [Guillain-Barre syndrome observed with adoptive transfer of lymphocytes genetically engineered with an NY-ESO-1 reactive T-cell receptor](https://doi.org/10.1186/s40425-019-0759-x) | clinical_safety | reviewed clinical attribution | used_in_selection | [Source](https://doi.org/10.1186/s40425-019-0759-x) |
| [NCT03709706 protocol: NY-ESO-1/LAGE-1A shared peptide and c259 receptor](https://cdn.clinicaltrials.gov/large-docs/06/NCT03709706/Prot_000.pdf) | clinical_safety | reviewed clinical attribution | used_in_selection | [Source](https://cdn.clinicaltrials.gov/large-docs/06/NCT03709706/Prot_000.pdf) |
| [NCT03132922 protocol: MAGE-A4 c1032 target peptide](https://cdn.clinicaltrials.gov/large-docs/22/NCT03132922/Prot_000.pdf) | clinical_safety | reviewed clinical attribution | used_in_selection | [Source](https://cdn.clinicaltrials.gov/large-docs/22/NCT03132922/Prot_000.pdf) |
| [First-in-human dose escalation trial to evaluate the clinical safety and efficacy of an anti-MAGEA1 autologous TCR-transgenic T cell therapy in relapsed and refractory solid tumors](https://doi.org/10.1136/jitc-2023-008668) | clinical_safety | reviewed clinical attribution | used_in_selection | [Source](https://doi.org/10.1136/jitc-2023-008668) |
| [Phase I clinical trial evaluating the safety and efficacy of ADP-A2M10 SPEAR T cells in patients with MAGE-A10+ advanced non-small cell lung cancer](https://doi.org/10.1136/jitc-2021-003581) | clinical_safety | reviewed clinical attribution | used_in_selection | [Source](https://doi.org/10.1136/jitc-2021-003581) |
| [NY-ESO-1-specific redirected T cells with endogenous TCR knockdown mediate tumor response and cytokine release syndrome](https://doi.org/10.1136/jitc-2021-003811) | clinical_safety | reviewed clinical attribution | used_in_selection | [Source](https://doi.org/10.1136/jitc-2021-003811) |
| [Cancer regression and neurological toxicity following anti-MAGE-A3 TCR gene therapy](https://pubmed.ncbi.nlm.nih.gov/23377668/) | clinical_safety | reviewed clinical attribution | used_in_selection | [Source](https://pubmed.ncbi.nlm.nih.gov/23377668/) |
| [Gong et al., Nature Communications (2021)](https://doi.org/10.1038/s41467-021-22695-y) | nomination | gene-list membership | annotation | [Source](https://doi.org/10.1038/s41467-021-22695-y) |
| [Gong et al., Nature Communications (2021)](https://doi.org/10.1038/s41467-021-22695-y) | nomination | gene-list membership | annotation | [Source](https://doi.org/10.1038/s41467-021-22695-y) |
| [Bradley et al., Nature Communications (2020)](https://doi.org/10.1038/s41467-020-19141-w) | nomination | gene-list membership | annotation | [Source](https://doi.org/10.1038/s41467-020-19141-w) |
| [Wang 2016](https://doi.org/10.1038/ncomms10499) | nomination | gene-list membership | annotation | [Source](https://doi.org/10.1038/ncomms10499) |
| [Bruggeman 2018](https://doi.org/10.1038/s41388-018-0357-2) | nomination | gene-list membership | annotation | [Source](https://doi.org/10.1038/s41388-018-0357-2) |
| [da Silva 2017 · testis-biased](https://doi.org/10.18632/oncotarget.21715) | nomination | gene-list membership | annotation | [Source](https://doi.org/10.18632/oncotarget.21715) |
| [da Silva 2017 · CT](https://doi.org/10.18632/oncotarget.21715) | nomination | gene-list membership | annotation | [Source](https://doi.org/10.18632/oncotarget.21715) |
| [Jamin 2021](https://doi.org/10.1002/1878-0261.12900) | nomination | gene-list membership | annotation | [Source](https://doi.org/10.1002/1878-0261.12900) |
| [Jamin 2021 · core](https://doi.org/10.1002/1878-0261.12900) | nomination | gene-list membership | annotation | [Source](https://doi.org/10.1002/1878-0261.12900) |
| [Chang 2019 · TGCT](https://doi.org/10.1002/cam4.2223) | nomination | gene-list membership | annotation | [Source](https://doi.org/10.1002/cam4.2223) |
| [Carter 2023](https://doi.org/10.1136/jitc-2023-007935) | nomination | gene-list membership | annotation | [Source](https://doi.org/10.1136/jitc-2023-007935) |
| [Seager 2024 · panel](https://doi.org/10.1186/s12967-024-04918-0) | nomination | gene-list membership | annotation | [Source](https://doi.org/10.1186/s12967-024-04918-0) |
| [Loriot 2025 · strict](https://doi.org/10.1371/journal.pgen.1011734) | nomination | gene-list membership | annotation | [Source](https://doi.org/10.1371/journal.pgen.1011734) |
| [Loriot 2025 · preferential](https://doi.org/10.1371/journal.pgen.1011734) | nomination | gene-list membership | annotation | [Source](https://doi.org/10.1371/journal.pgen.1011734) |
| [Gong 2021 · reproductive coding](https://doi.org/10.1038/s41467-021-22695-y) | nomination | gene-list membership | annotation | [Source](https://doi.org/10.1038/s41467-021-22695-y) |
| [Gong 2021 · reproductive noncoding](https://doi.org/10.1038/s41467-021-22695-y) | nomination | gene-list membership | annotation | [Source](https://doi.org/10.1038/s41467-021-22695-y) |
| [da Silva 2017 · tumor proteomics](https://doi.org/10.18632/oncotarget.21715) | nomination | gene-list membership | annotation | [Source](https://doi.org/10.18632/oncotarget.21715) |
| [Bai 2016](https://doi.org/10.1158/0008-5472.CAN-16-0225) | nomination | gene-list membership | annotation | [Source](https://doi.org/10.1158/0008-5472.CAN-16-0225) |
| [Boël et al., Immunity (1995)](https://doi.org/10.1016/s1074-7613(95)80053-0) | nomination | gene-list membership | annotation | [Source](https://doi.org/10.1016/s1074-7613(95)80053-0) |
| [Garg et al., Clinical Cancer Research (2007)](https://doi.org/10.1158/1078-0432.CCR-06-2340) | nomination | gene-list membership | annotation | [Source](https://doi.org/10.1158/1078-0432.CCR-06-2340) |
| [Sasaroli et al. (2011)](https://doi.org/10.4161/cbt.12.3.16260) | nomination | gene-list membership | annotation | [Source](https://doi.org/10.4161/cbt.12.3.16260) |
| [Hirbo et al. (2015)](https://doi.org/10.1371/journal.pone.0144155) | nomination | gene-list membership | annotation | [Source](https://doi.org/10.1371/journal.pone.0144155) |

| stage | rule |
| --- | --- |
| Reproductive or developmental expression | Full Reproductive bulk ≥5 nTPM OR supported reproductive-cell RNA ≥5 nTPM (≥20 reads; indirect sites additionally need complete outside-cluster coverage and every outside value <5, with no curated normal-induction conflict) OR reliable mapped positive reproductive IHC for every member without an antibody conflict OR supported developmental stage/lineage mean ≥5 TPM; ≥20 cells, ≥5 positive and ≥10% positive at ≥1 TPM. Presence evidence never overrides normal exclusions. |
| Cancer expression | At least one measured compatible tumor sample ≥2 in native TPM or HPA pTPM OR reviewed reliable tumor protein; do not add maxima from different samples. |
| HPA normal protein | Reliable non-Low positives in Vital/Protected supply Unsafe evidence; Manageable positives block CTA only; uncertain evidence remains unresolved. Low is flagged. |
| Required RNA coverage | Complete GTEx v11 bulk core (52 native contexts), required nonreproductive bulk anatomy across compatible sources, and distinct HPA tissue-cell/immune panels; every protein member measured. DICE and cornea coverage is optional. |
| Reproductive concentration | HPA v23 average normal/(normal+Reproductive maximum) above 10% blocks the adult/placental route when its reference is positive and the comparison is complete. Use 38 contexts/29 anatomy units, including Disposable and excluding Permissive thymus. Supported embryonic evidence remains independent; its contrast is unassessed. This is not safety evidence. |
| Normal bulk | One inclusive tissue/IHC ceiling: above it in Vital/Protected supplies Unsafe evidence; above it in Manageable blocks CTA. No second RNA boundary. |
| HPA tissue clusters | RNA requires ≥20 reads and uses the single ceiling. IHC credit requires an exact or declared compartment match; unrelated cells in the same tissue cannot inherit it. |
| HPA immune | Protected-tier native nTPM limits without matched negative-IHC credit; reproductive/placental location remains neutral. |
| Monaco immune | Protected-tier native TPM limits without matched negative-IHC credit. |
| DICE immune | Third-highest reported TPM value; native protected-tier limits. Unpaired group distributions retain bounds. Missing coverage alone does not exclude. |
| Adult cornea | Recurrent native CP10K: ≥5 positive cells and ≥20 counts in ≥2 sample groups. Nominal RNA limits /100; not a TPM conversion. Missing coverage alone does not exclude. |
| Accepted normal-protein evidence | Accepted positive protein findings exclude unless extracted Low intensity retains a flag; unquantified does not mean negative. |
| Clinical safety exclusion | Retain reviewed product-specific clinical injury counterevidence and causal-attribution limits; no admission shortcut from clinical use. |
| External validated normal protein | Use categorical/quantitative staining only within verified reagent, tissue and cell scope; no inherited HPA reliability label. |
| CTA and Other assignments | Other is the complete-normal-panel, resolved-evidence fallback when every normal gate passes but no candidate subtype qualifies. Tumor or developmental missingness alone does not block Other. Other:Tumor-supported requires low measured normal/developmental expression and two distinct tumor samples ≥10 TPM or reviewed reliable tumor protein; it is not target admission. Manageable holds remain unassigned, never cleared by breadth. |

![RNA cutoffs](figures/active-rna-cutoffs.png)

The negative-IHC columns require a usable profile control, complete group-member coverage and a supported tissue/cell match. Same-tissue IHC for an unrelated cell cannot increase a ceiling. Values retain their native units; pTPM denotes protein-coding **RNA** TPM. Tumor eligibility uses ≥2 in at least one compatible OncoRef TPM or HPA pTPM sample, without summing maxima from different samples.

## Filtering funnel

| layer | new_candidate_removals | new_Unsafe | CTA_eligible_so_far | Unsafe | Other | Unassigned |
| --- | --- | --- | --- | --- | --- | --- |
| Reproductive / embryonic association | 2173 | 0 | 17690 | 0 | 0 | 2173 |
| Cancer expression | 16 | 0 | 17674 | 0 | 0 | 2189 |
| HPA normal protein | 12044 | 10094 | 5630 | 10094 | 0 | 4139 |
| Required RNA coverage | 208 | 0 | 5422 | 10094 | 0 | 4347 |
| Reproductive concentration | 997 | 0 | 4425 | 10094 | 0 | 5344 |
| Normal bulk | 3444 | 7406 | 981 | 17500 | 0 | 1382 |
| HPA tissue clusters | 238 | 414 | 743 | 17914 | 0 | 1206 |
| HPA immune | 12 | 30 | 731 | 17944 | 0 | 1188 |
| Monaco immune | 7 | 17 | 724 | 17961 | 0 | 1178 |
| DICE immune | 19 | 55 | 705 | 18016 | 0 | 1142 |
| Adult cornea | 3 | 12 | 702 | 18028 | 0 | 1133 |
| Accepted normal-protein evidence | 3 | 2 | 699 | 18030 | 0 | 1134 |
| Clinical safety exclusion | 1 | 1 | 698 | 18031 | 0 | 1134 |
| External validated normal protein | 0 | 0 | 698 | 18031 | 0 | 1134 |
| CTA and Other assignments | 0 | 0 | 698 | 18031 | 600 | 534 |

![Human funnel](figures/human-funnel.png)

[Independent layer effects](../cta-selection-20261028/independent-layer-effects.csv) distinguish overlapping evidence from order-dependent first removals. All counts use proteins (identical selected reference sequences collapsed); gene-locus counts are explicitly labelled.

## Selected CTA expression

These plots compare adult/placental nTPM against the highest measured embryonic lineage mean TPM. **657 of 698 selected proteins (94%) have both coordinates; 41 (6%) lack an embryonic coordinate.** The plotted embryonic coordinate is the maximum across separate stage/lineage means from Tyser (gastrulation) and Petropoulos (days 3–7). It is not an average across studies. Petropoulos restores measured rows for PRAME, MAGE-A4, CT83, CTAG1A/B and XAGE1A/B. Missing rows remain unavailable; neither atlas establishes expression throughout later fetal development.

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

PRAME, CT83, CTAG1A/B and XAGE1B have zero reported reads across these 1,195 cells; MAGE-A4 has 14 reads in one cell. TBXT has 31,149 reported reads in 226 cells. XAGE1A has no separate column, so complete XAGE1A/B group coverage is not established in this raw-count source; the separate Petropoulos TPM matrix does measure both loci. These are observed nondetections or sparse signals in one early embryo, **not absence throughout fetal development**. The exact per-gene reason for exclusion from the atlas remains unavailable; the authors' counts support low detection but come from a separate quantification. They are not converted to TPM or inserted into the selection funnel. [Author data](https://github.com/ScialdoneLab/human-gastrula-shiny), [count evidence and gene IDs](embryonic-comparator-coverage.csv), [source hashes](embryonic-coverage-receipt.json).

Each figure now shows the 657 measured groups above and the 41 missing-embryonic groups in a separate categorical band below. All selected proteins are represented; the band has no embryonic-expression value. Gray points mean the healthy-expression color is unavailable. Overview labels highlight the widely studied antigens, including their gene-group names. Three companion panels label every selected protein ≥5 on both axes, with the same axes and one third of the names per panel. The latter is a descriptive expression criterion; developmental admission also requires cell support.

![Selected CTA reproductive versus embryonic expression](figures/selected-reproductive-embryonic.png)
![Selected CTA colored by healthy TPM](figures/selected-reproductive-embryonic-immune.png)
![Selected CTA colored by adult corneal CP10K](figures/selected-reproductive-embryonic-cornea.png)

Color uses the maximum supported value across non-disposable healthy compartments within each unit: HPA bulk/brain/cell/immune nTPM, direct GTEx/Monaco/DICE TPM, or corneal CP10K. No maximum mixes units. Reproductive/developmental locations are excluded from the color scale. Missing embryonic measurements are not plotted at zero. [Every selected protein's plot data](selected-expression-scatter.csv.gz).

### Groups high in both expression compartments

![Selected groups ≥5 on both axes, labels 1 of 3](figures/selected-high-both-1.png)
![Selected groups ≥5 on both axes, labels 2 of 3](figures/selected-high-both-2.png)
![Selected groups ≥5 on both axes, labels 3 of 3](figures/selected-high-both-3.png)

## Widely studied antigens

| protein_group_label | category | assessment_status | reproductive_ntpm | developmental_all_lineage_max_tpm | adult_bulk_max_ntpm | exclusion_reasons |
| --- | --- | --- | --- | --- | --- | --- |
| PAGE1 | CTA | screened | 16.7 | 26.69 | 0 |  |
| SMC1B | CTA | screened | 30.2 | 1.417 | 0.7 |  |
| CTAG2 | CTA | screened | 10.1 | 7.151 | 5.1 |  |
| MAGEA4 | CTA | screened | 21.5 | 203.6 | 0.6 |  |
| TBXT | CTA | screened | 0.2 | 245.1 | 2.5 |  |
| PRAME | CTA | screened | 102.2 | 22.89 | 3.6 |  |
| MAGEA1 | CTA | screened | 4.4 | 6.315 | 0 |  |
| CT83 | CTA | screened | 62.2 | 4.286 | 2.6 |  |
| MAGEA12 | Unsafe | screened | 2.7 | 102.1 | 0.5 | Absolute tissue/IHC RNA limits, qualifying normal protein and clinical safety evidence; required coverage is separate from absence of a detected exclusion.  /  Retain reviewed product-specific clinical injury counterevidence and causal-attribution limits; no admission shortcut from clinical use. |
| MAGEA3 | Unassigned | normal_evidence_hold | 7.8 | 71.75 | 0.2 | Absolute tissue/IHC RNA limits, qualifying normal protein and clinical safety evidence; required coverage is separate from absence of a detected exclusion.  /  colon: undifferentiated cells; 53.1 nTPM (ceiling 16) |
| NY-ESO-1 | CTA | screened | 15.5 | 0.9473 | 0.5 |  |
| SSX2/B | CTA | screened | 32.7 | 80.25 | 2.8 |  |
| XAGE1A/B | CTA | screened | 98.8 | 2.237 | 1.5 |  |

Each comparison has its own full-width figure. The embryonic comparison retains missing-measurement antigens in a separate band.

![Widely studied antigens: reproductive expression](figures/widely-studied-antigens-reproductive.png)
![Widely studied antigens: embryonic expression](figures/widely-studied-antigens-embryonic.png)

[Detailed tissue/cell/value explanations](human-explanations.csv.gz), [all comparator contexts](../cta-selection-20261028/widely-studied-contexts.csv.gz), and [every previously annotated CTA and its current category](../cta-selection-20261028/known-cta-decisions.csv.gz). Previous CTA labels are comparison annotations only. Historical-cohort counts, where available: {'CTA': 165, 'Unsafe': 48}.

## Serious adverse events from CTA-directed TCR-T

**Treatment toxicity and antigen-specific injury are different findings.** SAE seriousness (for example hospitalization or death) is distinct from CTCAE severity. This targeted review includes important severe events and ungraded serious cases; it is not an exhaustive SAE database or a pooled incidence estimate.

| Protein origins (current screen) | Receptor / product | Important events | RNA credit |
|---|---|---|---|
| CTAG1A: CTA; CTAG1B: CTA; CTAG2: CTA | [Letetresgene autoleucel / NY-ESO-1 c259](https://doi.org/10.1200/JCO-24-01466) | Fatal cardiac arrest, day 161. Guillain-Barre / inflammatory demyelinating neuropathy. | One step, capped |
| PRAME: CTA | [IMA203](https://doi.org/10.1038/s41591-025-03650-6) | Grade 4 transaminitis during a grade 3 CRS episode. | One step, capped |
| MAGEA4: CTA | [Afamitresgene autoleucel / c1032](https://doi.org/10.1038/s41591-022-02128-z) | Fatal aplastic anemia, day 55. Fatal ischemic cerebrovascular accident, day 17. | One step, capped |
| MAGEA3: Unassigned | [MAGE-A3 affinity-enhanced a3a receptor](https://doi.org/10.1126/scitranslmed.3006034) | Fatal cardiac injury with titin cross-recognition. | None |
| MAGEA1: CTA | [IMA202](https://doi.org/10.1136/jitc-2023-008668) | Fatal dyspnea. | One step, capped |
| MAGEA10: CTA | [ADP-A2M10 / MAGE-A10 c796](https://doi.org/10.1136/jitc-2021-003581) | Grade 4 CRS. Fatal marrow aplasia, sepsis and respiratory distress after repeat infusion. | One step, capped |
| CTAG1A: CTA; CTAG1B: CTA; CTAG2: CTA | [TBI-1301](https://doi.org/10.1136/jitc-2021-003811) | Steroid-responsive lung injury with TCR-T infiltration. | None |
| MAGEA3: Unassigned; MAGEA9: CTA; MAGEA12: Unsafe; MAGEA9B: CTA | [NCI HLA-A*02-restricted MAGE-A3/A9/A12-reactive TCR](https://pubmed.ncbi.nlm.nih.gov/23377668/) | Fatal neurological injury; MAGE-A12 cross-recognition implicated. | None |

The MAGE-A3 **a3a** cardiac cases implicated **titin (TTN)**, an unintended target; the separate NCI MAGE-A3/A9 receptor neurological cases implicated the cross-recognized **MAGE-A12** peptide. These findings do not establish that MAGE-A3 itself was present in heart or brain. TBI-1301's grade 3 lung injury retains an unresolved-specificity hold for that product. c259's late cardiac death and Guillain-Barre cases remain caveats without an established peptide-specific mechanism.

[Full attribution, epitope/HLA mappings, cohort denominators and limitations](../cta-tcr-trials-20261029/index.md). A favorable reviewed product supplies one capped RNA step, never protein or clinical-injury clearance. The exact tested peptide matters: NY-ESO-1 and **CTAG2/LAGE-1A** share SLLMWITQC; this does not clear LAGE-1B or other CTAG2 peptides.

## Evidence limits

Adult normal panels do not establish pediatric safety. Tyser contains one gastrulating embryo; Petropoulos adds preimplantation embryos. Later fetal and pediatric coverage remains incomplete. RNA expression does not measure peptide-HLA presentation; HPA negative staining has finite sensitivity. PRAME informed earlier parameter choices and is not independent validation. All gene IDs, equivalence-group members, source versions and native units remain in the exports.
