# TCR-T clinical evidence, serious events and RNA screening credit

This is a targeted, source-linked review of published CTA-directed TCR-T experience.
**A favorable trial is evidence for a particular peptide, HLA, receptor and regimen.**
The approved policy grants a modest protein-profile RNA credit, with those limitations
attached; it does not declare the whole protein clinically safe.

## How the evidence changes the screen

1. Require reviewed results from an actually treated cohort, with its denominator.
2. Identify the peptide, its protein origins, the HLA and the receptor/product.
3. Review whether any grade 4-5 event was suspected or confirmed to involve recognition
   of healthy tissue. Such a signal blocks credit for that receptor. Unassessed
   specificity and ungraded suspected specificity injury also block credit.
4. Otherwise grant **one RNA ladder step**, with the Enhanced-negative ceiling as the cap.
   Trials do not stack. The IHC grade is not changed or invented.
5. Keep independent protein and clinical-injury exclusions, required coverage,
   cancer expression, Reproductive evidence and diffuseness rules.

| Tier | Baseline: missing / Approved / Supported / Enhanced IHC | With qualifying clinical credit |
|---|---|---|
| Vital | 4 / 8 / 16 / 32 | 8 / 16 / 32 / 32 |
| Protected | 8 / 16 / 32 / 64 | 16 / 32 / 64 / 64 |
| Manageable | 16 / 32 / 64 / 128 | 32 / 64 / 128 / 128 |

Values are nominal TPM/nTPM screening parameters, not clinically measured toxicity limits.
Corneal CP10K retains its separate nominal /100 scale. An unfavorable product does not
automatically make every peptide of its intended protein Unsafe. For example, titin
cross-recognition implicates a receptor and a different protein. Accepted gene-level
normal-protein/clinical counterevidence is still handled by the existing safety layer.

CTAG2 was already CTA before this credit. Its HPA v23 heart value is 5.1 nTPM with
Enhanced-negative IHC and a 32 nTPM ceiling; the cap leaves that ceiling at 32.
The shared clinical peptide is in **LAGE-1A**, not every CTAG2 isoform. Reference
isoform IDs, peptide positions and sequence checksums are in the evidence ledger.
Gene-level RNA cannot identify which isoform supplied a normal-tissue measurement.

## Reading the safety findings

“No identified grade 4-5 specificity signal” is a reviewed attribution conclusion,
not an author-reported zero incidence or proof that risk is zero. Small cohorts,
HLA selection, follow-up, conditioning and incomplete mechanistic investigations
limit the inference. This initial review is not systematic or exhaustive. Unlisted
products remain unreviewed; registration, preclinical work and missing reporting do
not qualify. No minimum response rate is imposed on the independent safety credit.

The sentinel table includes deaths, significant organ injury and severe inflammatory
events. It is not a substitute for each paper's complete AE/SAE tables. Unknown CTCAE
grades remain null. Case reports do not supply new treated-cohort denominators.

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

[Full attribution, epitope/HLA mappings, cohort denominators and limitations](index.md). A favorable reviewed product supplies one capped RNA step, never protein or clinical-injury clearance. The exact tested peptide matters: NY-ESO-1 and **CTAG2/LAGE-1A** share SLLMWITQC; this does not clear LAGE-1B or other CTAG2 peptides.

## Measured impact

**0 protein category changes** versus the pinned oncoref.cta-selection.2026-10-29.3 snapshot. [Membership comparison](membership-changes.csv). RNA-ceiling changes and individual exclusion effects are in the linked context and assay tables; a changed ceiling need not change the final category.


## Letetresgene autoleucel / NY-ESO-1 c259

**RNA credit: one capped step.** Published MRCLS cohort and neurological case report; not a systematic review of every c259 trial.

HLA: HLA-A*02:01, HLA-A*02:05, HLA-A*02:06.

- **CTAG1A** (`ENSG00000268651`): `SLLMWITQC`. NY-ESO-1 157-165.

- **CTAG1B** (`ENSG00000184033`): `SLLMWITQC`. NY-ESO-1 157-165.

- **CTAG2** (`ENSG00000126890`): `SLLMWITQC`. LAGE-1A 157-165; no transfer to LAGE-1B or other CTAG2 peptides.


**20 treated**, Advanced/metastatic myxoid/round cell liposarcoma. No grade 4-5 event was attributed to peptide-specific recognition of normal tissue. Grade 4 cytopenias occurred; one late cardiac death remained potentially treatment-related without an established specificity mechanism. This is a curation conclusion about reported attribution, not an author-reported zero-event rate or proof of safety. Safety cutoff May 2022; survival follow-up January 2023. A fatal cardiac event occurred on day 161. [Letetresgene Autoleucel in Advanced/Metastatic Myxoid/Round Cell Liposarcoma](https://doi.org/10.1200/JCO-24-01466).


**Fatal cardiac arrest, day 161** (1 patient(s); grade 5). Treatment attribution: potentially related. Specificity: not established. Progressive disease, hypotension, renal insufficiency and intervening radiation/pazopanib complicate attribution. The report does not establish CTAG2-mediated heart injury. [Letetresgene Autoleucel in Advanced/Metastatic Myxoid/Round Cell Liposarcoma](https://doi.org/10.1200/JCO-24-01466).


**Guillain-Barre / inflammatory demyelinating neuropathy** (2 patient(s); grade unreported). Treatment attribution: suspected. Specificity: not established. Authors suspected treatment association and discussed indirect immune effects. No identified peptide cross-reaction explains the cases. Case reports have no independent treated-cohort denominator and confer no favorable credit. [Guillain-Barre syndrome observed with adoptive transfer of lymphocytes genetically engineered with an NY-ESO-1 reactive T-cell receptor](https://doi.org/10.1186/s40425-019-0759-x).


## IMA203

**RNA credit: one capped step.** Published 2025 phase 1a/1b monotherapy analysis; not every ongoing cohort.

HLA: HLA-A*02:01.

- **PRAME** (`ENSG00000185686`): `SLLQHLIGL`. PRAME-004; sequence positions 425-433 in matched reference proteins.


**40 treated**, PRAME-positive advanced solid tumors. No treatment-related deaths after IMA203 infusion were reported. Severe cytopenias and CRS occurred; no clinical renal toxicity was seen despite weak reactivity in one in-vitro renal donor. Normal peptide non-detection is assay-limited; this does not clear all PRAME peptides. Interim cutoff July 4, 2023; phase 1a and phase 1b monotherapy. [Autologous T cell therapy for PRAME+ advanced solid tumors in HLA-A*02+ patients: a phase 1 trial](https://doi.org/10.1038/s41591-025-03650-6).


**Grade 4 transaminitis during a grade 3 CRS episode** (1 patient(s); grade 4). Treatment attribution: related. Specificity: not established. Recorded as part of CRS, without evidence that the TCR recognized PRAME in hepatocytes. [Autologous T cell therapy for PRAME+ advanced solid tumors in HLA-A*02+ patients: a phase 1 trial](https://doi.org/10.1038/s41591-025-03650-6).


## Afamitresgene autoleucel / c1032

**RNA credit: one capped step.** Phase 1 and SPEARHEAD-1 cohort 1 publications; not a postmarketing surveillance census.

HLA: HLA-A*02:01, HLA-A*02:02, HLA-A*02:03, HLA-A*02:06.

- **MAGEA4** (`ENSG00000147381`): `GVYDGREHTV`. MAGE-A4 230-239; only peptide-containing isoforms.


**38 treated**, MAGE-A4-positive advanced solid tumors. Two possibly treatment-related deaths occurred: aplastic anemia and ischemic stroke. Marrow evaluation did not support MAGE-A4-directed injury; high-dose conditioning and inflammation were considered contributors. The fatal stroke has no demonstrated specificity mechanism; neither death is relabelled unrelated. Published phase 1 analysis; fatal events at days 17 and 55. [Autologous T cell therapy for MAGE-A4+ solid cancers in HLA-A*02+ patients: a phase 1 trial](https://doi.org/10.1038/s41591-022-02128-z).


**52 treated**, Synovial sarcoma or myxoid/round cell liposarcoma. No grade 5 adverse events or treatment-related deaths were reported; severe cytopenias and mostly low-grade CRS remained important regimen toxicities. No identified specificity event in this selected cohort does not establish a population-wide zero risk. Published cohort 1 analysis; later cohorts are outside this snapshot. [SPEARHEAD-1: a single-arm phase 2 trial of afamitresgene autoleucel in advanced synovial sarcoma and myxoid/round cell liposarcoma](https://doi.org/10.1016/S0140-6736(24)00319-2).


**Fatal aplastic anemia, day 55** (1 patient(s); grade 5). Treatment attribution: possibly related. Specificity: not supported. No marrow MAGE-A4 or enrichment of engineered T cells was detected; conditioning and systemic inflammation were considered contributors. [Autologous T cell therapy for MAGE-A4+ solid cancers in HLA-A*02+ patients: a phase 1 trial](https://doi.org/10.1038/s41591-022-02128-z).


**Fatal ischemic cerebrovascular accident, day 17** (1 patient(s); grade 5). Treatment attribution: possibly related. Specificity: not established. Occurred after grade 3 neurotoxicity; no demonstrated MAGE-A4-specific brain injury mechanism was reported. [Autologous T cell therapy for MAGE-A4+ solid cancers in HLA-A*02+ patients: a phase 1 trial](https://doi.org/10.1038/s41591-022-02128-z).


## MAGE-A3 affinity-enhanced a3a receptor

**RNA credit: none.** Two fatal cases and the mechanistic investigation; different MAGE-A3 receptors require separate review.

HLA: HLA-A*01:01.

- **MAGEA3** (`ENSG00000221867`): `EVDPIGHLY`. Intended MAGE-A3 target; adverse specificity is titin cross-recognition.


**2 treated**, One melanoma and one multiple myeloma patient. Fatal cardiac injury was attributed to T-cell cross-recognition, with titin identified as the likely unintended antigen. The finding blocks credit for this receptor; it does not prove that MAGE-A3 protein was present in heart or that every MAGE-A3 receptor is unsafe. Fatal cardiac injury within 4-5 days of infusion. [Identification of a Titin-Derived HLA-A1-Presented Peptide as a Cross-Reactive Target for Engineered MAGE A3-Directed T Cells](https://doi.org/10.1126/scitranslmed.3006034).


**Fatal cardiac injury with titin cross-recognition** (2 patient(s); grade 5). Treatment attribution: related. Specificity: suspected. Mechanistic work identified titin as the likely unintended target; MAGE antigens were not detected in patient heart samples. [Identification of a Titin-Derived HLA-A1-Presented Peptide as a Cross-Reactive Target for Engineered MAGE A3-Directed T Cells](https://doi.org/10.1126/scitranslmed.3006034).


## IMA202

**RNA credit: one capped step.** Published IMA202-101 cohort; no transfer from other MAGEA1 products.

HLA: HLA-A*02:01.

- **MAGEA1** (`ENSG00000198681`): `KVLEYVIKV`. MAGE-A1 peptide; other receptors are separate products.


**16 treated**, MAGE-A1-positive advanced solid tumors. Severe cytopenias and one grade 3 CRS occurred. No on-target/off-tumor toxicity was observed; fatal dyspnea was considered unrelated. The primary paper header lists other MAGEA1 trials; IMA202-101 identity is independently linked by the Columbia trial listing. No objective RECIST responses were required for this safety-only credit. Published cohort snapshot; later trials are not pooled. [First-in-human dose escalation trial to evaluate the clinical safety and efficacy of an anti-MAGEA1 autologous TCR-transgenic T cell therapy in relapsed and refractory solid tumors](https://doi.org/10.1136/jitc-2023-008668).


**Fatal dyspnea** (1 patient(s); grade 5). Treatment attribution: unrelated. Specificity: not supported. Investigators considered this death unrelated to treatment; severe conditioning cytopenias and grade 3 CRS occurred in the cohort. [First-in-human dose escalation trial to evaluate the clinical safety and efficacy of an anti-MAGEA1 autologous TCR-transgenic T cell therapy in relapsed and refractory solid tumors](https://doi.org/10.1136/jitc-2023-008668).


## ADP-A2M10 / MAGE-A10 c796

**RNA credit: one capped step.** Published NSCLC cohort; another solid-tumor cohort not yet adjudicated here.

HLA: HLA-A*02:01, HLA-A*02:06.

- **MAGEA10** (`ENSG00000124260`): `GLYDGMEHL`. MAGE-A10 peptide.


**11 treated**, Advanced non-small cell lung cancer. One grade 4 CRS and one treatment-related death after repeat infusion/intensified conditioning. Investigations favored conditioning over MAGE-A10-mediated marrow injury. Patient marrow antigen testing was limited by specimen quality; absence of engineered-cell enrichment does not prove no risk. Published cohort snapshot; later trials are not pooled. [Phase I clinical trial evaluating the safety and efficacy of ADP-A2M10 SPEAR T cells in patients with MAGE-A10+ advanced non-small cell lung cancer](https://doi.org/10.1136/jitc-2021-003581).


**Grade 4 CRS** (1 patient(s); grade 4). Treatment attribution: related. Specificity: not established. Required intensive supportive care; symptoms resolved. The subsequent death of this patient was attributed to disease progression. [Phase I clinical trial evaluating the safety and efficacy of ADP-A2M10 SPEAR T cells in patients with MAGE-A10+ advanced non-small cell lung cancer](https://doi.org/10.1136/jitc-2021-003581).


**Fatal marrow aplasia, sepsis and respiratory distress after repeat infusion** (1 patient(s); grade 5). Treatment attribution: related. Specificity: not supported. Authors favored intensified lymphodepletion; patient marrow antigen testing was limited. Fatality was 98 days after the second infusion. [Phase I clinical trial evaluating the safety and efficacy of ADP-A2M10 SPEAR T cells in patients with MAGE-A10+ advanced non-small cell lung cancer](https://doi.org/10.1136/jitc-2021-003581).


## TBI-1301

**RNA credit: none.** Published nine-person cohort; unresolved grade 3 lung-injury mechanism prevents favorable credit.

HLA: HLA-A*02:01, HLA-A*02:06.

- **CTAG1A** (`ENSG00000268651`): `SLLMWITQC`. NY-ESO-1 157-165.

- **CTAG1B** (`ENSG00000184033`): `SLLMWITQC`. NY-ESO-1 157-165.

- **CTAG2** (`ENSG00000126890`): `SLLMWITQC`. LAGE-1A 157-165; no transfer to LAGE-1B or other CTAG2 peptides.


**9 treated**, NY-ESO-1-positive advanced solid tumors. One grade 3 lung injury with TCR-T infiltration; resolved with steroids. Authors favored cytokine-mediated capillary leak but could not exclude normal-peptide recognition. Retain unresolved specificity as a review hold for this product; do not borrow favorable c259 evidence to clear TBI-1301. Published cohort snapshot; later trials are not pooled. [NY-ESO-1-specific redirected T cells with endogenous TCR knockdown mediate tumor response and cytokine release syndrome](https://doi.org/10.1136/jitc-2021-003811).


**Steroid-responsive lung injury with TCR-T infiltration** (1 patient(s); grade 3). Treatment attribution: related. Specificity: unresolved. Cytokine-mediated capillary leak was favored, but normal-peptide recognition remained a possibility. Below grade 4 does not mean clinically irrelevant. [NY-ESO-1-specific redirected T cells with endogenous TCR knockdown mediate tumor response and cytokine release syndrome](https://doi.org/10.1136/jitc-2021-003811).


## NCI HLA-A*02-restricted MAGE-A3/A9/A12-reactive TCR

**RNA credit: none.** Primary abstract and indexed full-text epitope/mechanism passages; normal neuronal antigen attribution remains proposed.

HLA: HLA-A*02:01.

- **MAGEA3** (`ENSG00000221867`): `KVAELVHFL`. Shared intended peptide.

- **MAGEA9** (`ENSG00000123584`): `KVAELVHFL`. Shared intended peptide.

- **MAGEA12** (`ENSG00000213401`): `KMAELVHFL`. Cross-recognized variant peptide; suspected normal brain source.

- **MAGEA9B** (`ENSG00000267978`): `KVAELVHFL`. Reference-verified shared peptide; MAGE-A9 protein-identical partner, not separately attributed clinical expression.


**9 treated**, Advanced cancers. Three neurological events, including two fatal cases. Authors implicated cross-recognized MAGE-A12 in brain. Brain staining recognized multiple MAGE-family proteins; MAGE-A12 attribution is supported but not uniquely proven. Published cohort snapshot; later trials are not pooled. [Cancer regression and neurological toxicity following anti-MAGE-A3 TCR gene therapy](https://pubmed.ncbi.nlm.nih.gov/23377668/).


**Fatal neurological injury; MAGE-A12 cross-recognition implicated** (2 patient(s); grade 5). Treatment attribution: related. Specificity: suspected. Two deaths among three neurological cases in nine treated patients; the third recovered. Attribution of normal brain antigen to MAGE-A12 is proposed rather than exclusive. [Cancer regression and neurological toxicity following anti-MAGE-A3 TCR gene therapy](https://pubmed.ncbi.nlm.nih.gov/23377668/).


## API, source storage and reproduction

```python
from oncoref import tcr_trials, reference_data

tcr_trials.products()       # product-level eligibility, exact peptide origins and HLA
tcr_trials.cohorts()        # treated denominators and explicit attribution review
tcr_trials.events()         # sentinel events, nullable severity and attribution
tcr_trials.protein_credit(["ENSG00000126890"])  # CTAG2; not isoform-wide clearance
reference_data.ensure("tcr_trial_lete_mrcls_2025", "snapshot-2026-10-29")
```

`oncoref tcr-trials --table events` exports JSON. `--table all` includes the complete
ledger. Original papers/protocols use the same checksum-pinned source API and
`<cache>/sources/<source>/<version>/<filename>` layout as other references;
`CANCERDATA_DATA_DIR` selects the cache root. Paper URLs can change upstream:
checksum mismatch fails explicitly instead of silently updating the reviewed snapshot.

Reproduce with `python scripts/build_tcr_trial_report.py` after the selection replay.
The curation lives in `oncoref/data/tcr-trial-evidence.json`; both the ledger and policy
code participate in the replay receipts. The report does not re-adjudicate literature
on each build. New evidence requires an explicit reviewed ledger update.

[Summary table](tcr-trial-summary.csv) · [Cohorts](cohorts.csv) ·
[Sentinel events](sentinel-events.csv) · [Complete evidence](evidence.json) ·
[Context-level RNA effects](clinical-credit-contexts.csv.gz) · [RNA effects by assay](clinical-credit-impact.csv)


[PDF report](cta-tcr-safety-1.24.0.pdf)
