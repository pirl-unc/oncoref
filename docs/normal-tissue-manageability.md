# Tissue injury and manageability: evidence for policy

This review concerns **significant immune injury**, its consequences, persistence
and control. Complete organ removal is one evidence source; it is not a model for
every pattern of immune injury. Reproductive, Disposable, Manageable, Protected and Vital are
project categories. A separate field records evidence for local treatment.

## Published classifications

**Bradley et al., Nature Communications (2020)** provides the closest explicit
precedent. Supplementary Table 2 divides human source tissues into four risk groups,
plus exceptions. The following is an anatomical summary of that table, not the
oncoref assignment:

| Published group | Tissue coverage, condensed from the source labels | Study maximum RNA |
|---|---|---:|
| Danger | Aorta, coronary and tibial arteries; sampled brain regions and cervical spinal cord; atrial appendage and left ventricle | 1 TPM |
| Hazard | Lung; tibial nerve; colon; esophageal junction/mucosa/muscularis; renal cortex; liver; small intestine; stomach | 3 TPM |
| Caution | Spleen; pancreas; thyroid; blood; adrenal; pituitary | 10 TPM |
| Non-essential | Mammary tissue; cervix; fallopian tube; ovary; prostate; suprapubic skin; uterus; vagina; subcutaneous adipose; bladder; minor salivary gland; skeletal muscle | 30 TPM |
| Other | Testis and placenta; transformed lymphocyte/fibroblast samples also appear in this source column | No TPM limit in that study |

Source: [article and Supplementary Table 2](https://pmc.ncbi.nlm.nih.gov/articles/PMC7577998/).
The table was checked visually in the supplementary PDF, page 25. Its GTEx TPM
values are **not adopted as HPA nTPM limits**. The authors expressly acknowledge
that recognition below a proposed threshold is possible. A source classification
of skeletal muscle, bladder or pituitary is not clinical validation of tissue loss.
Transformed cell lines are not normal-organ absence evidence.

**Kaufmann et al., Oncotarget (2019)** independently uses brain/heart/lung as an
essential group and scores expression in other tissues as predicted adverse
events. It provides a precedent for graded penalties rather than treating every
nonessential tissue as freely targetable. Its calculations use normalized counts,
not our nTPM policy. [Primary study](https://www.oncotarget.com/article/26808/text/)

**Blagosklonny (2003)** supplies a conceptual tissue-selective rationale, including
melanocyte and reproductive-organ targeting and thyroid-directed treatment. It
does not provide a universal safe-organ whitelist.
[Paper](https://www.nature.com/articles/6601256)

## What local steroids establish

| Evidence | What it supports | Boundary |
|---|---|---|
| Johnson et al. (2009), antigen-directed TCR therapy | Some eye/ear inflammation following melanocyte targeting was treated with local steroids | This is direct cell-therapy evidence, but not proof that all ocular, retinal or inner-ear cell killing is reversible |
| SITC (2021) guideline | Selected inflammatory skin manifestations can be treated locally; management depends on severity | Severe, blistering or mucosal disease can require escalation; checkpoint toxicity does not validate control of every engineered T-cell product |
| Warner et al. (2019), 20 patients with checkpoint-associated sicca | Symptoms often improved following intervention | Salivary output remained very low; supports a persistent-function-loss concern rather than a reliable local rescue assumption |

Sources: [Johnson clinical trial](https://pubmed.ncbi.nlm.nih.gov/19451549/),
[SITC guideline](https://pmc.ncbi.nlm.nih.gov/articles/PMC8237720/),
[Warner clinical series](https://pubmed.ncbi.nlm.nih.gov/30996010/).

**Project interpretation:** local access and steroid responsiveness are useful
attributes, but they must be attached to a specific compartment, manifestation,
severity and therapy modality. They cannot by themselves turn a whole organ into
Disposable. Suppressing inflammation and restoring cells already lost are
different outcomes. No paper above establishes that immune injury can always be
stopped locally while preserving antitumor activity.

## Current categories and the intermediate Manageable tier

The current policy separates **Manageable (lower protection priority)** from
**Protected (higher protection priority)**. Bradley's Caution/Hazard distinction
is a published anchor. Our exact assignments remain project inferences:

| Category | Assignment relevant to this review | Evidence and interpretation |
|---|---|---|
| Disposable | Cutaneous melanocytes, hair follicles, tonsil, appendix, gallbladder | Selected compartment-loss or organ-removal rationale; acute immune injury still needs review |
| Manageable | Whole skin, spleen, thymus, thyroid, adipose | Still screened; any positive IHC excludes. This does not imply easy local control or harmless loss |
| Protected | Salivary glands, cartilage, esophagus and other protected GI tissues, lymph nodes and remaining listed protected organs | Persistent dysfunction, structural damage or important mucosal/distributed function |
| Vital | Existing vital and critical sensory/neurological organs | Highest protection priority; local access does not remove irreversible-loss concerns |

Reproductive remains the separate anatomical/developmental union. See the
[complete root lists and IHC-dependent RNA grid](normal-tissue-scopes.md).

Two additional primary clinical studies argue against a default spleen/thymus
allowance. [Kristinsson et al. (2014)](https://haematologica.org/article/view/6942)
followed 8,149 splenectomized veterans and reported elevated infection and
thromboembolic risks. [Kooshesh et al. (2023)](https://pmc.ncbi.nlm.nih.gov/articles/PMC10557034/)
found associations between adult thymectomy and increased mortality/cancer.
These are observational surgical cohorts, not causal estimates of antigen-directed
immune injury. They support retaining protection and reviewing age/immune reserve;
they do not prove equivalence between thymus, spleen and skin.

For cartilage, [Heinemeier et al. (2016)](https://pubmed.ncbi.nlm.nih.gov/27384346/)
found minimal replacement of mature articular collagen. [Braaten et al. (2020)](https://pmc.ncbi.nlm.nih.gov/articles/PMC7031780/)
reported persistent inflammatory arthritis after checkpoint treatment cessation.
These support concern about lasting joint injury, but are indirect evidence for a
cartilage category: matrix turnover is not identical to repair capacity, arthritis
is not selective chondrocyte injury, and a generic cartilage sample does not locate
the joint or exclude an airway site.

The category grid is a research prioritization tool. It does not establish that
every esophageal injury is worse than every skin injury. Tissue compartment, extent,
therapy modality, persistence and treatment response remain separate considerations.
Checkpoint-toxicity experience cannot establish the controllability of a persistent
engineered T-cell response.

## Machine-readable representation

The taxonomy's `context_reviews` records:

- `allowance_basis`: removal tolerance, an explicit narrow exception, or no default allowance.
- `local_steroid_control`: the supported manifestation or absence of established control.
- `references`, `rationale` and `flags`: evidence, consequences and unresolved questions.

The RNA ceilings are explicit **project starting parameters**: Vital 2/4/6/12,
Protected 6/8/12/24, Manageable 10/20/30/60 nTPM for missing/uncertain, Approved,
Supported and Enhanced negative IHC. The user requested that every absolute RNA
ceiling be doubled, giving a minimum ceiling of 2 nTPM. The Manageable name
refers to potentially manageable or compensated injury, subject to context review. No GTEx-to-HPA conversion or
clinical outcome calibration is claimed. HPA reliability labels do not express a
probability that a negative stain rules out antigen or peptide-HLA presentation.

Review annotations do not themselves change thresholds. Complete routine negative
IHC is required to earn relaxation, and any positive IHC outside Reproductive or
Disposable excludes. Whole skin/salivary/thymus/spleen allowances from policy
2026-10-20.3 are retained only in immutable historical reports.

Human clinical rationale is not transferred to dogs or mice. The shared vocabulary
maps anatomy; tolerability evidence, numerical policies, age and patient context
remain separate. The historical [CTA report](audits/cta-tissue-priorities-20261020/index.md)
separates the effects of reclassification from the new RNA grid.

## Combined absolute and selectivity gates

The same local IHC evidence now sets both an absolute RNA ceiling and minimum
pairwise selectivity: highest Reproductive RNA divided by that reference plus the
screened tissue's RNA. Disposable is excluded from this comparison and cannot
supply the required reproductive-positive evidence. Missing or zero reproductive
RNA cannot establish a selectivity score.

The approved percentage ladders are Vital 75/70/65/60%, Protected 70/65/60/55%,
and Manageable 65/60/55/50%, in order of missing/uncertain, Approved, Supported,
and Enhanced negative IHC. Admission requires valid reproductive RNA >0 or positive
reproductive IHC. These are target-informed project parameters with a different
meaning from the historical global deflated fraction. They do not estimate clinical risk.
See the [complete matrix and formula](normal-tissue-scopes.md#tissue-by-tissue-selectivity)
and the [candidate comparison](audits/cta-recalibration-20261021/index.md).
