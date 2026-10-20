# Tissue injury and manageability: evidence for policy

This review concerns **significant immune injury**, its consequences, persistence
and control. Complete organ removal is one evidence source; it is not a model for
every pattern of immune injury. Reproductive, Disposable, Protected and Vital are
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

## Representation and current decisions

The taxonomy's `context_reviews` records:

- `allowance_basis`: removal tolerance, conditional partial injury or explicit
  project exception.
- `local_steroid_control`: the supported scope, or that control is not established.
- `references`, `rationale` and `flags`: evidence and unresolved questions.

These are annotations, not extra RNA thresholds. Appendix/gallbladder join the
conditional Disposable category with an in-situ inflammation review. Cutaneous
melanocytes, follicles and tonsil retain their existing allowances. Whole skin,
salivary glands, thymus and spleen retain the earlier requested research allowances
while this broader manageability policy is assessed. The salivary evidence
particularly argues against treating its allowance as established low-risk injury.
Retina/eye/inner ear remain Vital; local steroid case experience does not reverse
that assignment.

For a revised default, assess expected functional loss, reversible versus permanent
injury, feasibility and reliability of local control, need for systemic rescue,
immune reserve, age, pregnancy and treatment modality. Missing evidence remains
unknown. The human clinical rationale is not transferred to dogs or mice.
