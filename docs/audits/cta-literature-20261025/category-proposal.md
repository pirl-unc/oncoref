# Proposed general CTA assignment rules

Tracked in [#627](https://github.com/pirl-unc/oncoref/issues/627).

Status: proposal, 25 October 2026. No production memberships or thresholds have
changed. The examples below are review cases, not a gene-name whitelist or a
required retained count. Evidence adjudication can change an expected outcome.

Keep **CTA / Neither / Unsafe** as the three disjoint research buckets. Record
the biological evidence, uncertainty and therapeutic experience separately so
that “CTA” does not imply established clinical safety and “Unsafe” does not imply
that a treatment has injured a patient.

## 1. Classify findings before classifying genes

Represent every finding with the same fields whether it comes from an atlas,
an individual paper, proteomics or a reviewed clinical report:

| Field | Required meaning |
|---|---|
| Identity | Species, reference, gene/transcript/protein group; all plausible mappings and unresolved paralogs |
| Context | Tissue, cell type, developmental stage, sex, physiological state, donor/sample and healthy/nonmalignant/tumor status |
| Measurement | RNA, endogenous protein, presented peptide or clinical event; native units, assay limits, measured value or positive/negative/censored result |
| Support | Reagent/probe specificity, controls, orthogonal corroboration, sample count and independent-study lineage |
| Review | **Accepted**, **Unresolved** or **Rejected**, with reviewer, date, reason and exact source/figure/table |
| Scope | What the observation establishes and what it cannot establish; exact isoform, HLA, receptor and treatment context where relevant |

“Accepted” means fit for the stated inference, not an infallible experimental
result. A normal-protein claim needs a traceable primary observation, endogenous
nonmalignant cell identity, adequate gene/protein specificity, and controls or
corroboration sufficient for that context. One strong study can be enough;
several copies of the same study cannot. Antibody catalog ratings, a review's
CTA label, a search hit, and ectopic expression in a transfected cell are not
independent demonstrations of normal protein.

An abstract can establish a lead or a clearly reported observation. It cannot
certify unreported specificity controls. If the missing detail is material to
accepting a normal-protein veto, the finding stays Unresolved until checked.
Keep rejected results and their rejection reasons in the audit trail.

The same standards apply to HPA and publications. Do not automatically prefer
a curated negative over a credible specific positive, or automatically treat
every raw antibody stain as real protein. A positive and a negative in different
donors or cell states can both be correct. Reconcile their scope; do not average
them into a negative. Shared reagents and reused cohorts count once as evidence.

## 2. Keep the anatomy explicit

The tissue policy remains **Reproductive / Disposable / Manageable / Protected /
Vital**. Reproductive and Disposable are permitted by the research screen;
Manageable, Protected and Vital are not. Reproductive evidence is still required
for CTA admission. Disposable expression alone is not a CTA qualification.

Classify the actual cell compartment and developmental state. Adult hepatocytes
remain Vital; an explicitly fetal liver specimen can be Reproductive. Ordinary
epidermal keratinocytes do not inherit a hair-follicle allowance. A resident
macrophage does not automatically become a reproductive cell because its sample
came from testis or placenta. Protein made by an activated nonmalignant cell
remains normal-expression evidence. Therapeutic neutralization of a cytokine,
genetic loss of a protein and killing every expressing cell are different
interventions.

Use the common species fact/reference contract in [#599](https://github.com/pirl-unc/oncoref/issues/599)
and the existing tissue ontology. Orthology supplies a nomination, not measured
restriction. Human nTPM thresholds cannot be silently applied to mouse/dog or to
another assay. Preserve inapplicable anatomy, unresolved anatomy and absent data
as three different cases.

## 3. Apply one ordered decision policy

Evaluate all applicable gates and retain all reasons. Precedence determines the
single output bucket; it must not hide secondary findings.

| Order | General condition | Bucket |
|---|---|---|
| 1 | Accepted endogenous protein in a non-permitted normal compartment, a reviewed applicable safety veto, or a qualifying RNA exclusion | **Unsafe** |
| 2 | No decisive exclusion, but material unresolved normal-protein/identity conflict, qualifying donor variability, inadequate required coverage, failed selectivity or insufficient affirmative admission evidence | **Neither** |
| 3 | Supported reproductive expression, cancer nomination/evidence, adequate declared panel coverage, and all applicable normal gates pass with no unresolved material counterevidence | **CTA** |

Important consequences:

- **An accepted normal-protein positive is a veto even when RNA is low or zero.**
  Low staining still counts once accepted. Negative results elsewhere, a strong
  testis signal or a positive clinical trial cannot cancel that observation.
- **Uncertain protein is a hold, not accepted protein and not a negative.** A
  plausible conflicting normal stain yields Neither until identity/specificity
  is resolved. A specifically documented rejected off-target stain does not.
- **A reliable reproductive protein observation can satisfy admission below
  5 nTPM.** Reviewed reproductive evidence plus endogenous tumor-antigen
  recognition is another general admission route. Positive protein confined to
  Disposable tissue does not meet the reproductive requirement. Keep the
  5 nTPM reproductive floor for RNA-only admission.
- **A donor-variable exception is generalizable but narrow.** It requires
  reviewed vital-tissue donor data in a common cohort/version, reported zero or
  censoring in some donors and an applicable RNA-ceiling breach in others, plus
  usable negative local IHC. It yields Neither, never CTA; accepted normal
  protein, an applicable safety veto or normal RNA >30 still takes precedence.
- **Clinical experience is separately scoped.** A peptide/receptor-specific
  product can be promising even when unrestricted gene-level targeting is
  unsupported. A product cross-reacting with a different protein is not proof
  that the intended gene is expressed in that injured tissue. Keep reviewed
  precautionary holds explicit, with their scope and causal uncertainty.

Only Unsafe supplies the default forbidden-protein background. Neither remains
excluded from targets and from that default background. This is a screening
contract, not a claim that Neither proteins are absent or harmless.

## 4. Retain the approved numerical guardrails initially

The first implementation should fix evidence precedence and review handling
without retuning RNA numbers to preserve specific genes.

| Tissue | Missing/uncertain IHC | Approved negative | Supported negative | Enhanced negative |
|---|---:|---:|---:|---:|
| Vital ceiling, nTPM | 2 | 4 | 6 | 12 |
| Protected ceiling, nTPM | 6 | 8 | 12 | 24 |
| Manageable ceiling, nTPM | 10 | 20 | 30 | 30 effective cap |
| Vital selectivity | 80% | 75% | 70% | 65% |
| Protected selectivity | 75% | 70% | 65% | 60% |
| Manageable selectivity | 70% | 65% | 60% | 55% |

The underlying Enhanced Manageable ceiling is 60, but the approved hard >30
normal-expression exclusion limits it to 30. Keep the >30 TPM/nTPM check in each
source's native units and preserve its existing scope, including unsited contexts
that have not been established as permitted. Exactly 30 is not >30. Corneal
CP10K keeps its own rule. Invalidated measurements can be rejected with evidence;
an inconvenient valid value cannot be silently dropped or rounded down.

IHC relaxation requires an accepted negative covering the routine cells measured
in that tissue and the existing usable reproductive positive control. Missing,
unreliable, partial or disputed IHC gets no relaxation. A good antibody does not
relax unmeasured tissues. RNA and IHC from different donors do not establish a
paired RNA-to-protein calibration.

For the legacy bulk screen retain its documented reproductive maximum `R`,
normal value `N`, and selectivity `s`: absolute ceiling; `N < R`; and
`N <= R * (1-s)/s + 0.1`. The 0.1 allowance is a policy constant, not an inferred
confidence interval. Cell-state review triggers remain separate from those
bulk calculations and from the >30 cap. Failed enrichment alone yields Neither,
not Unsafe. These are research screening settings, not validated safety limits.

## 5. Add a continuous view of normal-expression breadth

The current reproductive maximum can make a diffuse gene appear selective.
Report the whole normal profile independently of that reproductive maximum.

For each comparable study and tissue category, define a descriptive burden:

`B = sum_j w_j * x_j / (x_j + L_j)`

Here `x_j` is a prespecified donor summary for normal compartment `j`, `L_j` is
the applicable RNA ceiling, and fixed anatomical weights sum to one. Report
median-based and upper-quantile-based versions separately where donor counts
support them. Do not substitute bulk ceilings for cell data without calibrating
that source. Partition weights by a fixed non-overlapping anatomy so adding
brain subregions or duplicating an atlas does not increase its weight.

This is a smooth, unitless expression summary relative to screening ceilings;
it is **not** a probability of injury, a reproductive RNA fraction or a measured
whole-body abundance. It rises with diffuse normal expression without requiring
a “number of tissues above 1 TPM” threshold. Raising reproductive RNA alone
cannot lower it. It never overrides a focal positive protein or a decisive RNA
exclusion. Calculate Vital, Protected and Manageable burdens separately.

Keep missing weight explicit. For the unmeasured portion of a fixed panel,
report a conservative possible range of contribution rather than renormalizing
over only favorable measured tissues. Sampling confidence intervals, when
estimable from independent donors, are separate from this missing-coverage range.
No donor replicates means no fabricated confidence interval. Analyze separate
studies in their native scales; do not pool TPM, nTPM and CP10K as one vector.

**Initially use breadth for review priority and sensitivity plots.** Calibrate
and lock any admission cutoff on a prespecified, diverse set of positive and
negative biological controls, then evaluate on held-out genes/families. Report
all gains/losses including familiar CTAs. If a validated breadth gate is later
added, a diffuse-but-below-ceiling failure belongs in Neither. No numerical
breadth threshold or additional retained count is claimed in this proposal.

## 6. Keep gene, protein and therapy conclusions distinct

At gene level, any accepted normal isoform blocks an **unrestricted gene-level**
CTA label. A reproductive-specific isoform can be separately reconsidered only
with source-bound expression and sequence identity. Different UTRs do not
necessarily create a different protein. Family-shared mass-spectrometry peptides
cannot establish every locus.

Preserve the existing complete-protein precedence `Unsafe > Neither > CTA` for
the default conservative export; it is not proof that every isoform was measured
in normal tissue. Any more precise isoform override needs explicit supporting
facts and its own versioned policy. Exact k-mer overlap remains a separate
operation, and absence of exact overlap does not establish TCR specificity.

Add normal and tumor peptide/HLA observations through the appropriate owning
data API and reference them in the shared evidence contract. The [HLA Ligand
Atlas](https://pmc.ncbi.nlm.nih.gov/articles/PMC8054196/) is a useful primary-source
starting point. Keep tissue, donor, HLA, sequence ambiguity and detection limits;
non-detection is not universal absence. The [IMA203 study](https://pubmed.ncbi.nlm.nih.gov/40205198/)
illustrates why a defined presented peptide and product can provide evidence
that cannot be transferred to every sequence from a gene.

## 7. Expected examples and release requirements

| Cases | Expected consequence of the general rules |
|---|---|
| FTMT, ABCC11, ACTRT1, CCDC158, YBX2 | Unsafe **if their normal-protein findings meet the acceptance criteria**; otherwise Neither pending material specificity/context questions |
| MAGEC2, LGALS13, PLCZ1 | Neither while the stated reagent, compartment or isoform conflict remains unresolved |
| CT55 | Low bulk RNA need not block admission when accepted reproductive/tumor-recognition evidence exists and normal gates pass |
| FOXR2 | Affirmative biology can satisfy admission; unresolved normal-antibody evidence still holds it in Neither |
| XAGE1A/B, CTAG1A/B | Preserve exact shared-protein attribution without pretending the data distinguish every locus |
| CTAG2 | Neither / Variable under the general reviewed donor rule, unless a stronger exclusion supersedes it |
| ROPN1 | >30 exclusion remains; proximity to the boundary prioritizes evidence review, not an exemption |
| SAGE1 and the other raw-antibody conflicts | Reject documented nonspecific signals, hold credible unresolved ones, exclude accepted normal positives; do not blanket-exclude all 180 |

The earlier 391/564/1,585 scenario assumes acceptance of five exclusions and
three holds. It is **not a promised result of the redesigned policy**. Applying
the same rules to all 2,540 nominations may uncover additional changes.

Implement as a versioned decision layer over source facts, with no gene-name
exceptions. Export every gate, all reasons, evidence status, context coverage,
policy/source identities and review task. Reuse [#620](https://github.com/pirl-unc/oncoref/issues/620)
for the active policy manifest and [#607](https://github.com/pirl-unc/oncoref/issues/607)
for consistent figure semantics.

Required checks include: renaming an otherwise identical synthetic gene does
not change its category; duplicate evidence adds no confidence; missing data
never becomes zero; improving a reproductive signal cannot erase accepted
normal counterevidence; unrelated negative IHC cannot relax a tissue; rejected
off-target staining cannot trigger exclusion; a stronger veto overrides Variable;
one normal isoform blocks unrestricted targeting without claiming all isoforms
were observed; wrong-reference/species/unit joins fail. Replay all 2,540 genes
and complete proteins, publish every category change, and update reports and
sequence manifests together before a release.

## Implementation and review issues

- [ ] [#628](https://github.com/pirl-unc/oncoref/issues/628): Adjudicate human neuronal FTMT protein and remove unrestricted CTA eligibility
- [ ] [#629](https://github.com/pirl-unc/oncoref/issues/629): Adjudicate normal neural ABCC11 protein before unrestricted CTA use
- [ ] [#630](https://github.com/pirl-unc/oncoref/issues/630): Resolve normal vascular and adrenal YBX2 protein against retained CTA status
- [ ] [#631](https://github.com/pirl-unc/oncoref/issues/631): Adjudicate renal CCDC158 protein and separate gene-level from isoform targeting
- [ ] [#632](https://github.com/pirl-unc/oncoref/issues/632): Curate normal keratinocyte ACTRT1 protein as non-reproductive evidence
- [ ] [#633](https://github.com/pirl-unc/oncoref/issues/633): Hold unrestricted MAGEC2 targeting while Purkinje-cell staining is unresolved
- [ ] [#634](https://github.com/pirl-unc/oncoref/issues/634): Resolve extra-placental PP13/LGALS13 protein and galectin specificity
- [ ] [#635](https://github.com/pirl-unc/oncoref/issues/635): Resolve pancreatic NYD-SP27/PLCZ1 evidence before unrestricted CTA use
- [ ] [#636](https://github.com/pirl-unc/oncoref/issues/636): Reconcile primary normal-antibody signals for 180 retained CTA assignments
- [ ] [#637](https://github.com/pirl-unc/oncoref/issues/637): Measure continuous normal-expression breadth without rescue by a reproductive maximum
- [ ] [#638](https://github.com/pirl-unc/oncoref/issues/638): Complete individual CTA literature adjudication beyond the 32-gene initial review
