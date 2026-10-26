# CTA definitions: a single decision hierarchy

**CTA means a reproductive-associated target candidate that passes this research screen.**
It does not mean clinically proven safe. The current policy is
`oncoref.target-sets.2026-10-25.2`: **213 CTA, 741 Neither and 1,586 Unsafe** among
2,540 mapped nominations. Across all 20,097 canonical human genes the counts are
213 / 2,179 / 17,705, respectively.

[Current report and changes](audits/cta-policy-20261025/index.md) ·
[Four-page PDF brief](https://github.com/pirl-unc/oncoref/releases/download/v1.20.0/cta-policy-1.20.0.pdf) ·
[All candidate decisions](audits/cta-policy-20261025/candidate-results.csv.gz) ·
[Reviewed findings and citations](audits/cta-policy-20261025/reviewed-findings.csv) ·
[Follow-up questions](audits/cta-policy-20261025/follow-up.csv)

## The three outcomes

| Outcome | Meaning | Target set | Forbidden background |
|---|---|---|---|
| **CTA** | Affirmative reproductive evidence; no accepted exclusion or recorded material hold; required measurement screen passes | Yes | No |
| **Unsafe** | Accepted normal-protein evidence, excessive normal RNA, or an applicable reviewed safety veto | No | Yes |
| **Neither** | Unresolved evidence, inadequate target evidence/coverage, or insufficient selectivity | No | No |

Unsafe means *excluded from unrestricted gene-level targeting*, not necessarily
observed clinical toxicity. Neither is not a safe reserve list or proof of absence.

## Read the rules in this order

1. **Identify the molecule and context.** Preserve species, reference gene/protein,
   tissue, cell type, developmental/activation state, assay, units and source
   version. An unresolved family name does not substitute for a mapped gene.
   Missing data are unknown, never zero.
2. **Exclude credible normal expression.** Accepted protein production in a
   non-permitted normal compartment, curated positive normal IHC (including Low),
   RNA exceeding its applicable absolute limit, or a reviewed gene-level safety
   veto gives **Unsafe**. Any audited non-permitted or unsited normal value **>30
   TPM or nTPM** also gives Unsafe. TPM and nTPM are checked separately.
3. **Hold material uncertainty.** An unresolved normal-expression/identity conflict,
   uncertain antibody evidence, a source-specific normal-cell RNA trigger, or a
   qualified donor-variable finding gives **Neither**, unless step 2 already
   supplies a decisive exclusion.
4. **Require affirmative reproductive evidence.** Admission needs **≥5 nTPM in
   a reproductive compartment**, **reliable reproductive protein**, or **accepted
   reproductive expression plus tumor expression/antigen recognition evidence**.
   Protein or reviewed biological support can satisfy admission below 5 nTPM.
   Tumor expression or a trial registration alone cannot.
5. **Require quantitative restriction and assessed coverage.** Each screened bulk
   normal tissue must pass both its absolute ceiling and selectivity comparison.
   The pinned candidate screen, routine RNA panel and identity/quality checks
   must pass. Failure of admission, coverage or selectivity alone gives **Neither**.
6. **Call the remainder CTA.** A source-backed candidate nomination is required;
   passing a normal-expression screen alone does not make an uncurated gene a CTA.

Every rule is evaluated. The report preserves **all applicable reasons** and
chooses one principal reason for navigation. Reproductive RNA, negative IHC and
clinical interest cannot cancel accepted normal-protein evidence. There are no
new gene-name branches or therapeutic-familiarity exemptions.

**Variable is a narrow reconsideration rule**, currently supported for CTAG2
only. Reviewed vital-tissue donors must include reported zero and a failing RNA
value in the same tissue/cohort/version/units, with appropriate negative local
IHC. It can turn the corresponding ordinary RNA exclusion into Neither; it does
not override >30, accepted normal protein, positive curated normal IHC or a safety
veto. RNA and IHC volunteers are not assumed to be the same people. An RNA-zero
volunteer is not thereby eligible for treatment.

## Which tissues count?

The anatomy policy is unchanged. **Reproductive and Disposable are permitted**;
Manageable, Protected and Vital are screened.

| Tissue class | Meaning and examples |
|---|---|
| Reproductive | Union of testicular, other male reproductive, ovarian, female tract, mammary, placental, other extraembryonic and explicitly embryonic/fetal compartments |
| Disposable | Narrow project allowances: cutaneous melanocytes, hair follicles, tonsil, appendix and gallbladder |
| Manageable | Whole skin, spleen, thymus, thyroid and adipose; still screened |
| Protected | Salivary glands, cartilage, lymph nodes and specified mucosal/other organs; still screened |
| Vital | Critical neurological, sensory and systemic organs, including heart, lung, brain, kidney, liver and adrenal; strictest RNA limits |

These are examples, not a replacement for the [complete anatomy lists and
source mappings](normal-tissue-scopes.md). A reproductive organ's resident immune
cells are not automatically reproductive. Fetal liver is developmental reproductive
context; adult maternal liver remains Vital. Disposable does not supply affirmative
reproductive evidence. These are general research sets, not pregnancy- or
fertility-preserving clearance.

The [tissue rationale review](normal-tissue-manageability.md) anchors the graded
approach in [Bradley et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC7577998/) and
organ-specific evidence; the project's assignments and RNA limits are not those
papers' validated safety boundaries. Shared anatomy mappings cover human, mouse
and dog. **This numerical CTA policy accepts human measurements only**; it does
not transfer human tolerability or thresholds to another species.

## Exact RNA and IHC rules

Each cell below gives **maximum nTPM / selectivity setting**. These are the
existing settings, unchanged by this literature update.

| Tissue | Missing, uncertain or unusable IHC | Approved negative | Supported negative | Enhanced negative |
|---|---:|---:|---:|---:|
| Vital | 2 / 80% | 4 / 75% | 6 / 70% | 12 / 65% |
| Protected | 6 / 75% | 8 / 70% | 12 / 65% | 24 / 60% |
| Manageable | 10 / 70% | 20 / 65% | 30 / 60% | 60 / 55%* |

\*The independent >30 rule makes the effective Manageable maximum **30**, even
with Enhanced negative IHC. Equality at an absolute limit passes; >30 fails.

Let **R** be the highest reproductive bulk nTPM, **N** the screened normal tissue's
nTPM, **C** its ceiling, and **s** its selectivity setting as a fraction. Require:

```text
N ≤ C
N ≤ R × (1 − s) / s + 0.1
N < R
```

The 0.1 nTPM allowance softens low-level ratios; it is **not** a fitted uncertainty
interval. The settings are not literal percentages of whole-body expression.
For R = 5 and baseline Vital s = 0.80, the selectivity limit is 1.35 nTPM:
N = 1.5 passes the absolute ceiling of 2 but fails selectivity. That isolated
selectivity failure leads to Neither, not Unsafe.

Relaxed columns require **complete negative routine local IHC**, using the
weakest usable confidence, plus a reliable reproductive-positive control in the
curated gene profile. Incomplete or uncertain staining, or no such positive
control, uses baseline limits. This is profile-level support, not proof of the
same reagent/batch/donor across tissues or a measured RNA-to-protein conversion.
Accepted positive normal staining wins over negative staining elsewhere.

## Cell-resolved evidence stays separate

Bulk tissue dilution can hide activated immune cells or rare ocular cells.
Their expression is therefore checked with source-specific rules, without mixing
native units or using a negative bulk stain to relax unpaired cell RNA.

| Source | Trigger for Neither, unless a stronger exclusion applies |
|---|---|
| HPA single-cell clusters | ≥5 nTPM and ≥20 reads in a non-permitted cluster, at least as high as the reproductive-cell maximum; missing reproductive coverage or unresolved anatomy remains a concern |
| Original DICE distributions | ≥5 TPM in ≥3 reported expression values for the same cell type/state |
| HPA/Monaco immune summaries | ≥5 in the source's native TPM or nTPM |
| Adult corneal atlas | Mean CP10K ≥1, ≥5 positive cells and ≥20 counts in ≥2 annotated sample groups |
| Primary antibodies | Existing uncertain-profile normal staining, or a specifically reviewed material protein/identity conflict |

The >30 TPM/nTPM rule overrides these holds. CP10K is not converted to TPM and
retains its separate trigger. DICE distribution columns are not stable donor
identities across genes; HPA's DICE summary is not independent replication.
Activated nonmalignant cells are normal-cell contexts. Rarity or cytokine-blocking
therapy does not automatically permit killing the producing cells.
[Sources, uniform cache/API and native-unit rules](normal-cell-evidence.md).

## What the literature changes

Findings have **Accepted**, **Unresolved** or **Rejected** status. Acceptance is a
curation judgement about a specific observation, not a paper-count score. Record
identity, endogenous versus engineered expression, normal versus malignant state,
anatomy, assay specificity and limitations. Negative results in another context do
not refute a positive automatically. Cross-reactivity or an unresolved isoform
can prevent a signal from becoming an accepted gene-level finding.

| Example | Current decision | Why |
|---|---|---|
| FTMT | Unsafe | Normal human neuronal protein with reported specificity controls. [Primary study](https://pubmed.ncbi.nlm.nih.gov/28386150/) |
| ABCC11 | Unsafe | Reported normal neural protein supported by multiple modalities; abstract-level review, full reagent details still a follow-up. [Study](https://pubmed.ncbi.nlm.nih.gov/16359813/) |
| YBX2 | Unsafe | Normal vascular/adrenal protein; specimen and polyclonal-control limitations retained. [Study](https://pubmed.ncbi.nlm.nih.gov/16479255/) |
| CCDC158 | Unsafe | Reported renal tubular protein/function; abstract-level review and isoform follow-up retained. [Study](https://pubmed.ncbi.nlm.nih.gov/39319391/) |
| ACTRT1 | Unsafe | Endogenous protein in differentiating normal keratinocytes; cultured state is explicit. [Study](https://pubmed.ncbi.nlm.nih.gov/33972689/) |
| MAGEC2 | Neither / Conflict | Purkinje staining requires reagent/family-specificity resolution. [Study](https://pubmed.ncbi.nlm.nih.gov/12920247/) |
| LGALS13 | Neither / Conflict | Historical extra-placental reactivity conflicts with placental restriction; identity/anatomy unresolved. [Early protein study](https://pubmed.ncbi.nlm.nih.gov/10527825/), [placental study](https://pubmed.ncbi.nlm.nih.gov/19497882/) |
| PLCZ1 | Neither / Conflict | Qualitative pancreatic NYD-SP27 RNA needs abundance/isoform reconciliation. The paper stains CFTR, not PLCZ1; the prior protein interpretation is rejected. [Study](https://pubmed.ncbi.nlm.nih.gov/14697805/) |
| CT55 | CTA | Reproductive plus leukemia-recognition evidence satisfies admission below 5 nTPM; normal gates pass. [Study](https://pmc.ncbi.nlm.nih.gov/articles/PMC5709753/) |
| FOXR2 | Neither / Antibody | Affirmative reproductive/cancer biology does not resolve normal-antibody concerns. [Study](https://pmc.ncbi.nlm.nih.gov/articles/PMC9437574/) |
| ABCB5, IL5, IL13, PAPPA | Unsafe, unchanged | Previously accepted normal-cell protein production; reproductive enrichment cannot override it. [Reviewed source records](audits/cta-policy-20261025/reviewed-findings.csv) |

**Clinical evidence is product-specific.** CTAG1A/B and XAGE1A/B remain CTA.
PRAME is Unsafe after accepted normal adrenal-protein evidence; MAGEA4 is
Neither because its raw normal-antibody signals remain unresolved. Their
clinical products retain their own evidence, summarized in the
[completed adjudication report](audits/cta-adjudication-20261025/index.md). A trial, response or lack of observed toxicity for a selected
peptide/HLA/receptor does not establish safety of every peptide or proteoform from
that gene. Conversely receptor-specific cross-reactivity is not automatically
proof of normal expression of the intended target gene.

**Review coverage is incomplete.** The literature search covered all 2,540
nominations, but only 32 genes received new individual primary-paper review;
46 have specific notes including prior work. The 178 cases remaining from the original 180-gene raw-antibody queue have
now all received completed adjudications: 177 Neither and 1 Unsafe; the other
two were addressed in the initial eight-gene review. No case from this queue
remains CTA. This resolves their classification, not every stain's molecular
identity. Other unmeasured contexts and the wider literature remain gaps. See the [review coverage and
open issues](audits/cta-literature-20261025/issues.md).

The maximum-based selectivity rule can still favor a gene with one high
reproductive compartment and diffuse low normal RNA. Continuous normal breadth,
empirical assay uncertainty and broad isoform adjudication are **not yet selection
rules**. These limitations remain visible rather than being hidden in an
unvalidated weighted score.

## API, reproduction and exact proteins

```python
from oncoref import target_sets, cta_policy, normal_rna_review

policy = cta_policy.manifest()      # current policy, thresholds, counts, hashes
findings = cta_policy.findings()   # accepted/unresolved/rejected cited facts
adjudications = cta_policy.adjudications()  # 178 completed antibody decisions
rows = target_sets.results()      # decisions, all reasons, measurement provenance
cta = target_sets.gene_ids("CTA")
neither = target_sets.gene_ids("Neither")
forbidden = normal_rna_review.forbidden_gene_ids()  # includes extra reference loci
followup = normal_rna_review.followup_reviews()
```

`cta_policy.classify(facts, findings)` is the pure decision function. It requires
explicit, reference-bound measurement-gate facts, not a previous final category.
`cta_policy.manifest()` and the result rows identify the **active** policy;
older classifier module constants identify their historical algorithms.

```sh
.venv/bin/python scripts/build_cta_policy.py \
  --output /tmp/cta-policy-report \
  --assets /tmp/cta-policy-proteins \
  --sequence-input /path/to/pinned-original-sequence-assets
```

The builder checks hashes, reconstructs facts from frozen measurement audits,
reproduces the previous categories with the original accepted findings, then
applies the reviewed additions. `--publish-dataset` updates the packaged data
after invariants pass. Outputs include all decisions, each reason, changed rows,
follow-up, exact proteins and checksums. Historical reports and source data remain
unchanged. [Current build inventory](audits/cta-policy-20261025/protein-assets.json).

Gene categories and complete exact-sequence protein categories are disjoint.
A protein encoded by multiple genes takes **Unsafe → Neither → CTA** precedence.
An unclassified encoding gene is Neither. Gene-level expression does not resolve
which isoform is expressed. XAGE1A/B and CTAG1A/B preserve shared protein identity.

The forbidden union is `unsafe.fasta.gz` plus `additional-forbidden.fasta.gz`
(incomplete/ambiguous translations from Unsafe genes only). Neither contributes
to neither target nor default forbidden sets. Full proteins can still share
k-mers; subtract the forbidden union without crossing ambiguous residues, stops
or record boundaries. Missing complete translations remain coverage gaps.

The 1,721 additional exact reference loci retain their RNA decisions and are
reissued under the current policy. The matching
[v1.20.0 release](https://github.com/pirl-unc/oncoref/releases/tag/v1.20.0) contains
**359 CTA / 3,230 Neither / 126,396 Unsafe** complete exact protein sequences.
Use matching policy/reference manifests, not an older FASTA with the current list.

**Migration:** the historical 624-gene CTA API is unchanged. It is a different
selection, not an alias for this policy. Downstream software must explicitly
adopt `target_sets.results()` or `target_sets.gene_ids("CTA")`, record
`cta_policy.manifest()`, and rebuild its peptide exclusion basis from matching
release assets. A version upgrade alone does not switch the legacy selection.

The current PDF brief is reproducible with `python scripts/build_cta_release_pdf.py`
using the `reports` extras. The separately saved `scripts/slides/` workflow
reproduces the historical 24 October slide deck; it does not relabel it as current.
