# OncoRef 1.25.0: Neither fallback and reproductive protein evidence

**Draft: further admission review is in progress.** The [Other-subtype audit](audits/cta-other-20261030/index.md) found a stale Reproductive bulk input and an omitted reproductive-cell RNA route. Its corrected counts are diagnostic; the packaged decisions below have not yet changed. PR #652 remains a draft and 1.25.0 has not been released.

**573 CTA proteins, 725 Neither, 18,031 Unsafe and 534 unassigned.** All 522
previous CTAs remain included, and Unsafe membership is unchanged. CTA represents
615 member loci after identical selected reference protein sequences are collapsed.

[Full report and funnel](audits/cta-source-walkthrough-20261028/index.md) ·
[Updated PDF slides](audits/cta-source-walkthrough-20261028/cta-selection-1.25.0.pdf) ·
[Changes and individual evidence](audits/cta-neither-20261030/update-1.25.0.md) ·
[Rules and API](cta-selection.md)

| Category | 1.24.0 | 1.25.0 |
|---|---:|---:|
| CTA | 522 | 573 |
| Tumor-restricted | 0 | 0 |
| Neither | 20 | 725 |
| Unsafe | 18,031 | 18,031 |
| Unassigned | 1,290 | 534 |

## What changed

**Neither now means the required normal screen passes but no candidate subtype
qualifies.** High or missing cancer expression no longer prevents Neither. Missing
required normal evidence or a normal-screen hold still prevents that assignment.
Seven old Neither proteins with Manageable-tissue holds return to unassigned.

**Reliable positive reproductive IHC again supplies an independent CTA association
route.** The selection rewrite had retained IHC as an antibody control but omitted
this admission route. The fix restores 51 of 52 reviewed proteins, including SAGE1
and GAGE10. CDY1/B remains Neither because it lacks qualifying cancer expression.
The rule requires reliable mapped positive staining for every member locus and
no unresolved antibody concern. All normal gates and cancer requirements still apply.

**Tumor-restricted is a separate candidate subtype.** It requires ≥10 native TPM
in two distinct tumor samples, or reviewed reliable tumor protein evidence, plus
low measured normal/developmental expression, resolved identity, required normal
coverage and a review of relevant normal cells and states. The existing CTA cancer
rule remains ≥2 TPM or HPA pTPM in one sample.

There are **66 numerical Tumor-restricted leads within Neither**, awaiting normal
coverage review. Forty are olfactory-receptor genes, illustrating why a bulk panel
alone cannot establish restriction from rare normal cells. No protein has yet been
admitted to the new subtype. The report provides their qualifying tumor samples
and coverage questions. These are research classifications, not clinical clearance.

## Downstream contract

```python
from oncoref import cta_selection

cta_selection.results(category="CTA")
cta_selection.results(category="Tumor-restricted")
cta_selection.candidates()  # CTA + Tumor-restricted
cta_selection.results(category="Neither")
cta_selection.audit().explain("SAGE1")
```

The accepted category vocabulary gains `Tumor-restricted`. Existing CTA queries
retain their scope. Unassigned normal holds now use `normal_evidence_hold`; missing
required normal coverage/resolution uses `incomplete_evidence`. Neither does not
imply absent Reproductive or cancer expression. Its profile labels are descriptive.
Pin policy **`oncoref.cta-selection.2026-10-30.1`** and the returned source snapshot.

This release does not silently change the independently versioned exact-protein
background, legacy CTA sets or k-mer exclusion APIs. Downstream code must choose
the intended policy and preserve the normal evidence for shared subsequences.
Full and cached selection replays produce byte-identical decisions. Regression
tests preserve all prior CTAs and Unsafe proteins, enforce IHC identity/quality,
and reject duplicate tumor observations as independent replication.
