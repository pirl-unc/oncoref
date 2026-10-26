# OncoRef 1.20.0: review and adoption

The [current policy report](audits/cta-policy-20261025/index.md),
[178 completed antibody adjudications](audits/cta-adjudication-20261025/index.md)
and [plain-language definitions](target-sets.md) describe the released research
selection. [PR #639](https://github.com/pirl-unc/oncoref/pull/639) records the
implementation, review, CI and merge trail.

[Four-page PDF](https://github.com/pirl-unc/oncoref/releases/download/v1.20.0/cta-policy-1.20.0.pdf) ·
[Release and exact sequence assets](https://github.com/pirl-unc/oncoref/releases/tag/v1.20.0)

## Reviewed outcome

| Scope | CTA | Neither | Unsafe |
|---|---:|---:|---:|
| 2,540 nominated genes | 213 | 741 | 1,586 |
| 20,097 canonical genes | 213 | 2,179 | 17,705 |
| 129,985 complete exact proteins | 359 | 3,230 | 126,396 |

The eight initial literature cases and all 178 remaining raw-antibody cases
have explicit curation dispositions. Accepted normal protein excludes; material
unresolved normal evidence holds. No antibody-discrepancy case remains CTA.
Raw gene assignments, shared family reagents and repeated source observations
are not counted as independent replication. No independent reread of every
histology image or experimental specificity resolution is claimed.

PRAME is Unsafe after accepted normal adrenal-protein evidence. MAGEA4 is
Neither because its raw normal stains remain unresolved. Their peptide/receptor
products retain their clinical evidence; these gene-level decisions do not claim
that those therapies cause injury. CTAG1A/B and XAGE1A/B remain CTA. CTAG2 remains
Neither/Variable. PLCZ1's historical protein interpretation was corrected to
qualitative RNA: the original pancreatic immunostaining detects CFTR.

## Verification performed

- Full local test suite: **2,489 passed, 6 skipped**, including historical checksum
  preservation and current package/report agreement.
- Focused policy, RNA/IHC, identity, provenance and sequence suite: **307 passed**.
- Independent `scripts/verify_cta_release.py`: all **129,985** complete sequence
  hashes and memberships, Unsafe > Neither > CTA precedence, XAGE/CTAG shared
  proteins, and **32,091** additional forbidden translation records verified.
- No audited non-permitted/unsited normal observation above 30 TPM or nTPM
  remains in CTA or Neither. Units are checked separately.
- Wheel import from an extracted artifact: 213 current CTA genes, 624 historical
  CTA genes, 178 adjudications and 196 reviewed findings; packaged hashes match.
- Source archive identity checks: **46 passed, 8 skipped**; repository-local
  `pip check`, Ruff and strict MkDocs builds passed.
- Four-page PDF rendered and visually inspected. Its generator is
  `scripts/build_cta_release_pdf.py`; historical slide reproduction remains in
  `scripts/slides/` with its original policy clearly labelled.

CI results for supported Python versions and Windows asset checks are recorded
on the PR. Release publication and documentation deployment are independently
visible in GitHub Actions. Source files remain pinned; no shared environment was
modified.

## Downstream migration is explicit

The historical `cta.cta_gene_ids()` selection remains **624 genes**. A package
upgrade alone does not replace it with the current research list. Choose the
intended policy deliberately:

```python
from oncoref import cta_policy, target_sets, normal_rna_review

manifest = cta_policy.manifest()
assert manifest["policy"] == "oncoref.target-sets.2026-10-25.2"
current_targets = target_sets.gene_ids("CTA")
forbidden_genes = normal_rna_review.forbidden_gene_ids()
findings = cta_policy.findings()
antibody_reviews = cta_policy.adjudications()
```

For exact proteins, use matching v1.20.0 release assets and SHA-256 values from
`manifest["sequence_artifacts"]`. The forbidden sequence union is
`unsafe.fasta.gz` plus `additional-forbidden.fasta.gz`; Neither is excluded from
both targets and this default forbidden background. Rebuild downstream caches
and k-mer exclusions with the matching reference and policy. Do not substitute
all non-CTA genes for Unsafe. Full-protein disjointness does not prevent shared
k-mers.

Coordinate exact OncoRef pins across consumers and validate their joint
installation before changing a shared environment. Tumor-marker panels and
product-specific therapy evidence have different purposes and should not be
silently reduced to the unrestricted gene-level target set.

## Remaining scientific work

[Continuous breadth (#637)](https://github.com/pirl-unc/oncoref/issues/637) and
[comprehensive individual literature assessment (#638)](https://github.com/pirl-unc/oncoref/issues/638)
remain open. The initial search covered all candidates but did not constitute
full-paper review of all 2,540 genes. The adjudication records specify the new
evidence needed for held candidates. These limitations do not make a held gene
an unfinished curation decision, a safe target, or proven absent from normal cells.
