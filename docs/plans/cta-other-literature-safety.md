# Other proteins: literature, normal evidence and sequence exclusions

Requested scope: review all 600 proteins classified Other in OncoRef 1.25.0,
search every member gene, preserve reusable literature findings, update the
selection rules and report, then review, merge and publish the release.

## Evidence and decisions

- Keep an immutable baseline and a search receipt for every gene, including
  zero hits, alias limitations and incomplete retrieval. Search results are
  leads, not accepted normal-expression findings.
- Store adjudicated findings with exact gene identity, species, normal cell
  type and activation state, evidence modality, source title/link, location,
  limitations and the resulting policy implication.
- Activated normal immune cells and rare normal neurons contribute to safety.
  Reproductive or cancer expression cannot override accepted normal evidence.
- Separate direct CTA eligibility from sequence-background treatment. Other
  does not automatically mean that every shared peptide is safe; unresolved
  normal evidence must remain explicit.
- Audit the tissues behind all 171 Other:Normal proteins, preserving native
  measurement units and raw versus supported observations. Skin/hair evidence
  must not be confused with neuronal or immune evidence.
- Check every integrated cancer source for CDY1/B and OR2T7. Permit positive
  cancer RNA below 2 TPM as cancer-presence evidence when the independent
  reproductive and normal-safety requirements pass; disclose low abundance.

## Review and release

Record evidence and implementation updates in the PR discussion. Recompute
the full funnel, explain each baseline Other decision and every changed CTA,
verify established CTA controls, exercise the public API and sequence exports,
run relevant tests and reproducibility checks, regenerate documentation and
figures, then merge and publish the tested artifacts. Report search coverage
separately from adjudication coverage and do not claim missing measurements
or a literature search with no hits demonstrates absence of normal expression.
