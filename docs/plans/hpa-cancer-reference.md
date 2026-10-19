# HPA cancer reference contract (#550)

- [x] Acquire fresh genome-wide HPA cancer IHC and sample RNA tables, with the
  provider's explicit release label, download date, raw hashes, and Ensembl
  version. Archive immutable copies. Preserve the old tsarina CSV hashes with
  their release marked unknown; do not assign them today's release.
- [x] Add independently downloadable reference sources and public IHC/RNA APIs.
  Keep TCGA and validation cohorts separate. IHC counts refer to patients;
  RNA retains native pTPM, finite measured denominators, threshold-positive
  counts, sums/means/maxima, and explicit missing-measurement counts.
- [x] Derive genome-wide summaries using bounded-memory processing. Reject
  duplicate identities, negative/invalid counts, and unrecognized cohort labels.
  Unknown IHC categories invalidate the denominator; zeros with a measured
  denominator remain genuine negatives. Preserve zero-denominator rows as gaps.
- [x] Ship an explicit 20-group crosswalk. Pool colorectal, lung, and renal RNA
  by measured sample counts only when every required member is present. Report
  glioma/GBM as a scope mismatch and carcinoid/lymphoma/skin as unmatched. Expose
  unpaired-cohort and antibody-resolution limitations separately from values.
- [x] Test missingness, count conservation, denominators, duplicate detection,
  source identity/checksums, weighted pooling, cohort separation, and all mapping
  outcomes. Rebuild the source summaries and six RNA/IHC comparison plots.
- [x] Document the API and migration from tsarina, open the oncoref PR, and make
  its artifacts reproducible from pinned sources. Downstream migration follows
  release of the new oncoref contract.

The new contract will not infer stained-cell percentages, paired RNA/protein
relationships, or per-antibody reliability from the aggregate cancer IHC table.
The existing normal-tissue reference version remains independent.

Implementation: [PR #552](https://github.com/pirl-unc/oncoref/pull/552).
Supporting sources: [HPA cancer v25.1-1](https://github.com/pirl-unc/oncoref/releases/tag/hpa-cancer-v25.1-1).

Validated all 17 uploaded assets against GitHub SHA-256 digests. The built wheel
loaded both runtime references through fresh public downloads and passed API
checks. Local suite: 1,450 passed; strict docs and wheel/sdist builds passed.
The tsarina migration remains a follow-up after the package release.
