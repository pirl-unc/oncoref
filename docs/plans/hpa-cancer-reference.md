# HPA cancer reference contract (#550)

- [ ] Acquire fresh genome-wide HPA cancer IHC and sample RNA tables, with the
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
- [ ] Document the API and migration from tsarina, open the oncoref PR, and make
  its artifacts reproducible from pinned sources. Downstream migration follows
  release of the new oncoref contract.

The new contract will not infer stained-cell percentages, paired RNA/protein
relationships, or per-antibody reliability from the aggregate cancer IHC table.
The existing normal-tissue reference version remains independent.
