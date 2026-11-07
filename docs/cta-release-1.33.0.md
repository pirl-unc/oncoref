# OncoRef 1.33.0: CTA paper gaps, concerns and localization

The [CTA report](audits/cta-literature-20261107/report.html) now evaluates all
33 proteins without qualifying paper coverage, summarizes the 29 flagged reviews,
and provides filterable localization evidence for all 602 CTAs.

The 33 gaps span duplicated regulatory families, germ-cell structure, reproductive
secretory families, maternal/oocyte programs and poorly characterized loci. All
pass the existing coding-reference rule. Eighteen enter through developmental
expression alone. Their median highest tumor RNA is 23.5 TPM, compared with
34.1 TPM for the 569 paper-covered proteins. Missing papers alone do not exclude
a candidate; five have specific normal-cell questions requiring adjudication.

The concern table separates direct normal protein, physiological compartments,
normal RNA, unresolved copies, plasma affinity assays, coding conflicts,
proteoform attribution and limited characterization. Five unnamed loci have
little literature; they do not have the experimental coding challenge reported
for FRG2/FRG2B.

Localization combines a checksum-pinned UniProt 2026_03 extract with scoped
primary-paper annotations. Fifty-six protein records have secretion or
extracellular-matrix annotations, 37 have cell-surface or integral
plasma-membrane annotations, and CRIPTO3 has extracellular detection with its
route unresolved. ZP4 occurs in the first two categories. Shared-entry and
partial-member annotations are separate; source evidence and model scope are
retained. These counts include inferred annotations, not only experiments.

New APIs are `cta_literature.gap_assessments()`, `cta_literature.concerns()` and
`cta_literature.localization()`. Each accepts an optional protein label or
Ensembl gene ID. CSV downloads and a comparison figure accompany the tables.
Candidate detail panels link to the localization evidence.

The classification snapshot remains **1.31.0: 602 CTA, 18,255 Unsafe and 1,006
Other**. This release adds evidence and interpretation without changing the
classification, expression datasets or paper-coverage calculation.

Validation covers all 33 gaps and 29 concerns, exact locus mapping, shared
annotations, organelle-versus-surface distinctions, evidence types, report
filters, offline reproduction, strict documentation builds and figure inspection.
