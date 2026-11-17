# Curation review — 17 November 2026

Policy `oncoref.target-selection.2026-11-17.1`, released in OncoRef 1.41.0.
[Category transitions from 1.40.1](category-changes.csv),
[data-only baseline](prior-receipt.json) and [independent checks](verification.json)
record the comparison.

Validation: 132 targeted tests passed, the strict documentation build and lint
checks passed, and the build verifier checked 204 generated outputs. Funnel,
Venn and threshold figures were visually reviewed. Wheel and source-distribution
metadata passed Twine validation.

## Small embryo populations

The twenty-cell floor excluded all sampled eight-cell embryos and many later
embryo lineages. The replacement has two routes: at least twenty cells, at least
five positive and at least 10% positive in one embryo; or at least five cells and
50% positive, replicated in two independent embryos of the same stage and lineage.
Positive means at least 1 native TPM. Each supporting embryo mean must reach the
selected fold. The smaller-population route also applies to sparsely sampled later
lineages; it is not a gene-specific rescue or an arbitrary day-three exception.

The supported expression ceiling uses the second independent embryo where
replication is required. It is distinct from the largest observed mean. This
restores PRAMEF1/2, TRIM43/49, CSAG1, SSX2/B and PRDM14 reproductive membership.
No normal ceiling was removed to obtain these assignments. The underlying human
preimplantation study spans the eight-cell stage through later blastocysts
([Petropoulos et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC4868821/)); the
build uses the separately pinned native-TPM reprocessing described in the source receipt.

## Normal-library and locus review

All 186 HPA individual tissue samples were screened within their tissue. Ten
normal testis libraries carry the expected spermatid signature; that is not a
normal-library failure. No actionable signature flags were found, including in
lung.167, thyroid gland.250, thymus.1599 and thymus.1603. The
[library ledger](../target-selection-20261110/rebuild/normal-library-qc.csv)
preserves every result. This tumor-derived screen is an audit aid, not a validated
normal-library quality model. Thymic expression remains in the denominator.

CAGE1/host-transcript overlap and CTAG2/paralog ambiguity remain attribution
hypotheses. Negative bulk RNA can dilute a sparse positive population, and negative
IHC has assay and sampling limits. Neither resolves those findings by itself.
Rare-cell summaries now distinguish replicated attribution uncertainty from
frequent single-study uncertainty and prioritize replicated findings at equal status.

## Names and consequences

RA retains its broad reproductive scope and CTA remains CA ∩ RA. Tissue-loss
consequences are a separate, cited axis: **Vital, Protected, Manageable, Disposable**.
The report shows acute/persistent consequences through the API, and persistent loss
and management in the table, including observed reproductive RNA below the tenfold
threshold. Whole-testis/ovary data retain Protected organ assessments; breast uses
its reviewed conditional assessment. Unknown embryo/cell scopes remain unresolved.
Anatomical absence is not inferred from a sex label, and no antigen is automatically
promoted on that basis.

## Verification and remaining work

Independent calculations reproduce all 918 exported patient-fold lists and all
918 assessed reproductive ceilings. The 18,253 Unsafe identities match 1.40.1
exactly. The changes are 43 CA-only → CTA, one Other → CTA and eleven Other →
RA-only. OR1N1 is the Other → CTA transition: its newly supported RA evidence
satisfies the existing home-tissue rule; the rule itself was not changed.

The new CLI and build runner share the existing owners. CI verifies their artifact
receipts. This update does not claim completion of a direct raw-normal-screen
replacement, wholesale legacy-module removal, a unified adjudication ledger, new
native-TPM germ-cell/oocyte panels or read-level attribution research. Those need
separate validation; [migration boundaries](../../curation.md#migration-boundaries)
are explicit. Existing historical receipts remain available without adding frozen
Python-code copies for this release.
