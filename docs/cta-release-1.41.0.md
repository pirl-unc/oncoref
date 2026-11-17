# OncoRef 1.41.0: supported embryo evidence and current antigen commands

Small embryo populations can now support reproductive-antigen membership without
requiring twenty sampled cells from an eight-cell embryo. A population with at
least five measured cells and at least 50% positive cells needs two independent
embryos at the same stage and lineage, each reaching the tenfold threshold.
The existing route for populations with at least twenty cells remains available.
Higher reproductive-fold filters preserve this replication requirement.

The resulting catalog has **290 CTA, 106 Cancer Antigen only, 113 Reproductive
Antigen only, 1,131 Other and 18,253 Unsafe proteins**. Relative to 1.40.1,
43 CA-only proteins become CTA, one Other protein becomes CTA and eleven Other
proteins become RA-only. Unsafe membership is unchanged. CTA remains CA ∩ RA;
the candidate union contains 509 proteins. These labels describe evidence and
membership, not demonstrated therapeutic safety.

The release adds:

- `oncoref antigens list`, `explain`, `funnel` and `subtract`, using the current
  catalog and existing exact, equivalence and weighted-edit peptide matching.
- Exact composition of tumor-fold and patient-count filters: two patients at
  30× means two known patients each reaching 30×.
- Tissue-loss consequences alongside reproductive expression, including exposure
  below the RA threshold. Bulk gonadal measurements retain organ-level consequences;
  reproductive membership and sex labels do not waive injury risk.
- An audit of all 186 individual HPA normal libraries. The signature screen finds
  no actionable flags; their measurements remain in the normal reference.
- One declared curation build graph, checksum verification in CI, content-addressed
  HPA/embryo inputs and data-only category comparisons.

CAGE1, CTAG2 and other unresolved normal-cell signals remain under review.
Host overlap, negative bulk RNA and negative staining do not by themselves establish
that a sparse signal is an artifact. The scientific `release_ready` flag remains
false while those attribution reviews are open; it does not block software publication.

Expression bundle **5.23.29** and raw source matrices **5.22.14** are unchanged.
The rebuild re-extracted individual normal and developmental observations and
regraded pinned rare-cell evidence, using the checksum-verified tumor cache.

[Review and verification](audits/curation-review-20261117/index.md) ·
[Current rules](target-selection.md) · [Commands and build workflow](curation.md) ·
[Funnel and Venn report](audits/target-selection-20261110/rebuild/summary.md) ·
[Filterable catalog](audits/target-selection-20261110/rebuild/report.html)
