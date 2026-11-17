# OncoRef 1.40.0: one enrichment rule for both antigen sets

**Cancer Antigens and Reproductive Antigens independently require at least 10×
the same measured non-reproductive normal reference.** CTA is their intersection;
their union is called candidate antigen proteins.

[Curation rules](target-selection.md) · [Funnel and Venn figures](audits/target-selection-20261110/rebuild/summary.md) ·
[Searchable catalog](audits/target-selection-20261110/rebuild/report.html) ·
[Optional filters](target-filters.md) · [Sequence APIs](antigen-sequences.md)

| Category | 1.39.0 | 1.40.0 |
|---|---:|---:|
| CTA | 437 | 246 |
| Cancer Antigen only | 34 | 149 |
| Reproductive Antigen only | 197 | 102 |
| Other | 972 | 1,143 |
| Unsafe | 18,253 | 18,253 |

There are **395 Cancer Antigens, 348 Reproductive Antigens and 497 candidates**.
Counts refer to identical selected reference protein groups, not gene loci.
[Every change from 1.39.0](audits/target-selection-20261110/rebuild/changes-from-1.39.csv).

## Shared comparison and independent normal evidence

Both positive compartments must reach `10 × max(N, D, 1 TPM)`: N is the supported
non-reproductive normal maximum and D the Disposable maximum. CA requires one
QC-passing tumor sample. RA requires one reproductive tissue sample or a supported
mean within one embryo and lineage. The latter requires at least 20 cells, five
cells at 1 TPM or more, and 10% positive cells; the mean itself must reach 10×.
Multiple cells or lineages from one embryo are not independent donors.

The reference now includes 186 individual HPA v23 tissue samples, including the
non-reproductive samples in the denominator. It also retains tissue/subregion
summaries, GTEx medians and supported immune summaries. This is not a maximum over
every healthy donor, and cross-study batch effects remain. Missing coverage is
not zero. Protein staining, UMI counts and nTPM cannot bypass the native-TPM fold.

Established normal-cell RNA and accepted normal protein evidence independently
exclude through the shared normal screen. Strong signals with unresolved locus
attribution or confirmation become Other:Review. A mapping-risk flag alone does
not establish artifact. This holds DAZL, DCST1 and NKX2-6 for review; PRDM14,
MAGEC2, PRAME and MAGEA4 remain CTAs. MAGEC2 and DAZ findings retain the shared
curation records and source-specific rationale introduced in 1.39.0.

Other distinguishes Review, Sequence, Tissue and Expression exclusions. The
existing Unsafe set is unchanged. Read-level adjudication and the hair-follicle
replay (#687) remain open; `release_ready: false` explicitly describes scientific
review completion, not software publication status.

## API and reporting changes

- Catalog rows expose actual 10× sample, known-patient, anonymous-sample and
  cancer-type counts, with every qualifying witness exported. The older
  above-normal counts remain explicitly weaker descriptors.
- `TargetFilter` adds `min_reproductive_fold` and `min_tumor_patients`. The latter
  counts known patients at the catalog's 10× threshold, even when a higher fold
  is selected separately. Publication counts are not independent replication.
- `antigen_sets` adds `Other` and `Unsafe|Other` backgrounds for exact or optional
  bounded-distance 8–15mer subtraction. The strict MAGE-prefix rule is preserved.
- Guides, tables, funnels, Venn diagrams and paired enrichment plots derive from
  the same policy. Historical 1.38 generators reject current-policy replay.

Package **1.40.0**, target policy `oncoref.target-selection.2026-11-16.2`.
The expression bundle remains **5.23.28** and raw source matrices **5.22.14**;
this release changes packaged curation and its API, not those expression artifacts.
