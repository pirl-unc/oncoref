# OncoRef 1.37.0: antigen sets and candidate filters

**This release adds a set-based antigen catalog and a filter layer on top of it. The
1.36.0 classification snapshot, sequence bundle and expression data (5.23.27) are
unchanged.** These are outputs of screening rules, not validated tumor antigens.

## Antigen sets

[Rules](target-selection.md) · [catalog and funnel](audits/target-selection-20261110/rebuild/index.md)

Policy `oncoref.target-selection.2026-11-10.1` assigns 19,893 proteins after one shared
normal-expression, coverage, identity, expression-bound and sequence screen.

| Category | Proteins |
|---|---:|
| CTA (both sets) | 626 |
| Cancer Antigen only | 245 |
| Reproductive Antigen only | 31 |
| Other | 750 |
| Unsafe | 18,241 |

The 902 proteins in either set are the Target candidates. `oncoref.target_selection`
exposes them. `oncoref.cta_selection` still returns the 1.36.0 catalog.

## Candidate filters

[Rules](target-filters.md) · [web summary](audits/target-filters-20261112/report.html) ·
[full report](audits/target-filters-20261112/index.md)

`oncoref.target_filters` narrows the 902 candidates. It never changes a category and
cannot add a protein. Six layers run in a fixed order: off-target RNA, warnings,
curation, set evidence, papers and exclusions by name.

| Filter mode | Kept | CTA | Cancer Antigen only | Reproductive Antigen only |
|---|---:|---:|---:|---:|
| `all` | 902 | 626 | 245 | 31 |
| `loose` | 473 | 442 | 25 | 6 |
| `standard` | 297 | 290 | 3 | 4 |
| `strict` | 145 | 145 | 0 | 0 |

PRAME, NY-ESO-1 and MAGEA4 pass every mode, and a test enforces it.

## Curation records used by the filters

[Literature sweep](audits/target-filters-20261112/normal-function-sweep.md)

`oncoref.target_curation` holds 256 cited normal-function findings on 207 candidates,
three acting family rules (olfactory receptors, bitter taste receptors, hair keratins)
and 21 read-through loci.

Limits:

- The findings came from an automated literature sweep. Each was checked by one second
  automated reader, not by a person, and mostly at abstract level.
- Seven tissue tiers (olfactory epithelium, taste buds, nasal mucosa, lacrimal gland,
  sensory ganglia, bone, dental pulp) are cited proposals in
  `oncoref/data/target-curation.yaml`. They are not in the hazard registry.
- These records are separate from the catalog's own literature review
  (`normal-expression-literature.json`) and act only in the filter layer. Folding them
  into the catalog would change categories and needs a new catalog policy.
- The `strict` MAGE rule is a judgment by family name, not a measurement.
