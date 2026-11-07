# CTA literature and source roles

[Open the filterable report](report.html). The review covers all **96 proteins**: 41 previously unmatched and 55 previously found only in a broad screening table.

76 have linked primary or related papers; 65 have an exact-human study credited under the declared coverage rule. Evidence types are shown separately: a DNA association is not a protein experiment.

The harmonized collection has 83 core papers and 67 supplementary sources. Its mapped coverage reaches 569 of the 602 released CTA records. 33 remain without qualifying paper coverage. The expression and reference datasets remain the basis of selection.

The 1.31.0 classifications are frozen. Normal-expression and identity findings in this review require adjudication before the next classification release.

- [Individual reviews](individual-reviews.csv)
- [Core and supplementary sources, roles and overlap](sources.csv)
- [Paper coverage for all 602 CTAs](cta-paper-coverage.csv)
- [All 33 paper gaps: biology, measurements and assessment](paper-gap-assessments.csv)
- [All 29 normal-expression and identity concerns](review-concerns.csv)
- [Localization of all 602 CTAs, with evidence and mapping scope](cta-localization.csv)
- [Complete registry, including source-row occurrences](registry.json.gz)
- [Search queries and retrieval receipts](search-ledger.json)
- [Figure: review outcomes](review-outcomes.svg)
- [Reproduction manifest](manifest.json)

Rebuild offline with `.venv/bin/python scripts/build_cta_literature_review.py`. Refresh indexed searches with `scripts/review_unlisted_cta_literature.py` and cited records with `scripts/fetch_cta_individual_review.py`; neither script automatically accepts a finding. The pinned workbook importer is `scripts/import_cta_program_tables.py`.
