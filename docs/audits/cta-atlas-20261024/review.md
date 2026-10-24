# Review record

This report preserves policy `oncoref.target-sets.2026-10-24.1` and the v1.18.0 protein assets. It adds evidence, source-linked flags and navigation. It does not silently clear or exclude genes on the basis of a new annotation.

The implementation review checked exact gene identity, source-native units, measured denominators, missing-data handling, shared UniProt entries, literature scope and PDF pagination. Clinical registry links remain separate from clinical outcomes. Cancer RNA and IHC are not assumed paired, and aggregate staining is not converted to an antibody confidence level.

Review found and corrected two report issues: positive fractions below 0.5% must not display as 0%, and gene headings must remain with their evidence block. The first PDF render revealed an orphaned heading; the renderer now keeps each heading and entry together. PDF bookmarks and the page index are checked against extracted text.

The literature review identified 37 retained-gene conflicts, including primary neural protein reports for ABCC11 and MAGEC2. These have prominent warnings, citations and follow-up actions. The frozen 399-gene count must not be read as an endorsement to target these genes before adjudication. Other flags distinguish isoform, low-RNA, shared-locus and disease-adjacent contexts.

This is an implementation-agent review, not an independent human approval or a systematic review of every gene's entire literature. Automated checks and the PR record provide the final validation status.

Local validation: 67 focused annotation/catalog/cache/dataset tests passed; changed-file Ruff checks passed; strict MkDocs build passed; all 2,540 gene bookmarks/page links matched their entry text; the browser verified counts, filters, source evidence and desktop/mobile layouts without JavaScript errors. Wheel and sdist builds passed Twine's strict checks. The repository-wide Ruff invocation also reports an existing unrelated unused variable in `scripts/regenerate_plots.py`; no new lint finding was introduced.
