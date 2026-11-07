# OncoRef 1.32.0: individual CTA literature and source roles

This release adds a centralized literature registry and a [filterable review](audits/cta-literature-20261107/report.html)
of the 41 previously unmatched proteins and 55 previously found only in a broad
screening table.

All 96 have individual assessments. Seventy-six have linked primary or related
papers; 65 have exact-human studies credited under the stated coverage rule.
The assessments distinguish gene-specific studies, genomic/RNA associations,
normal-expression findings, identity or translation conflicts, and related-species
or family evidence.

The source collection separates 14 core catalog/program papers and 68 focused
studies from 68 supplementary sources. Every source has a role, a contribution,
linked proteins and an overlap explanation. Twelve of the original 15 papers
preserve all 2,514 mapped proteins; eight preserve their 480 released CTAs.
Complete Chew and Smith program tables retain source occurrences, unmatched rows,
historical identities and upstream gene-set provenance. The 12,697-gene Chew
screening table does not count as a nomination list.

The expanded paper links cover 568 of the 602 released CTA records. The 34 remaining
coverage gaps are listed explicitly. Named genomic associations count as paper
coverage but do not establish protein production.

`oncoref.cta_literature` exposes the packaged registry, source roles, 96 reviews,
paper memberships and coverage for all 602 CTAs. The candidate explorer links
each reviewed protein to its new findings and publications.

The **classification snapshot remains 1.31.0**: 602 CTA, 18,255 Unsafe and 1,006
Other. Eighteen reviews identify normal-expression questions and eleven identify
identity or translation questions. Those findings require individual adjudication
and a classification rebuild; this evidence release does not mark them resolved.
Expression data bundle and source-matrix versions are unchanged.

Validation includes source-row and identity preservation, core coverage and
redundancy checks, complete review/search coverage, offline reproduction,
packaged-data access, report filters and figure inspection.
