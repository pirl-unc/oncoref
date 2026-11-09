# Scrutiny of the 250 Cancer Antigens only

**250 proteins pass the current expression rules. This is not a count of experimentally validated cancer antigens.** The category and its comparison methods require further review before release. PR #686 is draft.

## What is in the set?

| Family | Proteins |
|---|---:|
| Olfactory receptors | 141 |
| Other families | 50 |
| Keratin-associated proteins | 31 |
| Other sensory receptors | 11 |
| Unnamed Ensembl loci | 7 |
| SCYGR | 6 |
| USP17 | 4 |

The family counts use explicit symbol patterns, retained in the audit script. They are a descriptive breakdown, not a new annotation ontology.

## Tumor-expression concerns

- **22** have no measured tumor above 5 TPM; **69** have none at 10 TPM or higher.
- **53** qualify through the single-sample route without meeting the confirmed two-patient route. This does not mean only one measured tumor expresses them: patient identities are unresolved in several cohorts.
- **73** have their maximum in Hodgkin lymphoma (HL), including **69 in SRS3811196**. Another **64** peak in the Merkel cohort. These are concentrations of maxima, not proof of technical artifacts.
- Removing HL alone makes **21** lose qualification; removing Merkel alone makes **7** lose it. Removing both makes **31** lose it, leaving **219** under the otherwise unchanged rules. This is a sensitivity analysis, not a corrected final count.
- The same two-cohort removal affects **9 CTAs**. The concern applies to the shared cancer gate, not only the nonreproductive branch.

HL uses recount3 coverage-derived gene TPM; Merkel uses symbol-keyed gene counts divided by reference gene lengths. The source registry calls both `linear_rnaseq_tpm`, but that does not establish calibrated comparability with the normal reference. The comparison also takes tumor-sample maxima against tissue/cell summaries, including donor medians. It does not compare against the highest individual healthy donor.

The [HL study](https://pubmed.ncbi.nlm.nih.gov/30546079/) investigates LTR-driven transcription and includes laser-captured HRS cells. The [Merkel source](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE235092) is a bulk-RNA study. Bulk expression does not by itself identify which cell produced the RNA. None of these source descriptions proves that the suspicious values are false.

Olfactory receptors particularly need transcript-level scrutiny. A [primary RNA-seq study](https://pmc.ncbi.nlm.nih.gov/articles/PMC3566163/) reports cases, including OR5K2, where broad read coverage could be misinterpreted as receptor transcription. OR5K2 is in this set. That is a reason to inspect its tumor alignments, not evidence that this tumor measurement is itself an artifact. Exact gene-label mapping does not establish uniquely assigned reads or an intact coding transcript.

## Which are previously nominated?

The earlier **47** flag was inflated by assigning a source-wide nomination role to every member of a mixed paper intake. The da Silva paper has a 1,103-gene testis-biased background and a distinct 745-gene tumor-positive list. These are now separated at the exact source-table and gene level.

**24 of the 250 have exact catalog-nomination support** in the retained sources: 22 in [da Silva et al.](https://www.oncotarget.com/article/21715/text/) and five in [Wang et al.](https://www.nature.com/articles/ncomms10499), with three overlapping. These are computational expression-based nominations, not 24 demonstrations of tumor protein production or antigen presentation.

DBX1, BHLHE23, FUT5, PANX3, OR1K1, OR2D2, OR51Q1, TEX13B, OR2K2, OR5B12, HES3, OR11H4, OR52E2, OR51A7, TCHHL1, OR52L1, BPIFB3, OR51I2, OR13G1, OR5B21, OR13C3, OR1J4, LCE3C, OR5M11.

All **51** proteins with any retained paper link keep those links, including the 23 background-only records whose nomination flag was removed. A missing link in this curated collection does not prove that no relevant paper exists.

## Normal and developmental coverage

The 141 olfactory receptors, 31 keratin-associated proteins and other sensory receptors concentrate the question in specialized normal compartments. Ordinary tissue summaries can miss rare sensory cells. The reproductive/developmental gate records qualifying evidence in the imported sources; it does not establish absence from every fetal organ or developmental stage. Disposable allowances also deliberately permit some normal expression.

## Required follow-up

1. Review source-specific comparability, raw-count support, transcript coverage and mapping uniqueness for the shared cancer gate. Keep valid rare high expressors when their measurements are supported.
2. Audit rare normal sensory cells, skin appendages, and developmental compartments with exact-gene data. Apply the same normal-evidence rules to both cancer branches.
3. Resolve patient identities from original sample metadata, without treating aliquots as independent donors.
4. Keep source-table nominations separate from tumor protein, immune recognition and therapeutic evidence. Missing or conflicting required evidence belongs in Other:Review; a family name or a suspicious cohort alone does not establish Unsafe.

[All 250 support records](cancer-antigen-only-support.csv) · [Cohort-removal analysis for all Cancer Antigens](cohort-sensitivity-all-cancer-antigens.csv) · [Audit summary](summary.json) · [Catalog table](../rebuild/report.html).
