# Genome-wide normal-expression categories

Policy `oncoref.normal-expression-categories.2026-10-21.2`; 20,097 canonical human protein-coding genes, pinned HPA v23. This universe includes genes outside published CTA nominations.

| Category | Other subgroup | Genes |
|---|---|---:|
| Normal | — | 18,063 |
| Other | Disposable | 29 |
| Other | Silent | 372 |
| Reproductive | — | 1,023 |
| Unassessed | — | 610 |

## Definitions

- **Reproductive:** positive reproductive RNA >0 or IHC, complete routine RNA, and all approved absolute RNA/selectivity/IHC gates. Publication nomination and manual CTA review are independent: this is broader than the retained curated CTA set.
- **Normal:** measured expression fails an absolute ceiling, the positive-IHC gate, or (when reproductive-positive) the pairwise selectivity gate. It means normal-expression counterevidence under this policy, not demonstrated clinical toxicity.
- **Other:** no positive reproductive evidence, complete routine RNA, no disqualifying positive IHC, and all approved absolute category/IHC ceilings pass. No reproductive ratio is required. **Disposable** has a measured positive in Disposable tissue; **Silent** does not. Both mean **below-threshold expression**, not literal absence from every normal tissue.
- **Unassessed:** no decisive normal-expression counterexample, but missing routine RNA, unresolved anatomy, unusable/reviewed-unreliable evidence, or an undefined required reproductive ratio prevents classification into Reproductive/Other. It is a separate evidence status and never part of Silent.

Missing/partial IHC retains baseline RNA limits and is explicitly recorded; this does not establish protein absence. The 50 RNA/80 IHC routine contexts do not cover all developmental stages or anatomical compartments. Numerical limits apply only to human HPA v23. [Exact category/IHC grid](../../normal-tissue-scopes.md).

## Preparing for subsequence comparisons

The three main categories are disjoint research screening groups. Sequence membership is a separate question: a reproductive gene can share a subsequence with a normal paralog, and normal RNA below a screening limit is still observed expression.

`normal_expression_observed` preserves any valid RNA >0 or positive IHC in Manageable/Protected/Vital tissue, even below a cutoff. The conservative normal-reference export includes those genes, all Normal genes, and Unassessed genes. It can overlap Reproductive and Other. A later subsequence comparison must examine exact protein/transcript occurrences and their reference versions across both sets; this report does not establish sequence specificity. Disposable expression remains separately visible.

## Downloads

- [All genes, categories, coverage and reasons](gene-categories.csv.gz)
- [Reproductive](reproductive.csv.gz) · [Normal](normal.csv.gz) · [Other union](other.csv.gz)
- [Disposable subgroup](disposable.csv.gz) · [Silent subgroup](silent.csv.gz) · [Unassessed](unassessed.csv.gz)
- [Conservative normal sequence-reference candidates](conservative-normal-reference.csv.gz)
- [Counts](summary.csv) · [Input and implementation checksums](run-manifest.json)

The [CTA recalibration report](../cta-recalibration-20261021/index.md) applies publication, family and nomination reviews in addition to the Reproductive expression screen. The current doubled absolute ceilings preserve every previously passing Reproductive expression decision. `retained_recalibration_cta` records that pre-external-review stage; use the [current clinical curation decisions](../cta-curation-20261021-r2/index.md) for the final retained set and external safety holds.

Reproduce with `python scripts/normal_expression_categories_report.py --out <empty-directory> --workers 4`.
