# Curation summary

Policy `oncoref.target-selection.2026-11-16.1`, 19,893 human protein records. **437 CTA, 34 Cancer Antigen only, 197 Reproductive Antigen only.** These are outputs of screening rules, not validated tumor antigens. Every figure is drawn from the [rebuilt catalog](index.md); thresholds are read from the policy constants.

## Antigen membership

![Antigen sets](summary/antigen-sets.png)

CTA is the intersection. Reproductive-only proteins have no qualifying cancer evidence; the union name does not certify tumor targeting. Unsafe is a separate subtraction background.

## Rule-derived sets and published nominations

![Set overlap](summary/venn.png)

Of 437 CTA, 302 are published nominations and 135 are not. 1,820 published nominations are in neither set: 1,761 of all 2,233 nominations are Unsafe, 35 Other:Review and 7 Other:Sequence. A published nomination is an annotation; it never overrides an exclusion.

## Funnel by filtering idea

A protein leaves at the first idea that removes it, so later steps show only what each idea adds. The order is chosen so that every step removes something; the layers overlap heavily, as the next section shows.

![All proteins](summary/funnel-all.png)

| Filtering idea | Removed here | Remaining |
|---|---|---|
| All assessed proteins | 0 | 19,893 |
| Bulk tissue RNA above the tissue ceiling | 17,522 | 2,371 |
| Cell-type cluster RNA (HPA single cell) | 424 | 1,947 |
| Sorted immune-cell RNA | 96 | 1,851 |
| Rare-cell RNA established by donor | 107 | 1,744 |
| Normal-tissue protein (IHC) | 50 | 1,694 |
| Curated or published normal expression | 49 | 1,645 |
| Adult cornea RNA | 5 | 1,640 |
| Incomplete coverage, identity or translation | 660 | 980 |
| Nested or Unsafe-contained sequence | 62 | 918 |
| No qualifying expression or bounds fail | 167 | 751 |
| Home tissue no panel resolves | 83 | 668 |

![Published nominations](summary/funnel-published.png)

| Filtering idea | Removed here | Remaining |
|---|---|---|
| Published CTA / oncofetal nominations | 0 | 2,233 |
| Bulk tissue RNA above the tissue ceiling | 1,464 | 769 |
| Cell-type cluster RNA (HPA single cell) | 174 | 595 |
| Sorted immune-cell RNA | 25 | 570 |
| Rare-cell RNA established by donor | 68 | 502 |
| Normal-tissue protein (IHC) | 19 | 483 |
| Curated or published normal expression | 10 | 473 |
| Adult cornea RNA | 1 | 472 |
| Incomplete coverage, identity or translation | 35 | 437 |
| Nested or Unsafe-contained sequence | 7 | 430 |
| No qualifying expression or bounds fail | 7 | 423 |
| Home tissue no panel resolves | 10 | 413 |

Steps that remove no published nomination are omitted from the second funnel. [All proteins](summary/funnel-all.csv) · [Published](summary/funnel-published.csv)

## Normal-expression layers

![Unsafe layers](summary/unsafe-layers.png)

| Layer | Unsafe proteins flagged | Only reason |
|---|---|---|
| Normal bulk | 17,522 | 560 |
| HPA tissue clusters | 17,053 | 268 |
| Adult cornea | 12,095 | 5 |
| DICE immune | 11,632 | 45 |
| Monaco immune | 11,264 | 19 |
| HPA normal protein | 10,157 | 47 |
| HPA immune | 7,973 | 14 |
| Established rare normal-cell RNA | 129 | 106 |
| Reviewed normal-cell knowledge | 32 | 31 |
| Published normal-tissue expression | 19 | 18 |
| Accepted normal-protein evidence | 8 | 2 |

## Thresholds

![Thresholds](summary/thresholds.png)

Ceilings are tier base × 2 per negative-IHC reliability step (Vital 4, Protected 8, Manageable 16 without IHC credit). Exceeding a Vital or Protected ceiling is Unsafe; Manageable is Unsafe only above 128. These are screening conventions, not measured injury thresholds.

![Ceilings applied](summary/thresholds-applied.png)

The second table counts recorded bulk-RNA peak observations, not proteins: one protein can contribute a peak per source.

## Expression distributions and where the thresholds fall

![Distributions](summary/distributions.png)

One value per protein on each source's own scale; scales are never pooled. Counts use a compressed axis so the small candidate group stays visible.

## Cancer rule

![Tumor against normal](summary/tumor-normal.png)

## Rare-cell RNA

![Rare-cell grading](summary/rare-cell.png)

[Donor-level evidence](../../normal-cell-atlas-20261110/rare-cell/index.md)

## Inside the later stages

![Stage detail](summary/stages.png)

- Required normal coverage incomplete: 400
- Identity or translation unresolved: 102
- Normal signal blocks admission (not Unsafe): 158
- Patient identity unresolved in tumor samples: 0

A protein can use several routes, so route counts exceed set sizes.
