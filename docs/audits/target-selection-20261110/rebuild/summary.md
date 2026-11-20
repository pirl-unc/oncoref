# Curation summary

Policy `oncoref.target-selection.2026-11-20.1`, 19,893 human protein records. **327 CTA, 70 Cancer Antigen only, 128 Reproductive Antigen only.** These are outputs of screening rules, not validated tumor antigens. Every figure is drawn from the [rebuilt catalog](index.md); thresholds are read from the policy constants.

## Catalog evidence, biological sets and host selection

![Curation architecture](summary/source-to-tiers-architecture.png)

Milk-program labels are inferred from expression, not provider-confirmed lactation. Negative signature evidence does not establish universal resting-state absence.

## Host choices and screening outcomes

![Host choices and outcomes](summary/cta-tiers-and-profiles.png)

Reproductive protection is opt-in for every anatomy preset. Not protected, protect, accept loss and absent retain distinct meanings. Unresolved protected-tissue coverage fails the host screen, while biological CA/RA/CTA membership remains unchanged. [Rules and limitations](../../../host-selection.md).

| Host profile | Protected compartments | CTA | CA only | RA only | CA ∪ RA |
|---|---|---|---|---|---|
| all | none | 327 | 70 | 128 | 525 |
| female-patient | none | 327 | 70 | 128 | 525 |
| female-patient (protect ovary) | ovary | 277 | 67 | 113 | 457 |
| female-patient (lactating) | lactating_breast | 320 | 67 | 124 | 511 |
| female-patient (pregnant) | embryo, placenta, lactating_breast | 0 | 0 | 0 | 0 |
| male-patient | none | 327 | 70 | 128 | 525 |
| male-patient (protect testis) | testis | 54 | 29 | 24 | 107 |
| male-patient (protect adult tissues) | testis, epididymis, prostate, seminal vesicle, resting_breast | 40 | 23 | 11 | 74 |

[Profile table CSV](summary/cta-profiles.csv)

## Antigen membership

![Antigen sets](summary/antigen-sets.png)

CTA is the intersection. Reproductive-only proteins have no qualifying cancer evidence; the union name does not certify tumor targeting. Unsafe is the accepted-exclusion background; Unsafe|Other also subtracts unresolved and nonqualifying proteins.

## Rule-derived sets and published nominations

![Set overlap](summary/venn.png)

Of 327 CTA, 252 are published nominations and 75 are not. 1,877 published nominations are in neither set: 1,758 of all 2,233 nominations are Unsafe, 68 Other:Review and 7 Other:Sequence. A published nomination is an annotation; it never overrides an exclusion.

### Curated literature sources by antigen class

![Curated literature sources and dispositions](summary/literature-sources.png)

[All curated literature sources and disposition counts (CSV)](summary/literature-sources.csv)

## Funnel by filtering idea

A protein leaves at the first idea that removes it, so later steps show only what each idea adds. The order is chosen so that every step removes something; the layers overlap heavily, as the next section shows.

![All proteins](summary/funnel-all.png)

| Filtering idea | Removed here | Remaining |
|---|---|---|
| All assessed proteins | 0 | 19,893 |
| Bulk tissue RNA above the tissue ceiling | 17,507 | 2,386 |
| Cell-type cluster RNA (HPA single cell) | 436 | 1,950 |
| Sorted immune-cell RNA | 96 | 1,854 |
| Rare-cell RNA established by donor | 104 | 1,750 |
| Normal-tissue protein (IHC) | 50 | 1,700 |
| Curated or published normal expression | 49 | 1,651 |
| Adult cornea RNA | 5 | 1,646 |
| Coverage, identity, translation or unresolved evidence | 736 | 910 |
| Nested, Unsafe-contained or read-through sequence | 48 | 862 |
| No qualifying expression or bounds fail | 249 | 613 |
| Home tissue no panel resolves | 88 | 525 |

![Published nominations](summary/funnel-published.png)

| Filtering idea | Removed here | Remaining |
|---|---|---|
| Published CTA / oncofetal nominations | 0 | 2,233 |
| Bulk tissue RNA above the tissue ceiling | 1,454 | 779 |
| Cell-type cluster RNA (HPA single cell) | 182 | 597 |
| Sorted immune-cell RNA | 25 | 572 |
| Rare-cell RNA established by donor | 67 | 505 |
| Normal-tissue protein (IHC) | 19 | 486 |
| Curated or published normal expression | 10 | 476 |
| Adult cornea RNA | 1 | 475 |
| Coverage, identity, translation or unresolved evidence | 68 | 407 |
| Nested, Unsafe-contained or read-through sequence | 7 | 400 |
| No qualifying expression or bounds fail | 27 | 373 |
| Home tissue no panel resolves | 17 | 356 |

Steps that remove no published nomination are omitted from the second funnel. [All proteins](summary/funnel-all.csv) · [Published](summary/funnel-published.csv)

## Normal-expression layers

![Unsafe layers](summary/unsafe-layers.png)

| Layer | Unsafe proteins flagged | Only reason |
|---|---|---|
| Normal bulk | 17,507 | 558 |
| HPA tissue clusters | 17,053 | 271 |
| Adult cornea | 12,095 | 5 |
| DICE immune | 11,632 | 45 |
| Monaco immune | 11,264 | 19 |
| HPA normal protein | 10,152 | 47 |
| HPA immune | 7,973 | 14 |
| Established rare normal-cell RNA | 126 | 103 |
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

## Shared enrichment rule

Both CA and RA independently require 10× max(N, D, 1 TPM). Both normal compartments exclude thymus. The cancer panel below shows the tumor side.

![Tumor against normal](summary/tumor-normal.png)

## Rare-cell RNA

![Rare-cell grading](summary/rare-cell.png)

[Donor-level evidence](../../normal-cell-atlas-20261110/rare-cell/index.md)

## Inside the later stages

![Stage detail](summary/stages.png)

- Required normal coverage incomplete: 400
- Identity or translation unresolved: 103
- Normal signal blocks admission (not Unsafe): 233
- Other unresolved evidence or comparison: 0

A protein can use several routes, so route counts exceed set sizes.
