# CTA tissue-tier update — 20 October 2026

Policy `oncoref.cta-ihc-adaptive-screen.2026-10-20.2`. All 2,537 published candidates are reevaluated from the same pinned human HPA v23 observations.

## Tissue changes

- Tonsil and spleen become **acceptable (conditional)**, with explicit review flags.
- Lymph nodes stay **protected**; blood and marrow stay **vital**.
- Retina, eye and inner ear move to **vital / critical function**.
- Reproductive unions and other organ assignments are preserved.

The adaptive vital/protected ceilings remain **1/5, 2/5, 3/7.5 and 5/10 nTPM** for missing/uncertain, Approved, Supported and Enhanced negative IHC. The proposed 1/2/3/6 and 2/4/6/12 ladder is not adopted in this comparison.

## Results

| scenario | candidates | evidence_pass | retained | without_ihc | partial_ihc | complete_routine_ihc |
| --- | --- | --- | --- | --- | --- | --- |
| previous_adaptive | 2537 | 771 | 673 | 368 | 298 | 7 |
| current_fixed | 2537 | 711 | 613 | 359 | 248 | 6 |
| current_adaptive | 2537 | 752 | 654 | 359 | 288 | 7 |

![Retained candidates and IHC coverage](policy-comparison.png)

**Gained:** ADAMTS6, FAM9C, IL17F, LCE3C, NOX5, OTOA, TGM7.

**Lost:** C11orf42, C1orf146, C20orf141, CCDC172, CCDC183, CCDC62, CCDC87, DCST1, DNAH17, GAGE12J, GAGE13, GDPD4, LDHAL6A, LIN28B, NOXRED1, PIWIL1, SPATA31H1, SYCE2, TAF1L, TET1, TEX35, TKTL1, TPRX1, WDR88, ZNF560, ZPBP.

**228 retained candidates have positive observations requiring tissue review.** This is a research-screen result, not completed therapeutic acceptance. Review flags record RNA >0 nTPM or positive IHC in tonsil, spleen, thymus or lymph-node contexts; they do not independently change expression-pass or final-retention decisions.

## Evidence and interpretation

The full [tissue definitions and citations](../../normal-tissue-scopes.md) explain conditional organ allowances and their limits. Native tissue/cell labels and citations remain in [tissue-contexts.csv](tissue-contexts.csv). An organ allowance does not allow expression in the same immune-cell population elsewhere. Unsited immune cells remain unresolved. Mouse has no mapped tonsil; canine and mouse anatomical assignments do not transfer human thresholds or clinical tolerability.

Any positive IHC outside acceptable tissues, including Low, still excludes. RNA relaxation requires valid negative staining across every routine cell context of the same anatomical tissue. Missing/partial IHC is explicit. Complete routine coverage means 50 RNA contexts and 80 IHC tissue/cell pairs, not exhaustive anatomy. Whole eye and inner ear have no routine RNA measurements in this input; assigning their tier does not establish absence of expression in those organs. Existing family and manual curation gates remain separate. CTAG2's historical manual holdback is retained. Negative observations do not prove antigen absence.

## Downloads and reproduction

- [All candidate comparisons and reasons](candidate-comparison.csv)
- [Gained/lost candidates](membership-changes.csv)
- [Named targets, including the preceding sensitivity-analysis genes](notable-targets.csv)
- [Current fixed and adaptive decisions](candidate-decisions.csv)
- [Every RNA value, ceiling and supporting IHC](rna-threshold-audit.csv.gz)
- [Taxonomy, review conditions and bibliography](tissue-taxonomy.json)
- [Run manifest and checksums](run-manifest.json)

Run `python scripts/cta_tissue_tier_report.py --out <new-empty-directory>` at this report's implementation revision. The preceding [adaptive report](../cta-ihc-adaptive-20261020/index.md) remains an immutable baseline.
