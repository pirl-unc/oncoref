# Why CTA counts changed

**The current broad catalog contains 290 CTAs. The 175 and 111 counts are its
optional standard and strict selections. Neither PR #699 nor #700 removes CTAs
from the broad biological catalog. Host choices are another, separate selection.**

Counts are identical selected reference **protein groups**, not gene loci,
transcripts, cDNAs, all isoforms or peptide sequences. Identical protein sequences
do not require identical cDNA sequences because synonymous codons can differ.

## Catalog history

| Snapshot | Broad CTA | Main change |
|---|---:|---|
| 1.36 | 607 | Archived selected-protein snapshot; its later curation replay has 600 |
| 1.37 | 626 | New CA/RA set-based catalog; 74 prior CTAs leave and 93 enter |
| 1.38 | 453 | 153 become RA-only under the tumor 10× rule/TGCT exclusion; 11 become Other:Sequence and 9 Unsafe |
| 1.39 | 437 | 22 lose tumor support after library QC; 5 gain TGCT support and MAGEC2 is restored |
| 1.40 | 246 | Symmetric RA 10×, individual normal samples, stricter review holds and the small-embryo support problem |
| 1.41 | 290 | Small-embryo support fixed: 43 CA-only and 1 Other become CTA |
| PR #699 | 290 | Biological catalog membership unchanged; adds compartment/host views |
| PR #700 | 290 | Biological catalog membership unchanged; separates host decisions from membership |

The 1.39 → 1.40 transition accounts for the largest recent drop, **191**:

| New classification of former CTA | Count | Meaning |
|---|---:|---|
| CA-only | 121 | Cancer support remains; reproductive evidence does not meet the new gate |
| Other:Review | 35 | Unresolved normal-cell evidence/attribution, including CTAG2, CAGE1 and DAZL |
| Other:Tissue | 16 | Normal home-tissue evidence is unresolved by the available panels, including OR and KRTAP families |
| Other:Expression | 12 | Neither positive compartment qualifies after the shared comparison |
| RA-only | 7 | Reproductive support remains; tumor comparison no longer qualifies |

**None of those 191 became Unsafe.** The global Unsafe count stays **18,253**
from 1.38 onward. Changes in reference denominators, admission requirements and
review policy are not equivalent to newly demonstrated toxicity. The table gives
observed transitions, not a causal decomposition of overlapping gate changes.
The subsequent embryo fix restores 44 CTAs, including cleavage-stage candidates.

After that correction, the settled **437 → 290** difference is 147 proteins:
78 CA-only, 35 Other:Review, 15 Other:Tissue, 12 Other:Expression and 7 RA-only.
All 290 current CTAs were already among the 437; this net comparison has no new entrants.

## Optional evidence filters

The following attrition is sequential, starting from the same 290 current CTAs.
An earlier failing gate receives the count when a protein fails several gates.

| Mode | Tighter normal-RNA limit | Curated normal home | Publication count | Named MAGE exclusion | CTA retained |
|---|---:|---:|---:|---:|---:|
| all | 0 | 0 | 0 | 0 | 290 |
| loose | 0 | 0 | 53 | 0 | 237 |
| standard | 16 | 2 | 97 | 0 | 175 |
| strict | 23 | 4 | 141 | 11 | 111 |

Standard requires at least two counted papers; strict requires three. These are
publication-support choices, not independent tumor-patient counts or additional
measurements of normal-tissue toxicity. Strict also excludes MAGE-prefixed labels
except MAGEA4. The 111 strict candidates are all CTAs because strict requests the
intersection; the full CA/RA union still contains 509 proteins.

## Host selection is a third axis

PR #699's host presets implicitly allow gonadal exposure. PR #700 instead protects
present gonads unless their loss is explicitly accepted. It also independently
checks native cell-RNA measurements, unresolved protein findings and coverage.
These choices can yield fewer than 200 **host-eligible** CTAs without removing any
of the 290 biological CTAs. A whole-testis or whole-ovary restriction cannot be
interpreted as a selective fertility-only intervention from bulk RNA alone.

| Host context in corrected PR #700 | all mode | standard mode | strict mode |
|---|---:|---:|---:|
| No host screen | 290 | 175 | 111 |
| Female-anatomy preset, ovary protected | 239 | 144 | 90 |
| Male-anatomy preset, testis protected | 49 | 16 | 7 |

Accepting testis loss in the same male-anatomy preset raises its all-mode result
from **49 to 251**. That is a deliberate change in the permitted exposure, not new
cancer evidence or demonstrated tolerability. The corresponding female-anatomy
comparison is **239 to 272** when ovary loss is accepted. Presets retain explicit
anatomy and current-state assumptions; actual organs must be verified separately.

Report each result with its catalog version, set, evidence mode and host context.
For example, “290 broad CTAs; 175 standard CTAs; host context not applied” is
unambiguous. “175 CTAs” by itself conflates catalog and selection.

## Reproduction and review

[Catalog counts](catalog-counts.csv), [exact category transitions](transitions.json)
and [filter attrition, including independent failures](filter-attrition.json),
with [every protein entering or leaving CTA](cta-transitions.csv),
were computed from committed packaged tables. The CSV records the Git commit
for every snapshot. Tables were joined on stable protein entity IDs; the same
19,893 identities occur in every snapshot. For 1.36 the source is
`oncoref/data/cta-selection-human.csv.gz`; subsequent snapshots use
`oncoref/data/target-selection.csv.gz`. Filter attrition evaluates
`target_filters.apply(FILTER_MODES[mode], current_cta_rows)` in declared layer order.

Review of [#699](https://github.com/pirl-unc/oncoref/pull/699) identified missing-as-zero
host coverage, the protected-expression cap mismatch, implicit gonadal loss and
inconsistent sequence/background views (issues #693 through #698).
Those corrections are in [#700](https://github.com/pirl-unc/oncoref/pull/700).
Further review reproduced and corrected a convenience-API host override
([#701](https://github.com/pirl-unc/oncoref/issues/701)), false host exclusions after
tightening tumor evidence ([#702](https://github.com/pirl-unc/oncoref/issues/702)),
omitted HPA tissue aliases ([#703](https://github.com/pirl-unc/oncoref/issues/703)),
and two CI inventory/fixture failures ([#704](https://github.com/pirl-unc/oncoref/issues/704)).
The host summary's unescaped union-column separator was also corrected
([#705](https://github.com/pirl-unc/oncoref/issues/705)).
