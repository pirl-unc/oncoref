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

## What changed in 1.40

**Symmetric reproductive 10×.** CA and RA independently use the same native-TPM
reference, `H = max(non-reproductive normal N, Disposable D, 1 TPM)`. CA needs a
QC-passing tumor observation at least `10 × H`; RA needs a reproductive tissue
sample or supported embryo-lineage mean at least `10 × H`. Before 1.40 RA was
an association rule: detection at 1 TPM or 1 nTPM with the source's support checks,
or reliable reproductive protein staining, could qualify without a 10× fold.
A protein can
therefore remain CA but lose CTA membership if reproductive evidence is weaker.
For example, with `H = 2 TPM`, tumor 50 TPM qualifies CA, while reproductive
15 TPM does not qualify RA. The tumor need not exceed reproductive tissue by 10×.
Individual HPA normal biopsies were also added to N; that can raise both branches'
thresholds. nTPM and UMI counts cannot substitute for native TPM in this fold.

**Stronger review holds.** Consequential normal-cell signals whose attribution or
replication is unresolved now enter Other:Review instead of being dismissed just
because the locus has mapping risk or the signal comes from one study. The screen
uses donor-resolved multi-molecule support, cell frequency and study replication,
plus bulk RNA or gene-specific literature confirmation. A mapping flag says reads
could be misassigned; it does not prove they were. CAGE1, CTAG2 and DAZL illustrate
these holds. Resolution requires gene-specific evidence: unique exonic/junction
support, independent population replication, orthogonal localization, or a supported
artifact adjudication. Accepted normal expression can make a protein Unsafe;
demonstrated artifact can dismiss the signal and restore eligibility if other gates
pass. Negative bulk RNA alone does not settle a sparse-cell finding. These holds
are unrelated to the optional minimum-publication count.

**Overly restrictive embryo gate.** The 1.40 implementation moved from stage-pooled
cells to separate embryo × lineage means but kept a 20-measured-cell requirement.
An eight-cell embryo cannot supply 20 cells. The corrected 1.41 rule retains the
large-population route (at least 20 measured cells, at least 5 positive, at least
10% positive) and adds a small-population route (at least 5 measured cells, at
least 50% positive, with two independent embryos at the same stage/lineage reaching
10×). A positive cell has at least 1 TPM. Multiple cells from one embryo do not
count as independent embryos. This correction restored 44 CTAs; the smaller
population route still has an explicit replication requirement.

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

[All 97 proteins excluded at the standard paper step](paper-filter.md) are listed
with their source citations and tumor-patient counts. Of these, 49 have zero and
48 have one counted publication; 59 already have at least two known tumor patients
at 10×. Zero denotes the curated registry, not an exhaustive absence of literature.

## Host selection is a third axis

The default anatomy choice protects **no reproductive tissue**, for either men
or women. The same 290 broad CTAs remain available in both anatomy presets.
Protection is an explicit opt-in through `HostContext(protect=(...))` or
`protect_tissues=(...)`; it independently checks native cell-RNA measurements,
unresolved protein findings and coverage in the requested compartments.
Not protected, absent and accepted loss have distinct meanings. A whole-testis
or whole-ovary restriction cannot be interpreted as a selective fertility-only
intervention from bulk RNA alone.

| Host context in corrected PR #700 | all mode | standard mode | strict mode |
|---|---:|---:|---:|
| No host screen | 290 | 175 | 111 |
| Female-anatomy preset, default protection | 290 | 175 | 111 |
| Male-anatomy preset, default protection | 290 | 175 | 111 |

In all mode, protecting **only ovary retains 252 CTAs**; protecting **only testis
retains 62**. Both are explicit subsets of the 290, exposed through the same API.

Earlier review counts of 49 for the male preset and 239 for the female preset
included implicit protection of gonads and adult accessory tissues. They are no
longer defaults. Those smaller sets remain available by explicitly requesting
the corresponding tissues. Presets retain explicit anatomy and current-state
assumptions; actual organs must be verified separately.

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
