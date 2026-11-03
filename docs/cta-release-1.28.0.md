# OncoRef 1.28.0: every protein resolves, with a coding universe, cell-type tiers and sequence rules

**637 CTA / 18,278 Unsafe / 689 Other** representative proteins across
**19,604** proteins. No human protein is unassessed. Other includes
**601 Nonreproductive, 34 Held, 28 Nested, 20 Contained and 6 Diffuse**.

[Rules and API](cta-selection.md) ·
[Full source report](audits/cta-source-walkthrough-20261028/index.md) ·
[Category funnel](audits/cta-source-walkthrough-20261028/figures/category-funnel.png)

Downstream semantics are unchanged: CTA is targeted, Unsafe blacklists any matching
subsequence from CTAs, and everything else is ignored for targeting.

## What changed

| Change | Effect |
|---|---|
| **Universe** | Only proteins translated from a coding transcript or an established ncORF are classified (19,604). 259 groups are out of universe: 254 NMD-only or CDS-not-defined, 5 single-study Ribo-seq ORFs. They are listed with reasons in `out-of-universe.csv.gz`. |
| **No Unassessed class** | The 647 Unassessed proteins resolve to Unsafe or Other. `assessment_status` records whether coverage was incomplete or a hold applied. |
| **Coverage** | Only combined bulk (HPA consensus, HPA brain, GTEx v11) is required. The GTEx 52-context core and cell-cluster, immune, DICE and cornea panels never block a label. |
| **Unsafe, higher bar** | Manageable-tier RNA above 128 is Unsafe evidence (12 proteins, mostly skin-barrier genes). A reviewed unresolved Protected-tissue finding is an essential unknown and is Unsafe (113 proteins). |
| **Cell-type tiers** | One table (`oncoref.cta_cell_tiers`) feeds the RNA cascade and the rare-cell screen. Vital: neurons, glia, cardiomyocytes, alveolar types 1 and 2, endothelium, airway basal cells. Manageable: differentiated airway epithelium (ciliated, club, goblet, secretory, ionocytes), fibroblasts, bowel epithelium. |
| **Rare vital-tissue screen** | A protein is Unsafe when at least 1% of a Vital cell type (brain, spinal cord, heart, lung) expresses it in CELLxGENE normal summaries. Evaluated for every protein; 117 are Unsafe only because of it. No CTA is flagged. |
| **Contained** | A would-be CTA with at least 95% of its 9-mers inside Unsafe proteins is `Other:Contained`: its peptides are blacklisted. 20 proteins. |
| **Nested** | A reference-CDS microprotein (under 50 residues) nested in a same-strand protein-coding gene, without ncORF evidence, is `Other:Nested`. 28 proteins. |
| **Held** | Target gates pass, but a marginal or unresolved normal signal blocks CTA without establishing Unsafe: 34 proteins, including MAGEA3. |

Unsafe reconciles exactly: 18,062 previous, minus 14 that left the universe, minus 9
released because differentiated airway cells moved to Manageable, plus 239 new
(110 reviewed holds, 114 rare vital-tissue, 3 both, 12 Manageable above 128), gives
18,278.

## CTA membership: 701 to 637

**80 previous CTAs left** and 16 entered.
- **67 became Unsafe** through the rare vital-tissue screen alone, for example CSAG1
  (midbrain dopaminergic neurons 1.9%), PLAC1 (cerebellar granule cells 1.9%), DPPA2
  (airway basal cells 1.7%), SOX30, SMC1B, TDRD1, TPTE, GDF9, MORC1 and HORMAD2.
- **13 became Other:Contained**: CGB2, CSH2, CT45A3, DAZ2, DAZ3, GAGE2E, PRAMEF27,
  RIMBP3C, TAF11L9, TBC1D3G/H/K and ZC3H11C. Most multi-copy families (CT45, GAGE,
  PRAMEF, POTE) stay CTA at the 95% floor; their usable subsequences are left for review.
- PRAME, MAGEA4, MAGEA1, MAGEC2, NY-ESO-1, CTAG2, SSX2, CT83, XAGE1A/B, TBXT and
  GAGE10 remain CTA.

The 16 new CTAs include readthrough chimeras (ATF7-NPFF, EPPIN-WFDC6, SLCO1B3-SLCO1B7,
TGIF2-RAB5IF) and unnamed loci, plus MESP2, MSGN1 and TSSK2.

## Known limitations

- Several new CTAs rest on weak or embryonic-only reproductive RNA (for example
  ENSG00000258417, 0.0 nTPM adult reproductive) and readthrough products share RNA with
  their parent genes. They pass every gate and are reported as such; review is advised.
- The tumor gate remains one measured sample above 0, and the reproductive route
  accepts embryonic RNA of 5 TPM or more, so most proteins pass both.
- The rare-cell screen reports detection frequency among cells, not abundance, and can
  reflect ambient RNA, doublets or paralog multimapping; samples cannot be restricted to
  adults. It flags 15,888 proteins, most already Unsafe. Whether a cell type counts as
  Vital or Manageable is a policy judgement. `cta_rare_cell.UNSAFE = False` makes it a hold.
- MAGEA12 is Unsafe through its reviewed clinical exclusion and also now through the
  screen (midbrain dopaminergic neurons 1.5%). Research prioritization, not clinical safety.

## Downstream adoption

```python
from oncoref import cta_selection

results = cta_selection.results()        # category is never null for human
profiles = cta_selection.profiles()      # adds rare_cell_* and unsafe_sequence_overlap_pct
screen = cta_selection.rare_cell_screen()  # includes tier
roles = cta_selection.sequence_dispositions()
```

Update code that handles null human categories or the `Unassessed` label. New subtype
values are `Held`, `Contained` and `Nested`. Proteins outside the coding or
established-ncORF universe no longer appear in `results()`.
