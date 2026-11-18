# Catalog versions and optional selection losses

Counts are identical selected-reference protein groups. CA-only and RA-only exclude the CTA intersection. The either-set count includes CTA once.

| Version | CTA | CA-only | RA-only | Either CA or RA | Change in CTA | Main change |
|---|---:|---:|---:|---:|---:|---|
| 1.36 | 607 | — | — | — | — | Earlier catalog definition |
| 1.37 | 626 | 245 | 31 | 902 | +19 | CA/RA set algebra: 74 leave, 93 enter |
| 1.38 | 453 | 37 | 180 | 670 | -173 | 153 to RA-only; 11 to Other:Sequence; 9 to Unsafe |
| 1.39 | 437 | 34 | 197 | 668 | -16 | 22 to RA-only after tumor QC; 6 CTAs restored |
| 1.40 | 246 | 149 | 102 | 497 | -191 | RA 10×, individual normal samples, review holds and restrictive embryo gate |
| 1.41 | 290 | 106 | 113 | 509 | +44 | Embryo gate repaired; 44 CTAs restored |
| PR #699 | 290 | 106 | 113 | 509 | +0 | Compartment and host views; membership unchanged |
| PR #700 | 290 | 106 | 113 | 509 | +0 | Explicit host choices and gene-specific paper uses; membership unchanged |
| Current (oncoref.target-selection.2026-11-18.3) | 308 | 68 | 131 | 507 | +18 | GTEx individual reproductive witnesses and within-site normal donor P99; reproductive-program flags are advisory |

The frozen 1.39 → 1.41 difference is 147: 78 CA-only, 35 Other:Review, 15 Other:Tissue, 12 Other:Expression and 7 RA-only. None became Unsafe. These are observed destinations, not independent causal attributions for overlapping gates.

Each optional mode below starts from the same 308 CTAs; losses are sequential within that mode. Default anatomy protects no reproductive tissue and removes zero.

| Mode | Normal-RNA loss | Curated-home loss | Paper loss | MAGE-name loss | CTAs kept |
|---|---:|---:|---:|---:|---:|
| all | 0 | 0 | 0 | 0 | 308 |
| loose | 0 | 0 | 51 | 0 | 257 |
| standard | 31 | 2 | 101 | 0 | 174 |
| strict | 38 | 5 | 152 | 10 | 103 |

Standard requires two counted tumor-support papers; reproductive-only papers remain separate support. The [IHC counterfactual](index.md#how-much-would-reproductive-high-ihc-add) is not enabled.

Catalog `oncoref.target-selection.2026-11-18.3`; filter `oncoref.target-filters.2026-11-17.5`. Historical snapshots and exact transitions: [count history](../cta-count-history-20261117/index.md).
