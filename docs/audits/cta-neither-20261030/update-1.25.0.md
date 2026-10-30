# Neither, reproductive IHC and Tumor-restricted candidates

This update restores **51 CTA proteins**, expands Neither to **725**, and preserves every previous CTA and Unsafe assignment. Identical selected reference protein sequences remain collapsed; these are protein counts, not locus counts.

| Category | 1.24.0 | 1.25.0 |
| --- | --- | --- |
| CTA | 522 | 573 |
| Tumor-restricted | 0 | 0 |
| Neither | 20 | 725 |
| Unsafe | 18031 | 18031 |
| Unassigned | 1290 | 534 |

[Full current report and figures](../cta-source-walkthrough-20261028/index.md) · [1.25 logic slides](../cta-source-walkthrough-20261028/cta-selection-1.25.0.pdf) · [Original 20-Neither review](index.md) · [Every category change](changes-1.25.0.csv.gz) · [Transition table](transitions-1.25.0.csv)

## The corrected rules

1. **Normal risk comes first.** Existing Vital/Protected/clinical exclusions establish Unsafe. Missing required normal coverage or unresolved normal evidence remains unassigned. Manageable ceiling excess is a normal hold.
2. **CTA identity requires cancer and Reproductive evidence.** The existing tumor rule stays ≥2 native TPM or HPA pTPM in one sample. Reproductive RNA ≥5 nTPM, supported developmental RNA ≥5 TPM, or reliable positive reproductive IHC can establish association. The existing profile-concentration assessment still applies where RNA is comparable.
3. **Tumor-restricted is a separate candidate subtype.** It needs ≥10 native TPM in two distinct tumor samples, or reviewed reliable tumor protein evidence. It also needs low measured normal/developmental expression, no detected normal protein, resolved attribution, all normal gates passing and a review addressing relevant rare normal cells and physiological states.
4. **Neither is the fallback after a complete normal pass.** Cancer or Reproductive evidence may be high, low or missing. Failure to qualify as a candidate does not establish Unsafe.

All 725 Neither proteins remain below 30 in the measured non-permitted normal TPM/nTPM summaries: the maxima are 16.0 nTPM and 15.3163 TPM, kept as separate assays. Missing optional measurements remain missing.

Neither describes a research screen. It is not evidence that targeting a protein is safe. Existing exact-protein background/k-mer APIs keep their independent policy and are not silently rewritten by this selection update.

## Why the 52 IHC cases were overlooked

The newer selection-policy rewrite kept reproductive IHC as an antibody control but omitted its independent positive-association route. Admission used only reproductive or developmental RNA ≥5. The earlier screen had allowed protein evidence below that RNA floor. The migration tests did not cover that route.

The restored rule requires human HPA v23, a mapped reproductive compartment, Low/Medium/High staining with Approved/Supported/Enhanced reliability, positive evidence for every member locus, and no unresolved antibody concern. It cannot override normal exclusions, cancer admission or a measured concentration failure. A zero reproductive RNA reference leaves RNA concentration unassessed; an IHC category does not supply a fabricated RNA abundance.

**51 of the 52 become CTA. CDY1/B remains Neither because it lacks qualifying cancer expression.** All 522 earlier CTAs, including MAGE-A4, NY-ESO-1, PRAME, CTAG2 and XAGE1A/B, remain CTA. The full HPA observations retain cell type, staining, reliability, mapped anatomy and source version.

[Restored proteins and measurements](restored-CTA-1.25.0.csv) · [Individual positive IHC rows](restored-reproductive-IHC-1.25.0.csv)

| Protein | Reproductive nTPM | OncoRef tumor TPM | HPA tumor pTPM |
| --- | --- | --- | --- |
| ADAM7 | 2.2 | 146.341 | 194.141 |
| CDY2A/B | 2.0 | 8.641 | 0.692 |
| CGB2 | 1.5 | 61.58 | 194.544 |
| CRISP1 | 2.5 | 49.131 | 62.716 |
| CST11 | 0.4 | 25.74 | 48.847 |
| CT45A2/8/9 | 1.9 | 171.309 | 186.48 |
| DAZ2 | 2.9 | 72.118 | 179.182 |
| DAZ3 | 1.4 | 81.418 | 42.448 |
| DCAF8L1 | 0.9 | 85.503 | 35.932 |
| DEFB104A/B | 0.0 | 26.952 | 77.336 |
| DEFB106A/B | 0.0 | 1.51 | 58.895 |
| DEFB107A/B | 0.0 | 57.136 | 12.436 |
| DEFB118 | 0.0 | 22.7 | 32.997 |
| DEFB121 | 0.2 | 4.333 | 14.966 |
| DEFB129 | 0.0 | 30.75 | 15.191 |
| DEFB132 | 0.0 | 73.302 | 69.27 |
| DEFB136 | 0.0 | 33.089 | 20.23 |
| EDDM3A | 1.1 | 50.481 | 60.768 |
| EDDM3B | 0.0 | 177.785 | 37.622 |
| FNDC7 | 1.1 | 93.642 | 24.442 |
| GAGE10 | 1.6 | 200.769 | 82.793 |
| HHIPL2 | 4.7 | 143.777 | 138.104 |
| KRTAP24-1 | 0.0 | 37.99 | 14.782 |
| LALBA | 2.2 | 79700.437 | 94939.203 |
| LCN9 | 0.0 | 85.49 | 110.061 |
| MAGEB18 | 0.7 | 11.186 | 12.675 |
| NANOS2 | 3.5 | 13.42 | 9.804 |
| NLRP14 | 0.8 | 24.433 | 16.486 |
| OTUD6A | 4.6 | 86.213 | 92.992 |
| PATE1 | 0.0 | 163.018 | 146.834 |
| PKDREJ | 3.3 | 34.062 | 14.789 |
| PRDM7 | 1.0 | 183.776 | 25.788 |
| PSG7 | 4.0 | 55.42 | 19.25 |
| RBMXL3 | 2.6 | 12.91 | 17.158 |
| RBMY1B | 0.2 | 12.06 | 14.26 |
| RBMY1D | 1.3 | 4.21 | 1.375 |
| RBMY1E | 2.5 | 54.4 | 10.111 |
| RNASE12 | 1.3 | 13.93 | 27.584 |
| SAGE1 | 4.7 | 355.654 | 46.649 |
| SEMG2 | 0.2 | 3033.778 | 3560.2 |
| SMIM21 | 2.7 | 3.98 | 20.48 |
| SP7 | 0.8 | 590.49 | 170.822 |
| SPAG11A | 0.1 | 472.247 | 484.108 |
| SPAG11B | 0.1 | 957.414 | 557.161 |
| SPINT4 | 0.0 | 39.173 | 55.299 |
| TBC1D3H | 4.6 | 44.63 | 30.126 |
| TBC1D3K | 2.2 | 15.507 | 17.496 |
| TGM4 | 1.0 | 17629.448 | 10160.2 |
| TSPY8 | 1.2 | 1056.516 | 47.02 |
| USP26 | 0.4 | 9.921 | 29.997 |
| WFDC8 | 0.2 | 15.16 | 29.677 |

## Tumor-restricted leads need normal-context evidence

**66 proteins meet the numerical profile, but zero currently have an accepted normal-restriction review.** They remain Neither with a Tumor-restricted lead flag. The normal-low convention is <2 native TPM/nTPM, with a separate nominal CP10K floor; it is not a detection limit or a conversion between assays. Developmental coverage is required for this restrictive subtype. Missing optional sources remain explicitly missing.

Forty leads are olfactory-receptor genes. A general bulk normal panel cannot establish absence from rare sensory cells. [Durante et al., *Single cell analysis of olfactory neurogenesis and differentiation in adult humans* (2020)](https://pmc.ncbi.nlm.nih.gov/articles/PMC7065961/) demonstrates this coverage problem at the receptor-family level; it is not a gene-by-gene verdict on all forty leads. [Human dorsal-root-ganglion single-nucleus profiling](https://pmc.ncbi.nlm.nih.gov/articles/PMC8626086/) supplies relevant normal-cell evidence for SCN10A. [The GDF15–GFRAL receptor study (2017)](https://pubmed.ncbi.nlm.nih.gov/28846097/) reports GFRAL in the human area postrema and nucleus tractus solitarius. These findings motivate coverage review without converting missing bulk measurements into either absence or clinical toxicity.

[All 66 leads and coverage questions](tumor-restricted-leads-1.25.0.csv) · [Their two qualifying tumor observations](tumor-restricted-lead-observations-1.25.0.csv)

The tumor reader uses the existing checksum-pinned source-matrix cache and QC-selected samples. It collapses duplicate source views and TCGA participant/aliquot identities. Both observations must resolve to the protein's loci, and all identical-sequence member loci must be measured before summing within a sample. Counts stop at two and are lower bounds, not prevalence. HPA pTPM is not substituted for native TPM. The reliable-tumor-protein route is implemented, but its review ledger currently has no accepted entries; aggregate HPA cancer-IHC counts alone are insufficient reagent evidence.

USP17L3 and FAM90A9 are now Neither. Their low measured normals no longer become unassigned merely because cancer RNA exceeds the CTA cancer floor. Neither has two qualifying ≥10 TPM observations in this screen.

## Seven old Neither assignments now retain their normal holds

**CRYGN, KRT34, KRT85, KRT71, LEP, LCE2A, KRT25** previously entered a narrow multiple-Manageable-tissue profile despite a normal-screen hold. The fallback now requires every normal gate to pass. They remain unassigned with a normal-evidence hold, rather than gaining an unsupported normal-pass label. Unsafe membership is unchanged.

[Measurements and all holds for these seven proteins](previous-neither-normal-holds-1.25.0.csv)

## API and reproducibility

```python
from oncoref import cta_selection

cta_selection.results(category="CTA")
cta_selection.results(category="Tumor-restricted")  # empty until a review qualifies
cta_selection.candidates()  # union of both candidate subtypes
cta_selection.results(category="Neither")
cta_selection.audit().explain("SAGE1")
```

Policy: `oncoref.cta-selection.2026-10-30.1`. Run `scripts/build_tumor_restricted_support.py` with the pinned matrices cached, then `scripts/regenerate_cta_walkthrough.py --verify-replay` and `scripts/report_cta_125_changes.py`. `--download` on the tumor support builder fills missing cache artifacts through the existing uniform source API. Full and cached decision replay are byte-identical. The [current rules](../../cta-selection.md) document thresholds, exact species scope and clinical-credit limitations.
