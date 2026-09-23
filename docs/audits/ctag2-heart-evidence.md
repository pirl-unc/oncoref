# CTAG2 / LAGE-1 heart evidence review

Reviewed 2026-09-22 for [oncoref #544](https://github.com/pirl-unc/oncoref/issues/544).
CTAG2 remains outside the strict default and discoverable through the warning
and clinical tiers. The heart RNA signal is not corroborated by cardiomyocyte
IHC, but neither negative IHC nor the clinical record establishes cardiac safety.
The review answers the literature question; the biological mechanism remains open.

## Normal-tissue evidence

The pinned HPA v23 consensus estimates heart muscle at **5.1 nTPM**. Its matched
normal-tissue IHC table annotates CTAG2 in cardiomyocytes as **Not detected**,
with **Enhanced** reliability. The current HPA primary-data page also reports
non-detection with antibody HPA071467. These are assay-specific observations,
not a measurement of peptide-HLA display.
Sources: [v23 RNA](https://v23.proteinatlas.org/download/rna_tissue_consensus.tsv.zip),
[v23 IHC](https://v23.proteinatlas.org/download/normal_tissue.tsv.zip),
[HPA primary data](https://www.proteinatlas.org/ENSG00000126890-CTAG2/tissue/primary%2Bdata).

The existing `cta-reviewed-evidence` donor, single-cell, and isoform reviews
retain the other relevant distinctions: the bulk signal occurs in a minority
of donors; single-cell evidence does not reliably locate its cell of origin;
and transcript estimates include the epitope-containing LAGE-1A isoform in the
RNA-positive heart donors. RNA and IHC donors differ. Thus neither stromal
contamination nor expression solely of an epitope-lacking isoform is established
as an explanation. No reviewed evidence establishes cardiac presentation of
SLLMWITQC. The RNA/IHC disagreement should remain visible.

## Clinical evidence and target identity

**Afami-cel / TECELRA and SPEARHEAD-1 concern MAGE-A4**, not NY-ESO-1/CTAG2.
They cannot serve as clinical validation for this target. **Lete-cel** recognizes
the HLA-A*02-presented NY-ESO-1/LAGE-1A shared epitope SLLMWITQC.
Sources: [TECELRA label](https://dailymed.nlm.nih.gov/dailymed/drugInfo.cfm?setid=2ddd66e1-8036-4a4e-babe-4d673e660bf5),
[lete-cel primary report](https://aacrjournals.org/clincancerres/article-split/31/3/529/751208/Safety-and-Tolerability-of-Letetresgene-Autoleucel).

| Source reviewed | Finding | Interpretation |
| --- | --- | --- |
| [Robbins 2011](https://pmc.ncbi.nlm.nih.gov/articles/PMC3068063/) and [2015 follow-up](https://pmc.ncbi.nlm.nih.gov/articles/PMC4361810/) | No toxicity attributed to transferred NY-ESO-1-reactive cells. The follow-up reports a treatment-related death from E. coli septic shock during neutropenia. | Reassuring for those constructs and cohorts, but not proof of no rare cardiac events. The reports overlap and must not be counted as independent cohorts. |
| [Lete-cel MRCLS pilot, 2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC12084024/) | One potentially treatment-related fatal cardiac arrest on day 161 among 20 infused patients, after progression, with hypotension, renal insufficiency, and intervening radiation/pazopanib. | A cardiac event is documented; the confounded case does not establish CTAG2-mediated heart injury. |
| [NCT03709706 protocol, amendment 7, 2021-11-04](https://cdn.clinicaltrials.gov/large-docs/06/NCT03709706/Prot_000.pdf), Table 9, printed p. 42; program safety summary pp. 46–47 | Two unexpected fatal cardiac arrests, with late hypotension/renal confounding in one and early infection/pulmonary confounding in the other. The day-161 case also had early CRS with supraventricular tachycardia; investigator and sponsor attribution differed. | Program summaries overlap study reports. Do not add these as independent cases or describe the program as having no cardiac signal. |
| [Linette 2013](https://pmc.ncbi.nlm.nih.gov/articles/PMC3743463/) | Fatal cardiac toxicity with an affinity-enhanced MAGE-A3 receptor implicated titin-peptide cross-reactivity. | A receptor-specific off-target precedent, not evidence of native CTAG2 heart expression. |

This is a targeted review of the named reports and a relevant program safety
protocol, not a systematic estimate of event incidence across every NY-ESO-1
trial. Adverse-event attribution, native antigen expression, and receptor
cross-reactivity are separate questions. No reviewed report demonstrates
CTAG2-mediated cardiomyocyte injury, and the documented cardiac events preclude
using a supposed absence of such events to dismiss the RNA warning.

The clinical findings, source versions, and limitations are also available from
`cta_review.cta_reviewed_evidence()` and the `clinical_review` fields of
`cta_review.cta_evidence_summary()`. This curation does not promote CTAG2 or
establish safety for an individual receptor or patient.
