# OncoRef 1.34.0: extracellular sizes, cancer frequency and detailed concerns

The [37-protein membrane table](audits/cta-literature-20261107/report.html#surface)
now reports extracellular amino-acid intervals, total external residues and the
longest continuous external segment. It links every protein to UniProt, Ensembl
and HPA and retains each feature’s evidence code.

UniProt 2026_04 supplies exact external boundaries for 34 of the 37 reference
sequences. TEX38, XKR3 and XKR5 have transmembrane annotations without
outside-facing boundaries. Their sizes remain unknown.

Thirty-five reference sequences match the selected OncoRef sequence exactly.
SPAM1’s selected sequence has the same 455-aa mature chain but a different anchor
tail. SLCO1B3–SLCO1B7’s selected 748-aa sequence differs from the 687-aa UniProt
reference; its coordinates are not transferred. ZP4’s 487-aa external precursor
region is distinguished from its 444-aa annotated mature matrix chain.

Cancer frequency uses HPA 25.1 native pTPM. All 37 proteins have RNA measurements
across 21 TCGA cancer types; validation cohorts remain separate. The table shows
positive/measured sample counts at ≥1 and ≥5 pTPM, with filters for protein,
cohort and cancer. HPA IHC supplies separate staining counts for 14 proteins.
Missing measurements remain distinct from measured zero. This panel does not
cover every cancer, including NUTM, or measure selected proteoform/surface abundance.

The [29 concern explanations](audits/cta-literature-20261107/report.html#concerns)
now describe the assays, normal compartments, identity limits and next decisions.
FRG2/FRG2B experimental counterevidence remains separate from the five unnamed
coding loci. RLN1 nondetection retains the tested sample scope; the reported
REL2 serum detection limit is not assigned to REL1. An earlier SRY description is
corrected: SRY is the founding SOX-family member, but other SOX-gene studies do
not establish SRY expression.

New APIs: `cta_literature.surface_domains()`, `surface_cancer_rna()` and
`surface_cancer_ihc()`. Three CSV downloads, native feature evidence, source
receipts and a size/prevalence figure accompany the report. The builder runs
offline from pinned inputs.

The classification snapshot remains **1.31.0: 602 CTA, 18,255 Unsafe and 1,006
Other**. This release adds evidence without changing admission decisions.

Validation checks inclusive coordinates, GPI maturation, sequence mismatches,
missing versus zero expression, cohort separation, report filters, offline
reproduction and the exported figure.
