# Fetal and embryonic targets: source requirements

Assessment date: 26 October 2026. This is a source and evidence-gap assessment,
not an ingestion-complete developmental target list.

**Brachyury is TBXT (historical symbol T), ENSG00000164458.** The current
reference-bound expression screen places it in Neither / Antibody. It should
be evaluated as a developmental/lineage target with explicit adult-expression
review, rather than assuming that an embryonic label proves adult absence.

## Why brachyury needs additional adult coverage

A primary human study assessed gene expression in 60 and IHC in 56 nucleus
pulposus specimens and reported brachyury in adult cells across age and
degeneration groups. This is a concrete adult compartment often absent from
routine bulk-normal panels. Other studies report age-dependent findings, so
specimen state and antibody specificity must be reconciled rather than reducing
the evidence to a single presence/absence flag.
[Richardson et al., 2017](https://pubmed.ncbi.nlm.nih.gov/28473691/),
[Tang et al., 2016](https://pmc.ncbi.nlm.nih.gov/articles/PMC5321132/).

Chordoma clinical experience also needs careful interpretation. The randomized
GI-6301 vaccine/radiotherapy trial had one partial response in each arm and did
not establish the proposed benefit. A later single-arm BN-Brachyury/radiotherapy
study reported two partial responses among 26 evaluable patients and missed
its predefined response goal. These observations establish clinical experience,
not general efficacy or tolerability of all TBXT-targeting modalities.
[Randomized phase II study](https://pmc.ncbi.nlm.nih.gov/articles/PMC8100546/),
[Cote et al., 2024](https://acsjournals.onlinelibrary.wiley.com/doi/full/10.1002/cncr.35477).

## Data to add

| Evidence layer | Primary sources | What to retain and why |
|---|---|---|
| Early embryonic lineage | [Tyser human gastrulation atlas](https://www.nature.com/articles/s41586-021-04158-y) | Cell counts and annotations, exact embryo/stage and anatomical lineage. Early axial mesoderm/notochord programs can be missed by later fetal organ sampling. A single-embryo dataset cannot establish population prevalence. |
| Early organogenesis | [Human embryo atlas, 2023](https://www.nature.com/articles/s41556-023-01108-w); [post-gastrulation spatial atlas, 2026](https://doi.org/10.1038/s41586-026-10545-0) | Embryo, stage, dissected region, spatial coordinates and cell-state labels. Spatial context distinguishes lineage-localized expression from a broad whole-embryo signal. Confirm access and licensing for each actual matrix before registration. |
| Multi-organ fetal expression | [Cao et al. fetal atlas](https://pubmed.ncbi.nlm.nih.gov/33184181/) and [DESCARTES data portal](https://descartes.brotmanbaty.org/bbi/human-gene-expression-during-development) | Raw/processed expression and donor-level metadata across 15 organs. Useful for a general developmental nomination source; it does not cover every early transient lineage or every developmental stage. |
| Musculoskeletal development | [Human embryonic skeletal atlas](https://www.nature.com/articles/s41586-024-08189-z) | RNA and spatial evidence for joint/cartilage and related developmental states. Include adult disc/nucleus-pulposus studies separately; fetal skeletal coverage cannot stand in for adult absence. |
| Adult normal and rare cells | Existing pinned HPA/GTEx/cell evidence, plus the adult nucleus-pulposus primary studies above | Donor, age, degeneration/inflammation, precise disc compartment, assay identity and protein controls. Broaden other candidates to their relevant adult residual lineages and regenerative states. |
| Tumor recurrence and lineage | Patient chordoma bulk/single-cell cohorts, including [GSE239531](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE239531), linked to its [primary study](https://pmc.ncbi.nlm.nih.gov/articles/PMC10824073/) | Malignant-cell assignment, patient denominator, site, treatment status, tumor subtype and assay. Cell-line dependency data are useful but cannot substitute for patient occurrence or normal-cell coverage. |
| Protein and presented antigen | Gene-specific IHC/proteomics and [HLA Ligand Atlas](https://hla-ligand-atlas.org/welcome), with its [primary methods paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC8054196/) | Keep gene, isoform/protein, peptide and HLA identities separate. A missing mass-spectrometry observation is detection-limited; normal RNA and protein do not directly quantify peptide presentation. |
| Clinical evidence | Primary trial reports and trial registries | Modality, product, antigen identity, HLA where applicable, tumor response denominators, adverse events and attribution. Separate registration, immune responses, objective tumor response and toxicity findings. |

The existing normal-cell source registry already contains staged **fetal corneal
data from GSE155683**, alongside its adult data. That is a useful stage-aware
ingestion precedent, not general embryonic or fetal coverage. Its fetal rows
must not enter adult-normal exclusions by accident.

## Integration contract

Each downloaded artifact needs an accession/version, public source URL, license,
content hash and immutable source receipt. Keep native counts, normalized values
and their actual units; the current human HPA nTPM thresholds do not transfer
directly to sparse single-cell counts, CP10K, TPM or other species.

Every observation needs an exact gene ID and annotation release, transcript/
protein identity where measured, species, donor/embryo, independent specimen,
anatomy, cell type, assay, developmental stage and physiological state. Record
whether age means post-conception weeks, gestational weeks or another convention;
retain the original label and uncertainty. In-vitro embryo models, organoids,
adult stem cells and mouse atlases must be separate source classes, not human
in-vivo observations.

Summaries must expose the independent donor denominator, per-donor distribution,
detected-cell fractions and technical coverage. Do not treat thousands of cells
from one embryo as thousands of independent donors. Missing cell types, stages,
normal organs or non-detection remain explicit gaps. Maternal, placental and
fetal compartments require distinct labels.

Developmental nomination should be a separate evidence field beside adult
normal-expression results and tumor recurrence. It can nominate TBXT, AFP,
CLDN6 and other independently supported targets for review without silently
overriding their adult findings. In the current policy AFP has a normal-RNA
exclusion and CLDN6 a high-normal-RNA exclusion; a developmental label alone
does not revise those measurements.

## Recommended first implementation

Register and checksum the gastrulation, organogenesis and DESCARTES inputs;
retain their original donor/stage annotations. Build a reference-bound TBXT
case report joining embryonic lineage, adult disc, patient chordoma and clinical
records. Audit gene-specific reagent/isoform identity and measured adult normal
contexts before proposing any membership change. Then apply the same source
contract to other developmental nominees. No new developmental classification
cutoff or automatic fetal-expression admission is established by this assessment.
