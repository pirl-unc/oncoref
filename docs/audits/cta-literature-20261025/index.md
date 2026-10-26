# Initial CTA literature audit — 25 October 2026

**Historical review snapshot:** the proposals below preceded the completed
[adjudications](../cta-adjudication-20261025/index.md) and
[current 213-gene policy](../cta-policy-20261025/index.md). In particular, the
PLCZ1 pancreatic protein interpretation was corrected to qualitative RNA;
this frozen search table is not the current classification API.

**I would not endorse all 399 retained genes as established safe targets.** The
current categories reproduce the approved screen, but newly examined literature
supports five proposed exclusions and three review holds. These are proposed
research classifications; the released sets have not been changed.

[Search every gene](explore.html) · [All-gene table](all-genes.csv) ·
[Proposed changes](proposed-changes.csv) · [Paper leads](papers.csv) ·
[Exact searches](searches.csv) · [Reproduction and limitations](reproduction.md) ·
[General category proposal](category-proposal.md) · [Filed issues](issues.md)

## What the current decisions mean

| Current result | Genes | Interpretation |
|---|---:|---|
| CTA | 399 | Candidate passes the recorded screen; not established safe for every peptide, receptor, patient or tissue state |
| Neither | 561 | Neither a target nor default forbidden background; evidence is low, missing, variable or unresolved |
| Unsafe | 1,580 | Excluded under a normal-expression or safety rule; many have no demonstrated clinical toxicity |
| Total | 2,540 | Mapped source nominations, not all human genes or all historical CTA labels |

All 2,540 frozen category/detail assignments matched the 24 October API. No current
CTA or Neither row has an audited non-permitted normal maximum above 30 TPM or
nTPM. These are consistency checks, not independent validation of the biological
cutoffs. The two units are checked separately; corneal CP10K is not converted to
TPM.

The per-gene table explicitly distinguishes a policy-consistent decision from
completed primary-paper adjudication. **All 2,540 genes were searched in 6,013
queries; 32 genes received new primary-paper assessments and 46 have specific
assessment notes including carried-forward reviews.** The remaining 2,508 genes
have not received a new individual primary-paper adjudication in this pass.
There are 36,964 unique retrieved paper leads, 1,688 truncated searches and 160
genes with no hits in either topic query. No-hit searches do not establish safety.
Searching every gene does not mean every paper has been read or every result is
specific to that gene. `summary.json` records these coverage limits; unresolved
rows remain visibly provisional.

## Changes I would propose

| Gene | Current → proposed | Evidence and interpretation |
|---|---|---|
| **FTMT** | CTA → Unsafe | Human substantia nigra protein in neurologically unaffected postmortem donors, with a monoclonal antibody, absorption controls and an H-ferritin specificity check. This was absent from the previous 37-gene conflict list. [Primary study](https://pubmed.ncbi.nlm.nih.gov/28386150/) |
| **ABCC11** | CTA → Unsafe | Human brain immunoblotting and axonal immunofluorescence conflict with reproductive restriction. Abstract reviewed; full reagent assessment remains useful. [Primary study](https://pubmed.ncbi.nlm.nih.gov/16359813/) |
| **YBX2** | CTA → Unsafe | The original dbpC/contrin study reports vascular smooth-muscle and adrenal-cortical protein. Specimen context and polyclonal-reagent limitations are recorded; vascular staining is not cardiomyocyte staining. [Primary study](https://pmc.ncbi.nlm.nih.gov/articles/PMC2361212/) |
| **CCDC158** | CTA → Unsafe | Published renal tubular protein and endocytosis biology contradict unrestricted gene-level inclusion. A testis-specific isoform would require a separate sequence-level decision. [Primary study](https://pubmed.ncbi.nlm.nih.gov/39319391/) |
| **ACTRT1** | CTA → Unsafe | Endogenous protein in differentiating normal human epidermal keratinocytes, supported by multiple reagents and protein-band analysis. This is broader than permitted hair follicles. [Primary study](https://pmc.ncbi.nlm.nih.gov/articles/PMC8110579/) |
| **MAGEC2** | CTA → Neither | Historical Purkinje-cell staining conflicts with the curated normal-negative profile. Resolve polyclonal and paralog specificity before unrestricted use. [Primary study](https://pubmed.ncbi.nlm.nih.gov/12920247/) |
| **LGALS13** | CTA → Neither | Extra-placental antibody reactivity conflicts with placental RNA restriction and later galectin work. Resolve related-galectin specificity and adult versus fetal compartments. [Early study](https://pubmed.ncbi.nlm.nih.gov/10527825/), [later study](https://pubmed.ncbi.nlm.nih.gov/19497882/) |
| **PLCZ1** | CTA → Neither | Historical NYD-SP27 work reports endogenous human pancreatic-duct expression. Full-paper and isoform/reagent adjudication remain necessary. Another omission from the previous conflict list. [Primary study](https://pubmed.ncbi.nlm.nih.gov/14697805/) |

If all eight proposals were accepted, the mapped totals would become **391 CTA,
564 Neither and 1,585 Unsafe**. The public API, FASTA files and k-mer backgrounds
still use the released categories. None of these expression studies alone shows
that a particular therapy injured the corresponding tissue.

## Important findings that should not become automatic exclusions

| Finding | Why it matters | Appropriate response |
|---|---|---|
| **180 retained gene assignments / 1,707 primary-antibody normal-staining rows** | Curated IHC can be negative while individual primary records stain non-permitted compartments | Review specificity and the curator's rationale; do not treat all raw stains as accepted protein |
| **103 retained assignments have Medium/High primary staining in Vital compartments** | The discrepancy includes consequential cell types | Prioritize orthogonal evidence and reagent identity |
| **78 of the 180 have explicit off-target/cross-reactivity warnings** | Some positive stains were deliberately discounted | A blanket raw-IHC veto would discard legitimate targets and duplicate family-level observations |
| **PATE2, TPTE2, EQTN** | Historical normal-tissue expression is isoform-dependent | Resolve transcripts/proteins before gene-wide or k-mer claims |
| **ACTL7A and SPAM1** | Additional historical non-reproductive RNA not highlighted in the 37-gene list | Quantify and resolve identity; qualitative RT-PCR is not a current nTPM value or proof of protein |
| **SUN5, SLC22A14, FAM170A, HORMAD2 and other recorded conflicts** | Older tissue studies and modern measurements differ | Preserve both observations and distinguish species, RNA/protein and assay specificity |
| **PASD1** | RNA in histologically normal tissue from cancer patients | Relevant nonmalignant context, but not healthy-volunteer protein evidence |

The antibody counts are gene assignments, not independent proteins, donors or
clinical injuries. Shared family reagents can recur across loci. In particular,
SAGE1 profiles explicitly describe presumed off-target binding being disregarded.
An “Enhanced” antibody/profile label does not make every raw tissue stain
gene-specific, nor does it prove all unmeasured tissues are negative.

## Excluded genes worth revisiting

| Gene(s) | Current result | Assessment |
|---|---|---|
| **ROPN1** | Unsafe / High RNA | Salivary serous-cell RNA is 32.9 nTPM, close to the hard cutoff. Review source identity, ambient RNA and protein before any rescue. |
| **TDRD1, BRDT** | Neither / Cellular | Modest immune-RNA signals near review thresholds; worth targeted validation. Not demonstrated clinical toxicity. |
| **HORMAD1, DAZL** | Neither / Cellular | Repeated immune-RNA observations make a simple rounding rescue less persuasive; examine donors, mapping and protein. |
| **ACRBP, TEX101** | Unsafe / High RNA | Multiple immune observations exceed 30 in source-native units. Keep excluded under the current policy pending source-level adjudication. |
| **FOXR2** | Neither / Antibody | Strong biological nomination but unresolved antibody conflict. Valid affirmative protein/recognition evidence can address low bulk RNA; it cannot erase uncertain normal staining. |
| **CTAG2** | Neither / Variable | Keep the reviewed donor-dependent vital-tissue concern separate from a claim of proven myocardial toxicity. |
| **XAGE5** | Neither / Antibody | No gene-specific rescue follows from XAGE1 family interest; neural staining and low reproductive RNA remain unresolved. |
| **MAGEA3** | Unsafe / High RNA | Current colon-cell RNA exclusion is separate from historic receptor cross-reactivity with titin or MAGE-family proteins. |
| **MAGEA12** | Unsafe / Safety | Retain the independent neurological-safety hold; RNA-threshold relaxation is insufficient. |
| **SPAG9** | Unsafe / Protein | Published cancer antigenicity does not override broad normal protein and immune RNA. |

All threshold-boundary genes, including ones not named here, are searchable in
the explorer. A near-boundary result is a reason to review measurement quality,
not an automatic exemption.

## What the rules still miss

1. **Published normal protein has not been comprehensively curated.** FTMT is a
   concrete example of old tissue summaries missing later human protein data.
   The new decisions belong in the existing reviewed normal-protein mechanism,
   with assay, species, compartment, identity and contradictory evidence attached.
2. **Rare cells and biological states.** Bulk averages dilute neurons, ductal
   cells, keratinocyte differentiation states and activated immune populations.
   Current immune/corneal/HPA cell sources help, but no negative baseline atlas
   proves silence during inflammation, repair or treatment.
3. **A high reproductive maximum can mask diffuse normal expression.** Passing
   relative enrichment is not evidence of narrow distribution. Some retained
   genes have substantial normal cell RNA below the 30 cutoff. Bulk ceilings and
   cell-state review rules are distinct; a cell nTPM value must not silently be
   treated as a matched bulk measurement. Review the complete distribution and
   coverage, not only a single maximum or a thresholded count of tissues.
4. **Normal peptide presentation is missing from this gene-level audit.** For
   TCR targeting, endogenous peptide/HLA evidence is closer to the intervention
   than total RNA or IHC. The [HLA Ligand Atlas primary study](https://pmc.ncbi.nlm.nih.gov/articles/PMC8054196/)
   provides a benign-tissue reference; lack of detection remains limited by
   tissue/HLA coverage and assay sensitivity. This audit has not ingested that
   dataset or certified its absence calls for these candidates.
5. **Therapy evidence is product-specific.** PRAME's IMA203 study includes tumor
   and normal peptide-presentation work and clinical outcomes; it supports that
   defined peptide/receptor context, not every PRAME-derived sequence. [Primary trial](https://pubmed.ncbi.nlm.nih.gov/40205198/)
   MAGE-A4-directed treatment similarly supports a defined clinical use, not
   automatic acceptance of every disputed antibody stain or all MAGE paralogs.
   [Primary trial](https://pubmed.ncbi.nlm.nih.gov/38554725/)
6. **Cancer occurrence and actionable target abundance are separate.** Of the
   399 retained genes, the existing report shows at least one cancer RNA sample
   at ≥5 pTPM for 392; only 124 have Medium/High cancer IHC in the displayed
   matched groups. Another 107 have matched cancer IHC without such staining,
   and 168 lack exact-ID cancer IHC coverage. A single positive tumor is not
   broad indication coverage, clonality, peptide presentation or efficacy.
7. **Gene, proteoform and epitope decisions differ.** Identical XAGE1A/B or
   CTAG1A/B proteins do not provide locus-specific expression proof. Different
   untranslated transcripts can encode the same protein. Exact k-mer exclusions
   cannot by themselves rule out TCR recognition of a different peptide.

## Recommended next curation pass

Apply the five proposed normal-protein exclusions and three provisional holds
through the existing reviewed-evidence workflow after inspecting the recorded
limitations. Then prioritize the Vital antibody discrepancies, isoform conflicts
and boundary exclusions. Preserve a candidate-level status until normal protein,
presented peptide and product-specific safety have been adequately assessed.

Do not tune a universal threshold merely to preserve a familiar CTA name. Keep
the literature conflict, measurement uncertainty and any proposed rescue visible
alongside the unchanged released decision.
