# Brachyury: antibody agreement and adult-tissue evidence

**TBXT / ENSG00000164458**, reviewed 2026-10-27. Other antibodies do not
consistently reproduce the three HPA003322-positive normal compartments.
Low RNA expression in some B-cell populations and independent adult-tissue
reports prevent interpreting this disagreement as proof of universal absence.
The [current funnel](../cta-harmonized-20261027/index.md) selects TBXT under the
strict RNA fallback after removing the gene-wide antibody veto for every gene.
This review supplies context and flags, not a per-gene admission waiver.

## What generates the existing flag

The pinned HPA v23 primary data use HPA003322, a rabbit polyclonal reagent,
with **Uncertain** IHC validation and profile reliability. Medium staining is
reported in marrow hematopoietic cells, lung macrophages and lymph-node
non-germinal-centre cells. HPA reports very low agreement with RNA expression.
The separate curated normal-IHC table contains no accepted TBXT observations;
an empty table is missing accepted evidence, not a negative tissue panel.
The current HPA tissue page still lists HPA003322 and Uncertain reliability.
[HPA primary data](https://v23.proteinatlas.org/ENSG00000164458/tissue/primary+data),
[current HPA tissue page](https://www.proteinatlas.org/ENSG00000164458-TBXT/tissue).

## Compartment-by-compartment comparison

RNA values below come from the **pinned v23 source**, not a mixture of HPA
releases. Displayed zeros do not prove absolute biological absence. IHC and RNA
were not obtained from the same specimens, and the atlas clusters do not exhaust
every cell state within an organ.

| HPA003322-positive compartment | Other protein evidence | Pinned RNA evidence | Interpretation |
|---|---|---|---|
| Lung macrophages | 54-1 study reports no normal-lung staining and specifically describes HPA Prestige macrophage staining as background. HMV321 supplier panel is lung-negative. | Lung bulk 0.0 nTPM; all seven lung macrophage clusters 0 reads / 0.0 nTPM. | Strongest case for a reagent-specific HPA signal, although no same-donor multi-antibody adjudication was found. |
| Marrow hematopoietic cells | HMV321 supplier panel is marrow-negative. No verified second-antibody reproduction of the exact HPA compartment was found. | Bulk 0.0 nTPM; all ten sampled marrow clusters 0 reads / 0.0 nTPM. A separate 2016 qPCR study reports marrow expression. | Mixed RNA evidence; the medium IHC signal is not independently established. |
| Lymph-node non-germinal-centre cells | HMV321 supplier panel is node-negative; a trial protocol also summarizes negative normal-node IHC. | Bulk 0.2 nTPM; all seven sampled node clusters 0 reads / 0.0 nTPM. Separate blood B-cell populations have low measurable RNA. | Partial support for low-level lymphoid transcription, not confirmation of medium protein staining in this precise compartment. |

The direct 54-1 comparison is in
[Hamilton et al., 2015](https://www.oncotarget.com/article/3086/text/).
The whole-organ supplier results are in the
[HMV321 normal-tissue panel](https://ms-validatedantibodies.com/product/brachyury-hmv321/).
The independent marrow qPCR observation is described in
[Jezkova et al., 2016](https://www.oncotarget.com/article/7202/text/), Figure 1;
relative qPCR fold change is not a TPM measurement. The
[NCT03595228 protocol](https://cdn.clinicaltrials.gov/large-docs/28/NCT03595228/Prot_000.pdf)
reports weak RT-PCR amplification in 4/9 healthy-donor B-cell samples alongside
negative normal-node and spleen IHC. This is a trial sponsor's evidence summary,
not a new independently replicated tissue-IHC series.

The matched blood summaries are HPA naive B cells **0.8 nTPM** (native TPM 0.3),
Monaco exhausted-memory B cells **0.8 TPM**, and DICE naive B cells
**0.832306 TPM at the third-highest reported value**. These quantities have
different normalization and summary definitions. Blood populations cannot be
substituted for exact lymph-node subcompartments. All remain below the current
protected/missing-IHC candidate ceiling of 4 in their respective native units.

## What the alternative antibodies establish

| Reagent | Evidence found | Limit |
|---|---|---|
| MAb 54-1 | Published knockdown-dependent signal, testing against related T-box proteins, and normal-lung negative staining. A separate series studied 5,229 tumors and reported negative nuclear staining in its examined non-testicular normal tissues. | Normal panel is not exhaustive; nuclear-only scoring differs from studies accepting cytoplasmic staining. Tumor case count is not a normal-donor count. |
| HMV321 | Supplier normal-tissue images, reported RNA agreement, and a second-antibody comparison supporting chordoma staining. | Supplier-generated evidence; the page does not name its comparison reagent. Negative staining is limited by assay sensitivity and tissue coverage. |
| EPR18113 / ab209665 | Recombinant monoclonal with supplier chordoma IHC and documented University of Washington diagnostic-laboratory validation. | UW's reported diagnostic performance concerns distinguishing chordoma from other neoplasms; it does not validate negative marrow, macrophage or lymph-node subsets. |

Sources: [54-1 specificity study](https://www.oncotarget.com/article/3086/text/),
[large 54-1 series](https://doi.org/10.1097/PAS.0000000000000462),
[HMV321](https://ms-validatedantibodies.com/product/brachyury-hmv321/),
[EPR18113 supplier data](https://www.abcam.com/en-us/products/primary-antibodies/brachyury-bry-antibody-epr18113-ab209665),
[UW laboratory validation record](https://testguide.labmed.uw.edu/guideline/ihc_updates).
These are evidence descriptions, not newly assigned HPA Approved/Supported/Enhanced grades.

## Adult expression outside the disputed HPA compartments

- **Intervertebral-disc nucleus pulposus:** a 2017 human study reports RNA
  measurements (60 specimens), IHC (56), and flow cytometry (5). IHC used
  ab57480; flow used R&D IC20851A. The authors acknowledge conflicting earlier
  protein reports. This is independent of HPA003322 and concerns a compartment
  absent from the routine HPA tissue panel.
  [Richardson et al., 2017](https://www.nature.com/articles/s41598-017-01567-w).
- **Intestinal enteroendocrine cells:** the 2016 study reports cytoplasmic staining
  in a subset of chromogranin-A-positive normal crypt cells, with an ab140661
  knockdown control and similar results with two other antibodies reported as
  data not shown. This differs from nuclear-only diagnostic scoring and is not
  equivalent to broad bowel expression.
  [Jezkova et al., 2016](https://www.oncotarget.com/article/7202/text/).
- **Thyroid:** an earlier study reports positive staining in 4/6 normal thyroid
  specimens, whereas the 54-1 nuclear-scoring series and HMV321 supplier panel
  report negatives. This remains an assay-dependent disagreement.
  [Hamilton et al., 2012](https://aacrjournals.org/clincancerres/article/18/14/3868/76935/Brachyury-a-Driver-of-the-Epithelial-Mesenchymal).

The supported interpretation is **disputed HPA immune-compartment staining,
with separate adult-expression evidence still requiring reconciliation**.
Neither the uncertain HPA signal nor a cleaner alternative antibody alone
establishes clinical toxicity or safety. TBXT is selected in the current funnel,
with these evidence limitations preserved. The older Neither assignment resulted
from the additional gene-wide antibody veto, which has been removed uniformly;
no negative-IHC threshold boost has been assigned to these alternative reagents.

## Reproduction

Run `.venv/bin/python scripts/audit_tbxt_antibodies.py`.
[Source receipts and output hashes](manifest.json),
[HPA antibody observations](hpa_antibody_xml.csv),
[bulk RNA](hpa_rna_consensus.csv),
[cell-cluster RNA](hpa_single_cell_tissue.csv),
[HPA blood RNA](hpa_immune_cell.csv),
[Monaco blood RNA](hpa_immune_monaco.csv),
[funnel RNA peaks](funnel-rna-peaks.csv).
