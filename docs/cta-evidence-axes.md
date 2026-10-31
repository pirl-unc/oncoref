# Safety and CTA identity

**Only normal-tissue risk evidence can make a protein Unsafe.** CTA admission also
requires cancer expression, Reproductive evidence and adequate coverage.

![Independent CTA requirements](images/cta-evidence-axes.png)

[One-page PDF](images/cta-evidence-axes.pdf) · [Actual thresholds](cta-selection.md)

## Current decision logic

```text
Qualifying normal-tissue or applicable clinical safety exclusion -> Unsafe

Acceptable safety AND one qualifying cancer sample
AND a qualifying Reproductive route after concentration assessment
AND complete required coverage -> CTA

Otherwise: complete normal screen passes -> one named Other category
Remaining missing normal evidence or normal-screen holds -> Unassessed
```

Reproductive means reproductive, placental, extraembryonic, embryonic and fetal.
Named branches retain anatomical and developmental precision. The current measured
routes are full mapped HPA reproductive bulk RNA ≥5 nTPM, supported reproductive
cell RNA ≥5 nTPM with ≥20 reads, reliable positive reproductive IHC,
or supported Tyser/Petropoulos embryonic
lineage mean ≥5 TPM (≥20 cells, ≥5 positive cells, ≥10% positive at ≥1 TPM).
The current inputs do not cover every fetal stage.

Cancer requires one sample ≥2 TPM or HPA pTPM, or reviewed reliable tumor protein,
without a normal-expression ratio.
This establishes qualifying RNA, not necessarily translated protein, presentation
or expression in a specified fraction of patients.

Required coverage includes the fixed [GTEx v11](https://gtexportal.org/home/downloads/adult-gtex)
52-context bulk panel, every member locus, plus the existing normal bulk/cell/immune
panels. Optional DICE/cornea absence does not block admission; observed counterevidence
still counts. A missing measurement never becomes zero.

## Safety remains an independent assessment

Vital/Protected above-ceiling RNA or qualifying protein evidence can establish
Unsafe. Manageable excess blocks CTA without automatically establishing Unsafe.
Negative IHC adjusts the RNA ceiling only with the required antibody and local-cell
coverage. Missing or unrelated staining is not negative evidence. The exact numeric
rules remain in the [current policy table](cta-selection.md#uniform-human-rules).

Reviewed TCR-T clinical experience can add one RNA ladder step, capped at the
Enhanced-negative ceiling. It does not change the IHC grade or any cancer,
Reproductive, coverage or diffuseness criterion. The exact epitope/HLA/product
limitations and serious-event attribution remain attached to the protein-level
credit. [Clinical evidence and SAE review](audits/cta-tcr-trials-20261029/index.md).

## Profile concentration

The old rule requiring a Reproductive value above the highest normal measurement
has been removed. A single normal outlier should not define whether a protein is
diffuse. The score averages `N / (N + R)` over a fixed anatomy-weighted
normal panel, including Disposable tissues. N is that tissue's RNA; R is the source's
Reproductive reference. It is **not a percentage of tissues**.

The **HPA v23 >10% boundary is active** for completely measured adult/placental
profiles with a positive RNA reference and a qualifying RNA or IHC association.
RNA-only association still requires 5 nTPM. A zero RNA reference leaves concentration
unassessed even when IHC qualifies. Disposable tissues count;
**Permissive thymus is excluded from both the comparator and the reference**.
The fixed panel is 38 contexts / 29 anatomy units. A failed route with no other
qualifying Reproductive route and clear safety becomes Other:Diffuse. MTCP1 is one example
(24.53%). Independent Unsafe assignments do not change.

A weak adult reference does not decide embryonic identity. Supported embryonic
routes remain eligible with concentration explicitly unassessed until compatible
developmental/adult data exist. Missing HPA comparisons remain unassessed; the
GTEx core remains mandatory. [Method and measured impacts](audits/cta-unassigned-20261029/diffuseness-impact.md).

MAGE-A4, NY-ESO-1, PRAME, CTAG2 and XAGE1A/B retain CTA classification and GTEx
coverage. Removing the old maximum veto alone would admit TGIF2-RAB5IF and ATF7-NPFF;
both lack the required GTEx core and remain unassigned for that explicit reason.

All entries called **Proteins** collapse identical selected reference sequences;
originating loci remain listed. This is not a catalog of every alternative isoform
or modified proteoform. Human measurements and thresholds do not transfer to dog
or mouse through orthology.

## Other categories

| Category | Proteins | Definition |
|---|---:|---|
| Other:Diffuse | 8 | Complete HPA profile is dispersed (score >10% with a reference ≥2 nTPM), or the active concentration gate fails; no CTA route remains. |
| Other:Reproductive | 2 | Qualifying Reproductive evidence, but no qualifying cancer expression; the normal screen passes. |
| Other:Inducible | 5 | Cited expression in activated nonreproductive normal cells; a localized indirect atlas signal cannot establish Reproductive identity. |
| Other:Tumor-supported | 66 | Two distinct tumors ≥10 native TPM or reviewed reliable tumor protein, with low measured normal/developmental expression; not established tumor restriction. |
| Other:Disposable | 2 | Positive Disposable RNA/protein; measured RNA elsewhere below the descriptive 2-unit floor and no detected protein elsewhere, including low measured developmental RNA. |
| Other:Low | 254 | All measured normal RNA below the descriptive 2-unit floor and no detected normal protein. Cancer expression need not be low. |
| Other:Normal | 171 | Supported nonreproductive normal RNA reaches the descriptive 2-unit floor or reliable protein is detected; no preceding subtype or CTA definition qualifies. |
| Other:Weak | 92 | Reproductive signals are below admission requirements or lack support; no preceding subtype qualifies and no supported appreciable nonreproductive expression is established. |

Categories follow the table's priority order and are mutually exclusive. See the
[complete rules](cta-selection.md#categorical-assignments) for indirect reproductive
cell evidence, descriptive floors, source support and the categorical API.
