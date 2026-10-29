# Safety and CTA identity

**Only evidence of normal-tissue risk can make a target Unsafe.** A CTA must also
have cancer expression, adult restriction and reproductive/developmental
association. These are separate assessments with separate measurements and reasons.

![Four separate CTA requirements](images/cta-evidence-axes.png)

[One-page PDF](images/cta-evidence-axes.pdf)

## The decision logic

```text
Qualifying normal-tissue or applicable clinical safety exclusion -> Unsafe

Acceptable safety and adequate coverage
AND one qualifying cancer sample
AND reproductive or developmental association
AND adult restriction through a qualifying association route
-> CTA

Otherwise: measured Neither profile -> Neither; other cases -> unassigned
```

Failing cancer expression, association or restriction alone never supplies Unsafe
evidence. Safety is still evaluated for every target, including targets that
already failed another requirement. Missing evidence remains unknown.

## Safety: absolute tissue and protein evidence

Vital, Protected and Manageable describe the healthy compartment at risk. RNA
supports possible protein production; IHC adds protein evidence with antibody,
cell-coverage and specificity limitations. The [ceiling table](cta-selection.md#uniform-human-rules)
adjusts RNA tolerance using usable local negative IHC. Positive protein evidence
has its own intensity/confidence rules. Missing or unrelated IHC cannot become
negative evidence.

Above-ceiling Vital/Protected RNA or qualifying positive protein can establish
Unsafe. Manageable excess blocks CTA without automatically establishing Unsafe.
Low protein observations remain visible flags. These are research-policy
assignments, not measured probabilities of injury.

## Cancer: one qualifying sample

Require **at least one measured cancer sample ≥2 in native TPM or HPA pTPM**.
There is **no tumor-versus-normal ratio** in this gate. The existing 2-unit floor
establishes qualifying RNA expression, not necessarily high protein abundance.
Protein staining or antigen-presentation evidence must be identified separately.
One sample does not establish prevalence across cancers or patients, or malignant
cell origin in an admixed sample. Those remain descriptive evidence.

## Association: affirmative reproductive or developmental expression

Require **reproductive maximum ≥5 nTPM OR supported developmental lineage mean
≥5 TPM**. Developmental support requires ≥20 cells, ≥5 positive cells and ≥10%
positive at ≥1 TPM. Separate stage/lineage measurements retain their sources.
This requirement establishes presence in a qualifying context. Restriction is a
separate comparison.

## Restriction: relative to relevant adult normal expression

A qualifying association route must also satisfy its own native-unit comparison:

```text
reproductive nTPM > supported adult normal maximum nTPM
OR
supported developmental mean TPM > supported adult normal maximum TPM
```

Include supported, mapped **Vital, Protected and Manageable** compartments,
including bulk, cell and immune measurements. Exclude **Reproductive and
Disposable** from the adverse comparator. Comparing with every non-disposable
normal would incorrectly penalize expression in testis or placenta.

The normal maximum uses each adapter's supported statistic, such as a cluster
mean with enough reads, a GTEx tissue median, or DICE's third-highest reported
value. It is not the highest individual donor measurement in every source.
Required panel coverage still applies. Missing optional DICE/cornea data are
reported gaps, not zeros. All measured safety counterevidence still counts.

TPM compares with TPM and nTPM with nTPM. **CP10K retains its separate safety rule
and has no ratio to TPM/nTPM.** Native units avoid a direct mixed-unit ratio, but
source normalization and sampling still differ. This is a nominal enrichment
heuristic, not a calibrated fold change, tumor exclusivity or clinical safety.

Safety and restriction answer different questions. Reproductive 5 nTPM versus
normal 6 nTPM can pass a Protected absolute ceiling yet fail the reproductive
restriction route. A qualifying developmental route can still pass its own TPM
comparison. A restriction failure alone never establishes Unsafe.

## Named checks and boundaries

MAGE-A4, NY-ESO-1, PRAME, CTAG2 and XAGE1A/B pass the revised rules. Their clinical
products retain their own evidence and scope; this does not establish restriction
of every isoform, peptide or unmeasured normal compartment.

The independent restriction check removes TGIF2-RAB5IF and ATF7-NPFF from the
previous 548-group candidate set. Their supported developmental means are about
5.62 versus 5.80 normal TPM, and 5.12 versus 5.64 normal TPM, respectively. These
are close research boundaries, **not new Unsafe findings**. Their source values
and failed restriction reason remain in the audit for reconsideration.

The contract is species-independent. Numerical thresholds and measurement
coverage remain species-specific; human evidence does not transfer to dog or
mouse through orthology.

[Current results and actual thresholds](cta-selection.md) ·
[Full current audit](audits/cta-source-walkthrough-20261028/index.md)
