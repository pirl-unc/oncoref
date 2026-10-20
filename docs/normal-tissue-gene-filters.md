# Patient-context gene restriction filters

These are named **gene filtering policies**, separate from anatomical tissue
sets and from the general adaptive CTA screen. No disposable-tissue exception
is included. Both require positive evidence in their selected reproductive subset.

| Machine name | Permitted positive normal expression | Mandatory exclusion |
|---|---|---|
| `male_reproductive_maternal_fetal_excluded` | Explicitly postnatal male reproductive contexts in male samples | All female/maternal tissues, placenta, other extraembryonic contexts, embryo and fetus, including fetal testis |
| `female_developmental_male_reproductive_excluded` | Female reproductive contexts in female samples; embryonic/fetal and extraembryonic contexts | Male reproductive anatomy at every stage and every other postnatal context; known male mammary expression is excluded |

The first is the **pregnancy-oriented candidate filter**. The second is the
**male reproductive preservation candidate filter**. Neither name asserts that
a therapy is safe in pregnancy or preserves sperm/fertility. The second is
deliberately conservative about fetal male reproductive anatomy: its explicit
male-anatomy veto takes precedence over the larger fetal union.

The broad Reproductive set remains the union of all reproductive/developmental
contexts. Testis plus oocyte can fit that set while failing both narrower profiles.
Mammary tissue's anatomical grouping does not establish its sample sex.

## Positive compatibility versus supported filtering

`evaluate_normal_tissue_gene_filter()` distinguishes:

- `excluded`: a measured positive lies outside the selected policy.
- `insufficient_evidence`: no counterexample establishes exclusion, but a positive
  or required negative context/metadata is missing or unusable.
- `passes_declared_panel`: an allowed positive exists and every supplied required
  negative context is assessed and negative, with the minimum domains represented.

The result always reports `clinical_safety="not_assessed"`.
It also records the filter policy version, taxonomy fingerprint and a SHA-256
fingerprint of the exact required coverage panel.

The caller fixes a named coverage panel independently of the gene. The pregnancy
profile requires female somatic, female reproductive, pregnant maternal, embryonic,
fetal, fetal male gonadal, fetal female gonadal, placental and other extraembryonic
negative coverage. The male-preservation profile requires male reproductive,
male somatic and male mammary coverage. **These are minimum domain checks, not an
exhaustive panel**: the declared panel must specify relevant organs, routine cell
types, ages/stages and modalities. One sample per domain cannot demonstrate
biological absence throughout that domain.

Each negative is matched to the exact species, source/version, native
tissue/cell labels, source stage, sex, physiological state, context/sample ID,
assay and detection policy. Negative RNA cannot replace missing IHC; negative
adult tissue cannot clear a fetal sample; one dataset cannot clear another.
Every supplied positive overrides negatives, including positives outside the
declared panel. Unassessed observations prevent passing.

No validated pregnancy coverage panel or corresponding retained gene list is
bundled. The current adult HPA CTA report does not have the necessary maternal
and fetal coverage and must not be presented as such a list. Fetal expression of
MAGE-A1, GAGE and NY-ESO-1 has been observed in developing gonads; adult testis
restriction does not resolve this gap.
[Primary fetal-gonad study](https://pubmed.ncbi.nlm.nih.gov/17208940/)

## API

```python
from oncoref.normal_tissues import resolve_normal_tissue
from oncoref.normal_tissue_gene_filters import (
    NormalTissueEvidence,
    NormalTissueCoverageRequirement,
    normal_tissue_gene_filter_profiles,
    evaluate_normal_tissue_gene_filter,
)

normal_tissue_gene_filter_profiles()
annotation = resolve_normal_tissue(
    "testis", species="human", stage="adult", sex="male",
)
positive = NormalTissueEvidence(
    gene_id="ENSEMBL:example",
    annotation=annotation,
    context_id="reviewed-source:sample-1",
    assay="RNA",
    detection_policy_id="reviewed-RNA-rule-v1",
    result="detected",
    evidence_id="reviewed-source:record-1",
)
decision = evaluate_normal_tissue_gene_filter(
    [positive],
    profile="male_reproductive_maternal_fetal_excluded",
)
assert decision["positive_pattern_compatible"]
assert decision["status"] == "insufficient_evidence"
```

`NormalTissueCoverageRequirement` uses the annotation, context ID, assay and
detection-policy ID of each required negative context. Supply these as
`required_contexts` with a `coverage_panel_id`. Assay-specific code supplies
`detected`, `not_detected` or `unassessed`; absent rows, invalid values and unreliable
IHC must never be converted to negative calls. This API does not silently reuse
the general CTA screen's permissive RNA ceilings to assert absence.

## Datasets and species

Anatomical definitions are shared; source mappings and expression evidence remain
species-specific. Dog lacks the mapped seminal-vesicle/bulbourethral concepts;
mouse coagulating gland is a prostate descendant. Uterine horns occur in the
dog/mouse inventory. Mouse tonsil and dog/mouse vermiform appendix are inapplicable,
not measured negatives. Unmapped head/foot/soft-tissue contexts cannot be narrowed
by guessing. No mouse expression atlas is bundled.

The evaluator rejects mixed species, genes and taxonomy versions. Orthology does
not transfer expression absence. Canine expression sample annotations now retain
source sex where supplied; unknown source sex remains unspecified.
