# Patient anatomy and organ-loss filtering

**Biological origin, observed exposure, and tolerated tissue loss are separate axes.**
A CTA can fail a patient-profile screen without ceasing to be a biological CTA. These expression cutoffs support target screening; they are not clinical injury thresholds.

## Protection choices

| Choice | Meaning | What it does not establish |
|---|---|---|
| Not protected (default) | Apply no additional host restriction for this compartment | Does not mean absent, tolerable to lose, or clinically safe |
| Protect | Preserve the named tissue or physiological compartment | A low RNA value does not prove no injury |
| Accept loss | Explicitly accept loss of that compartment in the scenario | Acceptable lasting loss does not establish acceptable acute inflammation |
| Absent | Declare that the organ, pregnancy compartment or physiological state is absent/inactive | This does not erase expression elsewhere or sequence cross-reactivity |

A whole-organ label includes its support and endocrine cells. Bulk testis or ovary
RNA cannot justify assuming selective germ-cell injury. Review **acute loss, lasting
loss, repair and management** separately through
`tissue_tiers.consequences(...)` and `antigen_sets.loss_consequences(...)`.
Anatomy and tissue-loss choices do not override global Unsafe or Other findings.

## One explicit host object

```python
from oncoref import antigen_sets
from oncoref.host_context import HostContext

# Illustrative declarations, not a patient recommendation.
host = HostContext(
    protect=("ovary", "endometrium"),
    max_tpm=4,
    max_ntpm=4,
)
selection = antigen_sets.select("CTA", mode="standard", host=host)
proteins = selection.proteins()
decisions = selection.decisions()
assumptions = selection.metadata
```

Every unspecified compartment is **not protected**. `HostContext()` and both
anatomy presets retain the full candidate set at the chosen evidence mode.
`protect`, `absent` and `accept_loss` must be explicit, disjoint tuples.
Host protection can include vagina and either mammary state. `breast` expands to both. These conditional compartments
supplement the catalog's shared normal screen; they do not represent every organ.

Compatibility presets `male-patient` and `female-patient` encode **assumed anatomy**,
not anatomy inferred from sex or chromosomes. They protect no reproductive tissue by default, assume no
current pregnancy or lactation, and return those assumptions in metadata. Override
actual organs with `absent_tissues`, `protect_tissues` and `allow_tissues` (the legacy
name for accepting loss). Surgery must be represented by verified absence, not by
an automatic promise that no residual tissue exists. `all` means no host screen.

Use `reproductive_qualifying_compartments` for biological origin and `host_status` for the patient-profile decision.

## Evidence and requirements

Host exposure is built separately from RA admission. Finite, identity-resolved
observations inform exposure **including measured zero and
observations failing RA's positive-cell support rule**. Missing observations remain
missing. Flagged normal libraries retain their measured exposure and require QC
review; they cannot silently lower the host reference. A passed positive-evidence gate is not a negative-evidence coverage test.

For every protected compartment:

1. Native-TPM coverage and the declared developmental/state scope must be present.
2. TPM must be below `host.max_tpm`; available nTPM must separately be below
   `host.max_ntpm`. Both default to 4. No TPM/nTPM conversion occurs.
3. Positive reliable HPA protein findings require host-specific attribution review.
   They remain unresolved here; they do not silently rewrite global Unsafe.
4. For a CA meeting the requested cancer-evidence thresholds, cancer evidence must also satisfy the requested fold over
   **max(N, Disposable, protected-host TPM, 1)**. Known-patient replication uses
   the same enlarged denominator. Host rejection does not demote it to RA-only.

Tightening tumor fold or patient support can independently make a catalog CTA
RA-only. That is an evidence filter, not a host rejection, and does not add the
protein to the host subtraction background. A context protecting no compartments
leaves these biological decisions unchanged.

`max_off_target_tpm`, if requested, includes protected-host TPM. The catalog's
32 TPM Disposable cap never substitutes for this protection limit. Passing the
host screen is called **eligible**, not safe. Known failing values yield `excluded`;
otherwise incomplete evidence yields `unresolved`. All reasons remain visible,
even when a known failure takes precedence over uncertainty.

The packaged host receipt binds native-TPM observations, HPA v23 cell-cluster nTPM,
normal IHC, identity mapping and policy versions. Cell-cluster measurements do not
claim exhaustive cell coverage. Gene-level evidence does not prove translation of
every annotated isoform. Negative staining does not erase rare-cell signals.

## Mammary and developmental scope

The milk-program rule labels a sample `secretory_signature` when at least three of
LALBA, CSN1S1, CSN2 and CSN3 reach 100 TPM. Complete panels below that rule are
`no_high_secretory_signature`; missing markers are `unknown`. These are expression
inferences, **not provider-confirmed physiological states**. Legacy column names
`lactating_breast` and `resting_breast` remain for compatibility. An unknown-state
sample is applied to both states and prevents state-specific clearance.

The measured low values in three biopsies do not establish universal absence in
other donors, sexes or hormonal states. State-unspecified cell RNA and IHC apply to
both mammary states. Pregnancy protects the milk-secretory compartment too: mammary
differentiation changes during pregnancy, before postnatal feeding.
[HPA breast biology](https://www.proteinatlas.org/humanproteome/tissue/breast).

Petropoulos preimplantation and Tyser early gastrulation measurements do not cover
all fetal stages. Protecting the embryo compartment therefore remains unresolved,
even when its measured early-stage maximum is zero. A current nonpregnant assumption
also does not address future pregnancy during a persistent treatment.

Explicit `pregnant=True` and `lactating=True` are protection shortcuts for those
scenarios: pregnancy protects embryo, placenta and the milk-secretory compartment;
lactation protects the milk-secretory compartment. Other reproductive tissues
still require `protect_tissues` or `HostContext(protect=...)`.

Clinical experience with one target or product is not a universal host exemption.
For example, the MAGE-A4 product's label reports insufficient pregnancy data and
unknown fertility effects. [FDA TECELRA prescribing information](https://www.fda.gov/media/180565/download).

## One selection, consistent views and background

```python
sequences = selection.sequences()  # all annotated Ensembl 115 translations
background = selection.background_index("Unsafe|Other")
for sequence in sequences.sequence:
    retained_occurrences = background.subtract(sequence)
```

All views share one captured decision and its metadata. `proteins()` updates
category and membership booleans together for any tightened biological thresholds;
`catalog_*` retains the original values. `decisions()` retains all original catalog
rows, including host exclusions and their reasons.

The subtraction background combines the chosen global exclusions with **every
candidate excluded or unresolved by this host context**. Paper-count or label-only
filter losses are not automatically normal-tissue evidence. Missing gene translations
stop `background_index()`; inspect coverage rather than treating missing sequences
as absent. Exact, edit-distance and optional equivalence matching use the same
background. Sequence sharing, peptide presentation and receptor recognition remain
separate questions. [Peptide API and limits](antigen-sequences.md).

Selection assumptions are captured in DataFrame attributes and emitted as JSON to
stderr by the CLI. CSV alone does not preserve those attributes: retain the receipt
alongside exports. `selection.metadata` supplies the same receipt in Python.
