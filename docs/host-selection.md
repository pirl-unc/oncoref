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

## Plain-language sex protection vs. explicit `HostContext`

For standard patient queries, `antigens.table(..., sex="female" | "male")` (or `oncoref antigens list --sex female|male` on the CLI) automatically protects **all** non-pregnant, non-lactating adult reproductive and resting mammary organs of that biological sex in the chosen species and enforces `>=10×` tumor fold above protected-organ exposure (`258` of `270` `CTA-M + CTA-P` human CTAs for `sex="female"`; `58` of `60` `CTA-F + CTA-P` human CTAs for `sex="male"`), while letting you declare surgical absences (`absent=("ovary", "uterus")`) or accepted organ loss (`accept_loss=("prostate",)`):

```python
from oncoref import antigens

# Protects ovary, uterus, endometrium, fallopian tube, cervix, vagina, and resting_breast:
female_ctas = antigens.table("CTA", sex="female", ectopic_only=True)

# Maternal-fetal interface (scope_complete=True; 306/322 human CTAs pass):
maternal_interface = antigens.table("CTA", protect=("placenta", "lactating_breast"))

# Full pregnancy shortcut fails closed (0 rows) by default because embryo lacks 1st-3rd trimester coverage;
# pass allow_incomplete_scope=True to screen measured early-embryo + placenta + lactating RNA:
pregnant_default = antigens.table("CTA", sex="female", pregnant=True)
pregnant_screened = antigens.table("CTA", sex="female", pregnant=True, allow_incomplete_scope=True)
antigens.explain("LALBA", sex="female")  # excluded (resting_breast >= 8 nTPM)
antigens.explain("MAGEA4", sex="female")  # pass_host=True (CTA-M)
```

For custom RNA/nTPM thresholds (`max_tpm`, `max_ntpm`, `negative_ihc_multiplier`) or cross-species compartment inspection, use `HostContext`, `HostContext.for_host(...)`, `host_context.compartments(species)`, and `SEX_RESTRICTION_CODES`:

```python
from oncoref import antigen_sets, host_context
from oncoref.host_context import HostContext, SEX_RESTRICTION_CODES

dog_compartments = host_context.compartments("dog")
host_preset = HostContext.for_host("female-patient", lactating=True, species="dog")

# Illustrative declarations, not a patient recommendation.
host = HostContext(
    protect=("ovary", "endometrium"),
    max_tpm=8,
    max_ntpm=8,
)
selection = antigen_sets.select("CTA", mode="standard", host=host)
proteins = selection.proteins()
decisions = selection.decisions()
assumptions = selection.metadata
```

When constructing `HostContext` directly, every unspecified compartment is **not protected**.
`protect`, `absent` and `accept_loss` must be explicit, disjoint tuples.
Host protection can include vagina and either mammary state. `breast` expands to both. These conditional compartments
supplement the catalog's shared normal screen; they do not represent every organ.

Legacy compatibility presets `--profile male-patient` and `--profile female-patient` encode **assumed absent opposite-sex anatomy only**—unlike `--sex`, they protect **no** reproductive tissue by default. Prefer `--sex female` / `--sex male` (or `antigens.table(..., sex=...)`) when you want same-sex reproductive and mammary organs protected automatically. Surgery must be represented by verified absence (`absent=...`), not by
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
   `host.max_ntpm` (both default to 8, matching Stage 1 `Protected` organs, with
   `2×` = `16` headroom in human compartments where routine IHC is `not_detected`;
   non-human species without routine cell-resolved IHC validation default to `1×`).
   No TPM/nTPM conversion occurs.
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
