# Target filters and host selection

Catalog `$catalog_policy`; selection policy `$filter_policy`.
Filters narrow biological candidates. Optional host constraints are a separate axis
with their own versioned exposure evidence. [Full audit](audits/target-filters-20261112/index.md).

```python
from oncoref import antigen_sets, target_filters as tf
from oncoref.host_context import HostContext

selection = antigen_sets.select("CA", min_tumor_patients=2, min_tumor_fold=30)
# All unspecified compartments are protected; exceptions are explicit.
host = HostContext(absent=("embryo", "placenta"), accept_loss=("testis",))
with_host = antigen_sets.select("CTA", mode="strict", host=host)
decisions = with_host.decisions()
```

## Six ordered layers

| Layer | Input | Meaning |
|---|---|---|
| Off target | `max_off_target_tpm` | Maximum of N, Disposable and protected host TPM must be below the limit |
| Curation | `normal_function` | Optionally exclude recorded Vital, Protected or Manageable homes |
| Host | `HostContext` | Require coverage, bounded RNA in each native unit, resolved protein findings and cancer margin |
| Sets | `min_tumor_fold`, `min_reproductive_fold`, `min_tumor_patients`, `membership` | Tighten biological enrichment and known-patient replication |
| Papers | `min_papers` | Require scoped nominations or tumor reports |
| Exclusions | `exclude_prefixes`, `except_labels` | Explicit name-based judgments |

Every layer is evaluated; `filter_failed` records the first failure. Missing values
cannot pass an active requirement. `filtered_category` describes tightened biological
membership; `host_status` describes the independent host result.

Raising the tumor fold alone can move a CTA to RA-only in the filtered biological
view. Raising the reproductive fold alone can move it to CA-only. A **host rejection**
never performs that conversion: the protein is excluded from that host selection.
Two known patients at a 30× threshold means each patient's maximum reaches 30×,
including the enlarged protected-host denominator when one is used.

## Explicit host assumptions

`HostContext` distinguishes `absent` and `accept_loss`; every other compartment is
protected. Default host limits are **<4 native TPM and <4 nTPM**, independently.
These are research cutoffs, not injury thresholds. A stricter `max_off_target_tpm`
also applies. Missing coverage or unresolved protein findings prevents eligibility.

Compatibility profiles `male-patient` and `female-patient` assume anatomy and no
current pregnancy/lactation; they protect present gonads. Inspect their recorded
assumptions and override actual anatomy. `all` means no host screen. The former
`core` host shortcut is rejected; specify accepted loss explicitly.
[Host contract and limitations](host-selection.md).

| Scenario | Protected compartments | Retained CTA | CA only | RA only | Union |
|---|---|---:|---:|---:|---:|
$profile_rows

Counts are research-screen outputs. Pregnancy remains unresolved because early
embryo measurements do not cover all fetal/developmental stages. Milk-program
state labels are inferred and do not prove expression absence in other hosts.

## Modes

| Mode | Retained proteins without a host context | Changes |
|---|---:|---|
$mode_rows

Modes nest from all to strict. Baseline positive controls are regression checks;
any of them may fail a host requirement. Paper counts measure scoped published
attention, not independent experimental replication or evidence of tolerability.
Strict retains the explicit MAGE-prefix exclusion with MAGEA4 exception. Sequence
subtraction complements this judgment; it does not establish TCR cross-reactivity.

## Inspecting and subtracting

`selection.decisions()` returns every catalog row with membership, host status and
filter reasons. `selection.proteins()`, `.gene_ids()` and `.sequences()` share that
same selection. `.background_index()` includes global Unsafe/Other plus candidates
excluded or unresolved for the host. A host-rejected CTA receives no peptide exemption.

`tf.funnel(...)` reports sequential losses and `tf.describe(...)` states active
requirements. [Sequence API](antigen-sequences.md) · [CLI](curation.md).


## Normal nose and skin, which the expression panels do not sample

Recorded home tissues remain visible even when reproductive enrichment qualifies.
A host profile does not fill unmeasured normal-tissue gaps. See the
[home-tissue rule](target-selection.md#home-tissue-the-panels-do-not-resolve).
