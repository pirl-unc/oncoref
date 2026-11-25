# Patient profiles and candidate filters

Catalog `oncoref.target-selection.2026-11-25.1`; filter policy `oncoref.target-filters.2026-11-17.5`.
Filters narrow the biological candidate sets (`CA`, `RA`, `CTA`) across six explicit layers, while patient-anatomy constraints evaluate reproductive and mammary exposure independently without relabeling biological categories. [Antigen whitepaper](antigen-whitepaper.md) · [Interactive catalog](audits/target-selection-20261110/rebuild/report.html) · [Filter sweep tables](audits/target-filters-20261112/index.md).

```python
from oncoref import antigen_sets, target_filters as tf
from oncoref.host_context import HostContext

selection = antigen_sets.select("CA", min_tumor_patients=2, min_tumor_fold=30)
# No reproductive protection by default; request the smaller protected-tissue set.
host = HostContext(protect=("testis",))
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

`HostContext` distinguishes `protect`, `absent` and `accept_loss`; unspecified
compartments are not protected. Limits for explicitly protected tissues are
**<8 native TPM and <8 nTPM** independently (with **2×** = **<16** headroom in
human compartments where routine IHC is `not_detected`).
These are research cutoffs, not injury thresholds. A stricter `max_off_target_tpm`
also applies. Missing coverage or unresolved protein findings prevents eligibility.

Profiles `male-patient` and `female-patient` declare assumed anatomy with no current
pregnancy or lactation; neither preset protects reproductive tissue by default. Inspect recorded
assumptions and specify protected or absent organs explicitly. `all` applies no patient-anatomy filter.
[Patient-anatomy rules and limits](host-selection.md).

| Scenario | Protected compartments | Retained CTA | CA only | RA only | Union |
|---|---|---:|---:|---:|---:|
| `all` | none | **323** | 97 | 142 | 562 |
| `female-patient` | none | **323** | 97 | 142 | 562 |
| `female-patient (protect ovary)` | ovary | **293** | 95 | 137 | 525 |
| `female-patient (lactating)` | lactating_breast | **319** | 96 | 141 | 556 |
| `female-patient (pregnant)` | embryo, placenta, lactating_breast | **0** | 0 | 0 | 0 |
| `male-patient` | none | **323** | 97 | 142 | 562 |
| `male-patient (protect testis)` | testis | **74** | 73 | 29 | 176 |
| `male-patient (protect adult tissues)` | testis, epididymis, prostate, seminal vesicle, resting_breast | **58** | 63 | 16 | 137 |

Counts are research-screen outputs. Pregnancy remains unresolved because early
embryo measurements do not cover all fetal/developmental stages. Milk-program
state labels are inferred and do not prove expression absence in other hosts.

## Modes

| Mode | Retained proteins without a host context | Changes |
|---|---:|---|
| `all` | 562 | `{}` |
| `loose` | 396 | `{'min_papers': 1}` |
| `standard` | 253 | `{'max_off_target_tpm': 5.0, 'normal_function': 'no_vital', 'min_papers': 2}` |
| `strict` | 112 | `{'max_off_target_tpm': 4.0, 'normal_function': 'no_protected', 'membership': 'cta', 'min_papers': 3, 'exclude_prefixes': ('MAGE',), 'except_labels': ('MAGEA4',)}` |

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
