# Target filters

The catalog owns all evidence and category decisions. Filters only narrow its
candidate sets. Current catalog: `oncoref.target-selection.2026-11-17.1`; filter policy:
`oncoref.target-filters.2026-11-17.1`. [Full audit and sweep tables](audits/target-filters-20261112/index.md).

```python
from oncoref import target_filters as tf

all_candidates = tf.select("all")
replicated_ca = tf.select(tf.TargetFilter(
    membership="cancer_antigen", min_tumor_patients=2))
enriched_ctas = tf.select(tf.TargetFilter(
    membership="cta", min_tumor_fold=30, min_reproductive_fold=20))
decisions = tf.apply(tf.TargetFilter(min_papers=2))
```

## Five ordered layers

| Layer | Knob | Default | Meaning |
|---|---|---|---|
| Off target | `max_off_target_tpm` | None | Highest non-reproductive native TPM, including Disposable, must be below this value |
| Curation | `normal_function` | `any` | Optionally exclude recorded Vital, Protected or Manageable normal homes |
| Sets | `min_tumor_fold` | 10 | Tighten CA's tumor / max(N, Disposable, 1) |
| Sets | `min_reproductive_fold` | 10 | Independently tighten RA's reproductive observation / the same denominator |
| Sets | `min_tumor_patients` | 0 | Known patients with a qualifying sample at `min_tumor_fold`; anonymous samples do not count |
| Sets | `membership` | `target_candidates` | Keep union, CA, RA, CTA or either exclusive set after tightening |
| Papers | `min_papers` | 0 | Distinct scoped target nominations or tumor reports |
| Exclusions | `exclude_prefixes`, `except_labels` | Empty | Explicit name-based exclusions and exceptions |

Missing measurements never satisfy an active requirement. The fold knobs cannot
go below 10. A filter never promotes Other or Unsafe, changes the catalog, or adds
its own evidence. Every layer is evaluated; `filter_failed` names the first failure.
`filtered_category` reports membership after tightening.

Raising only the cancer fold can move a CTA to RA-only without removing it from the
union. Raising only the reproductive fold can move it to CA-only. Select `cta` when
both tightened memberships are required. The patient-count requirement tightens CA
only. Combining `min_tumor_patients=2` with `min_tumor_fold=30` requires two known
patients each reaching 30×. Multiple samples from one patient count once.

## Modes

| Mode | Kept proteins | Changes from the catalog defaults |
|---|---:|---|
| `all` | 509 | `{}` |
| `loose` | 370 | `{'min_papers': 1}` |
| `standard` | 247 | `{'max_off_target_tpm': 5.0, 'normal_function': 'no_vital', 'min_papers': 2}` |
| `strict` | 111 | `{'max_off_target_tpm': 4.0, 'normal_function': 'no_protected', 'membership': 'cta', 'min_papers': 3, 'exclude_prefixes': ('MAGE',), 'except_labels': ('MAGEA4',)}` |

Modes nest from all to strict. They are optional choices, not evidence grades.
PRAME, NY-ESO-1 and MAGEA4 are regression controls across these modes.

## What the measurements mean

`samples_at_tenfold`, `independent_patients_at_tenfold`,
`anonymous_samples_at_tenfold` and `cancer_types_at_tenfold` expose actual qualifying
tumor support. Unknown participant identity stays unknown. `*_above_normal` fields
describe a weaker historical threshold and must not be used as 10× witness counts.

`reproductive_fold` divides `reproductive_supported_tpm` by the same native-TPM
reference as CA. Small cell populations are limited by the second independent
embryo at the same stage and lineage; `reproductive_supported_context` identifies
that context. `reproductive_max_sample_tpm` retains the strongest cell-supported
individual observation, which can be higher. `reproductive_best` remains its
compatibility display alias. None of these fields pools TPM with nTPM.

<a id="normal-nose-and-skin-which-the-expression-panels-do-not-sample"></a>

The normal-function knob uses recorded tissues and tiers from shared curation.
[Unmeasured home tissues](target-selection.md#home-tissue-the-panels-do-not-resolve)
remain visible even when reproductive evidence qualifies.

## Papers and MAGE

Paper counts indicate published attention and scoped claims. Overlapping screens
are not independent experimental replication, and bibliographic support cannot
clear a normal-expression hold. Strong unresolved rare-cell signals now become
Other:Review in the catalog, before any paper filter.

Strict retains `exclude_prefixes=("MAGE",)` and `except_labels=("MAGEA4",)` as an
explicit judgment. It can remove MAGEC1 or MAGEC2 despite favorable expression.
Exact and bounded-distance [peptide subtraction](antigen-sequences.md) against
`Unsafe` or `Unsafe|Other` complements this rule. Sequence similarity alone does
not establish recognition by a particular TCR; the prefix rule is not silently replaced.

## Inspecting decisions

`tf.measures()` returns catalog evidence plus the quantities read by filters.
`tf.funnel("strict")` reports sequential losses; `tf.describe("strict")` names
every active knob. The [audit](audits/target-filters-20261112/index.md) provides
single-knob sweeps, every mode's members and overlapping failure reasons.
