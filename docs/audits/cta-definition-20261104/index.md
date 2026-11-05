# CTA decision and sequence-role contract

Current selection policy: `oncoref.cta-selection.2026-11-04.2`. Consumer-role policy:
`oncoref.sequence-roles.2026-11-04.1`. These identities are separate because changing
uncertainty handling need not change a biological classification.

**19,863 input protein groups → 19,604
classified groups + 259 out-of-universe groups.**
The classified results are **602 CTA,
18,243 Unsafe and 759 Other**.
These counts use protein-equivalence groups, not literature nominations or patients.

An OncoRef CTA candidate is a reference-resolved protein with supported cancer expression and qualifying reproductive or developmental association that passes the named normal-tissue evidence, coverage and selection rules. Sequence realization requires separately qualified translation records.

## Evidence status and consumer roles

The categories alone do not determine allowance. An Other:Nonreproductive record
can also have incomplete normal evidence. Check `assessment_status`, identity and
the explicit sequence role. Unsafe exclusion takes precedence; otherwise unresolved
assessments require review. Only assessed, identity-resolved records receive `allow`.
Allowance means no exclusion from that record under this policy; other occurrences
can still exclude shared sequence.

| Role | Protein groups |
|---|---:|
| Exclude | 18,243 |
| Review | 339 |
| Allow under this policy | 1,022 |

Compared with the pinned 1.29.0 snapshot, **225
classifications and 339 roles change**.
[Every changed record](changes-from-1.29.0.csv) retains baseline/current values.
The [full role table](sequence-roles.csv.gz) preserves identity, reasons and policy.

## Public API

```python
from oncoref import cta_selection

definition = cta_selection.contract()
gates = cta_selection.gates()
roles = cta_selection.sequence_dispositions()
coverage = cta_selection.coverage()
```

The [gate catalog](gate-catalog.csv) identifies independent questions; it does not
claim every assay is measured. The current [decision audit](../../cta-selection.md)
retains source observations and overlapping reasons. [Definitions](contract.json),
[coverage](coverage.json) and [reproduction receipt](manifest.json) are machine readable.

These are representative-protein decisions. The upcoming sequence qualification
and verified-record loader in the [implementation roadmap](../../plans/cta-definition-roadmap.md)
must establish usable translation records and all source occurrences. A gene mapping,
CTA label or exact sequence match does not validate every isoform or an approximate
patient-HLA similarity threshold. Existing historical FASTA bundles keep their own
policy and reference identities.
