#!/usr/bin/env python3
"""Generate the current decision contract and its separately counted role changes."""

from __future__ import annotations

import gzip
import hashlib
import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from oncoref import cta_selection  # noqa: E402

OUTPUT = ROOT / "docs/audits/cta-definition-20261104"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_table(frame, path):
    data = frame.copy()
    for column in data:
        if data[column].map(lambda value: isinstance(value, (tuple, list, dict))).any():
            data[column] = data[column].map(
                lambda value: (
                    json.dumps(value, sort_keys=True)
                    if isinstance(value, (tuple, list, dict))
                    else value
                )
            )
    content = data.to_csv(index=False, lineterminator="\n").encode()
    path.write_bytes(gzip.compress(content, mtime=0) if path.suffix == ".gz" else content)


def build():
    OUTPUT.mkdir(parents=True, exist_ok=True)
    baseline_receipt = json.loads((OUTPUT / "baseline.json").read_text())
    baseline_path = OUTPUT / baseline_receipt["snapshot"]
    if digest(baseline_path) != baseline_receipt["sha256"]:
        raise ValueError("The 1.29.0 baseline snapshot changed")
    baseline = pd.read_csv(baseline_path, keep_default_na=False)
    roles = cta_selection.sequence_dispositions()
    definition = cta_selection.contract()
    coverage = cta_selection.coverage()
    (OUTPUT / "contract.json").write_text(json.dumps(definition, indent=2) + "\n")
    (OUTPUT / "coverage.json").write_text(json.dumps(coverage, indent=2) + "\n")
    write_table(cta_selection.gates(), OUTPUT / "gate-catalog.csv")
    write_table(roles, OUTPUT / "sequence-roles.csv.gz")
    current = roles[
        ["entity_id", "label", "classification", "shared_sequence", "assessment_status"]
    ]
    delta = baseline[["entity_id", "classification", "shared_sequence"]].merge(
        current,
        on="entity_id",
        how="outer",
        suffixes=("_baseline", "_current"),
        validate="one_to_one",
        indicator=True,
    )
    delta["classification_changed"] = delta.classification_baseline.ne(
        delta.classification_current
    ) | delta._merge.ne("both")
    delta["role_changed"] = delta.shared_sequence_baseline.ne(
        delta.shared_sequence_current
    ) | delta._merge.ne("both")
    changed = delta[delta.classification_changed | delta.role_changed]
    write_table(changed, OUTPUT / "changes-from-1.29.0.csv")
    counts = roles.category.value_counts().to_dict()
    role_counts = roles.shared_sequence.value_counts().to_dict()
    universe = coverage["universe_counts"]
    text = f"""# CTA decision and sequence-role contract

[Read the illustrated, searchable report](report.html) · [Download the PDF](cta-definition.pdf)
· [Current report receipt](release-receipt.json) · [Full report data](report-data.json).

Current selection policy: `{coverage["policy_id"]}`. Consumer-role policy:
`{definition["role_policy_id"]}`. These identities are separate because changing
uncertainty handling need not change a biological classification.

**{universe["input_groups"]:,} candidate protein records → {universe["classified_groups"]:,}
classified proteins + {universe["out_of_universe_groups"]:,} out-of-universe records.**
The classified results are **{counts.get("CTA", 0):,} CTA,
{counts.get("Unsafe", 0):,} Unsafe and {counts.get("Other", 0):,} Other**.
Identical selected reference sequences count together; member gene loci remain
traceable. These are representative proteins, not a complete inventory of splice
isoforms or modified proteoforms. [Universe entry rules and excluded records](../../cta-selection.md#stage-0-universe).

{definition["definition"]}

## Evidence status and consumer roles

The categories alone do not determine allowance. An Other:Nonreproductive record
can also have incomplete normal evidence. Check `assessment_status`, identity and
the explicit sequence role. Unsafe exclusion takes precedence; otherwise unresolved
assessments require review. Only assessed, identity-resolved records receive `allow`.
Allowance means no exclusion from that record under this policy; other occurrences
can still exclude shared sequence.

| Role | Proteins |
|---|---:|
| Exclude | {role_counts.get("exclude", 0):,} |
| Review | {role_counts.get("review", 0):,} |
| Allow under this policy | {role_counts.get("allow", 0):,} |

Compared with the pinned 1.29.0 snapshot, **{int(delta.classification_changed.sum()):,}
classifications and {int(delta.role_changed.sum()):,} roles change**.
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

These are representative-protein decisions. [Translation qualification](../translation-qualification-20261104/index.md)
and the [verified-record loader](../../sequence-bundle.md)
retain usable translation records, original sequences and all source occurrences. A gene mapping,
CTA label or exact sequence match does not validate every isoform or an approximate
patient-HLA similarity threshold. Existing historical FASTA bundles keep their own
policy and reference identities.
"""
    (OUTPUT / "index.md").write_text(text)
    inputs = [
        ROOT / "oncoref/cta_contract.py",
        ROOT / "oncoref/cta_selection.py",
        ROOT / "oncoref/data/cta-selection.json",
        ROOT / "oncoref/data/cta-selection-human.csv.gz",
        Path(__file__),
        OUTPUT / "baseline.json",
        baseline_path,
    ]
    manifest = {
        "inputs": {str(p.relative_to(ROOT)): digest(p) for p in inputs},
        "outputs": {
            name: digest(OUTPUT / name)
            for name in [
                "contract.json",
                "coverage.json",
                "gate-catalog.csv",
                "sequence-roles.csv.gz",
                "changes-from-1.29.0.csv",
                "index.md",
            ]
        },
    }
    (OUTPUT / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(
        json.dumps(
            {
                "classification_changes": int(delta.classification_changed.sum()),
                "role_changes": int(delta.role_changed.sum()),
                "roles": role_counts,
            }
        )
    )


if __name__ == "__main__":
    build()
