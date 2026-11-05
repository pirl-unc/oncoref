#!/usr/bin/env python3
"""Generate the public sequence inventory and coverage explanation from its manifest."""

import json
from pathlib import Path

r = Path(__file__).resolve().parents[1]
out = r / "docs/audits/sequence-bundle-20261104"
m = json.loads((out / "manifest.json").read_text())
c = json.loads((out / "coverage.json").read_text())
text = f"""# Verified sequence bundle

The current bundle preserves **{m["record_count"]:,} original Ensembl 115 peptide records**
and **{m["occurrence_count"]:,} source occurrences**. It contains
**{m["record_roles"].get("exclude", 0):,} exclude**, **{m["record_roles"].get("review", 0):,} review**
and **{m["record_roles"].get("allow", 0):,} allow** records;
**{m["direct_target_records"]:,} records** are direct-target candidates before combining
all matching occurrences. Record counts are not protein-group counts.

[Loader API and examples](../../sequence-bundle.md) · [Manifest](manifest.json)
· [Coverage](coverage.json) · [Historical transition receipt](historical-transition.json)
· [Exact full-sequence comparison](historical-current-sequence-roles.csv.gz).

Selection policy: `{m["selection_policy_id"]}`.
Role policy: `{m["role_policy_id"]}`.
Translation policy: `{m["translation_policy_id"]}`.
Bundle identity: `{m["bundle_id"]}`.
Reference identity: `{m["reference"]["reference_id"]}`.
Native ORFdata identity: `{m["orfdata_identity"]}`.

## Coverage is explicit

The peptide inventory has **{len(c["source_gene_loci_without_current_decision"]):,} source loci
without a current decision**; these require review. **{len(c["classified_gene_loci_without_translation"]):,}
classified loci have no translation in this inventory**. The coverage payload retains
**{len(c["out_of_universe_groups"]):,} records outside the classification path** and
**{len(c["nc_orf_nominations_without_sequence_in_this_reference_bundle"]):,} ncORF nominations**
without claiming an empirical ncORF sequence inventory. These lists are different
coverage questions and must not be added together as a disjoint denominator.

Every original peptide, including incomplete/NMD/LoF and alternate-locus records,
remains present. Coordinate-bearing AA20 segments never bridge ambiguity or stops.
Exact raw-sequence grouping retains all records and occurrences. Exclude takes
precedence over review, then allow; an allowed CTA occurrence cannot cancel another
occurrence's exclusion. Missing records never supply implicit allowance.

## Reproduction

Build from the verified native ORFdata snapshot and exact reference manifest:

```sh
python scripts/build_cta_sequence_bundle.py SOURCE BUNDLE --release {m["release"]}
python scripts/build_cta_sequence_transition.py HISTORICAL BUNDLE
```

The builder checks source/reference identity, complete qualification-table coverage,
original sequence hashes, all payload hashes and counts, and an offline ORFdata
round trip. The final release rebuild binds the packaged data version to the
published ZIP pin. Historical assets retain their own explicit policy identity;
the comparison uses exact full-sequence hashes and never substitutes a partial
segment or approximate match for full equivalence.
"""
(out / "index.md").write_text(text)
