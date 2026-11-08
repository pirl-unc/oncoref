#!/usr/bin/env python3
"""Explain primary translation qualification and every adjacent-PR change."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from build_cta_contract_report import digest, write_table  # noqa: E402

from oncoref import cta_selection, translation_qualification  # noqa: E402

OUTPUT = ROOT / "docs/audits/translation-qualification-20261104"


def build():
    meta = json.loads((ROOT / "oncoref/data/translation-qualification.json").read_text())
    facts = pd.read_csv(
        ROOT / "docs/audits/cta-selection-20261028/protein-group-classifications.csv.gz",
        keep_default_na=False,
        low_memory=False,
    )
    columns = [
        "protein_equivalence_group_id",
        "protein_group_label",
        "member_gene_ids",
        "classification",
        "assessment_status",
        "translation_qualified",
        "translation_qualification_reasons",
        "representative_sequence_sha256s",
        "translation_qualification_policy",
        "representative_selection_policy",
        "representative_selection_approved",
        "representative_selection_reasons",
        "representative_candidate_sha256s",
        "representative_evidence_ids",
        "representative_transcripts",
    ]
    write_table(facts[columns], OUTPUT / "representatives.csv.gz")
    baseline = json.loads((OUTPUT / "baseline.json").read_text())
    path = OUTPUT / baseline["path"]
    if digest(path) != baseline["sha256"]:
        raise ValueError("Changed translation baseline")
    before = pd.read_csv(path, keep_default_na=False)
    after = cta_selection.sequence_dispositions()
    after["reasons"] = after.reasons.map(lambda x: json.dumps(list(x)))
    tracked = ["classification", "assessment_status", "shared_sequence", "reasons"]
    delta = before[["entity_id", "label", *tracked]].merge(
        after[["entity_id", *tracked]],
        on="entity_id",
        suffixes=("_before", "_after"),
        validate="one_to_one",
        how="outer",
    )
    flags = []
    for c in tracked:
        flag = c + "_changed"
        delta[flag] = delta[c + "_before"].ne(delta[c + "_after"])
        flags.append(flag)
    write_table(delta[delta[flags].any(axis=1)], OUTPUT / "changes-from-evidence-pr.csv")
    records = translation_qualification.records()
    cases = records[records.symbol.isin(["SERPINA2", "MROH5", "SELENOV", "BOLL", "PRAME"])]
    write_table(cases, OUTPUT / "worked-records.csv.gz")
    qualified = int(facts.translation_qualified.sum())
    text = f"""# Translation qualification

Policy `{meta["policy_id"]}` assesses **{meta["records"]:,} Ensembl 115 peptide records**;
**{meta["qualified_records"]:,}** pass the record-level target rule. Across the
{len(facts):,} classified proteins, **{qualified:,} selected representatives qualify**.
This is annotation-based admission, not experimental proof of production.

Every candidate record receives a classification. The former entry requirement
for coding annotation or established ncORF support is now a CTA admission gate;
{int((~facts.translation_support_eligible).sum()):,} records fail that requirement
and remain assessed for normal evidence. [Complete inventory audit](../cta-inventory-20261105/index.md).

## What qualifies

A target record must have an exact source gene/transcript/protein identity, a
primary chromosome locus, a coding transcript without LoF or NMD annotation,
annotated CDS and complete three-base start and stop codons, no incomplete tags,
and an uninterrupted AA20 sequence beginning with methionine. One terminal `*`
is permitted. MANE Select supplies the representative when present on the primary
locus; otherwise Ensembl canonical does. Every member must share the same approved
sequence. Conflicting designated sequences or unresolved coding evidence require
review. A shorter isoform or an alternate locus cannot silently replace a failed
representative. [Selection method and reviewed exceptions](../cta-transcript-selection-20261108/index.md).

Missing mappings, incomplete annotation, NMD, LoF and interrupted sequence create
**translation review**, never normal-tissue Unsafe evidence. Existing Unsafe
findings retain precedence. All normal holds and target-gate reasons remain in
the audit even when the displayed subtype is Translation or Nested.

These conservative record rules deliberately send selenocysteine-containing
sequences (for example SELENOV) to review under the AA20 consumer contract. That
is a sequence-support limitation, not a claim that selenoproteins are untranslated.
Likewise an absent annotated start codon is a reference-qualification limitation,
not proof that no protein can be made.

## Worked records

**SERPINA2:** primary `ENST00000553483` / `ENSP00000486888` encodes a 309-residue
LoF-annotated, incomplete reference record. The 421-residue coding record lies
on `HSCHR14_7_CTG1` and cannot rescue it. Its current result is
`{facts.set_index("protein_group_label").loc["SERPINA2", "classification"]}`.

**MROH5:** the primary 1,318-residue LoF record is interrupted; two short NMD
records remain separately visible. Its current result remains **Unsafe** through
normal evidence. Translation review cannot erase those exclusions.

The [worked transcript records](worked-records.csv.gz) preserve source versions,
chromosome, locus, biotype, tags, codon/CDS counts and sequence hashes. The full
record table is available offline through `translation_qualification.records()`.

## Changes and downstream use

Against the pinned evidence PR: **{int(delta.classification_changed.sum()):,}
classification-label changes**, **{int(delta.shared_sequence_changed.sum()):,}
sequence-role changes**, and **{int(delta.reasons_changed.sum()):,} reason changes**.
Label changes include Other subtype changes; inspect category transitions in the
[series-wide report](../cta-definition-20261104/index.md).
[Every changed record](changes-from-evidence-pr.csv) and
[all representative checks](representatives.csv.gz) are downloadable.

Qualification does not delete the original peptide. The verified sequence bundle
retains all records and coordinate-bearing AA20 segments, including exclusions
from incomplete or alternate records. Segments never bridge an X, internal stop,
selenocysteine or different record. An ncORF nomination is independent evidence;
it cannot validate an unrelated NMD reference peptide. This source is the Ensembl
peptide inventory and does not claim complete empirical ncORF sequence coverage.

Source: [Ensembl release 115](https://ftp.ensembl.org/pub/release-115/), with the exact
GTF, peptide, cDNA and ncRNA hashes in the [qualification receipt](qualification.json).
Reference identity: `{meta["reference"]["reference_id"]}`.

```python
from oncoref import translation_qualification, cta_selection
records = translation_qualification.records()
records[records.symbol.eq("SERPINA2")]
cta_selection.sequence_dispositions()
```
"""
    (OUTPUT / "index.md").write_text(text)
    (OUTPUT / "qualification.json").write_text(json.dumps(meta, indent=2) + "\n")
    inputs = [
        Path(__file__),
        ROOT / "oncoref/translation_qualification.py",
        ROOT / "oncoref/data/translation-qualification.json",
        ROOT / "oncoref/data/translation-qualification.csv.gz",
        ROOT / "oncoref/data/cta-selection-human.csv.gz",
    ]
    (OUTPUT / "manifest.json").write_text(
        json.dumps(
            {
                "inputs": {str(p.relative_to(ROOT)): digest(p) for p in inputs},
                "outputs": {
                    p.name: digest(p) for p in sorted(OUTPUT.iterdir()) if p.name != "manifest.json"
                },
            },
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    build()
