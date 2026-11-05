# Translation qualification

Policy `oncoref.translation-qualification.2026-11-04.1` assesses **245,535 Ensembl 115 peptide records**;
**185,993** pass the record-level target rule. Across the
19,863 classified proteins, **19,075 selected representatives qualify**.
This is annotation-based admission, not experimental proof of production.

Every candidate record receives a classification. The former entry requirement
for coding annotation or established ncORF support is now a CTA admission gate;
259 records fail that requirement
and remain assessed for normal evidence. [Complete inventory audit](../cta-inventory-20261105/index.md).

## What qualifies

A target record must have an exact source gene/transcript/protein identity, a
primary chromosome locus, a coding transcript without LoF or NMD annotation,
annotated CDS and complete three-base start and stop codons, no incomplete tags,
and an uninterrupted AA20 sequence beginning with methionine. One terminal `*`
is permitted. The common longest primary-locus sequence must qualify for every
member of the expression-equivalence group. Ties retain all common sequence hashes.
A shorter isoform or an alternate locus does not silently validate that representative.

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
`Other:Translation`.

**MROH5:** the primary 1,318-residue LoF record is interrupted; two short NMD
records remain separately visible. Its current result remains **Unsafe** through
normal evidence. Translation review cannot erase those exclusions.

The [worked transcript records](worked-records.csv.gz) preserve source versions,
chromosome, locus, biotype, tags, codon/CDS counts and sequence hashes. The full
record table is available offline through `translation_qualification.records()`.

## Changes and downstream use

Against the pinned evidence PR: **115
classification-label changes**, **93
sequence-role changes**, and **788 reason changes**.
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
Reference identity: `openvax:reference:v1:sha256:ba6f6186a301a3ef9217ea0199ec3b851f9da2eb43558b7ba7e0d6fcd9cd85e9`.

```python
from oncoref import translation_qualification, cta_selection
records = translation_qualification.records()
records[records.symbol.eq("SERPINA2")]
cta_selection.sequence_dispositions()
```
