# Verified sequence bundle

The current bundle preserves **245,535 original Ensembl 115 peptide records**
and **245,535 source occurrences**. It contains
**230,648 exclude**, **13,512 review**
and **1,375 allow** records;
**903 records** are direct-target candidates before combining
all matching occurrences. Record counts are not protein-group counts.

[Loader API and examples](../../sequence-bundle.md) · [Manifest](manifest.json)
· [Coverage](coverage.json) · [Historical transition receipt](historical-transition.json)
· [Exact full-sequence comparison](historical-current-sequence-roles.csv.gz).

Selection policy: `oncoref.cta-selection.2026-11-05.1`.
Role policy: `oncoref.sequence-roles.2026-11-04.1`.
Translation policy: `oncoref.translation-qualification.2026-11-04.1`.
Bundle identity: `oncoref:sequence-bundle:v1:sha256:607c4ede9af40e33657ce84b58eff11087e66be7942a9a156a81480e8d06335a`.
Reference identity: `openvax:reference:v1:sha256:ba6f6186a301a3ef9217ea0199ec3b851f9da2eb43558b7ba7e0d6fcd9cd85e9`.
Native ORFdata identity: `dbd1dd752e24ea08f4d787e48fe04e3464ebce248190572782d6e356e010822f`.

## Coverage is explicit

The peptide inventory has **3,770 source loci
without a current decision**; these require review. **0
classified loci have no translation in this inventory**. The coverage payload retains
**0 records outside the classification path** and
**1,063 ncORF nominations**
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
python scripts/build_cta_sequence_bundle.py SOURCE BUNDLE --release 1.31.0
python scripts/build_cta_sequence_transition.py HISTORICAL BUNDLE
```

The builder checks source/reference identity, complete qualification-table coverage,
original sequence hashes, all payload hashes and counts, and an offline ORFdata
round trip. The final release rebuild binds the packaged data version to the
published ZIP pin. Historical assets retain their own explicit policy identity;
the comparison uses exact full-sequence hashes and never substitutes a partial
segment or approximate match for full equivalence.
