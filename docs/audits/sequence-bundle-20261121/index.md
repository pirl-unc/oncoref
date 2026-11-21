# Verified sequence bundle (`target_selection` policy)

The current bundle preserves **245,535 original Ensembl 115 peptide records**
and **245,535 source occurrences**, bound to `oncoref.target_selection.results()`.
It contains **230,584 exclude**, **13,864 review** and **1,087 allow** records;
**570 records** are direct CTA candidates (`650` Cancer Antigen records,
`731` Reproductive Antigen records) before combining all matching occurrences.
Record counts are not protein-group counts.

[Loader API and examples](../../sequence-bundle.md) · [Manifest](manifest.json)
· [Coverage](coverage.json) · [2026-11-04 `cta_selection` bundle audit](../sequence-bundle-20261104/index.md).

Selection policy: `oncoref.target-selection.2026-11-20.2`.
Role policy: `oncoref.target-selection-sequence-roles.2026-11-21.1`.
Translation policy: `oncoref.translation-qualification.2026-11-04.1`.
Representative policy: `oncoref.transcript-selection.2026-11-08.1`.
Bundle identity: `oncoref:sequence-bundle:v1:sha256:938da7bdfa3903f582e086797288199d42d6f449b561fad7fc4780eab4badb9e`.
Reference identity: `openvax:reference:v1:sha256:ba6f6186a301a3ef9217ea0199ec3b851f9da2eb43558b7ba7e0d6fcd9cd85e9`.
Native ORFdata identity: `dbd1dd752e24ea08f4d787e48fe04e3464ebce248190572782d6e356e010822f`.

## Coverage and role invariants

- **Zero `Unsafe` allowances:** Every record with an `Unsafe` locus occurrence receives `shared_sequence == "exclude"`.
- **Complete candidate representative coverage:** All `327` `CTA`, `69` `Cancer Antigen only`, and `128` `Reproductive Antigen only` protein groups have at least one qualified `shared_sequence == "allow"` record matching `representative_sequence_sha256s`.
- **Unresolved and unscreened-home holds:** `Other:Review`, `Other:Sequence`, `Other:Unscreened tissue` (olfactory epithelium, nasal mucosa, taste bud, hair follicle, dental pulp, bone), nonrepresentative isoforms, and source loci outside the classified protein universe receive `shared_sequence == "review"`.

## Reproduction

```sh
python scripts/build_cta_sequence_bundle.py SOURCE BUNDLE --release 1.46.0
```
