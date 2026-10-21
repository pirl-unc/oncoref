# Complete normal-tissue registry — 20 October 2026

Policy `2026-10-20.4`. Every canonical concept is listed for human, dog and mouse, including explicitly inapplicable concepts. Every bundled native source/version/tissue/cell mapping is exported separately.

| species | reproductive | disposable | Manageable risk | High consequence | vital | unresolved | inapplicable | native_mapping_rows | native_sources | mapping_availability |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| human | 46 | 5 | 6 | 14 | 45 | 3 | 2 | 459 | 3 | bundled |
| dog | 45 | 4 | 6 | 14 | 45 | 3 | 4 | 126 | 6 | bundled |
| mouse | 48 | 3 | 6 | 14 | 45 | 3 | 2 | 0 | 0 | not_bundled |

Counts are concepts and mapping rows, **not samples or measured negative coverage**. The same source tissue can have several cell-type rows. Mouse source mappings are not bundled; no mouse expression data or absence calls are invented.

## Exact inventories

- [All species and canonical concepts](canonical-anatomy.csv)
- [Every native dataset mapping, category and citation](source-mappings.csv)
- [Taxonomy and source-artifact receipts](taxonomy.json)
- [Patient-context gene-filter definitions](gene-filter-profiles.json)
- [Checksums and implementation provenance](run-manifest.json)

## Interpretation

Reproductive, Disposable, Manageable risk, High consequence and Vital are the five resolved categories. Unresolved is a QC state. Disposable denotes the narrowed conditional project allowance, whose compatibility scope ID is `permitted_normal`.

Uterine horns are applicable to dog/mouse; coagulating gland is a mouse prostate descendant. Dog seminal vesicles/bulbourethral glands, mouse tonsil, and dog/mouse vermiform appendix are explicitly inapplicable. Applicability is checked through anatomical parents. Broad head, foot and soft-tissue labels stay unresolved.

Native mapping rows do not supply sex or developmental metadata that the source does not report. At sample annotation time, explicit fetal/embryonic context overrides adult organ categories, while pregnancy of the mother does not. Sex and physiological state are preserved separately. Gene filters require a declared negative panel and cannot pool species, sources, versions or assays.

See [definitions](../../normal-tissue-scopes.md), [manageability evidence](../../normal-tissue-manageability.md), and [patient-context gene filters](../../normal-tissue-gene-filters.md).

Reproduce with `python scripts/normal_tissue_registry_report.py --out <empty-directory>` at this report's implementation revision.
