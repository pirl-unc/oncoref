# TRIM64 nomination history and publication review

Reviewed 2026-09-24. Gene: TRIM64 / C11orf28 / TRIM64A, ENSG00000204450, Entrez 120146. The history and HPA values below refer to oncoref 1.8.204. In 1.8.205, the specificity audit demotes TRIM64 to candidate-only under [issue #560](https://github.com/pirl-unc/oncoref/issues/560).

**TRIM64 has primary publications, including a dedicated RNA/protein/function study. The unresolved questions are its original CTA nomination provenance and reproductive/tumor specificity, not whether anybody has published on the gene.**

## How it entered our panel

1. **2024-09-10 — OpenVax gene-lists.** Commit [736a19e](https://github.com/openvax/gene-lists/commit/736a19ee3342431638e85c4d7d3ec9353b8c07ca), “added CTAs,” introduces TRIM64 in both the CSV and a hard-coded `CTAs` list in `notebooks/Generate CTA list CSV.ipynb`. The README describes an intersection of CTpedia with HPA antibody staining restricted to testis/placenta. There is no TRIM64-specific publication or CTdatabase entry identifier in the notebook. Thus the intended method is documented, but TRIM64's individual CTpedia membership is not established by that code.
2. **2024-09-11 — pirlygenes.** Commit [512336f](https://github.com/pirl-unc/pirlygenes/commit/512336f9d05af612608993b55a1e7612ece5c4cd), “migrating data from OpenVax gene-lists,” copies the list and notebook, including TRIM64.
3. **2026-03-23 — explicit source tag.** pirlygenes commit [6f75f20](https://github.com/pirl-unc/pirlygenes/commit/6f75f20) adds the `source_databases` column and labels the existing TRIM64 row `CTpedia`. The row previously lacked a source field. This annotation is not evidence of a newly verified gene-specific CTpedia match.
4. **2026-03-27 — ctabase, later tsarina.** Initial release [acd69c2](https://github.com/pirl-unc/tsarina/commit/acd69c2059b116a1c9478838df6beee76d85d053) imports the HPA-filtered list with TRIM64 already present and tagged.
5. **2026-06-09 — current reference owner.** oncoref's predecessor imports the existing CTA definition in [e1fc174](https://github.com/pirl-unc/oncoref/commit/e1fc17457a413cc5eeff04cd1f3f560c95f4bfca). This is inherited membership, not a new TRIM64 nomination based on Gong, Bradley, or a recent literature scan.

The CTdatabase snapshot checked for the full-panel audit does not reproduce a TRIM64 entry under TRIM64, TRIM64A, or C11orf28. Historical differences remain possible; call its CTpedia provenance **unverified**, rather than claiming that the database certainly never contained it.

## Why the frozen default filter retains it

The oncoref 1.8.204 row records placenta-restricted HPA protein expression, `Approved` protein reliability, and placenta RNA of **1.2 nTPM**, versus **0** recorded testis, ovary, and maximum somatic RNA. It consequently has reproductive RNA fractions of 1.0, `passes_filters=True`, `never_expressed=False`, and `restriction=PLACENTAL`.

Its positive protein annotation keeps it from being treated as never expressed despite low RNA. These values explain the filter outcome; they do not independently establish robust tumor expression or validate the historical CTpedia label. A fraction of 1.0 here should be interpreted alongside the low absolute RNA abundance.

## What has actually been published

| Study | Verified finding | Interpretation |
|---|---|---|
| [Han, Lou & Sawyer, PLOS Genetics 2011](https://doi.org/10.1371/journal.pgen.1002388) | Table 1 explicitly maps **B2 to TRIM64**, RefSeq NM_001136486. Figure S2 places TRIM64 in the chromosome 11q14.3 cluster and marks its existing GenBank transcript with a black star. | Gene/sequence and evolutionary characterization. Crucially, B2 is **not** one of the study's testis RT-PCR products in S2B; the other TRIM paralogs' direct testis results cannot be transferred to TRIM64. The downloaded S2 PDF was visually inspected. |
| [Geng et al., Developmental Cell 2012](https://doi.org/10.1016/j.devcel.2011.11.013) | Table 2 lists TRIM64 in a DUX4-induced program in human myoblasts. | Experimental expression context, not a direct test of cancer-placenta/testis specificity. |
| [Zhu et al., Cell Biology and Toxicology 2023; online 2022](https://doi.org/10.1007/s10565-022-09768-4) | Figure 1E-F reports TRIM64 RNA and protein induction by oxidized LDL in THP-1-derived macrophages; subsequent knockdown/overexpression experiments examine NF-kappaB/I-kappaB-alpha signaling. | Dedicated expression and functional evidence. THP-1 is a leukemia-derived line used here as a macrophage model: neither direct healthy-donor macrophage validation nor a clinical tumor-specificity study. It warrants assay/paralog and immune-cell-expression review. |
| [Gramantieri et al., European Journal of Immunology 2024](https://doi.org/10.1002/eji.202350637) | Figure 3 reports TRIM64 RNA in patient PBMCs during an HCC immunotherapy study. | Nonmalignant blood-cell compartment; do not interpret the cancer diagnosis as assigning this signal to tumor cells. |

There is therefore ample reason to retain a bibliography for TRIM64, while keeping its CTA nomination **provisional/unverified**. The seven-paper nomination/expression cover and its contextual eighth-paper example are not an exhaustive TRIM64 bibliography. Normal-tissue specificity, probe/antibody discrimination among close TRIM paralogs, and the original source tag still need reconciliation. The 1.8.205 change preserves the raw candidate and HPA filter values, and excludes TRIM64 from default, filtered, and placental CTA helpers through the existing specificity audit. It is not promoted to clinical targets or the warning tier. The default panel changes only by TRIM64: 298 to 297 genes. Explicit unfiltered and relaxed HPA-only exploratory views remain available; they are not the reviewed default set.
