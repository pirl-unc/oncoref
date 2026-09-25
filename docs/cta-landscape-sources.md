# Complete CTA landscape intake

[Regenerated figures and counts](audits/cta-landscape-20260924/index.md) ·
[Combined vector PDF](audits/cta-landscape-20260924/oncoref-cta-landscape-figures.pdf)

The September 24, 2026 rebuild imports complete published candidate sets before
normal-tissue filtering. It expands the protein-coding candidate table from 439
to **2,499** loci. All 297 previously retained genes remain; 327 additional genes
pass the default rules, giving **624** retained loci. The 439 existing rows keep
their identity and HPA evidence values; only their source tags change.

These are source nominations and expression evidence, not a claim that every
gene has a uniquely measured protein, presented HLA peptide, or validated immune
response. Gene counts also do not count distinct protein sequences.

## Complete source sets

| Primary publication | Exact source and selection | Source entries |
|---|---|---:|
| [Wang 2016](https://doi.org/10.1038/ncomms10499) | Supplementary Data 3A, `ncomms10499-s4.xlsx`: every Ensembl-labelled coding CT row, including genes outside the narrower EECTG subset | 1,019 |
| [Bruggeman 2018](https://doi.org/10.1038/s41388-018-0357-2) | Supplementary Tables, `41388_2018_357_MOESM4_ESM.xlsx`, sheet `1D`: all GC genes; sheet `2` marks the TGCT-dependent subset | 756 |
| [da Silva 2017](https://doi.org/10.18632/oncotarget.21715) | Supplementary Table 1, `oncotarget-08-92966-s002.xls`, `1103 CTs predicted`: all normal-testis-biased nominations | 1,103 |
| da Silva 2017 | Supplementary Table 2, `oncotarget-08-92966-s003.xls`, `Sheet1`: expressed in **at least 10%** of a tumor type | 745 |
| [Jamin 2021](https://doi.org/10.1002/1878-0261.12900) | Table S6, `MOL2-15-3003-s007.xlsx`, `Expression data and classes`: SET/SEHET/PET/PEHET probes with a nonempty, non-dash UCNDH classification | 607 identities from 602 probes |
| Jamin 2021 | Table S2, `MOL2-15-3003-s003.xlsx`, `Annotation`: every Ensembl-labelled core CT row | 125 |
| [Chang 2019](https://doi.org/10.1002/cam4.2223) | Table S2, `CAM4-8-3511-s002.xlsx`, `Suppl Table 2`: C1, protein_coding, ratio1 at least 1% | 1,036 |
| [Carter 2023](https://doi.org/10.1136/jitc-2023-007935) | Supplementary Table 2, `jitc-2023-007935supp001.pdf`, PDF pages 17–27: CT-expression true and thymus-expression false | 103 |
| [Seager 2024](https://doi.org/10.1186/s12967-024-04918-0) | Methods, Data processing and statistical analysis: all 17 targets in the profiling panel | 17 |

The existing complete Gong and Bradley placental lists and historical resource
nominations are retained. Nested lists and overlapping papers must not be summed.
Noncoding-only CT-RNA screens are outside this protein-antigen intake; entries
in the selected sets that are currently noncoding or unmapped remain in the
source audit and are excluded at the corresponding funnel stage.

Jamin's 602 probes annotate 551 distinct named symbols and 56 unnamed probes,
giving 607 source identities. Every semicolon-separated member of a shared probe
is retained with its full probe annotation. This differs from the article's
headline 478 genes. Likewise, the core annotation table contains 125 Ensembl
rows, whereas the article reports 124 genes. Neither discrepancy is silently
trimmed. Mapped aliases collapse to distinct canonical loci in the plots.

Chang's evidence is specific to testicular germ-cell tumors. Bruggeman's
TGCT-dependent nominations, including NLRP9 and ZFP42, retain that qualification.
Seager is a targeted profiling study, not a genome-wide discovery screen; its
panel includes MLANA, a differentiation antigen that fails the normal-tissue
gate. Its historical BAGE and GAGE2 labels remain unmapped rather than being
assigned to guessed paralogs. SUN5 and ZFP42 are not in Carter's selected CT set.

## Mapping and filtering

Source Ensembl/Entrez identifiers use the shipped canonical reference maps.
An unknown identifier never falls back to a similarly named gene. Symbol-only
sources use unambiguous canonical symbols or shipped NCBI synonyms. Original
identities, exact spreadsheet rows/PDF pages, assay scope and mapping outcomes
are retained in `cta-publication-membership.csv`. A canonical locus lacking an
approved symbol uses its Ensembl ID as its display label.

| Stage | Remaining | Removed at this stage |
|---|---:|---:|
| Source union, including unresolved/noncoding identities | 2,842 | — |
| Canonical ID mapped | 2,694 | 148 |
| Protein-coding candidates | 2,499 | 195 |
| Family exclusions applied | 2,485 | 14 |
| HPA normal-tissue restriction | 850 | 1,635 |
| Specificity, expression and cancer-nomination requirements | 624 | 226 |

The existing HPA thresholds, low-expression rescue and reviewed exclusions are
preserved. Four newly nominated genes absent from HPA RNA retain missing values
and fail the restriction gate. TRIM64, CTAG2, CSAG1 and SPAG4 remain outside the
default set; SPAG4 remains in the separate candidate catalog with its pancreas
caveat. CSAG1 now also has a row in the HPA-assessed candidate table.

The expanded intake includes 160 genes whose sole source is da Silva's broad
normal-testis screen. These remain `candidate_only`: normal-testis restriction
alone does not establish cancer expression. This prevents 35 otherwise passing
genes from entering the default set. CT45A8 and CT47A8 retain their earlier
paralog-based nominations; their new broad-screen citations do not establish
locus-specific tumor expression. This rebuild does not resolve every evidence
qualification in issue #559.

## Reproduction and artifacts

Install `.[reports]`. Download the original supplementary archives from the
Europe PMC `source_url` values in `cta-publication-sources.csv`, extracting them
under `<input-dir>/PMC<id>-supp/`. Save Seager's full-text XML as
`<input-dir>/PMC10851610.xml`. The importer verifies all nine original-file
SHA-256 values in `oncoref.cta_landscape.INPUTS` before reading any selections.
Source metadata records the canonical-reference release/hash and retrieval date.
HPA candidate annotation uses the existing pinned v23 cache.

```sh
python scripts/import_cta_landscapes.py --input-dir /path/to/original-supplements
python scripts/render_cta_curation.py --out outputs/cta-landscape-current
```

The renderer writes all 11 figures as 300-dpi PNGs and vector PDFs, a combined
bookmarked PDF, exact source intersections, per-source and per-gene funnel
membership, candidate/default tables, and content hashes. The four Venn figure
files cover the three principal CT sets, historical resources, additional
landscapes/source groups, and placental sources. The all-source overlap matrix
includes every source. The pooled funnel starts before mapping and biotype
filtering; Venn diagrams count mapped coding candidates before HPA filtering.

This is a CTA source/curation rebuild. Tumor-cohort expression, survival,
coverage and clinical response analyses are separate reports and have not been
rerun for the expanded panel.
