# Complete CTA landscape intake

[Regenerated figures and counts](audits/cta-unified-20260924/index.md) ·
[Combined vector PDF](audits/cta-unified-20260924/oncoref-cta-landscape-figures.pdf)

The September 24, 2026 rebuild combines source expansion (#562), the citation
implementation for #559, and the downstream plotting contract for pirlygenes
#629. The starting pool is the **full union of ten selected papers**, with
**2,537 mapped coding candidates** and **624 retained genes**. All 297 genes
retained in oncoref 1.8.205 remain; complete source expansion contributes 327.
The raw historical evidence table has 2,546 rows: nine legacy-only candidates
are preserved for audit but are outside the active paper-defined intake.

An exhaustive set-cover search covers all 624 screened survivors, counting each
DOI once. Among the **eleven fully imported candidate-source papers**, the unique
minimum contains ten: Wang, Bruggeman, da Silva, Jamin, Chang, Carter, Seager,
Gong, Loriot and Bai. Bradley adds no uniquely required survivor. Every chosen
paper has an indispensable retained-gene witness; the machine-readable
[certificate](audits/cta-unified-20260924/cta-minimum-source-cover.json) records
all witnesses and alternatives. This is a bounded minimum, not a claim about all
published literature. The seven-paper result for the earlier 297-gene panel is
superseded by this expanded-panel calculation.

Selection does not truncate any paper to its survivors. All entries from every
selected candidate list are loaded **before** canonical mapping, biotype and HPA
filters. The selected union reproduces the full expanded 624-gene endpoint.
Targeted papers with incomplete retrieved lists remain additional citations,
not eligible complete-list sources in the optimization.

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

| da Silva 2017 | Supplementary Table 3, `oncotarget-08-92966-s004.xlsx`, `Supplementary_Table3`: complete tumor-proteomics gene-labelled list | 136 |
| [Gong 2021](https://doi.org/10.1038/s41467-021-22695-y) | Data 5 and 6: all Testis, Placenta and Ovary rows, preserving the source biotype and tissue | 744 coding-labelled + 1,067 noncoding-labelled |
| [Loriot 2025](https://doi.org/10.1371/journal.pgen.1011734) | S1 strict cancer-germline and S2 preferential cancer-germline lists, all Ensembl-labelled rows | 146 + 134 |
| [Bai 2016](https://doi.org/10.1158/0008-5472.CAN-16-0225) | Abstract: the study's EGFL6 target; cancer/vascular expression and function, not a claim of normal-tissue exclusivity | 1 |

The existing placenta-only Gong subsets and complete Bradley Figure 3 list are
retained as audit views. Historical resource tags do not substitute for primary
paper citations in the active intake. Nested lists and overlapping papers must not be summed.
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
| Selected paper union, including unresolved/noncoding identities | 3,895 | — |
| Canonical ID mapped | 3,654 | 241 |
| Protein-coding candidates | 2,537 | 1,117 |
| Family exclusions applied | 2,523 | 14 |
| HPA normal-tissue restriction | 880 | 1,643 |
| Default specificity, expression and new-nomination rules | 624 | 256 |

The existing HPA thresholds, low-expression rescue and reviewed exclusions are
preserved. Six nominated genes absent from HPA RNA retain missing values
and fail the restriction gate. TRIM64, CTAG2, CSAG1 and SPAG4 remain outside the
default set; SPAG4 remains in the separate candidate catalog with its pancreas
caveat. CSAG1 now also has a row in the HPA-assessed candidate table.

New normal-reproductive-only nominations from da Silva S1 and Gong remain
`candidate_only`, including genes in both lists. Normal reproductive restriction
alone does not establish cancer expression. CT45A8 and CT47A8 retain their prior
paralog nominations; their broad-screen citations do not establish locus-specific
tumor expression. Source coverage does not overrule existing specificity holds.
The provenance table flags **29 retained legacy genes** whose selected-paper
links are normal-reproductive expression only. Their pre-existing resource or
paralog nominations retain the prior policy; the new citations do not resolve
their tumor-expression or locus-specific antigen evidence. This caveat is
explicit in `nomination_scope` and `legacy_retention_evidence_caveat`.

`oncoref.cta_provenance.candidate_provenance()` supplies one DOI/anchor record for
every active coding candidate; `selected_membership()` supplies the complete
row-level source and assay detail. Both are exported beside the figures. `gene_citation_evidence_report()` separately
indexes published-list nomination, normal-reproductive expression, cancer
expression nomination, protein, HLA-peptide and antigen-specific T-cell links.
Empty evidence categories mean no imported link, not a biological negative.
`gene_publication_evidence()` now includes the targeted #559 findings (including
SUN5 Figure 1C–E, SSX4B shared-probe limitations, and GAGE10/GAGE12D publisher-preview
limitations). A checksum-labelled copy of the earlier audited mappings is in
`docs/audits/cta-targeted-provenance-20260924.json`; its hash describes audit JSON,
not original full-text bytes. None of these imported expression links is promoted
to locus-specific protein, HLA peptide or T-cell validation.

The historical-tag audit compares CTexploreR labels with Loriot's full lists and
da Silva protein tags with the complete S3. In particular, S3 FAM71E2 maps to
GARIN5B, not GARIN4; MEIOC is present in S3. These exact memberships drive the
new plots. Historical CTpedia tags have no primary-paper status. They remain
visible only in the explicitly labelled historical audit Venn.

RNA fraction regeneration uses `math.fsum`, eliminating Python-version-dependent
rounding without changing thresholds. Expanded intake includes RPL39L, which
fails the HPA filter; normalization's curated housekeeping categories remain
unchanged and all default CTAs remain outside its censored set.

## Reproduction and artifacts

Install `.[reports]`. Download the original supplementary archives from the
Europe PMC `source_url` values in `cta-publication-sources.csv`, extracting them
under `<input-dir>/PMC<id>-supp/`. Save Seager's full-text XML as
`<input-dir>/PMC10851610.xml`. Download Gong and Loriot supplementary files and
Bai BioC XML to the filenames in `oncoref.cta_landscape.INPUTS`, using the explicit
`INPUT_URLS` for those inputs. The importer verifies all fourteen original-file
SHA-256 values in `oncoref.cta_landscape.INPUTS` before reading any selections.
Source metadata records the canonical-reference release/hash and retrieval date.
HPA candidate annotation uses the existing pinned v23 cache.

```sh
python scripts/import_cta_landscapes.py --input-dir /path/to/original-supplements
python scripts/render_cta_curation.py --out outputs/cta-landscape-current
```

The renderer writes all 12 figures as 300-dpi PNGs and vector PDFs, a combined
bookmarked PDF, exact source intersections, per-source and per-gene funnel
membership, candidate/default tables, minimum-cover certificate, gene-level citations, historical-tag reconciliation,
legacy-only archive, and content hashes. The four Venn figure
files cover the three principal CT sets, historical resources, additional
selected paper lists, and placental audit sources. The all-source overlap matrix
includes every selected paper. The pooled funnel starts before mapping and biotype
filtering; Venn diagrams count mapped coding candidates before HPA filtering.

This is a CTA source/curation rebuild. Tumor-cohort expression, survival,
coverage and clinical response analyses are separate reports and have not been
rerun for the expanded panel.
