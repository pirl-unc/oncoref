# Published placental and cancer-placenta nominations

Oncoref retains the complete published memberships separately from the HPA CTA
candidate table and final default panel. Source membership never overrides the
normal-tissue restriction, expression or specificity gates.

| Publication and source | Published | Canonically mapped | Currently coding | Family eligible | HPA pass | Default |
|---|---:|---:|---:|---:|---:|---:|
| [Gong 2021](https://doi.org/10.1038/s41467-021-22695-y), Supplementary Data 5 | 71 | 70 | 69 | 68 | 31 | 26 |
| Gong 2021, Supplementary Data 6 | 74 | 65 | 1 | 1 | 1 | 1 |
| [Bradley 2020](https://doi.org/10.1038/s41467-020-19141-w), Figure 3 | 10 | 10 | 10 | 10 | 1 | 1 |

These source lists overlap and must not be summed. They contain 73 distinct
currently protein-coding genes, 31 already in the candidate table and 42 new.
The expanded candidate set has 439 rows; the public default has 298 genes,
including five new ones: INSL4, GCM1, CYP19A1, HTRA4 and KISS1. This is the
September 23, 2026 source-addition snapshot, for oncoref 1.8.202.

In oncoref 1.8.205, TRIM64 was moved to candidate-only status after its nomination
and specificity evidence were reviewed. The current public default has 297 genes;
the raw candidate table and the publication-source overlap counts above are
unchanged. See the [TRIM64 nomination audit](audits/trim64-nomination.md) and
[issue #560](https://github.com/pirl-unc/oncoref/issues/560).

## Evidence scope

Gong, *The RNA landscape of the human placenta in health and disease*, supplies
placental RNA enrichment, not tumor-antigen validation. The source annotations
are retained: ERVH48-1 was noncoding in S6 and is currently protein-coding;
DSCR4 was protein-coding in S5 and is currently lncRNA. TXNRD3NB from S5 and nine
S6 entries cannot be mapped by their original Ensembl IDs and remain explicitly
unmapped. Historical symbols are not used to guess replacement identities.

The Gong lists support 16/19 prior `placental_antigen` nominations and 8/9
retained genes. CGB1, CGB2 and CGB7 are not in these lists. The original internal
nomination tags are preserved alongside the new paper-specific tags.

Bradley, *Vestigial-like 1 is a shared targetable cancer-placenta antigen expressed
by pancreatic and basal-like breast cancers*, nominates VGLL1, PLAC1, CGB3, CGB5,
IGF2BP3, DEPDC1B, ADAM12, SLC38A9, CAPN6 and MMP11 in Figure 3. HLA-peptide and
antigen-specific T-cell validation from this study is recorded for VGLL1 only.
Only PLAC1 passes the current default gate; VGLL1 fails the HPA normal-tissue
restriction gate and remains available as a sourced candidate.

## Data and APIs

- `get_data("cta-publication-sources")`: exact publications, tables/figures,
  URLs, evidence scope, input hashes and canonical reference hash/release.
- `get_data("cta-publication-membership")`: all 155 source rows, including
  unmapped/noncoding entries, source identities and current canonical identities.
- `oncoref.cta_sources.intake_membership()` / `intake_counts()`: auditable
  mapping, coding, family, HPA and public-default stages.
- `get_data("cancer-testis-antigens")`: new coding candidates with source tags
  `Gong2021_placenta_PC`, `Gong2021_placenta_ncRNA` and `Bradley2020_CPA`.

Gong's input hash covers the original XLSX bytes. Bradley's hash covers the
transcribed Figure 3 gene-symbol list in the recorded order (UTF-8, one symbol
per line with a final newline); it is not a hash of the article PDF.

## Reproduce the import

Download the original Gong supplementary workbook using `GONG_URL` in
`scripts/import_cta_publications.py`, then run in an environment with oncoref,
openpyxl and the HPA v23 normal RNA/IHC inputs available:

```sh
python scripts/import_cta_publications.py --gong-xlsx /path/to/supplement.xlsx
```

The importer retains every source row, resolves against oncoref's canonical gene
space, and adds only mapped, currently protein-coding genes to the HPA table.
Existing evidence/filter decisions are unchanged; only new candidates have HPA
columns regenerated. Missing RNA stays missing and fails the gate. Unknown
canonical transcript IDs, full gene names and functional descriptions for new
candidates remain unassigned rather than being inferred. Re-running the import
is idempotent.

Pirlygenes renders the source intersections and intake/default funnels from
these owner tables, with CSV audits, 300 dpi PNGs and vector PDFs. The local
source-specific run records an explicit checkout import; it does not change
pirlygenes' released dependency pin or silently rewrite its earlier complete
figure batch.

PAGE4 now has a row in the HPA candidate table and is excluded from the default.
`cta_candidate_references()` therefore omits it from the pending watchlist;
`cta_candidate_references(include_in_table=True)` retains its historical
literature references and annotations for provenance and clinical/audit joins.
The historical watchlist reports prostate (166.6 nTPM); the canonical CTA
somatic-max scope excludes accessory reproductive tissues and reports smooth
muscle (23.3 nTPM). These are different tissue scopes, not a loss of the prostate
observation.

## Closing the CGB provenance gap

The two nomination lists leave CGB1, CGB2 and CGB7 uncovered. Targeted primary
papers now give all **19/19 prior placental nominations**, including **9/9
retained genes**, publication provenance. This statement means expression or
nomination support, not that every gene is a validated tumor antigen.

| Source | Genes used here | Evidence and resolution |
|---|---|---|
| [Rull and Laan 2005](https://doi.org/10.1093/humrep/dei261), Methods, Figures 2/4, Table II | CGB1, CGB2, CGB7 | Placental RT-PCR followed by gene-discriminating restriction digestion; CGB2 is distinguishable from CGB1 |
| [Rull et al. 2008](https://doi.org/10.1093/molehr/gam082), Figure 2/Table II | CGB1, CGB2 | Sensitive combined CGB1/CGB2 RNA assay in trophoblastic and normal tissues |
| [Kubiczak et al. 2013](https://doi.org/10.3390/ijms140612650), Figure 3/Table 1 | CGB1, CGB2 | Combined CGB1-2 RNA signal in ovarian cancers and term placentas; 13/32 cancer positives belongs to the shared assay |
| [Bialas et al. 2020](https://doi.org/10.3390/genes11091082), Figures 1-3/Tables 1-2 | CGB1, CGB2 | Separate TCGA gene-level RNA estimates; CGB2 reported in 26 of 33 cancer types; no gene-specific protein-product validation |
| [McKellar et al. 2025](https://doi.org/10.1101/2025.05.28.656535), Figures 1/2C, Table 1 | CGB7 | Preprint: tumor RNA and CGB7-specific qPCR in SCaBER; immune associations, not HLA/T-cell antigen validation |

`cta-gene-publication-evidence` records these ten **targeted evidence rows**;
it is not a transcription of every gene investigated in each article. Each row
retains a source anchor, assay gene scope, review status, limitations and NCBI
BioC input checksum. The raw full texts are retrieved through NCBI's public
BioC endpoint; they are not redistributed in the package.

`oncoref.cta_sources.gene_publication_evidence()` exposes these rows;
`placental_source_coverage()` reconciles the historical list to publications and
the default set. The importer attaches their tags without changing HPA values.
CGB2 stays in the default; CGB1/CGB7 stay excluded. No CGB-specific protein,
HLA-peptide or T-cell validation is asserted from this evidence. In particular,
combined assays and beta-hCG family protein measurements cannot establish a
unique CGB2 protein or epitope.
