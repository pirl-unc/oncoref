# Normal-cell data and CTA rules

The [normal-cell report](audits/cta-normal-cells-20261023/index.md) combines HPA's bulk
and curated IHC screen with cell-resolved RNA, primary antibody records and cited
human protein evidence. The [24 October review](audits/cta-neither-review-20261024/index.md)
added an absolute >30 TPM/nTPM ceiling outside explicitly permitted contexts.
The [current definitions](target-sets.md) apply the reviewed literature hierarchy
and retain 213 CTA genes. The source-specific RNA triggers below are unchanged.

## One download and cache contract

```python
from oncoref import catalog, reference_data, normal_cell_data as cells

catalog.datasets()                   # HPA, DICE, corneal atlas, other references
catalog.ensure("dice_cd4_stim")       # same datacache transfer wrapper as HPA
catalog.status("dice_cd4_stim")
reference_data.provenance("dice_cd4_stim")
reference_data.verify("dice_cd4_stim")  # full SHA256 verification
cells.sources()                      # versions, citations, hashes, cached status
```

The default path is `~/.cache/oncoref/sources/<source>/<version>/<filename>`.
`CANCERDATA_DATA_DIR` overrides the reference cache root, as for existing HPA
data. Downloads use locks, staging and atomic replacement. The packaged
`normal-cell-sources.json` pins both archive and extracted-content SHA256 plus
byte size. Source version is not inferred from a changing website. A corrupted
pinned download cannot replace a valid cached file.

```bash
oncoref data fetch dice_cd4_stim
oncoref data fetch cornea_collin_adult_counts
oncoref data status reference
oncoref data dir reference
```

All 44 new artifacts were downloaded for this audit. `catalog.fetch("reference")`
fetches the external DICE/corneal studies; HPA remains the `hpa` catalog kind.
`catalog.fetch("all")` covers both and the existing expression bundle. Fetching
large references is explicit, never an import-time side effect.

## Native tables and uniform observations

```python
cells.table("hpa_immune_cell", genes=["IL13"])
cells.table("hpa_antibody_xml", genes=["FOXR2"])
cells.observations("hpa_rna_consensus", genes=["CT55"])
cells.observations("hpa_single_cell_tissue", genes=["IL13"])
cells.observations("dice_cd4_stim", genes=["IL13", "ENSG00000169551"])
cells.observations("cornea_collin_adult_counts", genes=["ABCB5"])
cells.cornea_summary(genes=["ENSG00000004846"])
cells.publication_evidence()

# Native developmental/culture matrices stay separate from adult-normal evidence.
for frame in cells.iter_tables("cornea_collin_10pcw_counts", genes=["ABCB5"]):
    print(frame)
```

`publication_evidence()` preserves the original eight findings. For current
accepted, unresolved and rejected literature findings, use
`from oncoref import cta_policy; cta_policy.findings()`.

`table()` returns the source's native column names and units. `iter_tables()`
streams large native matrices. `observations()` has shared columns for source,
version, taxon, gene, tissue, cell type, state, developmental stage, donor/sample,
modality, unit, value, staining level, reliability and antibody. It keeps TPM,
nTPM and CP10K separate. A missing gene yields no measured observation; zeros
are values actually present in the source. Metadata tables use `table()`.

HPA convenience functions remain available alongside existing methods:

```python
from oncoref import human_protein_atlas as hpa
hpa.primary_ihc(genes=["FOXR2"])
hpa.immune_cells(study="monaco", genes=["IL13"])
hpa.single_cell_tissues(genes=["IL13"])
```

Corneal gene labels remain in `gene_name`; `gene_id`, `canonical_symbol` and
`gene_mapping` identify mapped or unresolved canonical identities. Ambiguous
many-to-one symbols remain unresolved rather than being summed after
normalization. No paralog substitutes for an absent gene. Native tables preserve
all original labels. HPA/DICE use their supplied Ensembl identities, with DICE
version suffixes removed in observations; source tables retain the suffixes.

## Assay and identity safeguards

- **Cornea:** adult author annotations must exactly match count-matrix barcodes.
  Counts are normalized by each cell's full library before averaging, including
  zero cells. Barcode suffixes are sample groups, not verified donor IDs.
  Developmental, diseased and cultured samples remain separate native matrices;
  they do not enter the adult-normal gate.
- **DICE:** numeric columns describe each gene's expression distribution. They
  are not stable donor IDs across genes/cell types. `cell_type` retains the study
  code and `cell_label` its author-provided name (for example, M2 means
  non-classical monocytes, not polarized M2 macrophages). HPA's Schmiedel/DICE summary
  is the same study, not another independent replication.
- **Antibodies:** primary staining, antibody validation, curated profile
  reliability, patient IDs, image links and notes are retained. Patient IDs do
  not turn an aggregate tissue/cell staining row into paired patient-level data.
  HPA primary observations can exist when curated `normal_tissue.tsv` has no row.
  Inherited result columns `IHC measurements` and `IHC coverage` describe curated
  IHC only; `Primary IHC measured` and the primary-antibody export are separate.
- **Anatomy:** immune cells sampled in testis/fallopian tube do not inherit the
  organ's reproductive permission. Cornea, limbus and conjunctiva are explicit
  eye compartments. All new measured mappings are human; the rules reject
  mouse/dog input instead of transferring human thresholds.

Large XML and corneal matrices get derived Parquet caches beside their raw
source. Cache keys include adapter version and input checksums. These caches
are atomic and bounded-memory builds; raw data remains reusable by other
analyses. Typical cached queries do not rewrite the source files.

## Decisions

Admission and normal counterevidence are separate. Reproductive RNA ≥5 nTPM,
reliable reproductive protein, or cited reproductive plus tumor/recognition
evidence can satisfy admission. Normal exclusions still take precedence.
CT55 passes the added-support rule; FOXR2 has an unresolved primary-antibody
conflict. Published nonmalignant protein production excludes ABCB5, IL5, IL13
and PAPPA from reproductive-restricted targeting.

New RNA/antibody signals produce **Neither**, not a clinical-toxicity claim.
Thresholds are source-specific review parameters: HPA clusters ≥5 nTPM/≥20
reads and normal ≥ reproductive-cell maximum; original DICE ≥5 TPM in ≥3
reported values; HPA/Monaco immune summaries ≥5 in native units; adult cornea
mean CP10K ≥1 with ≥5 positive cells and ≥20 counts in ≥2 sample groups.
Uncertain normal primary staining also prompts review. The [full report](audits/cta-normal-cells-20261023/index.md)
documents precedence, borderline cases, source dependence and limitations.

The pure functions in `normal_cell_screen` accept normalized observations and
return evidence flags or classifications. `classify()` requires the frozen
v1.16 baseline, a unique gene-level evidence summary and cited publication
facts. Its output records previous decisions, individual findings, current
policy and change. The [24 October review](audits/cta-neither-review-20261024/index.md)
additionally checks every previous CTA/Neither gene against the hard >30 ceiling,
including individual immune specimens and unsited cell-type summaries. High values
now enter Unsafe even with reproductive-positive IHC or a Variable designation.
Reproductive/Disposable contexts remain permitted; CP10K remains a separate rule.
The 23 October findings above remain available as the frozen baseline.
`target_sets.results()` exposes the published current table;
only Unsafe supplies forbidden proteins.

```python
from oncoref import normal_rna_review
normal_rna_review.followup_reviews()  # questions, evidence, maxima and current decisions
normal_rna_review.sequence_gene_reviews()  # reference loci outside the canonical table
normal_rna_review.forbidden_gene_ids()     # complete gene-ID union used by protein exports
```

## Reproduce

```bash
python scripts/build_normal_cell_mappings.py
python scripts/normal_cell_cta_report.py --previous-sequences /path/to/pinned/sequence-assets --assets /path/to/protein-assets
python scripts/render_normal_cell_report.py --assets /path/to/protein-assets --publish-dataset
```

The previous sequence directory must contain the checksum-pinned Ensembl 115
assets listed in the 21 October sequence manifest. The script verifies all of
them, regenerates gene/sequence groups and writes the new policy into the recipe.
The exact protein priority is Unsafe > Neither > CTA, including unresolved
encoding genes. The source inventory, mapping audit, native evidence tables,
benchmark decisions and implementation hashes accompany the report.
