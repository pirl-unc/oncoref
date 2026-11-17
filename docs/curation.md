# Antigen lookup and curation builds

The end-user commands read the current catalog. The maintainer command rebuilds
the catalog and its reports from the declared evidence inputs.

```sh
oncoref antigens list --set CTA --mode standard
oncoref antigens list --set CA --min-tumor-fold 30 --min-patients 2 --format csv
oncoref antigens explain CTAG2
oncoref antigens explain CTAG2 --loss --format json
oncoref antigens funnel --mode strict
oncoref antigens subtract protein.fa --background 'Unsafe|Other' --format csv
oncoref antigens subtract protein.fa --equiv ILV --max-distance 1
```

`--set` accepts `CTA`, `CA`, `RA`, `CA|RA`, `Unsafe`, `Other`, and `Unsafe|Other`.
Named modes and fold/patient constraints tighten the positive sets. Backgrounds
are exported without candidate filters. `--ids` exports member gene IDs.
Tables, CSV and JSON are supported. The old `oncoref cta` command retains its
compatibility snapshot and emits a notice pointing to the current commands.

Subtraction requires locally installed Ensembl 115 translations and stops if any
background gene lacks a translation. Defaults are exact 8–15mers against
`Unsafe|Other`. Use `--min-length` and `--max-length` to select lengths. Equivalence
groups are explicit, disjoint assumptions; they are not a universal recognition model.
Weighted edits use `--substitution-cost I:L:0.5`, `--insertion-cost`, and
`--deletion-cost`. Substitution costs are directed. JSON provenance is written to
stderr; retained peptide occurrences go to stdout. Sequence retention does not
demonstrate peptide presentation or therapeutic safety.

## One maintainer entry point

From a complete source checkout:

```sh
python -m oncoref.curation build --reuse-tumor-cache
python -m oncoref.curation check
# Equivalent verification-only spelling:
python -m oncoref.curation build --check
python -m oncoref.curation diff --against prior-categories.csv.gz
```

The declared dependency graph runs rare-cell grading, the catalog, plots, the rules
report, and the filter audit. Failures stop dependent steps. Omitting
`--reuse-tumor-cache` re-extracts tumors; supplying it uses the existing builder's
checksum and identity checks. `check` verifies implementation and output hashes
without downloading or reprocessing raw data. `diff` uses stable protein identities
and highlights named regression controls; its input is a data-only snapshot.

The graph and the individual HPA/Petropoulos/Tyser input bindings live in
`oncoref/data/curation-build.json`. Their provider URLs, versions, sizes and hashes
live in `reference_data`. `curation.reference_input(role)` resolves verified bytes
under a SHA-256 directory. An existing download can seed that cache with
`reference_input(role, import_path=path)`; it must match the registered checksum.

## Migration boundaries

The current build still consumes the pinned October normal-screen fact table and
its later adjudications, plus older identity and immune-context artifacts. The new
runner makes that dependency order explicit; it does not claim those inputs have
already been replaced by direct raw-assay evaluation. Removing that baseline needs
layer-by-layer parity and an explanation for every changed normal exclusion.

Historical receipts and public compatibility modules remain available. New category
comparisons use compact data-only snapshots. Existing historical fixtures are not
deleted or rewritten as part of a live policy update.

Unresolved locus attribution remains a research review, not a code failure. A
same-strand host or paralog provides a hypothesis; discriminating reads, junctions
or an independent gene-specific assay can adjudicate it. Negative bulk RNA or IHC
alone cannot establish that a sparse normal-cell finding is an artifact.

For toxicity consequences, use `tissue_tiers.consequences(tissue, cell_type,
resolution="cell")`. The result retains organ fallbacks, acute and persistent loss,
repair, management, evidence scope and linked sources. Reproductive membership does
not overwrite those assessments.
