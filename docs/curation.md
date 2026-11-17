# Antigen lookup and curation builds

The end-user commands read the current catalog (`oncoref.target_selection`). The
maintainer command rebuilds the catalog, figures, and reports from the declared
evidence inputs in topological order.

```sh
# Default largest CTA set (all 290 CTAs across all 12 compartments)
oncoref antigens list --set CTA

# Host anatomy & tissue-tolerance profiles
oncoref antigens list --set CTA --profile male-patient
oncoref antigens list --set CTA --profile male-patient --absent-tissue prostate --absent-tissue "seminal vesicle"
oncoref antigens list --set CTA --profile female-patient
oncoref antigens list --set CTA --profile female-patient --lactating
oncoref antigens list --set CTA --profile female-patient --protect-tissue ovary

# Combine with filter modes, fold cutoffs, or patient replication
oncoref antigens list --set CTA --profile female-patient --mode standard
oncoref antigens list --set CA --min-tumor-fold 30 --min-patients 2 --format csv

# Inspect every gate, mode, and host profile for one gene
oncoref antigens explain LALBA --loss
oncoref antigens explain CTAG2 --format json

# Sequential catalog + filter attrition funnel
oncoref antigens funnel --profile male-patient --mode standard

# Exact or bounded-distance peptide subtraction against Unsafe|Other
oncoref antigens subtract protein.fa --background 'Unsafe|Other' --format csv
oncoref antigens subtract protein.fa --equiv ILV --max-distance 1
```

`--set` accepts `CTA`, `CA`, `RA`, `CA|RA`, `Unsafe`, `Other`, and `Unsafe|Other`.
`--profile` accepts `all` (no host screen), `female-patient` and `male-patient`
(assumed anatomy, with present gonads protected). `--accept-loss` (legacy alias
`--allow-tissue`) explicitly accepts tissue loss; `--absent-tissue` records absent
anatomy/state; `--protect-tissue` overrides an exception. `breast` expands to both
mammary states. `--pregnant` also protects the milk-secretory state, and pregnancy
remains unresolved because developmental coverage is incomplete. The old `core`
host shortcut is rejected. [Concepts and exact rules](host-selection.md).

`--ids` exports gene IDs. Background sets remain global when listed directly;
`subtract --profile ...` adds candidates excluded or unresolved under that host
context to the chosen background. List, explain and subtraction emit selection
provenance to stderr. `explain --profile ...` reports all host reasons. The legacy
`oncoref cta` command retains its historical snapshot and migration notice.

Subtraction requires locally installed Ensembl 115 translations and stops if any
background gene lacks a translation. Defaults are exact 8–15mers against
`Unsafe|Other`. Equivalence groups are explicit, disjoint assumptions; they are not a
universal recognition model. Weighted edits use `--substitution-cost I:L:0.5`,
`--insertion-cost`, and `--deletion-cost`. JSON provenance is written to stderr;
retained peptide occurrences go to stdout. Sequence retention does not demonstrate
peptide presentation or therapeutic safety.

## One maintainer entry point

From a complete source checkout:

```sh
python -m oncoref.curation build --reuse-tumor-cache
python -m oncoref.curation check
# Equivalent verification-only spelling:
python -m oncoref.curation build --check
python -m oncoref.curation diff --against prior-categories.csv.gz
```

The declared dependency graph (`oncoref/data/curation-build.json`) runs seven steps in
topological order: `rare_cell` → `catalog` → `host_exposure` → `plots` → `report` → `filter_audit` →
`filter_web_report`. Failures stop dependent steps. Omitting `--reuse-tumor-cache`
re-extracts tumors; supplying it verifies the tumor cache receipt and input hashes first.
`check` verifies implementation and output hashes across all 210+ tracked files without
downloading or reprocessing raw data. `diff` uses stable protein identities and highlights
named regression controls (`PRAME`, `NY-ESO-1`, `MAGEA4`).

## Migration boundaries

The current build consumes the pinned normal-screen fact table and its replayed
adjudications, plus checksum-verified HPA v23 individual tissue biopsies and
Petropoulos/Tyser human embryo single-cell datasets via `curation.reference_input(role)`.
Historical receipts and public compatibility modules (`cta_selection`, `target_sets`,
`cta`) remain frozen snapshots and are never overwritten by live policy updates.

Unresolved locus attribution remains an explicit research hold (`Other:Review`), not a
code failure. For reviewed tissue-loss consequences alongside observed reproductive RNA,
use `antigen_sets.loss_consequences(gene)` or `tissue_tiers.consequences(tissue,
cell_type, resolution="cell")`.
