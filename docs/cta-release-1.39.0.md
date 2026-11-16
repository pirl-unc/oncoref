# OncoRef 1.39.0: antigen sets, library QC and peptide subtraction

**The current catalog contains 437 CTAs and 668 candidate antigen proteins.**
CTA remains Cancer Antigens ∩ Reproductive Antigens; cancer-only proteins retain
their own category. The union is available without redefining CTA.

[Curation rules](target-selection.md) · [Funnel and Venn figures](audits/target-selection-20261110/rebuild/summary.md) ·
[Searchable catalog](audits/target-selection-20261110/rebuild/report.html) ·
[Set and peptide APIs](antigen-sequences.md)

| Category | 1.38.0 | 1.39.0 |
|---|---:|---:|
| CTA | 453 | 437 |
| Cancer Antigen only | 37 | 34 |
| Reproductive Antigen only | 180 | 197 |
| Other | 970 | 972 |
| Unsafe | 18,253 | 18,253 |

[Every classification change from 1.38.0](audits/target-selection-20261110/rebuild/changes-from-1.38.csv).

There are 471 Cancer Antigens, 634 Reproductive Antigens and 668 members of their
union. Counts refer to groups of identical selected reference proteins. Membership
does not certify every isoform or peptide as a tumor antigen.

## Curation and library evidence

The rule remains `T ≥ 10 × max(N, D, 1)` in at least one passing tumor sample.
N and D are supported normal-context summaries, not the maximum across every
normal donor. TPM comparisons are a screening convention across studies, not
batch-corrected differential expression. Independent normal RNA, protein,
rare-cell, identity, coverage and sequence gates still apply.

The library screen removes one Hodgkin library with broad olfactory-receptor
background and seven TGCT libraries with a coordinated late-spermatid program.
The latter may reflect normal testis admixture; RNA alone does not establish
the producing cell type. A further rule for four unrelated tissue programs
adds no exclusions. The measurements and reasons remain in the
[QC table](audits/target-selection-20261110/rebuild/library-qc.csv.gz).

Passing TGCT libraries now count. DAZL, PRDM14, ZNF705D, TMEM202 and GOLGA6D become
CTA. MAGEC2 also becomes CTA (681.6×), after a source-specific review using
HPA062230 and GTEx resolves the historical TC-1 hold as unsupported. Its MAGE
family filter remains in force. Removing the Hodgkin library removes cancer
support from 22 prior CTAs and three cancer-only proteins; all 22 remain
reproductive antigens.

DAZ1/DAZ4 remain Unsafe on normal gastric RNA and DAZ2/DAZ3 remain Other:Sequence.
Their shared findings preserve family-level RNA support, uncertain individual
paralog assignment and the lack of protein confirmation. Reviewed findings now
support explicit supersession: historical observations remain visible, while
only active findings impose holds. The web catalog exposes those records for
every gene through the same mechanism.

## APIs

- `antigen_sets.proteins()` and `gene_ids()` accept `CTA`, `CA`, `RA`, `CA|RA` and
  `Unsafe`. CA and RA include CTA. Historical APIs retain their pinned meanings.
- `antigen_sets.sequences()` returns all locally available Ensembl 115 protein
  translations, their gene/transcript occurrences, transcript biotypes, reference
  identity and gaps, including loss-of-function and decay-associated transcripts.
- `PeptideIndex.matches()` maps a query to every occurrence in either set.
  `subtract()` returns surviving 8–15mer occurrences with coordinates.
- Exact Unsafe subtraction is the default. Optional Levenshtein distance,
  directed weighted edit costs and disjoint reduced alphabets support sensitivity
  screens, including I/L and I/L/V. Stops, ambiguous residues and record boundaries
  are never crossed. Similarity is not a validated TCR cross-reactivity prediction.

## Versions and reproduction

Package **1.39.0**, target policy `oncoref.target-selection.2026-11-16.1`, expression
bundle **5.23.28**, unchanged raw matrix version **5.22.14**. Bundle 5.23.28 rebuilds
HL and TGCT summaries, percentiles, representatives and within-sample artifacts.
Per-cohort metadata distinguishes these v4 QC rebuilds from unchanged v3 artifacts;
the bundle-level metadata reports mixed policy versions. The old 1.36 sequence
bundle remains a separate historical snapshot.

```sh
python scripts/build_tissue_marker_panels.py
python scripts/rebuild_target_catalog.py
python scripts/build_target_selection_report.py
python scripts/render_target_catalog.py
python scripts/build_target_filter_audit.py
python scripts/build_target_filter_web_report.py
```

The full rebuild verifies pinned source hashes. `--reuse-tumor-cache` can reuse a
verified tumor matrix after changes outside extraction, checking member identities,
raw matrices, sample metadata, mapping code and library-QC inputs. The existing
catalog limitations remain explicit, including uncalibrated cross-study maxima,
unmeasured normal contexts and the separate hair-follicle replay tracked in #687.
