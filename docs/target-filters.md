# Filtering the Target candidates, loosely or strictly

**The catalog decides what is a candidate. These knobs only narrow that list, and
every knob reads a column of the catalog: the filters hold no evidence of their own.
Nothing here changes a category, an Unsafe call or the catalog policy
`oncoref.target-selection.2026-11-15.1`.**

[Web summary](audits/target-filters-20261112/report.html) ·
[report with every count, filter-mode funnel and familiar CTA](audits/target-filters-20261112/index.md) ·
[literature sweep](audits/target-filters-20261112/normal-function-sweep.md)

```python
from oncoref import target_filters as tf

tf.select("strict")                                   # a filter mode
tf.select(tf.TargetFilter(min_papers=2, min_tumor_fold=50))      # your own knobs
tf.funnel("strict")                                   # what each layer removes
tf.apply(tf.TargetFilter(min_papers=2))               # every protein, with why it failed
print(tf.describe("strict"))                          # the rule in words
tf.measures()                                         # every value a knob reads
```

## How it fits the catalog

The catalog runs one shared screen, then asks two separate questions, then names a
category. [Full rules](target-selection.md).

| Step | Catalog rule | Filter layer that tightens it |
|---|---|---|
| Shared screen | IHC-dependent RNA ceilings of 4 / 8 / 16 by tissue tier, immune and single-cell cluster rules, established rare-cell RNA, coverage, identity, sequence (read-through loci included), Disposable at most 32 | 1 `off_target` |
| Published normal expression and home tissue | Every candidate reviewed. A quoted, gene-specific report in a Vital or Protected tissue is Unsafe; weaker findings, family rules and the nose and skin reference record a home tissue. A home no panel resolves needs reproductive evidence (Other:Tissue) | 2 `curation` |
| Cancer evidence | One QC-passing tumor sample, including TGCT, at 10 × max(N, Disposable, 1) | 3 `sets` (`min_tumor_fold`) |
| Reproductive evidence | A reproductive or developmental RNA route at 1 or more, or reliable protein | none |
| Category | CTA = both. Cancer Antigen only. Reproductive Antigen only | 3 `sets` (`membership`) |
| Published support | Recorded, never required | 4 `papers` |
| | | 5 `exclusions` (your judgment, by name) |

## The rules

1. **A protein is kept only if it passes all five layers.** `filter_failed` names the
   first layer it fails.
2. **Every default is the catalog.** `TargetFilter()` returns the plain 668 candidates.
3. **`max_*` means below the value. `min_*` means at least the value.**
4. **A missing measurement never passes.** Unmeasured is not low.
5. **Knobs cannot add a protein.** Other and Unsafe proteins stay out, and
   `min_tumor_fold` cannot go below the catalog's 10.

## Every knob

| Layer | Knob | Default | Rule |
|---|---|---|---|
| 1 `off_target` | `max_off_target_tpm` | none | The larger of the catalog's normal reference N and the Disposable maximum, in native TPM, must be below this. N is the highest value across non-reproductive tissues, brain subregions and immune-cell summaries |
| 2 `curation` | `normal_function` | `any` | Drops proteins whose recorded normal home is at or above a tier: `no_vital`, `no_protected` (Vital and Protected) or `no_manageable`. Disposable tissues never remove a protein |
| 3 `sets` | `min_tumor_fold` | 10 | To stay a **Cancer Antigen**: highest tumor ÷ max(N, Disposable, 1). The catalog requires 10 |
| | `membership` | `target_candidates` | Which set to keep: `target_candidates` (either), `cancer_antigen`, `reproductive_antigen`, `cta` (both), `cancer_antigen_only`, `reproductive_antigen_only` |
| 4 `papers` | `min_papers` | 0 | Distinct papers, by DOI, that nominate the protein or report it in tumors |
| 5 `exclusions` | `exclude_prefixes` | none | Drop every protein whose label starts with one of these |
| | `except_labels` | none | Labels the prefix rule spares |

Earlier knobs that acted on flags the catalog does not exclude on (HPA nTPM, dismissed
rare-cell signals, literature conflicts, mapping risk), on patient counts or on paper
families are gone. The flags stay on the catalog rows for reading.

### Layer 2 acts on the recorded normal home

**The catalog and layer 2 read the same records.** The review either excludes a
protein (Unsafe in the catalog, so it never reaches a filter) or records a finding
with its tissue and tier. Layer 2 acts on the recorded ones.

| Record | Where it lives | Catalog effect | Layer 2 |
|---|---|---|---|
| Quoted, gene-specific report in a Vital or Protected tissue | `normal-expression-literature.json` | Unsafe | Not reached |
| Confirmed finding below the rubric | Same file, `recorded_tissue` and `recorded_tier` | Other:Tissue when the tissue is one no panel resolves and no reproductive evidence qualifies | Acts by tier |
| Family rule by protein name | Same file, `family_rules` | Other:Tissue when the family's home is a tissue no panel resolves and no reproductive evidence qualifies | Acts by tier when `filter_acts` |
| Tissue tier | `normal-tissue-tiers.yaml` and its dated updates | Decides R2 | Decides the tier |
| Measured in normal nose or skin at medium or high grade | `nose-skin-reference.csv`, columns `reference_recorded_*` | Same, for nose tissues | Acts by tier |
| Read-through locus | `readthrough-loci.csv` | Other:Sequence | Not reached |

A protein's `curated_tier` is the most protective tier among its recorded finding, its
nose and skin reference record and an acting family rule; `curated_sources` keeps the
links. Bitter taste receptors are Protected by family rule and olfactory receptors
Manageable. Hair keratins and PRAME (skin melanocytes) are Disposable, so they are
recorded and never removed. The catalog rows carry the
columns `literature_recorded_tissue`, `literature_recorded_tier`, `literature_family`
and `literature_family_tier`.

```python
from oncoref import target_curation as tc

tc.entries("LRRC52")      # verdict, effect, recorded tissue and tier
tc.findings("LRRC52")     # every evidence item: source, quote, second reader's note
tc.families()             # family rules and their sources
tc.compartments()         # tissues added to the tissue-tier registry, with references
tc.readthrough()          # read-through loci and parent genes
tc.nose_skin("GLYATL3")   # olfactory bulk RNA and nose, skin, scalp, tongue single cells
tc.rna_grade(12.3)        # "medium": an RNA value read on the staining calibration
```

`target_curation` owns no data; it reads the catalog's stores. To change a record, edit
the literature file or the tissue-tier registry and rebuild the catalog
([steps](audits/target-selection-20261110/rebuild/index.md#evidence-and-reproduction)).

### Normal nose and skin, which the expression panels do not sample

**Four candidates are recorded in normal nose or skin: GLYATL3 and PABPC3 in olfactory
mucosa, FKBP6 in nasal ciliated cells and EFCAB3 in skin cells.** Nasal mucosa and
skin are Protected, so `strict` removes FKBP6 and EFCAB3; olfactory epithelium is
Manageable. All four have reproductive evidence, so none changes category under the
[home-tissue rule](target-selection.md#home-tissue-the-panels-do-not-resolve).

The catalog's panels hold no nasal tissue, and bulk skin dilutes single cell types. Two
references fill that gap for the 868 proteins that were candidates before the
home-tissue rule (`oncoref/data/nose-skin-reference.csv`):

| Reference | What it covers | Recorded when |
|---|---|---|
| Olfactory mucosa bulk RNA, four people, plus one nasal respiratory sample ([Olender 2016](https://doi.org/10.1186/s12864-016-2960-3), Table S2) | 643 of the 868 match a gene symbol | 1 FPKM or more in three of four people, and the median grades medium or high |
| Adult normal single cells in the CELLxGENE Census: nose (110,033 cells, 26 donors, 2 usable datasets), skin (188,971 cells, 95 donors, 7 datasets), scalp and tongue (1 dataset each) | All 868 | 5% of a cell type positive, 1% with two or more molecules, in two datasets, at medium or high grade |

Ten more candidates are in olfactory mucosa at low grade (ZC3H11B, PGAM4, OR10G3, RTP2,
CDCP2, FKBP6, TRIM43, RGPD3, OR11A1, OR8B3). They carry the value and no recorded tissue.

**Grades come from antibody staining.** HPA tissue RNA was paired with the strongest
staining level of the same gene in the same tissue, using only antibodies of Enhanced
reliability: 153,739 pairs, 3,951 genes, 44 tissues
(`oncoref/data/rna-ihc-calibration.json`).

| Tissue RNA (nTPM) | Pairs | Stained | Medium or high staining |
|---|---:|---:|---:|
| 0 to 0.5 | 40,183 | 6% | 2% |
| 0.5 to 1 | 9,624 | 15% | 7% |
| 1 to 2 | 11,215 | 24% | 12% |
| 2 to 5 | 18,448 | 38% | 23% |
| 5 to 10 | 17,138 | 56% | 39% |
| 10 to 20 | 18,414 | 67% | 50% |
| 20 to 50 | 19,776 | 75% | 59% |
| 50 to 100 | 9,084 | 80% | 66% |
| 100 or more | 9,857 | 81% | 70% |

Low starts at 1, medium at 5 and high at 20. Staining levels overlap widely in RNA, so
the grade says how likely the protein is to be seen, not how much there is. FPKM and
single-cell counts per million are read on the same scale as an approximation.

Limits:

- Olfactory biopsies mix olfactory and respiratory epithelium, and 225 candidates have
  no matching symbol in the 2016 table. Those are unmeasured, not absent.
- PABPC3 is a retrocopy of PABPC1. Its olfactory value may be misassigned reads.
- The nasal single cells are respiratory brushings. They hold no olfactory neurons, so
  olfactory receptors are measured only by the bulk table.
- Scalp and tongue have one dataset each, so nothing is recorded from them. The scalp
  cells are follicular keratinocytes, not hair-shaft cortex, where hair keratins and
  keratin-associated proteins are made.
- One skin dataset shows 86 candidates in 5% or more of some cell type, at similar rates
  across unrelated cell types and in no other dataset. Requiring two datasets keeps it
  from recording a tissue. EFCAB3 is the one that a second dataset (single nuclei)
  supports.
- Melanocytes are left out.

Rebuild with `scripts/extract_nose_skin_cell_counts.py`,
`scripts/build_rna_ihc_calibration.py` and `scripts/build_nose_skin_reference.py`.

### Layer 3 tightens Cancer Antigens, then picks a set

The fold knob does not delete proteins directly. It decides whether a protein is still
a Cancer Antigen, and `filtered_category` records the result.

| Protein | `min_tumor_fold=50` | Result |
|---|---|---|
| CTA with its best tumor at 30x | Loses Cancer Antigen, keeps Reproductive Antigen | Kept as Reproductive Antigen only |
| Cancer Antigen only at 30x | Loses its only membership | Removed |
| Reproductive Antigen only | Not affected by a tumor knob | Kept |
| Any of these with `membership="cta"` | Must hold both after tightening | Only CTAs at 50x or more are kept |

### What counts as a paper

| Family | Counts for a protein when | Sources |
|---|---|---|
| `cta` | The paper **nominates** the protein in a cancer/testis, cancer-germline or cancer-placenta list. Appearing only in a background table does not count | da Silva 2017, Wang 2016, Chang 2019, Loriot 2025, Carter 2023, Jamin 2021, Bruggeman 2018, Seager 2024 |
| `oncofetal` | The protein is in an embryonic gene program that the paper finds re-expressed in cancer | Smith 2023, Chew 2019 |
| `tumor` | The paper reports the protein's tumor expression | Single-gene reports |
| `reproductive` | The paper reports normal reproductive enrichment. Not counted: it makes no tumor claim | Gong 2021, Hirbo 2015 and others |

A DOI counts once however many of its tables list the gene. Overlapping lists are not
independent experiments. [Every source with its recorded scope](audits/target-filters-20261112/paper-families.csv).

## Filter modes

`all` → `loose` → `standard` → `strict`. Each mode keeps a subset of the one before,
and a test enforces it. Modes are starting points, not policy.
[Counts and funnels](audits/target-filters-20261112/index.md#4-filter-modes).

| Filter mode | Off target | Curation | Sets | Papers | Exclusions |
|---|---|---|---|---|---|
| `all` | | | | | |
| `loose` | | | | ≥1 | |
| `standard` | TPM < 5 | No Vital home | | ≥2 | |
| `strict` | TPM < 4 | No Vital or Protected home | CTA only | ≥3 | MAGE except MAGEA4 |

| Filter mode | Kept | CTA | Cancer Antigen only | Reproductive Antigen only |
|---|---:|---:|---:|---:|
| `all` | 668 | 437 | 34 | 197 |
| `loose` | 453 | 333 | 6 | 114 |
| `standard` | 296 | 235 | 1 | 60 |
| `strict` | 141 | 141 | 0 | 0 |

**Every mode must keep PRAME, NY-ESO-1 and MAGEA4** (`POSITIVE_CONTROLS`); a test fails
otherwise. No mode raises the fold above the catalog's 10: doing so relabels CTAs as
Reproductive Antigen only rather than removing anything, so it is left as a knob.

## The MAGE rule

`strict` drops every protein whose label starts with `MAGE` except
MAGEA4. **This is a judgment by family name, not a measurement.** Eleven of the 18 it drops
would pass `strict` on their own data. [Per-protein table](audits/target-filters-20261112/index.md#6-the-mage-family).

The reasons sit outside this repo and are not verified against sources here:
MAGE-A3/A12 T-cell receptor therapies are reported to have caused severe toxicity
through normal-tissue expression, and the near-identical MAGE-A paralogs are hard to
tell apart by RNA. `exclude_prefixes=()` turns the rule off. The catalog rows carry
`paralog_ambiguous` and `host_overlapped` for a data-driven reading.

### MAGEA3 and MAGEB2: what drives their Unsafe call

Both are Unsafe on "HPA tissue clusters" alone, but not for bulk tissue. The cause is
one cluster: **colon "undifferentiated cells" (cluster c-8)**.

| Gene | Colon c-8 nTPM (reads) | Highest bulk |
|---|---|---|
| MAGEA3 | 53.1 (251) | 0.2 TPM |
| MAGEB2 | 118.4 (560) | 0.1 TPM |
| MAGEA6 | 0.0 (0) | 0.2 TPM |

I tested whether c-8 is an artifact cluster and **it is not obviously one**: about 25% of
genes reach 20 nTPM with at least 20 reads in every colon, bronchus, eye and adrenal
cluster I checked, c-8 included.

What the data does show: MAGEA3 and MAGEA6 are near-identical, behave alike in testis
(62 and 54 nTPM), and differ completely in colon c-8 (53 against 0 reads). That is what
either a real MAGEA3-specific signal or a read-assignment artifact would look like, and
counts alone cannot separate them. MAGEA3 and MAGEB2 stay Unsafe under the catalog's
rule; nothing here overrides it. Settling it needs the reads, or an assay specific to
MAGEA3 in colon.

These results describe screening rules. They are not validated tumor antigens.

Regenerate the report with `PYTHONPATH=. python scripts/build_target_filter_audit.py`.
