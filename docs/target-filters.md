# Filtering the Target candidates, loosely or strictly

**The catalog decides what is a candidate. These knobs only narrow that list.
Nothing here changes a category, an Unsafe call or the catalog policy
`oncoref.target-selection.2026-11-10.1`.**

[Web summary](audits/target-filters-20261112/report.html) ·
[report with every count, filter-mode funnel and familiar CTA](audits/target-filters-20261112/index.md) ·
[literature sweep](audits/target-filters-20261112/normal-function-sweep.md)

```python
from oncoref import target_filters as tf

tf.select("strict")                                   # a filter mode
tf.select(tf.TargetFilter(min_papers=2, max_off_target_tpm=1.0))   # your own knobs
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
| Shared screen | IHC-dependent RNA ceilings of 4 / 8 / 16 by tissue tier, immune and single-cell cluster rules, established rare-cell RNA, published normal expression, coverage, identity, sequence, Disposable at most 32 | 1 `off_target`, 2 `warnings` |
| Cited literature | Reviewed for about 99 proteins | 3 `curation` |
| Cancer evidence | Two independent tumors above the normal reference N, or one tumor at 10 × max(N, 1) | 4 `sets` |
| Reproductive evidence | A reproductive or developmental RNA route at 1 or more, or reliable protein | 4 `sets` |
| Category | CTA = both. Cancer Antigen only. Reproductive Antigen only | 4 `sets` (`membership`) |
| Published support | Recorded, never required | 5 `papers` |
| | | 6 `exclusions` (your judgment, by name) |

## The rules

1. **A protein is kept only if it passes all six layers.** `filter_failed` names the
   first layer it fails.
2. **Every default is off.** `TargetFilter()` returns the plain 902 candidates.
3. **`max_*` means below the value. `min_*` means at least the value.**
4. **A missing measurement never passes.** Unmeasured is not low.
5. **TPM and nTPM are never converted.** Each has its own knob.
6. **Knobs cannot add a protein.** Other and Unsafe proteins stay out.

## Every knob

| Layer | Knob | Default | Rule |
|---|---|---|---|
| 1 `off_target` | `max_off_target_tpm` | none | **Tissue level.** The larger of the catalog's normal reference N and the Disposable maximum, in native TPM, must be below this. N is the highest value across non-reproductive tissues, brain subregions and immune-cell summaries |
| | `max_off_target_ntpm` | none | **HPA nTPM, single-cell clusters included.** The larger of the HPA healthy maximum and the Disposable maximum must be below this. The healthy maximum takes in single-cell clusters with at least 20 reads, so one cluster can set it |
| | `max_disposable` | none | Appendix, gallbladder and tonsil alone must be below this on both scales. The catalog already caps them at 32 |
| 2 `warnings` | `rare_cell` | `any` | `no_dismissed` drops candidates that carry a dismissed rare-cell signal |
| | `literature` | `any` | `no_conflict` drops candidates whose literature review is `conflicting_or_nonspecific` or `specialized_tissue_expression` |
| | `mapping` | `any` | `no_risk` drops candidates flagged `paralog_ambiguous` or `host_overlapped`: reads could belong to another locus |
| 3 `curation` | `normal_function` | `any` | Drops proteins with a cited normal function or expression in a tissue at or above a tier: `no_vital`, `no_protected` (Vital and Protected) or `no_manageable`. Disposable tissues never remove a protein |
| | `readthrough` | `any` | `exclude` drops loci that share exons with two or more other protein-coding genes |
| 4 `sets` | `min_patients` | 0 | To stay a **Cancer Antigen**: independent tumors above N |
| | `min_tumor_fold` | 0 | To stay a **Cancer Antigen**: highest tumor ÷ max(N, 1) |
| | `min_reproductive` | 0 | To stay a **Reproductive Antigen**: the highest of adult reproductive nTPM, reproductive cell nTPM and developmental mean TPM, each in its own unit |
| | `membership` | `target_candidates` | Which set to keep after the three knobs above: `target_candidates` (either), `cancer_antigen`, `reproductive_antigen`, `cta` (both), `cancer_antigen_only`, `reproductive_antigen_only` |
| 5 `papers` | `min_papers` | 0 | Distinct papers, by DOI, in the counted families |
| | `families` | `cta`, `oncofetal`, `tumor` | Which families count (below) |
| | `min_families` | 0 | The counted papers must span this many families |
| | `exact_members_only` | off | Drops the six proteins whose nomination names only some of their member genes |
| 6 `exclusions` | `exclude_prefixes` | none | Drop every protein whose label starts with one of these |
| | `except_labels` | none | Labels the prefix rule spares |

### Layer 3 acts on cited literature

`oncoref.target_curation` holds the records. Nothing in it changes the catalog.

| Record | What it holds | Acts when |
|---|---|---|
| Finding | One protein, one normal tissue, the claim, its source and a second reader's verdict | The verdict is confirmed: human, adult, normal tissue, gene-specific assay |
| Family rule | A name pattern, a tissue and sources | Marked `acts`: olfactory receptors, bitter taste receptors, hair keratins |
| Compartment | A tissue or cell type and its tier | Tier from `tissue_hazards`, or a cited proposal for the seven tissues it had no call for |
| Read-through locus | The locus and the genes whose exons it shares | `readthrough="exclude"` |

A protein's `curated_tier` is the most protective tier among its acting records.
`curated_sources` keeps the citations. Proposed tiers: olfactory epithelium, taste buds,
nasal mucosa and lacrimal gland are Protected; sensory ganglia and bone are Vital; dental
pulp is Manageable. Hair follicle and cutaneous melanocytes are Disposable in the
registry, so hair keratins and PRAME are recorded and never removed.

```python
from oncoref import target_curation as tc

tc.findings("LRRC52")     # claim, source, verdict, tier
tc.compartments()         # every tissue, its tier and where the tier comes from
tc.families()             # family rules and their sources
tc.readthrough()          # read-through loci and parent genes
```

Edit `oncoref/data/target-curation.yaml`, then run
`PYTHONPATH=. python scripts/build_target_curation.py`.

### Layer 4 tightens each set, then picks one

The evidence knobs do not delete proteins directly. They decide which set a protein
still belongs to, and `filtered_category` records the result.

| Protein | `min_patients=5` | Result |
|---|---|---|
| CTA with 3 tumors above N | Loses Cancer Antigen, keeps Reproductive Antigen | Kept as Reproductive Antigen only |
| Cancer Antigen only with 3 tumors | Loses its only membership | Removed |
| Reproductive Antigen only | Not affected by a tumor knob | Kept |
| Any of these with `membership="cta"` | Must hold both after tightening | Only CTAs with 5 or more tumors are kept |

### What counts as a paper

| Family | Counts for a protein when | Sources |
|---|---|---|
| `cta` | The paper **nominates** the protein in a cancer/testis, cancer-germline or cancer-placenta list. Appearing only in a background table does not count | da Silva 2017, Wang 2016, Chang 2019, Loriot 2025, Carter 2023, Jamin 2021, Bruggeman 2018, Seager 2024 |
| `oncofetal` | The protein is in an embryonic gene program that the paper finds re-expressed in cancer | Smith 2023, Chew 2019 |
| `tumor` | The paper reports the protein's tumor expression | Single-gene reports |
| `reproductive` | The paper reports normal reproductive enrichment. Not counted by default: it makes no tumor claim | Gong 2021, Hirbo 2015 and others |

A DOI counts once however many of its tables list the gene. Overlapping lists are not
independent experiments. [Every source with its recorded scope](audits/target-filters-20261112/paper-families.csv).

## Filter modes

`all` → `loose` → `standard` → `strict`. Each mode keeps a subset of the one before,
and a test enforces it. Modes are starting points, not policy.
[Counts and funnels](audits/target-filters-20261112/index.md#4-filter-modes).

| Filter mode | Off target | Curation | Sets | Papers | Exclusions |
|---|---|---|---|---|---|
| `all` | | | | | |
| `loose` | | No read-throughs | | ≥1 | |
| `standard` | TPM < 5 | No read-throughs; no Vital-tissue function | Cancer Antigen needs ≥2 tumors | ≥2 | |
| `strict` | TPM < 4 | No read-throughs; no Vital or Protected-tissue function | Cancer Antigen needs ≥5 tumors and ≥10×; CTA only | ≥3, exact members | MAGE except MAGEA4 |

**Every mode must keep PRAME, NY-ESO-1 and MAGEA4** (`POSITIVE_CONTROLS`); a test fails
otherwise. That is why no mode uses `max_off_target_ntpm`, `rare_cell` or `literature`:
each of those removes PRAME. They remain knobs you can set yourself. For "below 1
everywhere outside tumors and reproductive tissue", use
`TargetFilter(max_off_target_tpm=1.0, max_off_target_ntpm=1.0)`.

## The MAGE rule

`strict` drops every protein whose label starts with `MAGE` except
MAGEA4. **This is a judgment by family name, not a measurement.** Ten of the 17 it drops
would pass `strict` on their own data. [Per-protein table](audits/target-filters-20261112/index.md#6-the-mage-family).

The reasons sit outside this repo and are not verified against sources here:
MAGE-A3/A12 T-cell receptor therapies are reported to have caused severe toxicity
through normal-tissue expression, and the near-identical MAGE-A paralogs are hard to
tell apart by RNA. `exclude_prefixes=()` turns the rule off. `mapping="no_risk"` is the
data-driven alternative, but it removes half of all candidates.

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
