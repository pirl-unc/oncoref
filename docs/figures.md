# Figures

The plot CLI generates reference figures from packaged data. Analysis scripts
also generate reports from pinned cached source data. Dated documentation audits
include their figures, tables and provenance manifests; ordinary `figures/` and
`outputs/` runs remain local. Rebuild with the recorded inputs to update a figure.

## Where figures go

```
figures/
  run_20260915-104233/
    cta-curation/            # multi-figure families get a directory
    expression-provenance/
    cta-covering-set.png     # single figures get a file
  latest -> run_20260915-104233
```

Omitting `--out` opens a fresh `figures/run_<YYYYMMDD-HHMMSS>/` directory and
repoints `figures/latest` at it. A run never overwrites an older one, so a figure
in a draft can always be traced back to the run that made it. Passing `--out`
explicitly overrides this and writes exactly where you say.

`figures/` is gitignored. Set `ONCOREF_FIGURES_DIR` to relocate the root.

## Presets

Every command takes `--preset`:

| | `print` (default) | `slide` |
|---|---|---|
| Type | 10pt base | 17pt base |
| Page | up to 20in | 16:9, up to 12 x 6.75in |
| Rows | all of them | top 12, disclosed on the axis |
| Scatter labels | every point | top 12 points; all points still drawn |
| Category labels | `non_hodgkin_lymphoma` | `non hodgkin lymphoma` |

```bash
oncoref plot burden-category-bars --preset slide
python scripts/regenerate_plots.py --preset slide
```

A slide is not a small manuscript page. A figure in a paper is read at 30cm by
one person who can stop and study it; a figure on a slide is read at 5m by a room
that gets it for thirty seconds. Bigger type alone does not bridge that — the
figure also has to carry fewer things, which is why `slide` trims rows and rations
labels rather than only scaling the font.

Trimming is never silent: a trimmed axis says `— top 12 of 37`. Scatter plots keep
**every** point and ration only the labels, because dropping points would misstate
the distribution.

Slide runs land in their own `run_<stamp>-slide/` directory, so the two presets
never overwrite each other.

## Highlighting one cancer type

```bash
oncoref plot apd1-vs-tmb --highlight BRCA_TNBC
```

Picks one cancer type out of a landscape figure: it gets the highlight colour, a
larger mark, a dark ring and a bold label, while the rest fade toward the page but
keep their lineage hue so the landscape still reads as a landscape.

Nothing is filtered, re-ranked or re-scaled — the highlighted figure and its plain
twin show the same points in the same places. A highlighted point is labelled even
when it falls outside the `slide` preset's label budget, because it is the point of
the figure. Ranked slide charts retain the highlighted row alongside the highest
ranked rows. An unknown code, or a code absent from that plot, is an error rather
than producing an all-muted figure.

Only plots whose marks are per-cancer-code can carry a highlight. A burden-category
chart has no single row to light up, so those are excluded rather than approximated.

`scripts/regenerate_plots.py` also writes a `tnbc/` subdirectory of highlighted
landscape variants into every run, and those pages appear in the combined PDF
alongside the plain figures.

The highlighted type is `BRCA_TNBC`, the registry's receptor-defined TNBC cohort.
Its molecular measurements use primary TCGA tumors, while checkpoint-response
estimates come from metastatic TNBC trials. Highlighted response figures carry
that population note below the plotting axes. It remains attached to the saved
artifact for provenance and sits outside the core panel so the panel can be
cropped cleanly when reused.

## Style

All figures render publication-ready by default: opaque white ground (figure,
axes and saved raster), 300 dpi, top and right spines removed, a colour-blind-safe
categorical palette, and Type-42 fonts so text stays text in vector exports.

In-figure titles are stripped — a caption belongs in the manuscript, not burned
into the raster. The information a title carried lives in the axis labels. For
interactive work where a figure travels without a caption, call
`oncoref.figure_style.set_minimal(False)` to keep titles.

## Commands

### Reproductive and disposable CTA filtering

The [current report](audits/cta-tissue-scopes-20261019/index.md) follows all
registered published nominations through mapping, coding/family eligibility,
and reproductive or reproductive + disposable filtering. It shows the existing
adaptive thresholds and a separate strict positive-evidence sensitivity screen.

```bash
python scripts/cta_scope_report.py --out outputs/cta-tissue-scopes-current
```

Four PNG/vector-PDF figures cover source intake, filtering stages, retained
candidates by paper, and observed RNA/IHC compatibility and coverage. The combined
PDF includes methods and citations. CSVs identify gains, outside-tissue signals,
missing evidence and review holdbacks for every candidate. This report requires
the pinned HPA v23 bulk RNA and normal-tissue IHC files in addition to packaged data.

### CTA definition

How the cancer-testis-antigen panel is defined and filtered. Reads packaged data
only — no downloads, no cached matrices.

```bash
oncoref plot cta-curation
oncoref plot cta-curation --tissue-scope both --out outputs/cta-curation-both
```

From a source checkout, use `python -m oncoref.cli plot ...` at the repository
root to ensure the command uses the checkout rather than an older installation.

| Figure | Shows |
|---|---|
| `cta-source-venn.png` | Complete Wang, Bruggeman and da Silva coding nominations before HPA filtering |
| `cta-stage-funnel.png` | Sequential attrition, source union → shipped set |
| `cta-filter-funnel.png` | Final default panel vs outside default, per source |
| `cta-filter-outcome.png` | Default / HPA pass outside default / HPA fail, per source |
| `cta-deflated-frac-dist.png` | Where genes sit against the RNA thresholds |
| `cta-protein-vs-rna.png` | The tiered protein-reliability × RNA-fraction rule |
| `cta-publication-funnel.png` | All ten selected papers: source identities → mapped → coding → family → HPA → default |
| `cta-legacy-source-venn.png` | Historical CTpedia/CTexploreR/proteomics tags, labelled as an audit |
| `cta-landscape-source-venn.png` | Additional complete imported paper lists |
| `cta-source-overlap.png` | All-paper candidate and final-panel intersections |
| `cta-placental-source-overlap.png` | Gong, Bradley and historical placental nominations |
| `cta-placental-evidence-coverage.png` | Publication coverage of all 19 historical placental nominations |

`--tissue-scope both` writes all twelve figures under each of `core/` and `extended/`,
plus `cta-panel-comparison.png` and a gene-level `cta-panel-comparison.csv`.
Use `--tissue-scope extended` for that panel alone; core remains the default.
Every PNG has a vector PDF companion. CSVs record the plotted evidence, source
memberships/counts, publication intake, stages and thresholds. Source lists
overlap and cannot be summed. The source Venn and source bars cover coding
candidates; the publication intake also retains noncoding and unmapped rows.
Dashed/solid RNA thresholds come directly from the filter (0.80/0.90/0.95/0.97).
The [panel audit](audits/cta-panel-curation.md) explains the two definitions.

### Expression provenance

What reference expression data oncoref ships, where it came from, and how much
of the cancer space it covers.

```bash
oncoref plot expression-provenance            # US burden weighting
oncoref plot expression-provenance --region world
```

| Figure | Shows |
|---|---|
| `expr-source-samples.png` | Samples contributed by each upstream source |
| `expr-samples-per-code.png` | Depth per cancer type, in samples |
| `expr-source-quality.png` | Exact linear TPM vs proxy quantification scales |
| `expr-family-coverage.png` | Cancer types with and without a source, by family |
| `expr-burden-coverage.png` | Coverage weighted by incidence share |

### CTA panel design

```bash
oncoref plot cta-covering-set                      # q3 > 30 TPM
oncoref plot cta-covering-set --stat median --threshold-tpm 10
```

Greedy weighted set cover: at each step, the CTA adding the most still-uncovered
US incidence. Two curves — coverage by cancer type and by patients. Defaults to
`q3` rather than `median` because CTAs are subset antigens, so "at least 25% of
patients" is the meaningful bar; a median cut hides real targets.

Needs the expression bundle (percentile artifacts).

### Everything else

`oncoref plot --help` lists the full set, including the aPD-1/ICI response
scatters, incidence-vs-mortality, the CTA expression heatmaps, and the
per-patient coverage panels (which need cached per-sample matrices).

### All figures at once

```bash
python scripts/regenerate_plots.py
```

Writes every figure into one `figures/run_<timestamp>/` snapshot, organised by
family, plus an `index.md` listing what was produced and what was skipped, and a
combined `oncoref-all-figures.pdf` contact sheet. Figures that cannot be drawn (missing
per-sample matrix, empty data) are reported and skipped rather than aborting the
batch.
