# Figures

Every figure oncoref ships is generated from packaged data by `oncoref plot`.
Nothing is hand-drawn, and nothing is committed to the repo — regenerate instead
of copying, so a figure always matches the data version that produced it.

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

## Style

All figures render publication-ready by default: opaque white ground (figure,
axes and saved raster), 300 dpi, top and right spines removed, a colour-blind-safe
categorical palette, and Type-42 fonts so text stays text in vector exports.

In-figure titles are stripped — a caption belongs in the manuscript, not burned
into the raster. The information a title carried lives in the axis labels. For
interactive work where a figure travels without a caption, call
`oncoref.figure_style.set_minimal(False)` to keep titles.

## Commands

### CTA definition

How the cancer-testis-antigen panel is defined and filtered. Reads packaged data
only — no downloads, no cached matrices.

```bash
oncoref plot cta-curation
```

| Figure | Shows |
|---|---|
| `cta-source-venn.png` | Overlap of the three primary source databases |
| `cta-stage-funnel.png` | Sequential attrition, source union → shipped set |
| `cta-filter-funnel.png` | Kept vs excluded, per source |
| `cta-filter-outcome.png` | Kept / kept-but-weak / excluded, per source |
| `cta-deflated-frac-dist.png` | Where genes sit against the RNA thresholds |
| `cta-protein-vs-rna.png` | The tiered protein-reliability × RNA-fraction rule |

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
