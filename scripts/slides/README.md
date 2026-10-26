# Reproducible report slides

These scripts reproduce `output/pdf/cta-screen-logic-20261024.pdf`, including
all 26 slides, the two funnel diagrams, threshold tables, reason-count charts,
references and PDF navigation. They do not change the gene classifications.

## Run

Use Node 22+ with `@oai/artifact-tool` available and Python with OncoRef's
`reports` extras (`reportlab`, `pypdf`, `Pillow`). The original renderer was
`@oai/artifact-tool` 2.8.99. This is a separate rendering dependency, not part of
OncoRef's Python runtime. If supplied by an artifact-tool runtime, point
`ONCOREF_ARTIFACT_TOOL_MODULE` to its `dist/artifact_tool.mjs`; alternatively use
`--artifact-tool /path/to/dist/artifact_tool.mjs`. No machine-specific paths are
embedded in the scripts.

From the repository root:

```bash
node scripts/slides/cta_logic.mjs --out-dir tmp/pdfs/cta-logic-slides
python scripts/slides/assemble_pdf.py \
  tmp/pdfs/cta-logic-slides output/pdf/cta-screen-logic-20261024.pdf
```

`--help` works without loading the rendering dependency. Paths are resolved
relative to the scripts/repository, so the builder also works from another
working directory. The output directory is configurable.

The builder writes high-resolution slide PNGs, layout JSON, a content/link
manifest, and an editable `source-deck.pptx`. The assembler packages the rendered
slides as an **image-based PDF**, with bookmarks and clickable source links.
This is the original export method: native PDF export had a missing standard-font
resource in the available runtime. The PDF does not contain searchable slide text;
the manifest and native presentation preserve that text.

## Reuse and update

- `report_deck.mjs`: reusable native slide helpers for titles, citations, tables,
  bar charts, diagrams, footers and manifests. Category colors and page geometry
  live here. Charts use explicit fonts and category axes.
- `cta_logic.mjs`: the reviewed narrative, figures and policy snapshot. It reads
  the checked-in release/review summaries and checks their counts before rendering.
- `assemble_pdf.py`: generic packaging, bookmark/link checks and atomic output
  replacement; it derives page counts and dimensions from the manifest.

The CTA narrative and thresholds are deliberately frozen to policy
`oncoref.target-sets.2026-10-24.1`, source commit `2e40a0f`. It is **not** a live
classification engine. For a new policy, update the report sources, numerical
assertions, figure values, narrative, citations and snapshot labels together.
Do not simply remove failed assertions or relabel old slides as current.

To create another report, use `createReportDeck(Presentation, metadata)` from
`report_deck.mjs`; author its content separately and export the same manifest
and numbered slide images. The PDF assembler accepts any nonempty deck using
that contract.

## Verify

Check counts against the source tables, then render and inspect every PDF page:

```bash
pdftoppm -scale-to 1280 -png output/pdf/cta-screen-logic-20261024.pdf \
  tmp/pdfs/cta-logic-slides/check
```

Pay particular attention to chart labels, two-line titles, table wrapping,
formula symbols, arrows, denominators and citations. Source updates can alter
layout even when the export succeeds. Generated artifacts remain ignored by Git.
