# Rebuild the CTA analysis reports

The scripts in `scripts/cta_*` are a reproducible analysis pipeline. Install the
report extra in a checkout (`pip install -e '.[dev,reports]'`). They read cached
source matrices and write to an explicit output directory. Historical local
`outputs/cta_*_20260916` and `20260917` directories predate the corrected analysis:
their thresholds, selections and narrative must not be used as current results.
`outputs/` is ignored by Git; generated reports are not shipped in the wheel.

The corrected analysis uses the intersection of canonical biological gene IDs
across all registered cohorts as its fixed percentile background. It writes the
IDs to `common_background.csv`. The proteoform analysis collapses exactly those
loci for its background, while candidate sums retain all registered member loci.
Missing background measurements stop analysis; missing candidate measurements
cannot qualify. The gene screen requires at least 10 patients; the proteoform
shortlists require at least 20. The gene p30/p50 screens and the details script's
`all` cohort stratum are exploratory, including degenerate zero cutoffs. Gene
exports count loci, which can encode identical proteins; use the proteoform
shortlists to count distinct protein identities.

TCGA donor groups retain one specimen class, preferring the primary tumor when
both primary and metastatic specimens are present. Within a class, repeat
specimens are averaged in linear space. The sample audit records exclusions.
MTC requires an explicit GSE32662 series matrix with donor titles; its contents
are hashed, and missing donor mappings stop the run.

Run the following in order, using new output folders. Supply the actual local
cache and donor-matrix paths; do not reuse a historical output directory.

```sh
export CANCERDATA_PER_SAMPLE_CACHE=1
export MPLCONFIGDIR="$PWD/tmp/cta-matplotlib"
SOURCE_CACHE="$PWD/tmp/cta_threshold_report/source-matrices"
GENES="$PWD/outputs/cta_current/genes"
PROTEINS="$PWD/outputs/cta_current/proteins"
MTC_MATRIX="/path/to/GSE32662_series_matrix.txt.gz"

python scripts/cta_threshold_report.py --out "$GENES" \
  --source-cache "$SOURCE_CACHE" --mtc-donor-matrix "$MTC_MATRIX" --no-plots
python scripts/cta_proteoform_report.py --base "$GENES" --out "$PROTEINS" \
  --source-cache "$SOURCE_CACHE"
python scripts/plot_cta_proteoform_report.py --out "$PROTEINS"
python scripts/cta_hpa_tissue_report.py --out "$PROTEINS"
python scripts/cta_mortality_coverage_report.py --out "$PROTEINS"
python scripts/cta_primary_panel_report.py prepare --out "$PROTEINS"
python scripts/cta_burden_coverage.py --out "$PROTEINS"
python scripts/render_cta_proteoform_report.py --out "$PROTEINS"
python scripts/cta_decision_guide.py --out "$PROTEINS"
```

The protein FASTA must match the pinned Ensembl 112 checksum; the default path
is the Ensembl 112 pyensembl cache. A different release is rejected before any
annotation is written. The two analysis scripts accept `--resume` and verify
source, implementation, policy and checkpoint-output hashes before reuse.

For supplementary gene-level cohort and patient views, run
`cta_threshold_report_details.py --base "$GENES" --out <new-details-directory>
--source-cache "$SOURCE_CACHE"`, followed by
`render_cta_analysis_all_figures.py --out <new-details-directory>`.

Each analysis writes a validation result, a run manifest and a content-hash
receipt. Downstream stages verify the relevant inputs and outputs; a failed or
stale run is rejected. PDF text, list counts and bookmarks come from current
CSV data and the page index. Every indexed figure category is included.
PDFs, gzip tables and ZIP metadata are deterministic. The ZIP preserves analysis
checkpoints and provenance for inspection. Layout and image-library versions
can still affect rendering across environments.

Mortality comparisons use broad histology cohorts (PAAD for pancreas; HNSC for
head/neck). They do not use pancreatic neuroendocrine tumors or thymoma as
proxies. These are observed cohort estimates, not worldwide patient prevalence.
Burden represented is a separate sum of categories touched by a qualifying
cohort; even a subtype hit represents the parent category and is not a
population coverage estimate.

Normal-tissue RNA uses within-tissue sums of the same loci as tumor RNA. Each
HPA member measurement is rounded: a reported zero is below 0.05 nTPM, and a
sum of two reported zeros is below 0.1 nTPM. All tied maximum tissues are retained.
IHC observations, unavailable measurements and RNA estimates remain distinct;
none establishes peptide presentation or clinical safety.
