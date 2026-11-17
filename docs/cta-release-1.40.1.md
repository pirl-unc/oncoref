# OncoRef 1.40.1: preserve cohorts in shared expression summaries

A targeted TGCT rebuild in data bundle 5.23.28 replaced a shared TCGA summary
file with TGCT-only rows. The other 31 cohorts were still advertised by the
source registry but returned no rows through the summary-based expression API.
This affected downstream heatmaps and response comparisons using those summaries.

Data bundle **5.23.29** restores **1,071,701 rows across 31 cohorts** from 5.23.27
as code-specific source shards. Every restored cell was compared with the prior
bundle; the corrected TGCT file is preserved byte for byte. The other 5.23.28
artifacts, including HL, are unchanged. The public API was checked for all 32
cohorts, including TGCT.

The incremental merge now uses the existing cohort-aware summary merge and
unlinks hardlinked destinations before writing. A regression test verifies that
a focused rebuild preserves unrelated cohorts and every base-file checksum.

The release also fixes an HPA cache race found by parallel CI: readers could
open a Parquet file while another process was still writing it. Cache files
are now published atomically; damaged optional caches fall back to the source
TSV. A synchronized concurrent-reader regression covers the failure. This
changes cache I/O only, not the HPA measurements or their interpretation.

The [1.40.0 curation policy](cta-release-1.40.0.md) and antigen memberships are
unchanged: **246 CTA, 149 CA-only, 102 RA-only, 1,143 Other and 18,253 Unsafe**.
Raw source matrices remain **5.22.14**. This is a derived-summary restoration,
not a raw-data reanalysis.

[Repair audit](audits/summary-overlay-20261116/repair.json) ·
[Public API verification](audits/summary-overlay-20261116/owner-api-verification.csv) ·
[Bundle manifest](audits/summary-overlay-20261116/manifest.json) ·
[Current curation report](audits/target-selection-20261110/rebuild/report.html)
