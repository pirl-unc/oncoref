# Reproduce the literature audit

This directory is a review overlay on the 25 October 2026 checkout. It does not
change `oncoref.target_sets`, the released gene sets, protein sets or exclusions.

```bash
# Public Europe PMC/PubMed searches; resumable, rate-limited and cached.
.venv/bin/python scripts/audit_cta_literature.py

# Optional targeted primary papers; metadata and available open-access full text.
.venv/bin/python scripts/fetch_cta_review_papers.py 28386150 39319391

# Join the frozen atlas, source flags, search ledger and explicit decisions.
.venv/bin/python scripts/build_cta_literature_audit.py
```

The query script requires network access. Rebuilding the review after retrieval
is offline. Dependencies are the repository's pandas/requests/PyYAML environment.
The builder verifies cached response checksums, refuses incomplete retrieval by
default and preserves the original category in every row. `--allow-incomplete`
creates a progress report explicitly labelled incomplete.

## Files

| File | Meaning |
|---|---|
| `explore.html` | Searchable, offline report with one entry per mapped gene |
| `all-genes.csv` | Original evidence plus assessment, review depth and proposed category |
| `proposed-changes.csv` | Recommendations that differ from the released category |
| `adjudications.yml` | Explicit gene-specific reasoning and primary references; manually maintained |
| `searches.csv` | Exact queries, timestamps, hit counts, truncation, response hashes and cache names |
| `papers.csv` | Deduplicated gene–paper leads with titles and links; not automatically accepted evidence |
| `retained-primary-normal-signals.csv` | Discordant primary-antibody rows for retained genes |
| `summary.json` | Counts and input checksums |
| `search-summary.json` | Search run scope and status |
| `category-proposal.md` | Proposed general evidence and category rules; not an implemented policy |
| `issues.md`, `issue-links.json` | Filed GitHub issues and their durable links |
| `issue-drafts.json` | Exact self-contained child issue bodies filed from this audit |

Full abstracts and public full texts stay in the local cache and are not copied
into the public report. Default caches are
`tmp/cta-literature-20261025/raw/` and `tmp/cta-literature-20261025/targeted/`.
These are ignored working files, not versioned source assets. Keep/copy them for
an exact offline rebuild, or use `--cache`/`--out` to put them in durable storage.
The versioned ledger contains the exact public queries and hashes; a later fresh
download can change as Europe PMC updates its records and is not guaranteed to
reproduce an earlier checksum. It must be reported as a new retrieval.

## Scope and limitations

Every mapped nomination receives two title/abstract searches (normal expression
and therapy/safety), using its current symbol and selected annotation-linked
aliases. Existing tissue-specific UniProt PubMed references are fetched in a
separate query so a relevance/date-ranked search cannot bury them. New searches
are limited to PubMed-indexed records dated on or before 25 October 2026 and the
first 20 results per topic. `Hits`, `Retrieved` and `Truncated` expose that limit.
The pinned annotation-citation query retrieves its explicit PMID list separately.

This is **not a systematic review**: no complete pagination, citation-network
search, independent second reviewer or universal full-text adjudication is
claimed. Searches can miss papers that do not name a gene in the title/abstract;
older aliases are incomplete; short names and family labels can match unrelated
papers. In particular, ABCC11's historical MRP8 alias also names S100A8. A search
hit is never sufficient to change a category. Zero hits never establish absence.

`Assessment basis` distinguishes newly inspected primary full text/abstracts,
prior targeted reviews carried forward, existing source evidence, and unresolved
rows. Full text means the relevant methods/results were examined; it does not
mean the experiments were independently replicated. An abstract-only assessment
cannot certify controls or figure details that the abstract does not report.

Primary-antibody flags use the already mapped `normal_signal` column from the
pinned normal-cell audit. They do not undo its explicitly ignored reagents, and
do not equate raw staining with an accepted HPA protein annotation. Profile
warnings, family-shared reagents and the cell compartment remain visible.

The 2,540-gene table covers mapped nominations, not the whole human genome and
not unmapped historical nominations such as BAGE/GAGE2B. Those intake gaps remain
outside this denominator; no paralog is substituted.

## Checks

```bash
.venv/bin/python -m pytest tests/test_cta_literature_audit.py -q
.venv/bin/ruff check scripts/audit_cta_literature.py \
  scripts/fetch_cta_review_papers.py scripts/build_cta_literature_audit.py \
  tests/test_cta_literature_audit.py
```

The tests cover duplicate-query cache consistency, cache integrity, keeping
search leads separate from accepted findings, complete per-gene coverage,
preserving released categories, and rejecting hidden >30 normal-RNA failures in
Neither. They do not claim biological or clinical validation of the thresholds.

The separate reusable PDF slide scripts are documented in
`scripts/slides/README.md` in the repository.
