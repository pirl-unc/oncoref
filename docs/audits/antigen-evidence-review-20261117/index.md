# Reproductive IHC and broader publication coverage

Catalog policy `oncoref.target-selection.2026-11-17.1`; publication-count correction
`oncoref.target-filters.2026-11-17.5`. No reproductive compartment is protected by
default. The IHC analysis below is a counterfactual; it changes no catalog memberships.

## How much would reproductive High IHC add?

Using pinned HPA v23 normal-tissue IHC, require High in at least one reproductive
cell compartment and Not detected in every reported non-reproductive cell/tissue
row. Both positives and negatives must have Approved, Supported or Enhanced
reliability. Every member locus of a protein group must have non-reproductive
observations. Tissue/cell mappings use the existing curation registry; shared
normal, coverage, translation, review, expression-bound and sequence gates remain.

| Set | Current native-TPM rule | With the additional IHC route | Change |
|---|---:|---:|---:|
| CTA | 290 | 299 | +9 (3.1%) |
| RA, including CTA | 403 | 423 | +20 (5.0%) |
| CA or RA | 509 | 520 | +11 |

The nine CA-only groups that would become CTAs are **DEFB136, LCN9, MAGEB18,
OTUD6A, PKDREJ, RBMY1B, RBMY1E, SP7 and SSX5**. All nine already pass tumor 10×;
their reproductive RNA folds are 1.20–8.92. Their IHC controls are Enhanced or
Supported; including Approved does not add further CTAs.

The eleven other potential RA additions are CDY1, CDY1B, CDY2A/B, DEFB132,
GOLGA6B, LYZL6, NANOS2, RBMXL3, RBMY1D, RSPH6A and SMIM21. These are shared-screen
eligible but do not meet the current cancer criterion. CTAG2 and DAZL have a
matching IHC pattern but retain their independent review holds.

**“All reported tissues” is not “all tissues.”** The nine potential CTAs have
negative observations in 35–37 non-reproductive source tissue labels. None has
complete coverage of every non-reproductive tissue represented anywhere in this
HPA export. Missing tissues are listed explicitly in the downloadable audit.
The consensus export does not resolve matched antibody/run/donor pairs, and
qualitative High/Not detected is not a calibrated 10× protein ratio. A named
protein-evidence RA route would require individual antibody and cell-type review.

SP7 illustrates the coverage issue: human variants implicate it in bone biology
([2010 human genetics](https://pubmed.ncbi.nlm.nih.gov/20579626/),
[2022 mechanistic study](https://pmc.ncbi.nlm.nih.gov/articles/PMC8816926/)).
These findings do not themselves establish normal adult protein abundance or
therapeutic toxicity. [Issue #707](https://github.com/pirl-unc/oncoref/issues/707)
tracks adult bone-cell coverage before using IHC to admit it as an RA.

[Every reproductive-High group, negatives, missing tissues and remaining gates](ihc-counterfactual.csv.gz)
· [Source hashes and counterfactual totals](ihc-counterfactual.json).

## Zero counted tumor papers does not mean no literature

Of the original 49 CTAs with zero counted tumor-support publications:

- **43 have curated literature links** after joining follow-ups by exact gene ID.
- **41 have at least one coverage-eligible record**, which may describe development,
  normal expression or function. That is not tumor validation.
- **27 have eligible reproductive/developmental publications.**
- Six have no linked record in this registry: **KRTAP20-4, SPADH, SSU72L2, SSU72L3,
  SSU72L4 and SSU72L5**. This is a curation gap, not a claim that no paper exists.

| Protein(s) | What the literature contributes | Tumor-support credit |
|---|---|---|
| DPRX | [Human embryo transcripts and regulatory experiments](https://pmc.ncbi.nlm.nih.gov/articles/PMC4944136/); another developmental study is also curated | 0 |
| TPRX2 | Two existing developmental papers plus newly added [human embryo translation and TPRX perturbation evidence](https://pubmed.ncbi.nlm.nih.gov/36074823/) | 0 |
| NANOGNB | [Human preimplantation expression and transcriptional function](https://pmc.ncbi.nlm.nih.gov/articles/PMC5413911/); a separate mutation-specific immune record remains ineligible for reference-protein coverage | 0 |
| FOXL3; SSU72L6 | **Testis-enriched RNA**, in the comparison table of the [placental transcriptome study](https://pmc.ncbi.nlm.nih.gov/articles/PMC8113443/); no demonstrated reproductive function | 0 |
| SCGB1D4 | Normal lymphoid and [nasal expression/cytokine-response evidence](https://link.springer.com/article/10.1186/1465-9921-12-28); retained as counterevidence with resting-pattern and assay-specificity limitations | 0 |
| SSU72L1 | A [cancer-cell RNA perturbation table](https://pmc.ncbi.nlm.nih.gov/articles/PMC7868764/) names its previous symbol SSU72P7; no copy-specific protein validation | 1, after fixing the join |

The broader curation remains available independently of prioritization filters.
The base catalog requires zero papers. `min_papers` counts qualifying nominations,
tumor-expression or cancer/embryo-program records; reproductive papers and normal
counterevidence have their own roles. Counting every biological paper as tumor
validation would erase the claim distinction. Distinct publications may reuse the
same cohort and are not independent patient replication.

### Central paper → gene → claim records

The central current ledger, `oncoref/data/cta-literature-followup.json` schema 2,
contains one record per reviewed paper with its exact gene list, native gene names,
reported context, assay, claim, limitations and allowed uses. It currently has
13 reviewed papers and 16 gene-specific records. Original nomination tables remain
immutable source inputs; imported scope is labeled separately from a reviewed
gene claim. New reviews belong in this ledger, including corrections to old papers.

The original Gong workbook was retrieved again and its SHA-256 matched the pinned
`b91d126e9e0a270309cb3f6ad5253f1153db3a21f0a8c19d8ef55080410a2007`.
FOXL3 appears as **WI2-2373I1.1**, row **1274**, and SSU72L6 as **SSU72P8**,
row **1242**, on `Data 5 - tissue-enriched PC`. Both rows say **Testis**.
The earlier description of these as placental evidence was incorrect. A paper's
title cannot substitute for the tissue and gene on the source row.

The catalog builder, API and publication filters consume resolved per-gene uses;
source-wide roles cannot override them. The old source records remain available
for provenance. A normal-expression concern cannot earn reproductive or tumor
credit, and an expression study cannot silently become a functional validation.
The TPRX knockdown cited by the 2023 LEUTX paper is recorded under its original
2022 paper, not counted as a second independent functional experiment.

[Reviewed papers, gene lists and their uses](reviewed-paper-claims.md)
· [Structured claims](reviewed-paper-claims.json)
· [Version and filter comparison](version-comparison.md).

### Corrected joins and assay scope

The catalog already retained SSU72L1's follow-up, but the filter joined on an absent
historical protein-group ID. `current_attestations()` now joins exact member loci
to current groups while preserving source records and historical identity. It
does not transfer a nomination across a historical group split. The filter also
honors `eligible_for_coverage`: an alternative GALP proteoform and a YAP1::LEUTX
fusion paper no longer count for their selected reference proteins. Their records
remain visible. `current_coverage()` now uses the current target catalog, replacing
its stale legacy-catalog input. [Issue #706](https://github.com/pirl-unc/oncoref/issues/706).

The same 97 CTAs still fail the standard publication step, now **48 with zero and
49 with one counted source**. Current CTA selections are **290 all, 238 loose,
175 standard and 110 strict**. SSU72L1 gains loose membership; LEUTX loses strict
membership; their biological CTA classification is unchanged.

```python
from oncoref import cta_literature, target_filters

cta_literature.current_attestations("SSU72L1")  # Sources, scopes, exact gene IDs.
cta_literature.current_coverage()              # Broader coverage of current CTAs.
cta_literature.papers("DPRX")                  # Paper records with exact gene lists.
cta_literature.paper_gene_claims("DPRX", use="reproductive_support")
target_filters.paper_counts()                 # Separate families and counted tumor support.
```

`oncoref antigens explain DPRX` displays both its two curated sources and its zero
counted tumor-support publications. The new TPRX translation paper is recorded in
the shared follow-up ledger for TPRX1 and TPRX2; it adds no inferred TPM value or
automatic tumor credit.

`oncoref antigens explain FOXL3 --papers --format json` includes the original
name, Testis context, supplementary row, allowed uses and limits. The DAG rebuilds
these paper views and the IHC audit with the catalog; they are generated views of
the central curation, not additional editable curation stores.

[All original 49, with broader counts](original-49-literature-coverage.csv)
· [Every linked paper and claim scope](original-49-paper-links.csv)
· [Summary](literature-summary.json).

Reproduce from the pinned checkout with
`python scripts/audit_reproductive_ihc_route.py` and
`python scripts/audit_antigen_paper_counts.py`. These audit scripts do not edit
catalog membership. The current catalog and reports are rebuilt separately by
`python -m oncoref.curation build --reuse-tumor-cache`.
