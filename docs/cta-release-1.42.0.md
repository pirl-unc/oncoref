# OncoRef 1.42.0: antigen sets, explicit host protection and gene-specific literature

This release combines PRs #699 and #700. The catalog contains **290 CTA,
106 Cancer Antigen only, 113 Reproductive Antigen only, 1,131 Other and
18,253 Unsafe proteins**, covering **19,893** selected-reference protein groups.
CA includes 396 proteins, RA includes 403, and their union includes 509.

## Broad membership and optional protection

CA and RA independently require 10× enrichment over the same non-reproductive
reference; CTA is their intersection. Reproductive scope includes testis, ovary,
embryo, placenta, prostate, epididymis, seminal vesicle, fallopian tube,
endometrium, cervix and both mammary states. Prostate and mammary expression do
not enter the default non-reproductive denominator.

**No reproductive compartment is protected by default**, including for the
`male-patient` and `female-patient` anatomy presets: each retains all 290 CTAs
before optional evidence filters. Protection is an explicit host choice and
does not rewrite biological catalog membership or the global Unsafe set.

```python
from oncoref import antigen_sets
from oncoref.host_context import HostContext

all_ctas = antigen_sets.select("CTA")
protected = antigen_sets.select("CTA", host=HostContext(protect=("ovary",)))
proteins = protected.proteins()
background = protected.background_proteins("Unsafe|Other")
```

Host choices distinguish **Not protected**, **Protect**, **Accept loss** and
**Absent**. Explicit protection requires measured coverage, resolved evidence
and RNA below separately applied native-TPM and nTPM limits (default <4 in each
unit). Protected TPM also enters tumor fold, independent-patient support and
active absolute off-target limits. Missing observations never become zero.
Protecting only ovary retains 252 CTAs; protecting only testis retains 62.
Pregnancy selection remains unresolved because the embryo panel does not cover
all developmental stages. These research screens do not establish injury
thresholds or clinical safety.

## Evidence, API and reports

- `antigen_sets.select()` carries one selection through protein, gene, annotated
  sequence and peptide-background views. Host-excluded or unresolved candidates
  enter that selection's subtraction background; original catalog labels remain
  available as provenance.
- Mammary sample labels retain an expression-inferred milk-program state,
  marker measurements and missingness. They are not provider-confirmed
  lactation or resting-state diagnoses. The older `cta_tier` annotation remains
  a compatibility description of biological origin, not a toxicity ranking;
  `core` is no longer a host-protection shortcut.
- The central publication ledger records exact gene lists, native names,
  assays, context, limitations and permitted uses. `cta_literature.papers()`,
  `paper_gene_claims()` and `oncoref antigens explain GENE --papers` expose the
  same records consumed by the catalog and publication filters. Reproductive
  RNA, translation, functional experiments and normal-tissue concerns retain
  distinct claims. SSU72L1's SSU72P7 cancer-cell record now receives credit;
  FOXL3 and SSU72L6 correctly retain Testis RNA evidence from Gong's table.
- The broad CTA catalog remains 290. Optional modes retain 238 loose,
  175 standard and 110 strict. Standard's paper step removes 97 proteins for
  fewer than two counted tumor-support publications, not for lacking all
  reproductive or functional literature.
- The curation DAG rebuilds the catalog, exposure evidence, figures, filter
  report, paper views and IHC counterfactual. The release receipt verifies
  **221 outputs**. The IHC alternative is an audit only and is not enabled.

Expression bundle **5.23.29** and raw source matrices **5.22.14** are unchanged.

[Current rules](target-selection.md) ·
[Host selection](host-selection.md) ·
[Filters](target-filters.md) ·
[Version comparison](audits/antigen-evidence-review-20261117/version-comparison.md) ·
[Reviewed papers](audits/antigen-evidence-review-20261117/reviewed-paper-claims.md) ·
[Commands and build workflow](curation.md) ·
[Filterable catalog](audits/target-selection-20261110/rebuild/report.html)
