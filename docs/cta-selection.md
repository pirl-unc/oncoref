# Current CTA and oncofetal selection

The current policy supplies **522 selected human proteins, 18,031 Unsafe, 20 Neither and 1,290 unassigned** across 19,863 proteins. These are research-policy assignments, not demonstrated clinical safety. The [current report](audits/cta-source-walkthrough-20261028/index.md) includes sources, figures, the funnel and tissue-specific explanations. The [1.21.0 report](audits/cta-release-1.21.0/index.md) remains frozen.

## Safety, cancer expression and Reproductive evidence

**Only safety can establish Unsafe.** Cancer expression and Reproductive evidence
are separate admission requirements. Missing required coverage remains unassigned.

**Reproductive** is the umbrella for adult reproductive, placental, extraembryonic,
embryonic and fetal contexts. Branches describe the actual tissue and stage;
“embryonic-only evidence” does not prove absence from unmeasured branches.

| Assessment | Question | Current rule |
|---|---|---|
| Safety | Is normal-tissue expression concerning under the tissue/IHC policy? | Absolute RNA ceilings, separate positive-protein rules and reviewed clinical exclusions. Vital/Protected exclusions can establish Unsafe. Manageable excess blocks CTA. |
| Cancer expression | Is there expression in at least one measured cancer? | One sample ≥2 TPM or HPA pTPM. No tumor-versus-normal ratio. |
| Reproductive evidence | Is there positive expression in a qualifying branch? | HPA adult/placental maximum ≥5 nTPM OR supported embryonic lineage mean ≥5 TPM. |
| Coverage | Is the comparison adequately observed? | Complete GTEx v11 core: 52 bulk contexts for every member locus; required normal bulk/cell/immune panels also apply. DICE/cornea coverage remains optional. |

The former maximum-versus-maximum restriction rule is removed. The
[continuous diffuseness gate](audits/cta-unassigned-20261029/diffuseness-impact.md)
uses **HPA v23 consensus**: average `N / (N + R)` over **29 equally weighted
anatomical units / 38 contexts**. R is the highest adult/placental Reproductive
nTPM in the same matrix. Disposable tissues count; **Permissive thymus** belongs
to neither the comparator nor the Reproductive reference. Its safety tier is unchanged.

A score **above 10%** blocks the adult/placental route when R meets the existing
**5 nTPM** floor and the full HPA comparison is measured. Missing comparisons
remain explicitly unassessed. An independently supported embryonic route remains
eligible; its concentration needs a compatible developmental/adult panel. GTEx
remains the required coverage core and a descriptive comparison, not a second veto.

A complete diffuse profile with no remaining qualifying route and clear normal
safety becomes **Neither**. This moves **MTCP1 (24.53%)** from CTA to Neither,
with no change to Unsafe membership. The result is stable from 7.5% through 20%.
[Logic figure and examples](cta-evidence-axes.md).

## One API for each species

```python
from oncoref import cta_selection

selected = cta_selection.results(species="human", category="CTA")
all_decisions = cta_selection.results(species=9606)
scope = cta_selection.coverage(species="dog")
sources = cta_selection.sources(species="mouse", category="developmental_expression")
```

Equivalent bound accessors are `human.cta_selection()`, `canine.cta_selection()` and `mouse.cta_selection()`, with `cta_selection_coverage()` and `cta_sources()` on each module. They also live on `oncoref.for_species("human").antigens`. Species names, aliases and taxonomy IDs use one resolver. Unknown species and categories fail explicitly; nothing silently substitutes human data.

| Species | Returned entities | Selection status | Meaning |
|---|---:|---|---|
| Human | 19,863 proteins (identical selected reference sequences collapsed) | Available; evidence partial | Complete evaluation of the named reference universe under `oncoref.cta-selection.2026-10-29.3` |
| Dog | 14 curated native nominations | Unavailable; evidence partial | Unassigned / unassessed; these are assay or family nominations, not 14 proven CTA genes |
| Mouse | 0 | Unavailable | No integrated nomination/selection snapshot; this does not mean mice have no CTAs |

The contract is identical even when the evidence is not. Human thresholds, antibody grades, developmental measurements and clinical observations are never inherited through orthology. Canine matrices still need unit/reference reconciliation and a species policy; the registered tumor source has unresolved unit/reference provenance. Mouse sources are identified, not ingested.

## Result fields and identities

`results()` returns a detached pandas DataFrame. `taxon_id` is integer; scalar text uses pandas string dtype; the four collection fields are tuples. Empty results have the same columns and types. Reads do not download source matrices.

| Field | Meaning |
|---|---|
| `taxon_id`, `species` | Explicit organism boundary |
| `entity_id`, `entity_kind`, `label` | Stable result key; human protein-equivalence group versus native nomination; display label |
| `gene_ids`, `gene_identities` | Every supported member/alternative; IDs, symbols, namespaces, historical/alternate IDs, synonyms and identity status |
| `identity_status`, `reference_id` | Resolution completeness and exact reference receipt |
| `category`, `assessment_status`, `policy_id` | CTA / Unsafe / Neither or null; screened, incomplete evidence, outside category definitions or unassessed; versioned policy |
| `reasons`, `evidence_flags` | All contributing exclusions and retained observations, including Low protein, antibody review and stem/progenitor expression |
| `source_snapshot_id` | Content identity of the decision snapshot |

Unknown history is `None`, not an empty list pretending the search was exhaustive. Ambiguous canine identities remain alternatives and do not become protein-equivalence groups. Human reference identity pins the Ensembl-115 primary gene-space projection; it does not claim exhaustive sequence coverage. A representative-protein equivalence group is also not proof that every alternative isoform is identical. The complete locus crosswalk retains identities and explicitly labels earlier decisions as historical/baseline.

TPM/nTPM is summed only for matched group-member measurements in the same source context. Missing members do not become zeros. Unpaired DICE distributions and corneal detection bounds retain uncertainty. A single measured tumor member above threshold establishes a group lower bound; maxima from different samples are not added.

## Uniform human rules

| Gate | Rule |
|---|---|
| Adult/placental Reproductive evidence | HPA maximum ≥5 nTPM; no normal maximum comparison |
| Embryonic Reproductive evidence | Supported Tyser or Petropoulos lineage mean ≥5 TPM; ≥20 cells, ≥5 positive cells and ≥10% positive at ≥1 TPM |
| Tumor expression | ≥2 in one measured sample, in its native TPM or HPA pTPM; no normal comparison |
| Required coverage | GTEx v11 core (52 native bulk contexts) and specified bulk/cell/immune panels must be complete; DICE and cornea coverage is optional, but their measured signals count |
| Normal RNA | One inclusive tissue/IHC ceiling; above it supplies Unsafe evidence in Vital/Protected and blocks CTA in Manageable |
| Normal protein | Low remains a positive flag; stronger HPA or accepted unquantified literature positives retain their hold/exclusion; validated external normal positives and clinical exclusions remain active |

| Candidate ceiling | Missing / unusable IHC | Approved negative | Supported negative | Enhanced negative |
|---|---:|---:|---:|---:|
| Vital | 4 | 8 | 16 | 32 |
| Protected | 8 | 16 | 32 | 64 |
| Manageable | 16 | 32 | 64 | 128 |

These are policy parameters, not measured toxicity limits. TPM/nTPM assays share nominal limits while retaining their units; corneal CP10K uses nominal limits divided by 100, **not a TPM conversion**. HPA clusters require ≥20 reads; DICE requires three reported values; cornea requires ≥5 positive cells and ≥20 counts in at least two sample groups. The source-specific support statistics are preserved rather than calling every measurement a comparable sample maximum.

Exact adult tissue/cell matching or a declared source-specific histological-compartment correspondence, reliable negative staining, a positive profile control and no unresolved antibody conflict are required for IHC credit. Every expected member needs a positive profile control for a group negative. Low positive staining grants no negative-IHC credit. Author-described Low and HPA Low retain different scales, with unknown fractions/H-scores left unknown.

Mapped HPA adult fibroblasts use Manageable; bulk heart and cardiomyocytes retain their organ tiers. Colon/rectum distal enterocytes, goblet cells, Paneth cells and undifferentiated cells use Manageable in both RNA and HPA protein assessment. Undifferentiated expression retains a stem/progenitor flag. This permissive mucosal assumption and the numeric cutoffs have not been clinically validated. Reproductive/placental immune cells receive neutral anatomical treatment, not an assertion that all immune cells are disposable. Adult normal panels do not establish pediatric safety.

Every layer evaluates every group. **Unsafe** takes precedence when Vital/Protected or clinical exclusion evidence is present. **Neither** requires complete normal RNA coverage, resolved evidence and one of two measured profiles: all measured RNA below the shared 2-unit expression floor with no detected protein, or no reproductive/developmental association and RNA ≥2 in at least two Manageable tissue concepts with no Vital/Protected exclusion. A third route is a complete diffuse adult/placental identity profile with no remaining qualifying reproductive route and no normal-safety hold; it does not require cancer expression. The low-expression route additionally requires mapped tumor and developmental measurements. The current counts are reported in the generated funnel; all measured-profile requirements must hold.

Other failures are **unassigned**, with `category` null and an explanatory `assessment_status`; they are not automatically Neither. Only groups passing all CTA gates are CTA. The final funnel partitions CTA, Unsafe, Neither and unassigned counts. Intermediate rows display candidate eligibility, Unsafe and unassigned holds until the final Neither-profile evaluation. Final assignments do not depend on layer order.

## Evidence categories across species

`sources()` returns one row per source/category with native unit, scope, URL and integration status. [Paper titles and database versions](audits/cta-source-walkthrough-20261028/references.md) provide readable links to explanatory pages, with artifact URLs retained separately. A source can serve several roles; source-category rows must not be summed as independent studies or donors.

| Category | Human selection | Dog resources | Mouse resources identified |
|---|---|---|---|
| Nomination | Published CTA/placental lists as annotations | Two native CTA studies | No curated list |
| Reproductive expression | HPA bulk | Bgee matrices, BarkBase, native study assays | Requires reproductive panel integration |
| Developmental expression | Tyser gastrulation | BarkBase matrix; DoGA metadata | MOCA, gastrulation atlas, GXD |
| Tumor expression | 132 compatible OncoRef cohorts + HPA cancer | Native studies; cohort provenance unresolved | Native tumor cohort needed |
| Normal bulk RNA | HPA consensus | Bgee, BarkBase; DoGA/Pfizer RNA metadata | GXD assay resource |
| Normal cell RNA | HPA tissue clusters, adult cornea | No equivalent integrated atlas | Tabula Muris/Senis |
| Immune RNA | HPA, Monaco, DICE | No equivalent integrated panel | Tabula Muris/Senis |
| Normal protein | HPA IHC, curated external IHC/protein | Pfizer iBAQ matrix, native IHC study | GXD assay resource |
| Clinical safety | Curated product-specific findings | No transferred human evidence | No transferred human evidence |

`coverage()` distinguishes sources used in selection, available observations/matrices, and resources merely identified. DoGA's 5′ expression assay, BarkBase TMM-CPM, Pfizer iBAQ and mouse UMIs cannot be silently labeled TPM. Source inventories and species limitations are checksum-pinned in the package.

## Versioning and reproduction

Use this API for the new funnel. Older `cta_genes`, nomination snapshots and exact-protein forbidden-background APIs remain independently versioned compatibility/reference products. Their counts must not be presented as this selection; a CTA admission does not automatically remove a protein from a near-self background. See [the exact-protein policy](target-sets.md).

Run from a checkout with the pinned reference inputs available:

```sh
python scripts/regenerate_cta_walkthrough.py
python scripts/audit_cta_profile_concentration.py
node scripts/build_cta_source_slides.mjs
python scripts/regenerate_cta_walkthrough.py --finalize-receipts
```

The slide script requires the runtime paths documented in its header. `--reuse-contexts` is permitted only when the complete-context replay inputs and output hashes verify. Frozen release files retain their original policy and code bytes. Current selection and audit exports are checksum-pinned separately from presentation artifacts.
