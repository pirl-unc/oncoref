# Current CTA and oncofetal selection

The current policy supplies **701 CTA proteins, 18,062 Unsafe, 454 Other and 646 unassessed** across 19,863 proteins. These are research-policy assignments, not demonstrated clinical safety. The [current report](audits/cta-source-walkthrough-20261028/index.md) includes sources, figures, the funnel and tissue-specific explanations. The [1.21.0 report](audits/cta-release-1.21.0/index.md) remains frozen.

## Safety, cancer expression and Reproductive evidence

**Only safety can establish Unsafe.** Cancer expression and Reproductive evidence
are separate admission requirements. Missing required coverage remains unassigned.

**Reproductive** is the umbrella for adult reproductive, placental, extraembryonic,
embryonic and fetal contexts. Branches describe the actual tissue and stage;
“embryonic-only evidence” does not prove absence from unmeasured branches.

| Assessment | Question | Current rule |
|---|---|---|
| Safety | Is normal-tissue expression concerning under the tissue/IHC policy? | Absolute RNA ceilings, separate positive-protein rules and reviewed clinical exclusions. Vital/Protected exclusions can establish Unsafe. Manageable excess blocks CTA. |
| Cancer expression | Is there expression in at least one measured cancer? | One sample >0 TPM or HPA pTPM, or reviewed reliable tumor protein. No tumor-versus-normal ratio. |
| Reproductive evidence | Is there positive expression in a qualifying branch? | Full HPA reproductive bulk maximum ≥5 nTPM OR supported reproductive cell RNA ≥5 nTPM with ≥20 reads OR reliable positive reproductive IHC or reviewed reproductive protein OR supported embryonic lineage mean ≥5 TPM. |
| Coverage | Is the comparison adequately observed? | Complete GTEx v11 core: 52 bulk contexts for every member locus; required normal bulk/cell/immune panels also apply. DICE/cornea coverage remains optional. |

The former maximum-versus-maximum restriction rule is removed. The
[continuous diffuseness gate](audits/cta-unassigned-20261029/diffuseness-impact.md)
uses **HPA v23 consensus**: average `N / (N + R)` over **29 equally weighted
anatomical units / 38 contexts**. R is the highest adult/placental Reproductive
nTPM in the same matrix. Disposable tissues count; **Permissive thymus** belongs
to neither the comparator nor the Reproductive reference. Its safety tier is unchanged.

A score **above 10%** blocks a qualifying adult/placental route when R is positive
and the full HPA comparison is measured. RNA-only admission still requires **5 nTPM**;
reliable positive IHC can establish association below that floor. Missing comparisons
remain explicitly unassessed. An independently supported embryonic route remains
eligible; its concentration needs a compatible developmental/adult panel. GTEx
remains the required coverage core and a descriptive comparison, not a second veto.

A complete diffuse profile with no remaining qualifying route and clear normal
safety becomes **Other:Diffuse**. MTCP1 (24.53%) is one example. Diffuse identity alone never changes Unsafe membership. [Logic figure and examples](cta-evidence-axes.md).

## Reviewed normal-cell knowledge

The [600-protein literature audit](audits/cta-other-literature-20261030/index.md) adds
source-linked findings through `gene_knowledge.findings(gene)`. Established adult human
normal protein, or localized RNA with orthogonal normal-cell function, in Vital/Protected
can establish Unsafe. Rare cells and activation states count. Sparse RNA and unresolved
species, stage or isoform evidence receive a review hold, with reconsideration requirements.
A search hit or annotation alone is not accepted protein evidence. Normal findings take
precedence over reproductive protein or cancer presence.

Cancer RNA now needs only **>0** in one compatible measured sample. This admits CDY1/B
and OR2T7 without claiming high tumor abundance. Reviewed germ-cell protein supports
BPY2/B/C independently of its RNA floor. All 698 prior CTAs remain retained.

`cta_selection.sequence_dispositions()` exposes `direct_target_candidate` and
`shared_sequence`: `allow` for CTA/Other, `exclude` for Unsafe, `review` for Unassessed.
These roles describe this research screen. They do not establish peptide presentation,
cover all isoforms or change the independently versioned legacy forbidden backgrounds.
Downstream systems must preserve the policy/reference receipt and handle review explicitly.

## One API for each species

```python
from oncoref import cta_selection

selected = cta_selection.results(species="human", category="CTA")
all_decisions = cta_selection.results(species=9606)
other_low = cta_selection.results(category="Other", subtype="Low")
labels = all_decisions[["label", "classification"]]
scope = cta_selection.coverage(species="dog")
sources = cta_selection.sources(species="mouse", category="developmental_expression")
```

Equivalent bound accessors are `human.cta_selection()`, `canine.cta_selection()` and `mouse.cta_selection()`, with `cta_selection_coverage()` and `cta_sources()` on each module. They also live on `oncoref.for_species("human").antigens`. Species names, aliases and taxonomy IDs use one resolver. Unknown species and categories fail explicitly; nothing silently substitutes human data.

| Species | Returned entities | Selection status | Meaning |
|---|---:|---|---|
| Human | 19,863 proteins (identical selected reference sequences collapsed) | Available; evidence partial | Complete evaluation of the named reference universe under `oncoref.cta-selection.2026-10-30.3` |
| Dog | 14 curated native nominations | Unavailable; evidence partial | Unassigned / unassessed; these are assay or family nominations, not 14 proven CTA genes |
| Mouse | 1 P1A/Trap1a nomination | Unavailable; evidence partial | Unassigned / unassessed; qualitative normal and tumor RNA, with strain and reference limits |

The contract is identical even when the evidence is not. Human thresholds, antibody grades, developmental measurements and clinical observations are never inherited through orthology. Canine matrices still need unit/reference reconciliation and a species policy; the registered tumor source has unresolved unit/reference provenance. Mouse P1A literature evidence is imported; the listed mouse atlases remain unintegrated.

## Result fields and identities

`results()` returns a detached pandas DataFrame. `taxon_id` is integer; scalar text uses pandas string dtype; the four collection fields are tuples. Empty results have the same columns and types. Reads do not download source matrices.

| Field | Meaning |
|---|---|
| `taxon_id`, `species` | Explicit organism boundary |
| `entity_id`, `entity_kind`, `label` | Stable result key; human protein-equivalence group versus native nomination; display label |
| `gene_ids`, `gene_identities` | Every supported member/alternative; IDs, symbols, namespaces, historical/alternate IDs, synonyms and identity status |
| `identity_status`, `reference_id` | Resolution completeness and exact reference receipt |
| `category`, `subtype`, `classification` | CTA / Unsafe / Other or null; mutually exclusive Other subtype; flat label such as Other:Low or Unassessed |
| `assessment_status`, `policy_id` | Screened, incomplete evidence, normal-evidence hold or unassessed; versioned policy |
| `reasons`, `evidence_flags` | All contributing exclusions and retained observations, including Low protein, antibody review and stem/progenitor expression |
| `source_snapshot_id` | Content identity of the decision snapshot |

Unknown history is `None`, not an empty list pretending the search was exhaustive. Ambiguous canine identities remain alternatives and do not become protein-equivalence groups. Human reference identity pins the Ensembl-115 primary gene-space projection; it does not claim exhaustive sequence coverage. A representative-protein equivalence group is also not proof that every alternative isoform is identical. The complete locus crosswalk retains identities and explicitly labels earlier decisions as historical/baseline.

TPM/nTPM is summed only for matched group-member measurements in the same source context. Missing members do not become zeros. Unpaired DICE distributions and corneal detection bounds retain uncertainty. A single measured tumor member above threshold establishes a group lower bound; maxima from different samples are not added.

## Uniform human rules

| Gate | Rule |
|---|---|
| Adult/placental Reproductive evidence | Full mapped HPA bulk maximum ≥5 nTPM OR supported reproductive cell RNA ≥5 nTPM with ≥20 reads OR reliable positive reproductive IHC or reviewed reproductive protein; no normal maximum comparison |
| Embryonic Reproductive evidence | Supported Tyser or Petropoulos lineage mean ≥5 TPM; ≥20 cells, ≥5 positive cells and ≥10% positive at ≥1 TPM |
| Tumor expression | >0 in one measured sample, in native TPM or HPA pTPM, or reviewed reliable tumor protein; no normal comparison |
| Required coverage | GTEx v11 core (52 native bulk contexts) and specified bulk/cell/immune panels must be complete; DICE and cornea coverage is optional, but their measured signals count |
| Normal RNA | One inclusive tissue/IHC ceiling; above it supplies Unsafe evidence in Vital/Protected and blocks CTA in Manageable |
| Normal protein | Low remains a positive flag; stronger HPA or accepted unquantified literature positives retain their hold/exclusion; validated external normal positives and clinical exclusions remain active |

| Candidate ceiling | Missing / unusable IHC | Approved negative | Supported negative | Enhanced negative |
|---|---:|---:|---:|---:|
| Vital | 4 | 8 | 16 | 32 |
| Protected | 8 | 16 | 32 | 64 |
| Manageable | 16 | 32 | 64 | 128 |

A reviewed TCR-T cohort with no identified grade 4-5 specificity-related toxicity
can supply **one additional RNA step, capped at Enhanced**. The resulting ladders
are Vital **8/16/32/32**, Protected **16/32/64/64**, and Manageable **32/64/128/128**.
Trials do not stack or change the observed IHC grade. A severe suspected specificity
signal, ungraded suspected specificity injury or unassessed event attribution blocks
credit for that receptor. Protein and clinical-injury exclusions still apply.

This is a protein-profile screening inference with explicit peptide, isoform, HLA,
receptor and population limitations. It does not clear every epitope of a protein.
**CTAG2/LAGE-1A and NY-ESO-1 share the c259 target peptide**; CTAG2 already passed
before the credit. Its Enhanced-negative heart ceiling is capped and stays at 32.
[Serious adverse events, attribution, clinical cohorts and exact epitope mappings](audits/cta-tcr-trials-20261029/index.md).

`oncoref.tcr_trials` exposes the ledger, cohorts, sentinel events and clinical-credit
assessment. `oncoref tcr-trials --table all` exports it as JSON. Original source
artifacts use the shared versioned reference cache with checksums. Unreviewed or
unreported outcomes never become evidence of no toxicity.

These are policy parameters, not measured toxicity limits. TPM/nTPM assays share nominal limits while retaining their units; corneal CP10K uses nominal limits divided by 100, **not a TPM conversion**. HPA clusters require ≥20 reads; DICE requires three reported values; cornea requires ≥5 positive cells and ≥20 counts in at least two sample groups. The source-specific support statistics are preserved rather than calling every measurement a comparable sample maximum.

Exact adult tissue/cell matching or a declared source-specific histological-compartment correspondence, reliable negative staining, a positive profile control and no unresolved antibody conflict are required for IHC credit. Every expected member needs a positive profile control for a group negative. Low positive staining grants no negative-IHC credit. Author-described Low and HPA Low retain different scales, with unknown fractions/H-scores left unknown.

Mapped HPA adult fibroblasts use Manageable; bulk heart and cardiomyocytes retain their organ tiers. Colon/rectum distal enterocytes, goblet cells, Paneth cells and undifferentiated cells use Manageable in both RNA and HPA protein assessment. Undifferentiated expression retains a stem/progenitor flag. This permissive mucosal assumption and the numeric cutoffs have not been clinically validated. Reproductive/placental immune cells receive neutral anatomical treatment, not an assertion that all immune cells are disposable. Adult normal panels do not establish pediatric safety.

## Categorical assignments

Every layer evaluates every protein. **Unsafe** takes precedence when Vital/Protected or clinical exclusion evidence is present. **CTA** requires all safety, coverage, cancer and Reproductive identity gates. **Other** requires complete normal coverage, resolved evidence and no normal-screen hold, but fails CTA identity. High, low or missing tumor RNA does not alone prevent Other.

Each Other protein receives exactly one subtype, in the following priority order. The flat `classification` field is the categorical result. Secondary evidence remains in the audit; there is no unnamed Other remainder.

| Category | Proteins | Definition |
|---|---:|---|
| Other:Diffuse | 6 | Complete HPA profile is dispersed (score >10% with a reference ≥2 nTPM), or the active concentration gate fails; no CTA route remains. |
| Other:Reproductive | 0 | Qualifying Reproductive evidence, but no qualifying cancer expression; the normal screen passes. |
| Other:Inducible | 0 | Cited expression in activated nonreproductive normal cells; a localized indirect atlas signal cannot establish Reproductive identity. |
| Other:Tumor-supported | 46 | Two distinct tumors ≥10 native TPM or reviewed reliable tumor protein, with low measured normal/developmental expression; not established tumor restriction. |
| Other:Disposable | 1 | Positive Disposable RNA/protein; measured RNA elsewhere below the descriptive 2-unit floor and no detected protein elsewhere, including low measured developmental RNA. |
| Other:Low | 210 | All measured normal RNA below the descriptive 2-unit floor and no detected normal protein. Cancer expression need not be low. |
| Other:Normal | 124 | Supported nonreproductive normal RNA reaches the descriptive 2-unit floor or reliable protein is detected; no preceding subtype or CTA definition qualifies. |
| Other:Weak | 67 | Reproductive signals are below admission requirements or lack support; no preceding subtype qualifies and no supported appreciable nonreproductive expression is established. |

These descriptive floors mean 2 native TPM/nTPM, or the separate nominal 0.02 CP10K floor. They are not safety ceilings and the units are not converted or pooled. Low does not mean never expressed. The Weak category records insufficiently supported signals, including high raw cluster values with too few reads. Its normal coverage still passes. All current Other proteins have supported non-permitted normal maxima below 30 TPM/nTPM. Raw values that fail source-specific read/cell support remain in the audit. The three baseline Other proteins with raw non-permitted values >30 TPM/nTPM now have explicit review holds.

Reliable **positive reproductive IHC** is an independent association route. It requires a mapped human HPA v23 reproductive compartment, Low/Medium/High staining with Approved/Supported/Enhanced reliability, positive evidence for every member locus and no unresolved antibody concern.

**Reproductive cell RNA** uses the same pinned HPA v23 atlas, source-specific tissue/cell labels, at least 5 nTPM and 20 reads. Every same-protein member must be measured before summing a cluster. Generic immune/endothelial or mixed cells can provide indirect sampling-site support when the outside-reproductive cluster panel is complete and every outside cluster is below 5 nTPM. This establishes association with the sampled organ, not the producing cell. [Macrophage passenger transcripts](https://pmc.ncbi.nlm.nih.gov/articles/PMC7459098/) motivate that distinction. Published normal induction invalidates an indirect-only localization claim: IFNA2/14 are now Unsafe after the normal-protein review, extending [human plasmacytoid dendritic-cell expression](https://pmc.ncbi.nlm.nih.gov/articles/PMC6457222/). Curated counterevidence is identified by gene, species, context and source; it is not copied to unreviewed family members.

No presence route overrides normal ceilings, cancer requirements or a measured concentration failure. A missing/zero reproductive bulk reference leaves RNA concentration explicitly unassessed; cell RNA or IHC does not invent a bulk abundance.

**Other:Tumor-supported** identifies leads with two distinct tumors ≥10 native TPM or reviewed reliable tumor protein, low measured normal/developmental RNA and no detected normal protein. Source views and TCGA participant/aliquot duplicates do not count as replication. Aggregate cancer-IHC counts lack the reagent validation required for protein admission. The 46 current leads await assessment of relevant rare normal cells and physiological states. They are not declared tumor-exclusive or added to `candidates()`.

Missing required normal coverage, unresolved evidence or an outstanding normal-screen hold yields **Unassessed**, with a null biological category and an explanatory status. Manageable excess cannot be rescued into Other by expression in a second Manageable tissue. `cta_selection.candidates()` returns CTA only. The deprecated `category="Neither"` query returns Other with a warning. Final classifications do not depend on layer order.

## Evidence categories across species

`sources()` returns one row per source/category with native unit, scope, URL and integration status. [Paper titles and database versions](audits/cta-source-walkthrough-20261028/references.md) provide readable links to explanatory pages, with artifact URLs retained separately. A source can serve several roles; source-category rows must not be summed as independent studies or donors.

| Category | Human selection | Dog resources | Mouse resources identified |
|---|---|---|---|
| Nomination | Published CTA/placental lists as annotations | Two native CTA studies | One P1A/Trap1a literature nomination |
| Reproductive expression | HPA bulk, supported tissue-cell RNA and reliable reproductive IHC | Bgee matrices, BarkBase, native study assays | Qualitative P1A testis RNA; no complete reproductive panel |
| Developmental expression | Tyser gastrulation | BarkBase matrix; DoGA metadata | MOCA, gastrulation atlas, GXD |
| Tumor expression | 132 compatible OncoRef cohorts + HPA cancer | Native studies; cohort provenance unresolved | Qualitative J558 cell-line RNA; no cohort |
| Normal bulk RNA | HPA consensus | Bgee, BarkBase; DoGA/Pfizer RNA metadata | Five BALB/c organs for P1A; GXD remains unintegrated |
| Normal cell RNA | HPA tissue clusters, adult cornea | No equivalent integrated atlas | Tabula Muris/Senis |
| Immune RNA | HPA, Monaco, DICE | No equivalent integrated panel | Tabula Muris/Senis |
| Normal protein | HPA IHC, curated external IHC/protein | Pfizer iBAQ matrix, native IHC study | GXD assay resource |
| Clinical safety | Curated product-specific findings | No transferred human evidence | No transferred human evidence |

`coverage()` distinguishes sources used in selection, available observations/matrices, and resources merely identified. DoGA's 5′ expression assay, BarkBase TMM-CPM, Pfizer iBAQ and mouse UMIs cannot be silently labeled TPM. Source inventories and species limitations are checksum-pinned in the package.

## Versioning and reproduction

Use this API for the new funnel. Older `cta_genes`, nomination snapshots and exact-protein forbidden-background APIs remain independently versioned compatibility/reference products. Their counts must not be presented as this selection; a CTA admission does not automatically remove a protein from a near-self background. See [the exact-protein policy](target-sets.md).

The older `normal_expression_categories.classify_normal_expression()` API still
provides **Other → Disposable/Silent** under
`oncoref.normal-expression-categories.2026-10-21.2`. Those are below-threshold
expression sets using that older RNA/IHC policy. They are not aliases for this
release's Other categories, which now includes qualifying permitted-normal profiles through its normal-pass fallback.
Downstream analyses must preserve the policy ID and must not silently transfer
these classifications into the current selection or k-mer exclusions.
The migration is tracked in [#647](https://github.com/pirl-unc/oncoref/issues/647).

Run from a checkout with the pinned reference inputs available:

```sh
python scripts/regenerate_cta_walkthrough.py
python scripts/audit_cta_profile_concentration.py
node scripts/build_cta_source_slides.mjs
python scripts/regenerate_cta_walkthrough.py --finalize-receipts
```

The slide script requires the runtime paths documented in its header. `--reuse-contexts` is permitted only when the complete-context replay inputs and output hashes verify. Frozen release files retain their original policy and code bytes. Current selection and audit exports are checksum-pinned separately from presentation artifacts.
