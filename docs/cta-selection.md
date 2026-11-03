# Current CTA and oncofetal selection

The current policy supplies **625 CTA proteins, 18,351 Unsafe and 628 Other (550 Nonreproductive, 27 Held, 28 Nested, 19 Contained, 4 Diffuse)** across 19,604 proteins in the coding or established-ncORF universe. No protein is unassessed. These are research-policy assignments, not demonstrated clinical safety. The [current report](audits/cta-source-walkthrough-20261028/index.md) includes sources, figures, the funnel and tissue-specific explanations. The [1.21.0 report](audits/cta-release-1.21.0/index.md) remains frozen.

## Safety, cancer expression and Reproductive evidence

**Only safety can establish Unsafe.** Cancer expression and Reproductive evidence
are separate admission requirements. Every protein receives CTA, Unsafe or Other; missing coverage is disclosed in `assessment_status`, never as a separate class.

**Reproductive** is the umbrella for adult reproductive, placental, extraembryonic,
embryonic and fetal contexts. Branches describe the actual tissue and stage;
“embryonic-only evidence” does not prove absence from unmeasured branches.

| Assessment | Question | Current rule |
|---|---|---|
| Safety | Is normal-tissue expression concerning under the tissue/IHC policy? | Absolute RNA ceilings, separate positive-protein rules and reviewed clinical exclusions. Vital/Protected exclusions and Manageable excess above 128 establish Unsafe; lower Manageable excess blocks CTA and lands in Other:Held. |
| Cancer expression | Is there expression in at least one measured cancer? | One sample >0 TPM or HPA pTPM, or reviewed reliable tumor protein. No tumor-versus-normal ratio. |
| Reproductive evidence | Is there positive expression in a qualifying branch? | Full HPA reproductive bulk maximum ≥5 nTPM OR supported reproductive cell RNA ≥5 nTPM with ≥20 reads OR reliable positive reproductive IHC or reviewed reproductive protein OR supported embryonic lineage mean ≥5 TPM. |
| Coverage | Is the comparison adequately observed? | Complete combined bulk (HPA consensus, HPA brain and GTEx v11) for every member locus. Cell-cluster, immune, DICE and cornea panels count when present but never block a label. |

The former maximum-versus-maximum restriction rule is removed. The
[continuous diffuseness gate](audits/cta-unassigned-20261029/diffuseness-impact.md)
uses **HPA v23 consensus**: average `N / (N + R)` over **29 equally weighted
anatomical units / 38 contexts**. R is the highest adult/placental Reproductive
nTPM in the same matrix. Disposable tissues count; **Permissive thymus** belongs
to neither the comparator nor the Reproductive reference. Its safety tier is unchanged.

A score **above 10%** blocks a qualifying adult/placental route when R is positive
and the full HPA comparison is measured. RNA-only admission still requires **5 nTPM**;
reliable positive IHC can establish association below that floor. Missing comparisons
are disclosed in `assessment_status`, never as a separate class. An independently supported embryonic route remains
eligible; its concentration needs a compatible developmental/adult panel. GTEx
remains the required coverage core and a descriptive comparison, not a second veto.

A complete diffuse profile with no remaining qualifying route and clear normal
safety becomes **Other:Diffuse**. MTCP1 (24.53%) is one example. Diffuse identity alone never changes Unsafe membership. [Logic figure and examples](cta-evidence-axes.md).

## The selection funnel

Selection is a funnel. A protein leaves at the first stage it fails, and every later
stage sees only what remains.

| Stage | Question | Outcome |
|---|---|---|
| 0. Universe | Is the protein translated from a coding transcript or an established ncORF? | Otherwise it is **out of universe** and not classified. |
| 1. Safety | Does normal-tissue, reviewed or clinical evidence exclude it? | **Unsafe** |
| 2. Target space | Does it have reproductive evidence, cancer expression and concentration? | Otherwise **Other** with the unmet requirement as its subtype |
| 2a. Nested microprotein | Is it a reference-CDS protein under 50 residues inside a same-strand protein-coding gene, without ncORF evidence? | **Other:Nested**: its RNA belongs to the host gene |
| 3. Sequence | Is at least 95% of its 9-mers inside Unsafe proteins? | **Other:Contained**: Unsafe sequence is blacklisted, so almost no peptide survives |
| 4. Outcome | Did everything pass with no unresolved normal signal? | **CTA**; a pass blocked only by a normal-tissue hold is **Other:Held** |

**Stage 0 universe.** The universe is 19,604 protein groups: Ensembl 115 reference CDS
from a coding-biotype transcript, plus non-canonical ORFs that orfdata records with
mass-spectrometry peptides or at least two independent Ribo-seq studies. The 259
groups left out are 254 whose only reference CDS comes from NMD or CDS-not-defined
transcripts (mostly readthrough genes) and 5 supported by a single Ribo-seq study.
They are listed with their reason in `out-of-universe.csv.gz`, so none disappears
silently. Rebuild the gene table with `scripts/build_cta_protein_universe.py`.

**Downstream semantics.** CTA is targeted. Unsafe blacklists any matching subsequence
from a CTA. Other is ignored for targeting. Unsafe therefore needs positive evidence,
or an essential unknown, and never rests on missing optional panels.

| What puts a protein in Unsafe | Example |
|---|---|
| Vital/Protected expression above the tissue/IHC ceiling | Established adult protein or localized RNA with function |
| **Manageable-tier excess above 128** (the top Manageable ceiling) | Resolved Manageable compartments above the 128 ceiling |
| **Reviewed unresolved Protected-tissue finding** (an essential unknown) | Findings marked “resolve before shared-sequence clearance” |
| Reviewed clinical safety exclusion | MAGEA12 |

Below 128, Manageable excess and unresolved HPA protein evidence block CTA only;
the protein becomes Other and keeps its reasons. Clinical history stays
receptor-specific: the titin and MAGE-A12 cross-reactivities of anti-MAGE-A3 TCRs
are peptide liabilities, so MAGEA3 is assessed on its own expression.

## Cell-type tiers and the rare vital-tissue screen

The [central tissue-hazard ontology](tissue-hazards.md) supplies organ and cell calls
for RNA, IHC and the rare-cell screen. Its complete report provides the rationale,
evidence designs and limitations for every rule.

| Tier | Examples |
|---|---|
| Vital | Neurons, cardiomyocytes, AT1/AT2, club and basal airway cells, intestinal crypt progenitors |
| Protected | Ciliated/goblet airway cells, generic fibroblasts, differentiated bowel epithelium |
| Manageable | Resolved mature B cells, thyroid follicular cells |
| Disposable | Cutaneous melanocytes and hair-follicle compartments |

These are consequence-based policy calls, not measured probabilities of toxicity.
Bulk tissue retains its broader organ concern.

Pooled panels average away rare normal cells, so a protein can read near zero in
every bulk, cluster and immune measurement while a small population expresses it.
**MAGEA12** is the reference case: its HPA bulk, GTEx and single-nucleus cluster
values are near zero, but 293 of 19,760 midbrain dopaminergic neurons (1.5%)
express it in the CELLxGENE summaries, matching the rare positive neurons reported
after anti-MAGE-A3 TCR neurotoxicity.

The available genes are screened (`scripts/build_cta_rare_cell_screen.py`) against
CELLxGENE normal-sample summaries in brain, spinal cord, heart and selected lung
populations. A cell type needs at least 3,000 cells and 10 expressing cells. The
existing exclusions of broad parent terms, immune lineages, fibroblasts/stroma
and developmental labels are coverage limits, not evidence of tolerance.

At least 1% expression in a Vital cell type supplies Unsafe evidence. At least 1%
in a Protected population creates a candidate review hold. The complete funnel
and packaged profiles report the current counts; the central ontology extends
beyond the organs measured by this existing rare-cell panel.

It reports detection frequency among cells, not abundance: expressing cells typically
carry about one UMI, and ambient RNA, doublets and close-paralog multimapping can
create detections. Samples cannot be restricted to adults, so developmental cell types
are removed by name only. `cta_rare_cell.UNSAFE = False` turns the layer into a
review hold instead.

`cta_selection.profiles()` carries `rare_cell_screened`, `rare_cell_flag`,
`rare_cell_hold`, `rare_cell_max_pct` and `rare_cell_type`, and
`cta_selection.rare_cell_screen()` returns every supporting cell type with its `tier`
and `review_status = screen_unreviewed`.

```python
from oncoref import cta_selection

profiles = cta_selection.profiles()
rare_only = profiles[profiles.rare_cell_flag & profiles.category.eq("Unsafe")]
detail = cta_selection.rare_cell_screen()   # tissue, cell type, n_expressing, n_cells, tier
```

Every filter is evaluated for every protein, including those already Unsafe, so the
order of filtering changes only the funnel's attribution, never the final category.
The [filter-status table](audits/cta-selection-20261028/filter-status-summary.csv) lists each
protein against each filter.

## Contained: too similar to Unsafe

Unsafe proteins blacklist any matching subsequence from a CTA, so a protein whose
sequence lies almost wholly inside Unsafe proteins has no usable peptide. A would-be
CTA with **at least 95% of its 9-mers also present in an Unsafe protein** is assigned
**Other:Contained** (the boundary is inclusive; an undefined overlap, meaning no
9-mer, is never contained). It is a sequence rule, not expression evidence, and it
never creates Unsafe. The floor is deliberately strict: copies that keep a few
distinguishing 9-mers stay CTA, and their usable subsequences are for later review.

The overlap uses every protein-coding transcript of Ensembl 115
(`oncoref.cta_sequence_overlap`) and is reported as `unsafe_sequence_overlap_pct` in
`cta_selection.profiles()`. **20 proteins are Contained** (median overlap 99%), among
them DAZ2/3, PRAMEF27, CT45A3, GAGE2E and TBC1D3G/H/K. Most multi-copy families such as
the CT45, GAGE, PRAMEF and POTE copies stay CTA.

## Nested microproteins

A reference-CDS protein under 50 residues nested inside a same-strand protein-coding
gene, without ncORF evidence (two Ribo-seq studies or any mass-spectrometry peptide),
has RNA that cannot be told from its host's. Its expression says nothing about the
peptide, so it can never be a CTA: **28 proteins are Other:Nested** (48 such loci in the
universe table, tagged `nested_microprotein`).

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
`shared_sequence`: `allow` for CTA/Other and `exclude` for Unsafe. Every classified protein resolves; `review` remains only for sequences outside the universe or the nonhuman nominations.
These roles describe this research screen. They do not establish peptide presentation,
cover all isoforms or change the independently versioned legacy forbidden backgrounds.
Downstream systems must preserve the policy/reference receipt and handle review explicitly.

## One API for each species

```python
from oncoref import cta_selection

selected = cta_selection.results(species="human", category="CTA")
all_decisions = cta_selection.results(species=9606)
other = cta_selection.results(category="Other")
profiles = cta_selection.profiles()
labels = all_decisions[["label", "classification"]]
scope = cta_selection.coverage(species="dog")
sources = cta_selection.sources(species="mouse", category="developmental_expression")
```

Equivalent bound accessors are `human.cta_selection()`, `canine.cta_selection()` and `mouse.cta_selection()`, with `cta_selection_coverage()` and `cta_sources()` on each module. They also live on `oncoref.for_species("human").antigens`. Species names, aliases and taxonomy IDs use one resolver. Unknown species and categories fail explicitly; nothing silently substitutes human data.

| Species | Returned entities | Selection status | Meaning |
|---|---:|---|---|
| Human | 19,604 proteins (identical selected reference sequences collapsed; 259 out-of-universe groups listed separately) | Available; evidence partial | Complete evaluation of the coding or established-ncORF universe under `oncoref.cta-selection.2026-10-31.1` |
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
| `category`, `subtype`, `classification` | CTA / Unsafe / Other or null; mutually exclusive Other subtype; flat label such as Other:Nonreproductive or Other:Held; human results never have a null category |
| `assessment_status`, `policy_id` | Screened, incomplete evidence or normal-evidence hold (why an Other was not cleared); versioned policy |
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
| Required coverage | Combined bulk (HPA consensus, HPA brain, GTEx v11) complete for every member locus; cell, immune, DICE and cornea panels are optional and never block, but their measured signals count |
| Normal RNA | One inclusive tissue/IHC ceiling; above it supplies Unsafe evidence in Vital/Protected and above 128 in Manageable; lower Manageable excess blocks CTA |
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

Adult human fibroblasts and differentiated bowel epithelium use Protected; intestinal crypt/undifferentiated populations use Vital. Bulk organs retain their organ tiers. Sorted immune and corneal populations use resolved cellular calls from the same registry. Reproductive membership remains independent, and adult evidence does not establish pediatric safety.

## Categorical assignments

The decision follows these questions, in this order:

1. **Is the protein in the universe?** Otherwise it is not classified.
2. **Does normal-tissue or applicable clinical evidence trigger exclusion?** Assign **Unsafe**. This includes Manageable excess above 128 and reviewed unresolved Protected findings.
3. **Is there qualifying reproductive evidence?** Otherwise assign **Other:Nonreproductive**. Reproductive includes placenta, embryo and fetus; the label records an unmet evidence requirement, not demonstrated absence.
4. **Is there measured cancer expression?** Otherwise assign **Other:Reproductive**.
5. **Does reproductive concentration pass?** Otherwise assign **Other:Diffuse**.
6. **Is at least 95% of its sequence inside Unsafe proteins?** Then assign **Other:Contained**: its peptides are blacklisted.
7. **Is the normal screen complete and unresolved-free?** Passing everything assigns **CTA**. If target gates pass but a marginal or unresolved normal signal remains, assign **Other:Held**, and `assessment_status` records whether coverage was incomplete or a hold applied.

All safety layers evaluate every protein, even when an identity requirement fails.
The RNA/IHC safety ceilings and the reproductive/cancer presence requirements
are unchanged by this presentation. A high tumor value does not erase normal
counterevidence or establish reproductive identity.

**Other has one primary reason, not an abundance or confidence class.** The former
Low, Weak, Normal, Tumor-supported and Disposable categories mixed different
questions. Their useful evidence is now available as independent attributes:

```python
from oncoref import cta_selection

results = cta_selection.results()
profiles = cta_selection.profiles()
# Example columns: tumor_max_tpm, replicated_tumor_rna, normal_max_ntpm,
# adult_rna_below_floor, unconfirmed_normal_signal, adult_disposable_only,
# reproductive_evidence, developmental_measured
```

| Attribute | Meaning |
|---|---|
| Cancer maxima and replicated tumor RNA | Native TPM and HPA pTPM remain separate. At least two identified tumors ≥10 TPM is a support flag, not a CTA admission requirement. |
| Normal RNA maxima | Source-supported Vital/Protected/Manageable summaries in nTPM, TPM and CP10K, reported separately. Permitted-compartment expression is reported independently. A supported value can be below the safety ceiling and still represent real normal expression. |
| Adult RNA below the descriptive floor | All recorded adult raw RNA is below 2 native TPM/nTPM or the separately specified 0.02 CP10K floor. This is an abundance attribute, not evidence of complete absence. |
| Unconfirmed normal signal | Raw adult RNA reaches that descriptive floor while the source-supported summary does not. The exact read/cell/recurrence evidence remains in the audit. |
| Adult Disposable expression | An appreciable RNA or reliable protein observation in an adult compartment designated Disposable. |
| Adult Disposable-only pattern | Positive Disposable evidence, other measured adult RNA below the descriptive floor, no detected protein elsewhere in the adult panel, and complete required normal coverage. |
| Reproductive and developmental evidence | Reported separately from adult tissue disposition, with native units, measured stages and missing coverage. |

**Disposable designates an adult tissue compartment permitted by the research
policy.** It never requires low or measured embryonic/fetal/placental expression.
A protein can have an adult Disposable-only pattern and also be expressed during
development. These are separate observations. The pattern does not override an
Unsafe decision or unresolved normal evidence. Whole skin is Manageable; mapped
hair follicles and cutaneous melanocytes have conditional Disposable allowances.

The 2-unit descriptive floor answers whether a signal is appreciable for reporting;
it is not a second safety cutoff. Missing optional assays remain missing. Sparse
RNA, protein evidence and clinical injury retain distinct strengths. Native human
thresholds and findings do not transfer to mouse or dog through orthology.

[Simple funnel, exact counts, four example proteins and migration](audits/cta-simplification-20261031/index.md).

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
