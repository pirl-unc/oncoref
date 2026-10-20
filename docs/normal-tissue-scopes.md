# Normal tissue scopes

oncoref uses one cited anatomical vocabulary for human, dog and mouse. A **scope
is a set of permitted normal anatomical/developmental contexts**. A detected
protein or gene is compatible with that scope only if **every detection** belongs
to it. Tissue sets are unioned before evaluating expression patterns. Positive
observations alone do not establish absence outside a scope.

Definitions and citations ship in `oncoref/data/normal-tissue-taxonomy.yaml`;
exact dataset mappings ship in `oncoref/data/normal-tissue-mappings.csv`.
Policy version `2026-10-20.3` and a SHA-256 identity over both files accompany
each annotation. Source-artifact receipts record the bytes used to review labels.
All access is offline once expression data are loaded.

## Read the hierarchy

Each row below shows inclusion into successively larger tissue sets. Sibling
branches are combined at their shared parent.

| Starting context | Successive unions |
|---|---|
| Testis | `testis` → `male_reproductive` → `reproductive` → `reproductive_or_disposable` |
| Oocyte | `oocyte` → `ovary` → `female_reproductive` → `reproductive` → `reproductive_or_disposable` |
| Mammary glandular cells | `mammary_glandular` → `mammary` → `female_reproductive` → `reproductive` → `reproductive_or_disposable` |
| Placenta | `placental` → `extraembryonic` → `reproductive` → `reproductive_or_disposable` |
| Explicit embryo/fetus context | `embryonic` or `fetal` → `reproductive` → `reproductive_or_disposable` |
| Cutaneous melanocytes or hair follicles | `disposable_normal` → `reproductive_or_disposable` |
| Whole skin, salivary glands, thymus, tonsil, spleen, appendix or gallbladder | `permitted_normal` → `reproductive_or_permitted_normal` |

`reproductive_or_disposable` is a subset of the expanded
`reproductive_or_permitted_normal` union. `disposable_normal` is a subset of
`permitted_normal`; the wider name records a project screening allowance.
`acceptable_normal` and `reproductive_or_acceptable_normal` are readable aliases
for those existing scope identifiers. They resolve to the same sets.

For a gene detected in both testis and oocyte, **the entire observed pattern**
fits `reproductive_core` and `reproductive`, while neither sex-specific scope
contains both detections. The tissue union is evaluated before selecting genes.
The aliases `male_strict` and `female_strict` name anatomical scopes; a claim of
strict gene restriction additionally needs adequate negative tissue coverage.

The [original CTA source and filtering report](audits/cta-tissue-scopes-20261019/index.md)
compares reproductive with reproductive + disposable using both the existing
adaptive CTA policy and a strict positive-evidence sensitivity screen. It includes
figures, exact gene decisions, source citations and coverage limitations.
The [expanded permitted-tissue report](audits/cta-permitted-tissues-20261019/index.md)
adds the requested organs, vital/protected RNA ceilings, and an independent IHC
exclusion gate, with included/excluded decisions for published targets.
The [IHC-adaptive absolute RNA report](audits/cta-ihc-adaptive-20261020/index.md)
compares those fixed ceilings with the approved tissue-specific negative-IHC ladder.
These dated reports preserve their original policies and results. The
[tissue-tier update](audits/cta-tissue-tiers-20261020/index.md) compares the current
anatomy policy with the preceding adaptive screen, holding RNA thresholds fixed.
The [context and removable-organ update](audits/cta-tissue-contexts-20261020/index.md)
adds appendix/gallbladder and records four tissue categories and manageability reviews.

## Reproductive scopes

The following are source contexts, not assertions of immune privilege. Anatomical
subcompartments inherit their organ scope. For example, an oocyte belongs to both
the oocyte and ovarian scopes, but an observation of bulk ovary does not identify
oocytes.

| Scope | Exact definition |
|---|---|
| `testis` (`male_strict`) | Testis, including explicitly mapped seminiferous tubules, germ cells, Sertoli, Leydig and peritubular cells. Whole testis does not imply germ-cell-only expression. |
| `oocyte` (`female_strict`) | Oocytes/ova. Whole ovary, ovarian stroma and follicle-cell observations do not qualify for this narrower scope. |
| `ovary` | Whole ovary and mapped ovarian compartments: oocytes, follicles, granulosa, theca, stroma and corpus luteum. |
| `mammary_glandular` | Explicit mammary glandular epithelium, including luminal and myoepithelial cells. |
| `mammary` | Whole breast/mammary tissue and its subcompartments. Bulk tissue remains broader than glandular-cell evidence; lactating state is retained when reported. |
| `male_reproductive` (`male_extended`) | `testis` ∪ spermatozoa ∪ epididymis ∪ vas deferens ∪ prostate ∪ seminal vesicle ∪ bulbourethral gland ∪ penis, filtered for species applicability. Spermatozoa alone do not imply a testis sampling site. Mouse coagulating gland is included under prostate. |
| `female_reproductive` (`female_extended`) | `ovary` ∪ fallopian tube/oviduct ∪ uterus (including endometrium, myometrium, cervix and applicable uterine horns) ∪ vagina ∪ vulva ∪ clitoris ∪ `mammary`. The grouping does not infer the sample's sex. |
| `placental` | Placenta and its annotated subcompartments. A placenta-specific observation qualifies without resolving maternal versus fetal compartment. |
| `extraembryonic` | `placental` ∪ amnion ∪ chorion ∪ yolk sac ∪ allantois ∪ umbilical cord. |
| `embryonic` | Whole embryo or an anatomical sample explicitly labeled embryonic. Original organ, stage and age metadata remain available. |
| `fetal` | Whole fetus or an anatomical sample explicitly labeled fetal/foetal. Original organ, stage and age metadata remain available. |
| `reproductive_core` (`core`) | `testis` ∪ `ovary` ∪ `placental`, retaining the historical core organ union. |
| `reproductive` (`extended`, `permissive`) | `male_reproductive` ∪ `female_reproductive` ∪ `extraembryonic` ∪ `embryonic` ∪ `fetal`. |

The organ groupings follow the [male](https://openstax.org/books/anatomy-and-physiology/pages/27-1-anatomy-and-physiology-of-the-male-reproductive-system)
and [female reproductive anatomy](https://openstax.org/books/anatomy-and-physiology/pages/27-2-anatomy-and-physiology-of-the-female-reproductive-system),
which includes mammary glands as accessory reproductive organs. The broader
embryonic/fetal inclusion is an explicit **project definition**, supported by
retaining the source distinctions used in [BarkBase](https://doi.org/10.3390/genes10060433)
and the [human fetal atlas](https://pubmed.ncbi.nlm.nih.gov/33184181/).

An embryonic liver sample can belong to `embryonic` and `reproductive`; adult
liver does not. Gestational-day cutoffs are not copied between species, and
pregnancy of the donor does not turn maternal organs into fetal tissue. A source's
embryonic/fetal terminology is preserved rather than reinterpreted from age alone.

## Conditional disposable normal scope

`disposable_normal` contains **cutaneous melanocytes and hair-follicle
compartments**. `reproductive_or_disposable` is its union with `reproductive`.
These are research policy names, not declarations that destruction is harmless
or acceptable for every patient or treatment.

| Included compartment | Source and precise boundary |
|---|---|
| Cutaneous melanocytes | The explicit example list in [Blagosklonny 2003](https://doi.org/10.1038/sj.bjc.6601256) names melanocytes, prostate, testis, breast and ovary. The four reproductive organs are assigned above. Here melanocytes are limited to skin: clinical antigen-directed T-cell therapy also caused [ocular and auditory toxicity](https://pmc.ncbi.nlm.nih.gov/articles/PMC2929689/), so an unsited melanocyte label, ocular melanocytes, inner-ear cells and retinal pigment epithelium do not inherit this exception. |
| Hair follicles | Included as the separately requested conditional alopecia exception. [NCI](https://www.cancer.gov/about-cancer/treatment/side-effects/hair-loss) describes treatment-associated alopecia and variable regrowth. [Purba et al. 2019](https://pubmed.ncbi.nlm.nih.gov/31512803/) distinguishes follicular matrix and stem-cell injury. A bulk follicle observation does not establish matrix-only expression or reversible injury. This evidence does not validate immune ablation of follicles. |

Whole skin, epidermis, generic keratinocytes, sebaceous/sweat glands, fibroblasts,
retina, marrow, thyroid, lymphoid organs and unsited immune-cell labels do not
automatically qualify. Thymus is a somatic tissue, including when expression
reflects immune tolerance. No unlisted normal tissue is admitted by default.

## Vital, protected and acceptable tissues

The four user-facing categories are **Reproductive, Disposable, Protected and
Vital**. `normal_tissue_category()` returns their lowercase IDs.
Reproductive and Disposable share the existing numeric screening tier `permitted`.
Here **Disposable means the expanded conditional non-reproductive allowance**:
its canonical scope remains `permitted_normal`. The earlier `disposable_normal`
scope keeps its narrow cutaneous-melanocyte/hair-follicle meaning for compatibility
and historical comparisons. Do not confuse that scope ID with the four-way category.

There are three resolved numerical screening tiers. **Acceptable means conditionally
allowed in this research screen**, with consequences and evidence limitations
recorded separately. **Vital / critical function** covers both potentially fatal
injury and major irreversible neurological or sensory disability. **Protected**
covers other functionally important tissues without a project allowance.
Unknown or mixed anatomy is **unresolved**, a QC status outside the four resolved categories.

The `permitted_normal` scope is `disposable_normal` ∪ whole skin ∪ salivary glands
∪ thymus ∪ tonsil ∪ spleen ∪ appendix ∪ gallbladder. `reproductive_or_permitted_normal` unions this with
`reproductive`. Skin permission includes explicitly skin-sited keratinocytes,
glands, fibroblasts, immune cells and other mapped skin compartments. It does not
extend the exception to those cell types elsewhere, unsited cells, eye, retina or
inner ear. Parotid inherits salivary gland. Thymic cortex/medulla inherit thymus.
Anatomical mappings and organ labels remain available at their original resolution.

### Full assignment

| Tier | Complete root list (mapped descendants inherit) |
|---|---|
| **Reproductive** | All reproductive contexts defined above, including mammary, placental, embryonic and fetal contexts. |
| **Disposable (conditional)** | Cutaneous melanocytes; hair follicles; whole skin; salivary glands; thymus; tonsil; spleen; human appendix; gallbladder. Mapped descendants include nail bed and parotid. |
| **Vital / critical function** | Heart; lung; bronchus; brain; spinal cord; liver; kidney; bone marrow; blood; adrenal gland; pituitary gland; pancreas; small intestine; skeletal muscle; smooth muscle; peripheral nerve; vascular endothelium; retina; eye; inner ear. |
| **Protected** | Adipose tissue; thyroid gland; parathyroid gland; lymph node; esophagus; stomach; colon; rectum; urinary bladder; ureter; cartilage; oral mucosa; nasopharynx; tongue. |
| **Unresolved** | Head; foot; generic soft tissue; unrecognized, conflicting or unsited anatomy. |

Retinal pigment epithelium inherits retina. The eye label is conservative because
a bulk eye sample can include critical structures. The full vital label replaces
an interpretation based only on mortality; protected tissue injury can also be
serious. These are project assignments, not a clinical ranking of every possible
pattern of organ injury.

Appendix and gallbladder move out of Protected because their removal is compatible
with normal life ([NHS](https://www.nhs.uk/conditions/appendicitis/),
[NIDDK](https://www.niddk.nih.gov/health-information/digestive-diseases/gallstones/treatment)).
In-situ immune injury still needs review; adjacent intestine, bile ducts and liver
do not inherit the allowance.

The [manageability evidence review](normal-tissue-manageability.md) anchors the
ongoing policy discussion in published tissue classifications and clinical
immune-toxicity observations. Existing skin, salivary and lymphoid allowances
remain explicit project choices. Local-steroid control is recorded independently;
it does not automatically relax RNA or IHC gates.

### Lymphoid organs and conditional acceptance

| Source context | Assignment and explicit review conditions |
|---|---|
| **Tonsil** | Acceptable with immune-cell distribution and age/immune-reserve review. Palatine tonsillectomy is the human rationale; an unspecified dataset `tonsil` label retains that resolution. This does not map adenoids, all mucosal lymphoid tissue or mouse NALT into tonsil. |
| **Spleen** | Acceptable with immune-cell distribution, age/immune-reserve and asplenia infection-risk review. Surviving organ removal does not imply low-consequence systemic injury. |
| **Thymus** | Retains the requested acceptable assignment, with immune-cell distribution and age/immune-reserve review. |
| **Lymph node** | Protected, with lymphatic/immune-function review. Regional removal does not justify allowing injury across all nodes or lymphatic stroma. |
| **Blood / bone marrow** | Vital. No general immune-lineage exception is introduced. Unsited B-cell and T-cell source labels remain unresolved. |

The [NHS tonsillectomy decision aid](https://www.england.nhs.uk/wp-content/uploads/2023/10/decision-support-tool-making-a-decision-about-recurrent-tonsillitis-in-children-and-adults.pdf)
documents clinical removal and overlapping immune functions.
[CDC asplenia guidance](https://www.cdc.gov/meningococcal/risk-factors/index.html)
documents increased infection risk.
[NCI lymph-node guidance](https://www.cancer.gov/about-cancer/diagnosis-staging/staging/sentinel-node-biopsy-fact-sheet)
documents immune-cell traffic and increasing lymphedema risk with more extensive
removal. Our inference is that these support **conditional organ allowances and
cell-distribution review**, not equivalence between surgery and antigen-directed
systemic injury.
[NEI retinal guidance](https://www.nei.nih.gov/sites/default/files/nei-pdfs/nationalplan1.pdf)
and the [melanocyte-directed TCR study](https://pmc.ncbi.nlm.nih.gov/articles/PMC2929689/)
support distinguishing ocular/inner-ear injury from skin or lymphoid-organ loss.
None establishes the numerical thresholds.

Review conditions are machine-readable in `context_reviews` and returned by
`normal_tissue_screening_details()`. Candidate reports flag measured positive
RNA (>0 nTPM) or positive IHC in these contexts. A flag remains visible when the
candidate passes; it is not a new numerical tier or a completed clinical review.
Measurements in every other organ retain that organ's independent RNA/IHC gates.

### Expression rules

The **fixed tiered CTA screen** uses the following exact research policy:

| Tier | Normal RNA rule | Normal IHC rule |
|---|---|---|
| Acceptable (machine ID `permitted`) | No upper ceiling | Positive staining allowed; conditional review annotations retained |
| Vital / critical function (machine ID `vital`) | Every measurement ≤1 nTPM | Any positive staining excludes, including Low |
| Protected normal | Every measurement ≤5 nTPM | Any positive staining excludes, including Low |
| Unresolved | Positive evidence prevents retention | Positive evidence prevents retention |

**Head, foot and generic soft tissue remain unresolved** because they can span
tiers. A new or unknown tissue never silently receives the 5 nTPM allowance.
Explicit embryonic/fetal source context retains the requested reproductive
permission; adult organs retain their adult tier.

The numerical rules use **raw human HPA v23 nTPM**, with inclusive comparisons
and no deflation or averaging across tissues. IHC Low, Medium, High, Ascending
and Descending are positive. A gene additionally needs ≥1 nTPM or positive IHC
inside the permitted set and complete coverage of all 50 routine bulk RNA
contexts. Unresolved positive anatomy, unusable measurements or reviewed
unreliable IHC prevents retention. Existing family exclusions and
nomination/manual-review decisions remain separate, visible gates.

Missing IHC remains **unassessed**, not negative. The report separates RNA-only,
partial-IHC and complete-routine-IHC results; 80 routine tissue/cell pairs define
the last category. Negative measurements and coverage are distinct. The gate
does not establish adequate evidence of absence or clinical safety.

This tier list and the 1/5 nTPM values are **project screening definitions**.
[NCI organ-inflammation guidance](https://www.cancer.gov/about-cancer/treatment/side-effects/organ-inflammation)
supports organ-specific review of immune injury, and
[NIDDK adrenal guidance](https://www.niddk.nih.gov/health-information/endocrine-diseases/adrenal-insufficiency-addisons-disease/definition-facts)
supports treating the adrenal/pituitary axis conservatively. Neither source
validates the numerical ceilings or a claim that permitted-organ injury is
tolerable. The skin allowance is whole-organ by request; the
[HPA skin anatomy](https://www.proteinatlas.org/learn/dictionary/normal/skin)
documents its multiple compartments and functions.

Anatomical tiers apply to human, dog and mouse with species filtering. Applying
the human HPA thresholds to canine or mouse assays requires a separate policy;
cross-species anatomical mapping does not transfer assay units or tolerability.

```python
from oncoref.normal_tissues import (
    normal_tissue_screening_policy,
    normal_tissue_screening_details,
    normal_tissue_screening_tier,
    resolve_normal_tissue,
)

policy = normal_tissue_screening_policy()
assert policy["vital"]["maximum_rna"] == 1
assert policy["protected"]["maximum_rna"] == 5
assert normal_tissue_screening_tier(resolve_normal_tissue("pancreas", species="human")) == "vital"
assert normal_tissue_screening_tier(resolve_normal_tissue("retina", species="human")) == "vital"
assert normal_tissue_screening_tier(resolve_normal_tissue("lymph node", species="human")) == "protected"
assert normal_tissue_screening_tier(resolve_normal_tissue("skin", species="dog")) == "permitted"
details = normal_tissue_screening_details(resolve_normal_tissue("spleen", species="human"))
assert details["screening_tier_label"] == "Acceptable (conditional)"
assert "asplenia_infection_risk_review" in details["review_flags"]
```

## IHC-adaptive absolute RNA ceilings

The `ihc_adaptive_cta` policy uses the same permitted union and anatomical tiers,
with these **inclusive maximum human HPA v23 nTPM values**:

| IHC evidence in the RNA tissue | Vital / critical function | Protected |
|---|---:|---:|
| Missing, incomplete, uncertain or unusable | 1 | 5 |
| Approved, all routine contexts Not detected | 2 | 5 |
| Supported, all routine contexts Not detected | 3 | 7.5 |
| Enhanced, all routine contexts Not detected | 5 | 10 |

The routine tissue/cell inventory is independent of the candidate being assessed:
the report freezes contexts measured for at least half of the source's measured
genes (80 IHC tissue/cell pairs). This inventory is not exhaustive biological coverage.
Every expected context for that anatomical tissue must be measured and negative;
the weakest reliability among all its IHC observations determines the ceiling.
Additional positive or invalid observations also prevent relaxation. Kidney RNA
requires both routine glomerular and tubular IHC contexts. Stomach RNA requires
both source labels, `stomach 1` and `stomach 2`. Caudate-only IHC cannot relax
broader basal-ganglia RNA. Negative staining in another organ supplies no support.

Positive IHC anywhere outside the permitted set still excludes, including Low.
Complete routine RNA, positive evidence inside the permitted set, source resolution,
and the existing family/manual-review gates still apply. Missing IHC remains a
coverage limitation. The ladder is an approved **project screening policy**:
HPA's Enhanced label does not establish a quantitative detection limit or absence
of peptide-HLA presentation, and these numbers are not validated clinical limits.

```python
from oncoref.cta_normal_screen import adaptive_normal_rna_limits
from oncoref.normal_tissues import normal_tissue_screening_policy

policy = normal_tissue_screening_policy("ihc_adaptive_cta")
assert policy["negative_ihc_rna_limits"]["Enhanced"] == {"vital": 5, "protected": 10}

# One gene's annotated HPA v23 observation table; an independently specified
# inventory of expected source tissue/cell pairs is required.
audit = adaptive_normal_rna_limits(
    observations,
    routine_ihc_contexts=routine_ihc_contexts,
)
```

The [comparative report](audits/cta-ihc-adaptive-20261020/index.md) preserves each
RNA value, ceiling, supporting IHC and exclusion. It also explains CTAG2's
**pinned consensus heart value of 5.1 nTPM**: HPA's own cohort and the separately
displayed GTEx estimates are different datasets. Current portal values do not
replace the pinned release. This absolute screen is separate from the older
IHC-conditioned **RNA fraction** policy used by the core/extended CTA panels.
Numerical adaptation is currently implemented only for human HPA v23;
cross-species anatomical definitions do not transfer these assay thresholds.

## Dataset mappings and coverage

Mappings are keyed by **taxon + source ID + source version + tissue + cell type**.
The original labels are always returned. Matching normalizes case and whitespace;
it does not guess using substrings. Unknown labels or versions remain unresolved.
A known organ can support an organ-level annotation when the cell subtype is
unresolved; that mapping cannot support a narrower cell-type scope.

| Data source | Mapping behavior |
|---|---|
| HPA `hpa_rna_consensus`, `v23` | Whole-organ labels, including breast, testis, ovary and placenta. Skin RNA cannot resolve follicles or melanocytes. |
| HPA `hpa_normal_tissue`, `v23` | Preserves tissue and cell type together. `breast / glandular cells` and `breast / myoepithelial cells` map to glandular mammary tissue; `breast / adipocytes` retains whole-organ context. `skin 1`/`skin 2` with `melanocytes` or `hair follicles` map to the corresponding disposable compartment. The `hair` tissue's annotated epithelial compartments also map to hair follicle. `endometrium 1` and `endometrium 2` map to endometrium while preserving the original label. |
| HPA `hpa_single_cell`, `v23` | Oocytes, testicular cell types, mammary glandular cells, trophoblasts and other anatomically resolved cell labels map individually. HPA's skin-derived `Melanocytes` label maps to cutaneous melanocytes. Generic fibroblasts, ductal cells and other unsited cell types remain unresolved. |
| Canine `bgee-srp009687`, `bgee-srp181218`, `Bgee 15.2` | Uses the bundled curated tissue fields, including testis and ovary. Original anatomical IDs, source labels and any corrections remain in sample metadata. `life cycle` does not imply adult or fetal. |
| Canine `barkbase`, `retrieved-2026-10-08` | Native labels such as `frontalcortex` and `kidneycortex` retain their labels and map to frontal cortex and renal cortex. Explicit embryonic versus adult sample stages are preserved in annotation and tissue summaries. |
| Canine `pfizer-dog-protein`, `MSV000099744-update-2026-03-23` | `Breast/Mammary` → whole mammary tissue; `Fallopian_Tube` → fallopian tube; `Uterus_Cervix` → cervix; `Uterus_Endometrium` → endometrium; ovary, prostate and testis map directly. Bulk `Skin` does not qualify for the disposable scope. |
| Canine `hindriks2025`, `2025-05-06-publisher-version`; `nemec2019`, `2019-article-abstract-europepmc-2026-10-13` | Curated normal observation labels map at their reported resolution. Placenta qualifies directly. Tumor observations are excluded from normal-scope annotation. |
| Mouse | Canonical anatomical categories are available now. No native mouse expression/CTA atlas is bundled yet; no source mapping or measured coverage is fabricated. |

Every CSV row cites its source. The HPA vocabulary is checked against the pinned
v23 source files; [HPA histology](https://www.proteinatlas.org/learn/dictionary/normal/skin)
and [pigment-cell documentation](https://www.proteinatlas.org/humanproteome/single%2Bcell/single%2Bcell%2Btype/pigment%2Bcells)
support the compartment distinctions. Canine source evidence and checksum receipts
remain in the existing species metadata.

Species filtering excludes seminal vesicles and bulbourethral glands in dogs,
consistent with [comparative veterinary anatomy](https://veterinaryanatomy.com/atlas/hindlimb-pelvis/splanchnology/male-reproductive-pelvic-and-external-genitalia/).
Mouse reproductive anatomy follows the [Jackson Laboratory anatomy resource](https://www.informatics.jax.org/cookbook/chapters/contents2.shtml).
Anatomical transfer does not transfer evidence of expression or clinical
tolerability: animal annotations explicitly report that tolerability is unvalidated.

The [complete registry audit](audits/normal-tissue-registry-20261020/index.md)
exports every canonical concept for each species and every bundled native mapping.
`normal_tissue_inventory(species=...)` includes **inapplicable anatomy** explicitly;
`normal_tissue_mapping_inventory(species=..., source_id=..., source_version=...)`
reports native labels, mapping precision, citations and four-way categories.
No mapping inventory for a source means **not bundled**, never expression absent.
Species applicability is validated through every anatomical parent.

Source sex and physiological state can be supplied separately from anatomy.
Unknown sex stays unknown; the grouping of mammary tissue with female reproductive
anatomy never proves the donor is female. Unknown stage labels now remain unresolved
instead of inheriting an adult allowance. No gestational cutoffs or ortholog
expression calls are copied between species.

The [patient-context gene filters](normal-tissue-gene-filters.md) evaluate actual
sex/stage/source-specific detections and declared negative coverage. Their names
describe expression restrictions, not pregnancy or fertility safety.

## Use

```python
from oncoref import human, canine, mouse
from oncoref.normal_tissues import classify_normal_tissue_pattern

mouse.normal_tissue_scopes()
observations = [
    mouse.resolve_normal_tissue("testis"),
    mouse.resolve_normal_tissue("mammary gland"),
]
pattern = classify_normal_tissue_pattern(observations)
assert pattern["narrowest_scopes"] == ("reproductive",)
assert pattern["absence_evidence"] == "not_assessed"

follicles = human.resolve_normal_tissue(
    "skin 1", cell_type="hair follicles",
    source_id="hpa_normal_tissue", source_version="v23",
)
assert "disposable_normal" in follicles.scopes

canine.cta_tissue_annotations(candidate_id="CTCFL")
# dataset = canine.expression.read("barkbase", offline=True)
# dataset.tissue_annotations()  # current taxonomy, indexed by original sample ID
# dataset.tissue_summary()      # groups tissue AND source stage
```

`matching_scopes` includes all ancestors, and `narrowest_scopes` identifies the
minimal compatible scopes (more than one is possible for overlapping anatomy and
stage). Testis + oocyte qualifies for the reproductive core and broader reproductive
unions, but neither male-only nor female-only. Testis + skin melanocytes qualifies
with the combined reproductive-or-disposable union as its narrowest scope. Testis + adult heart
qualifies for none of these scopes. Empty detections never qualify vacuously;
unmapped detections block a supported pattern classification.

Detection thresholds remain assay-specific. The caller supplies detected
observations; no shared numerical threshold is applied to RNA, IHC or PCR.
Coverage of negative observations must be assessed separately before claiming
strict restriction. The `cta_tissue_annotations()` view explicitly applies the
current taxonomy; saved CTA v1/v2 facts, roles and source values are not rewritten.
Human summaries carry an explicit `HPA:v23` lineage marker for source mapping;
an older archive without a declared HPA version remains unresolved in this view.

## Existing human CTA tables

Both HPA-derived tables are regenerated with the corrected anatomy. The extended
RNA numerator now includes breast. `permissive` is an alias of extended in the
shared taxonomy, and the legacy tier audit's two columns now agree. Thymus is
counted in somatic maxima/detection counts and cannot establish reproductive
protein restriction. Explicit `*_and_thymus_frac` summaries retain their name and
meaning.

Single-organ labels now require the complete detected-tissue set to match that
organ. Testis + prostate or testis + breast yields `REPRODUCTIVE`, not `TESTIS`.
RNA restriction considers all detected normal tissues for either numerator scope.
The synthesized label unions RNA/protein scope evidence; somatic counterevidence
cannot retain a testis-only label merely by lowering confidence. These legacy
labels summarize observations and are not substitutes for adequate negative
tissue coverage. The richer shared scopes are available through the new API.

The core panel remains 624 genes. The extended panel is now 695 genes, gaining
ABCC11, CSN2, CSN3, CST9, LALBA, ROPN1 and ROPN1B. The low-expression rescue
requires positive evidence inside the selected scope, so an accessory-only trace
cannot enter the core panel through an all-zero deflated denominator.
