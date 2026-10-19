# Normal tissue scopes

oncoref uses one cited anatomical vocabulary for human, dog and mouse. A **scope
is a set of permitted normal anatomical/developmental contexts**. A detected
protein or gene is compatible with that scope only if **every detection** belongs
to it. Tissue sets are unioned before evaluating expression patterns. Positive
observations alone do not establish absence outside a scope.

Definitions and citations ship in `oncoref/data/normal-tissue-taxonomy.yaml`;
exact dataset mappings ship in `oncoref/data/normal-tissue-mappings.csv`.
Policy version `2026-10-19.1` and a SHA-256 identity over both files accompany
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

For a gene detected in both testis and oocyte, **the entire observed pattern**
fits `reproductive_core` and `reproductive`, while neither sex-specific scope
contains both detections. The tissue union is evaluated before selecting genes.
The aliases `male_strict` and `female_strict` name anatomical scopes; a claim of
strict gene restriction additionally needs adequate negative tissue coverage.

The [current CTA source and filtering report](audits/cta-tissue-scopes-20261019/index.md)
compares reproductive with reproductive + disposable using both the existing
adaptive CTA policy and a strict positive-evidence sensitivity screen. It includes
figures, exact gene decisions, source citations and coverage limitations.

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
only for the combined reproductive-or-disposable union. Testis + adult heart
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
