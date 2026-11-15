# Normal-tissue tiers

OncoRef 1.29.0 centralizes **Vital / Protected / Manageable / Disposable** calls for adult human organs and resolved cell populations. The [complete report](audits/tissue-hazards-1.29.0/index.md) includes 107 assessments, 50 linked references, four generated figures, every canonical human anatomy concept and the native source mappings. A [searchable HTML report](audits/tissue-hazards-1.29.0/report.html), [PDF](audits/tissue-hazards-1.29.0/tissue-hazard-ontology.pdf) and CSV/JSON exports accompany it.

The tiers describe **consequences of sustained injury**, including loss of repair capacity. They are policy judgments, not measured fatality probabilities. The evidence ledger distinguishes clinical cell-therapy injury, selective animal depletion, regeneration experiments, functional disease analogies, surgical cohorts and anatomical references. An anatomy-only citation supports identity, not tolerability. This is a curated evidence review, not an exhaustive systematic review.

| Tier | Interpretation | Examples and boundaries |
|---|---|---|
| Vital | Highest protection priority for life-threatening failure, major irreversible sensory/neurological injury or indispensable repair | Heart, neurons, podocytes, AT1/AT2, airway club/basal and intestinal crypt compartments. Compensation after finite depletion remains relevant. |
| Protected | Persistent injury is consequential and compensation or control is insufficiently established | Ciliated/goblet airway cells, salivary tissue, generic fibroblasts, whole skin, spleen and thymus. This is not a finding of clinical safety. |
| Manageable | Defined loss has a credible compensation or treatment precedent | Thyroid and mature B cells; replacement/support does not control every inflammatory manifestation. Adipose and mammary calls retain explicit inference and extent limitations. |
| Disposable | Narrow conditional compartment-loss allowance | Cutaneous melanocytes, hair follicles, appendix, gallbladder and tonsil. Shared cells elsewhere and in-situ inflammation require separate assessment. |

## How a tier is chosen

**A tier answers one question: if this tissue were lost and stayed lost, what would
follow, and is there a precedent for living with it?** Each call weighs acute and
persistent loss, repair, compensation and management, with references, in the
registry record. Where cancer treatment already causes the same loss, that record is
the most direct evidence there is, so it is used the same way for every tissue.

| Oncology precedent | Tier it supports | Examples |
|---|---|---|
| The loss is graded at most Grade 2 by CTCAE v5, standard cancer treatments that cause it stay in use, and it is not a recorded reason to stop or reduce them | Manageable, or Disposable when the loss has no lasting consequence | Olfactory epithelium (anosmia is Grade 1 only); hair follicle (alopecia) |
| Standard treatment causes it, but it limits dose, causes discontinuation or is graded up to Grade 3 or 4 | No lower than Protected | Taste buds (taste loss stops vismodegib in 6%); inner ear (cisplatin hearing loss is dose-limiting); sensory ganglia (dose-limiting neuropathy) |
| No precedent | Decided on repair, compensation and disease analogues | Most organs |

A precedent outside cancer counts against tolerance: intranasal zinc was withdrawn
after anosmia reports, because the loss was not worth it for a cold. A tier is a
policy judgment about tolerating loss in cancer treatment, not a finding of safety.

### Sources for the precedent

These are recorded with the tissue they support in
`oncoref/data/normal-tissue-tier-updates.yaml` (`tissue_tiers.update_references()`).
Abstracts and summaries were read, not full texts.

| Source | Finding | Used for |
|---|---|---|
| [CTCAE v5.0, Anosmia](https://www.cismef.org/page/detail/en/NCI_CO_C143230); [Olfactory nerve disorder, Grade 2](https://www.cismef.org/page/detail/en/NCI_CO_C144871); [trial protocol table](https://cdn.clinicaltrials.gov/large-docs/06/NCT03954106/SAP_001.pdf) | Anosmia has one grade, Grade 1 "Present" | Olfactory epithelium |
| [Taste and smell disorders in cancer treatment, 2023 review](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9916934/) | Smell disorders in 8–45% and taste disorders in 17–86% of patients on standard chemotherapy | Olfactory epithelium, taste buds |
| [Treatments for taste and smell disturbance in cancer, 2021 scoping review](https://pubmed.ncbi.nlm.nih.gov/32734392/) | Few management tools; no firm conclusions | Olfactory epithelium |
| [Radiotherapy and smell in head and neck cancer, 2025 trial](https://www.epistemonikos.org/en/documents/389d1578f48ebdba37920f7c8460c14bd6e3dc1c) | Dose-dependent smell loss above about 22 Gy, lasting up to a year; radiotherapy given | Olfactory epithelium |
| [Olfactory bulb after nasopharyngeal radiotherapy, 2014](https://pubmed.ncbi.nlm.nih.gov/24746667/) | Smaller bulbs and lower smell scores a year or more after curative treatment | Olfactory epithelium |
| [FDA consumer update on intranasal zinc, 2009](https://webarchive.library.unt.edu/web/20090807210616/http://www.fda.gov/ForConsumers/ConsumerUpdates/ucm166931.htm); [news report](https://abcnews.com/Health/ColdandFluNews/story?id=7853178) | More than 130 anosmia reports; products recalled | Olfactory epithelium: the limit of tolerance outside cancer |
| [Croy, Nordin and Hummel 2014](https://pubmed.ncbi.nlm.nih.gov/24429163/); [depression, olfaction and quality of life](https://www.doi.org/10.3390/BRAINSCI8050080) | A subset report lower quality of life and more depression | Olfactory epithelium: why not Disposable |
| [Pinto et al. 2014](https://knowledge.uchicago.edu/record/10483/files/journal.pone.0107541.pdf) | Anosmia with 3.4 times the odds of death in five years in older adults; an association | Olfactory epithelium: why not Disposable |
| [Vismodegib adverse events, including STEVIE](https://academic.oup.com/oncolo/article-abstract/21/10/1218/6401021); [managing Hedgehog inhibitor effects](https://practicaldermatology.com/articles/2017-jan/practical-tips-for-managing-hedgehog-pathway-inhibitor-side-effects); [product information](https://www.cancercareontario.ca/en/node/44341); [sonidegib taste-bud recovery](https://pmc.ncbi.nlm.nih.gov/articles/PMC5715770); [sonidegib and vismodegib discontinuation](https://www.analysisgroup.com/Insights/publishing/real-world-treatment-patterns-of-patients-with-basal-cell-carcinoma-using-sonidegib-and-vismodegib-discontinuation-rates-and-clinical-conditions-during-treatment/) | Taste disturbance in 54%, ageusia in 22%; taste loss causes 6% of discontinuations and weight loss 5%; some cases last a year after stopping | Taste buds stay Protected |
| [FDA approval of sodium thiosulfate, 2022](https://www.fda.gov/drugs/resources-information-approved-drugs/fda-approves-sodium-thiosulfate-reduce-risk-ototoxicity-associated-cisplatin-pediatric-patients); [ASCO Post](https://ascopost.com/news/september-2022/fda-approves-sodium-thiosulfate-to-reduce-the-risk-of-ototoxicity-associated-with-cisplatin) | Cisplatin hearing loss in 63% of children; cisplatin kept and prevention added | Lasting sensory loss can be accepted, but hearing loss is dose-limiting, so the inner ear is not lowered |

Bulk observations retain their organ-level concern. A cell-specific exception requires a recognized cell identity in a cellular assay. Unknown cells inherit the organ rule; mixed populations retain the more protective matched concern. Ambiguous cell rules, unknown anatomy and contexts outside adult human scope remain unresolved in the public assessment API. The screening cascade uses the most protective matched tier for ambiguous cell rules and retains incoming concern when anatomy is unresolved. Reproductive/developmental membership is independent: an ovary can be Reproductive for CTA association and Protected for functional injury.

```python
from oncoref import tissue_tiers

bulk = tissue_tiers.assess("lung")
assert bulk["tier"] == "vital"
cell = tissue_tiers.assess("lung", "ciliated cells", resolution="cell")
assert cell["tier"] == "protected"
assert cell["assessment_id"] == "airway_ciliated"

rules = tissue_tiers.assessments()  # rationale, acute/persistent loss, repair, references
sources = tissue_tiers.references()  # design, model, supported finding, limitation
inventory = tissue_tiers.anatomy_inventory()  # selection scope, tier, localization gaps
policy = tissue_tiers.registry()  # detached data plus SHA-256 policy identity
```

The current inventory's `selection_category` groups placental, extraembryonic,
embryonic and fetal contexts under **Reproductive**. `tier` separately records
the adult injury assessment where one exists. **Soft tissue, head and foot**
require localization; they do not receive a permissive tier. The
[current inventory and coverage audit](audits/target-selection-20261109/expression-maxima/index.md)
separates these broad labels from specific compartments needing new curation.

## Registry updates

**Seven tissues were added on 2026-11-13 as dated updates; no base call changed.
Olfactory epithelium moved from Protected to Manageable on 2026-11-15, on the oncology
precedent above.**
The base registry is `oncoref/data/normal-tissue-tiers.yaml`. Additions live in
`oncoref/data/normal-tissue-tier-updates.yaml` and are replayed over it, as
gene-knowledge updates are. `tissue_tiers.assess()` resolves both;
`tissue_tiers.registry_assess()` answers from the base registry alone. An update can
only add a tissue the base registry does not resolve.

| Tissue | Tier | Confidence |
|---|---|---|
| Sensory ganglia (dorsal root, trigeminal) | Vital | high |
| Bone, osteoblasts and osteocytes | Vital | medium |
| Olfactory epithelium | Manageable | medium |
| Taste buds | Protected | medium |
| Lacrimal gland | Protected | medium |
| Nasal respiratory mucosa | Protected | low |
| Tooth and dental pulp | Manageable | medium |

They were written from abstracts of clinical and ablation analogues; full texts were
not read. They tier the tissues that published normal-expression findings point to
([literature rule](target-selection.md#published-normal-expression)). Fold them into
the base file at the next full data rebuild.

```python
from oncoref import tissue_tiers

tissue_tiers.assess("olfactory epithelium")  # tier "manageable", reason "registry_update"
tissue_tiers.assess("lung")                  # unchanged base call
tissue_tiers.updates()                       # the added rules, with rationale and references
tissue_tiers.update_references()             # their evidence ledger
tissue_tiers.policy_id()                     # base registry hash + updates hash
```

RNA and IHC use the same registry, with cellular resolution for HPA clusters/staining, sorted immune panels, DICE and the adult cornea panel. Bulk RNA cannot borrow a cell exception. Existing numeric thresholds remain research parameters; this review does not recalibrate them. Dog and mouse screening roots retain their previous compatibility assignments; the human tier API returns unresolved for nonhuman species.

The existing CELLxGENE screen still covers brain, spinal cord, heart and selected lung populations. Qualifying Vital expression supplies the existing Unsafe screen; qualifying Protected expression creates `rare_cell_hold`, which blocks CTA admission without itself establishing Unsafe. Other organs' ontology entries do not imply a new expression panel. Excluded lineages and unresolved donor ages remain coverage limitations. Detection frequency is distinct from abundance, receptor recognition and tissue injury. The [2026-11-10 target policy](target-selection.md#rare-normal-cell-rna) regrades this panel per adult donor, ignores single-molecule cells and ends every signal as established or dismissed; the released 1.36.0 catalog keeps the pooled rule described here.

Edit `oncoref/data/normal-tissue-tiers.yaml` to change a call. Each record requires organ scope, tier, mechanism, acute/persistent consequences, repair, compensation, rationale and references. The reference ledger requires a biological model, evidence design, supported finding and limitation. Keep receptor-specific evidence separate from source-gene RNA. Conflicting rules must be resolved explicitly; registry order cannot choose the more permissive rule.

Regenerate the ontology artifacts with `python scripts/build_tissue_tier_report.py` (Matplotlib and ReportLab). After a policy change, rebuild the rare-cell screen from its recorded raw cache and run `python scripts/regenerate_cta_walkthrough.py --verify-replay`, then rebuild the ontology report. The [manifest](audits/tissue-hazards-1.29.0/manifest.json) fingerprints inputs and outputs; [classification changes](audits/tissue-hazards-1.29.0/classification-changes.csv) compare against the recorded 1.28.0 baseline. Historical reports retain their original rules.
