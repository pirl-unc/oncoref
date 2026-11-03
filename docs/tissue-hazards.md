# Normal-tissue hazard ontology

OncoRef 1.29.0 centralizes **Vital / Protected / Manageable / Disposable** calls for adult human organs and resolved cell populations. The [complete report](audits/tissue-hazards-1.29.0/index.md) includes 107 assessments, 50 linked references, four generated figures, every canonical human anatomy concept and the native source mappings. A [searchable HTML report](audits/tissue-hazards-1.29.0/report.html), [PDF](audits/tissue-hazards-1.29.0/tissue-hazard-ontology.pdf) and CSV/JSON exports accompany it.

The tiers describe **consequences of sustained injury**, including loss of repair capacity. They are policy judgments, not measured fatality probabilities. The evidence ledger distinguishes clinical cell-therapy injury, selective animal depletion, regeneration experiments, functional disease analogies, surgical cohorts and anatomical references. An anatomy-only citation supports identity, not tolerability. This is a curated evidence review, not an exhaustive systematic review.

| Tier | Interpretation | Examples and boundaries |
|---|---|---|
| Vital | Highest protection priority for life-threatening failure, major irreversible sensory/neurological injury or indispensable repair | Heart, neurons, podocytes, AT1/AT2, airway club/basal and intestinal crypt compartments. Compensation after finite depletion remains relevant. |
| Protected | Persistent injury is consequential and compensation or control is insufficiently established | Ciliated/goblet airway cells, salivary tissue, generic fibroblasts, whole skin, spleen and thymus. This is not a finding of clinical safety. |
| Manageable | Defined loss has a credible compensation or treatment precedent | Thyroid and mature B cells; replacement/support does not control every inflammatory manifestation. Adipose and mammary calls retain explicit inference and extent limitations. |
| Disposable | Narrow conditional compartment-loss allowance | Cutaneous melanocytes, hair follicles, appendix, gallbladder and tonsil. Shared cells elsewhere and in-situ inflammation require separate assessment. |

Bulk observations retain their organ-level concern. A cell-specific exception requires a recognized cell identity in a cellular assay. Unknown cells inherit the organ rule; mixed populations retain the more protective matched concern. Ambiguous cell rules, unknown anatomy and contexts outside adult human scope remain unresolved in the public assessment API. The screening cascade uses the most protective matched tier for ambiguous cell rules and retains incoming concern when anatomy is unresolved. Reproductive/developmental membership is independent: an ovary can be Reproductive for CTA association and Protected for functional injury.

```python
from oncoref import tissue_hazards

bulk = tissue_hazards.assess("lung")
assert bulk["tier"] == "vital"
cell = tissue_hazards.assess("lung", "ciliated cells", resolution="cell")
assert cell["tier"] == "protected"
assert cell["assessment_id"] == "airway_ciliated"

rules = tissue_hazards.assessments()  # rationale, acute/persistent loss, repair, references
sources = tissue_hazards.references()  # design, model, supported finding, limitation
inventory = tissue_hazards.anatomy_inventory()  # includes explicit unresolved concepts
policy = tissue_hazards.registry()  # detached data plus SHA-256 policy identity
```

RNA and IHC use the same registry, with cellular resolution for HPA clusters/staining, sorted immune panels, DICE and the adult cornea panel. Bulk RNA cannot borrow a cell exception. Existing numeric thresholds remain research parameters; this review does not recalibrate them. Dog and mouse screening roots retain their previous compatibility assignments; the new human hazard API returns unresolved for nonhuman species.

The existing CELLxGENE screen still covers brain, spinal cord, heart and selected lung populations. Qualifying Vital expression supplies the existing Unsafe screen; qualifying Protected expression creates `rare_cell_hold`, which blocks CTA admission without itself establishing Unsafe. Other organs' ontology entries do not imply a new expression panel. Excluded lineages and unresolved donor ages remain coverage limitations. Detection frequency is distinct from abundance, receptor recognition and tissue injury.

Edit `oncoref/data/normal-tissue-hazards.yaml` to change a call. Each record requires organ scope, tier, mechanism, acute/persistent consequences, repair, compensation, rationale and references. The reference ledger requires a biological model, evidence design, supported finding and limitation. Keep receptor-specific evidence separate from source-gene RNA. Conflicting rules must be resolved explicitly; registry order cannot choose the more permissive rule.

Regenerate the ontology artifacts with `python scripts/build_tissue_hazard_report.py` (Matplotlib and ReportLab). After a policy change, rebuild the rare-cell screen from its recorded raw cache and run `python scripts/regenerate_cta_walkthrough.py --verify-replay`, then rebuild the ontology report. The [manifest](audits/tissue-hazards-1.29.0/manifest.json) fingerprints inputs and outputs; [classification changes](audits/tissue-hazards-1.29.0/classification-changes.csv) compare against the recorded 1.28.0 baseline. Historical reports retain their original rules.
