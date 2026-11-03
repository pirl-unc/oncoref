# OncoRef 1.29.0: central normal-tissue hazard ontology

**625 CTA / 18,351 Unsafe / 628 Other** representative proteins across **19,604** proteins. Other includes **550 Nonreproductive, 28 Nested, 27 Held, 19 Contained and 4 Diffuse**. No human protein has an unassigned category.

[Ontology and API](tissue-hazards.md) · [Complete evidence report](audits/tissue-hazards-1.29.0/index.md) · [Searchable report](audits/tissue-hazards-1.29.0/report.html) · [PDF](audits/tissue-hazards-1.29.0/tissue-hazard-ontology.pdf) · [Classification changes](audits/tissue-hazards-1.29.0/classification-changes.csv)

## Central calls and their evidence

One versioned registry supplies **39 bulk assessments and 68 cell-population assessments**, backed by **50 references**. Every record separates immediate consequences, persistent injury, repair, compensation, rationale and evidence scope. The bibliography distinguishes human cell-therapy injury, animal depletion/repair experiments, functional analogies, surgical cohorts and anatomical references. Four generated figures explain the ontology, organ coverage, lung evidence and bulk-versus-cell distinctions. CSV/JSON inventories cover all 123 canonical human anatomy concepts and 1,072 native mapping rows, including unresolved contexts.

Vital, Protected, Manageable and Disposable are qualitative protection priorities. They do not estimate fatality rates. Expression, accessible antigen or peptide-HLA, receptor recognition, inflammatory amplification and final organ failure remain distinct questions. Reproductive membership stays separate from functional hazard.

| Compartment | Current call | Main reason |
|---|---|---|
| Ciliated and goblet lung cells | Protected | Clearance/defense loss and persistent injury warrant protection; direct selective immune-killing evidence is limited. |
| AT1 and AT2 | Vital | Alveolar function and repair. AT2 has stronger lineage-specific experimental lethality evidence; finite AT1 loss can be compensated. |
| Club and basal airway cells | Vital | Sustained regenerative-compartment injury can defeat repair. |
| Generic fibroblasts and differentiated bowel epithelium | Protected | Matrix, barrier and regenerative support cannot be presumed manageable from turnover alone. |
| Intestinal crypt/undifferentiated populations | Vital | Broad loss threatens epithelial renewal; selective Lgr5 depletion demonstrates alternative repair, not general progenitor dispensability. |
| Adult whole skin, spleen and thymus | Protected | Barrier injury or lasting immune consequences require stronger protection than the previous Manageable assignment. |
| Mature B cells | Manageable | Clinical depletion and immunoglobulin-support precedent, with explicit limitations. |
| Plasma cells/plasmablasts | Protected | BCMA CAR-T evidence documents prolonged humoral deficiency; this compartment is distinguished from mature circulating B cells. |

Bulk observations cannot borrow a narrow cell exception. Unknown cell labels retain organ concern; mixed and ambiguous rules cannot silently choose a more permissive match. Known immune aliases resolve explicitly. Shared immune cells within a removable organ do not inherit Disposable status. The public hazard API leaves nonadult and nonhuman contexts unresolved; dog/mouse screening roots preserve their previous compatibility policy.

## Effect on the current selection

**104 classifications change** relative to the recorded 1.28.0 snapshot:

| Previous | Current | Proteins |
|---|---|---:|
| CTA | Unsafe | 12 |
| CTA | Other:Held | 5 |
| Other | Unsafe | 74 |
| Unsafe | CTA | 5 |
| Unsafe | Other | 8 |

The five new holds are **FKBP6, PRSS54, SPMIP3, TXNDC8 and WFDC6**, through Protected rare-cell expression. The five entries are **DAZL, FSCN3, C12orf71, PSG3 and PAGE2**, after resolved cell rules replace broader tissue assumptions. Entry is a research-policy result, not demonstrated tolerability.

**MAGEA3 moves from Held to Unsafe**: its colon undifferentiated-cell signal is 53.1 nTPM against the Vital 4 nTPM ceiling. **MAGEB2 moves from CTA to Unsafe** under the same regenerative-compartment rule, with 118.4 nTPM against its 32 nTPM ceiling. The [changed-protein RNA contexts](audits/tissue-hazards-1.29.0/changed-protein-rna-contexts.csv) preserve native units and source-specific ceilings. Prior and current exclusion reasons accompany every changed classification. PRAME, MAGEA4, NY-ESO-1, CTAG2, CT83, XAGE1A/B, SSX2/B and TBXT retain CTA status.

The CELLxGENE panel retains its existing brain/spinal-cord/heart/lung scope and source snapshot. Qualifying Vital detection still supplies Unsafe evidence; Protected detection now creates a review hold. Profiles expose `rare_cell_hold` independently of `rare_cell_flag`. The screen records 15,874 proteins with a Vital flag and 13,399 with a Protected hold flag; these overlap heavily with other exclusions. Neither a new tissue expression panel nor clinically calibrated thresholds is implied.

## API and reproduction

```python
from oncoref import tissue_hazards, cta_selection

assert tissue_hazards.assess("lung")["tier"] == "vital"
assert tissue_hazards.assess(
    "lung", "ciliated cells", resolution="cell"
)["tier"] == "protected"

rules = tissue_hazards.assessments()
references = tissue_hazards.references()
profiles = cta_selection.profiles()  # includes rare_cell_hold
```

`oncoref/data/normal-tissue-hazards.yaml` is the authoritative registry. RNA, IHC and rare-cell rules consume it. `scripts/build_tissue_hazard_report.py` produces the report and inventories with content receipts. The full selection replay compares complete and cached outputs. Historical decisions retain their original hashes and archived source inputs. The wheel now explicitly includes rare-cell JSON provenance alongside the registry and screen data.

This release is an adult-human evidence curation, not an exhaustive systematic review. Some cell calls rely on physiology or adjacent-compartment evidence, labeled explicitly. Unreviewed developmental, gestational and nonspecific anatomical contexts remain visible. Persistent immune recognition, inflammatory amplification, target distribution and disease/repair states can defeat compensation seen after finite injury. The existing RNA/IHC numeric parameters and 1% rare-cell screening line remain uncalibrated research thresholds.
