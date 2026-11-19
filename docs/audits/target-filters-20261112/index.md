# Target filters: categories, knobs and where proteins land

**The catalog assigns five categories and they are not changed here. Filters only
remove proteins from the 523 Target candidates, through six layers that
follow the catalog's own shape.** Filter policy `oncoref.target-filters.2026-11-17.5`, applied to
catalog `oncoref.target-selection.2026-11-19.2`. Every number below is computed by
`scripts/build_target_filter_audit.py`. [Rules reference](../../target-filters.md).

These are outputs of screening rules. They are not validated tumor antigens.

## 1. The categories

| Category | Proteins | Rule |
|---|---|---|
| CTA | 327 | Passes the shared screen; cancer evidence and reproductive evidence both qualify |
| Cancer Antigen only | 69 | Passes the shared screen; only cancer evidence qualifies |
| Reproductive Antigen only | 127 | Passes the shared screen; only reproductive evidence qualifies |
| Other | 1,117 | Not excluded, but fails coverage, identity, expression bounds, sequence, the home-tissue rule or both evidence routes |
| Unsafe | 18,253 | An accepted normal-tissue or clinical exclusion |

Target candidates are the first three rows: 523 proteins. A filter never moves
a protein into or out of Other or Unsafe.

**One cancer rule.** A protein is a Cancer Antigen when one tumor sample, from any
QC-passing comparable TPM cohort including germ cell tumors, is at least 10 times its highest
non-reproductive normal value (normal tissue, immune cells and Disposable tissue; values
below 1 TPM count as 1). The 396 Cancer Antigens: 115 at 10 to 20x, 116 at 20 to 50x, 46 at 50 to 100x, 106 at 100 to 1000x, 13 at 1000x or more.
97 have exactly one qualifying sample;
244 have qualifying samples from
at least two known patients. Anonymous samples are reported separately.
RA uses the same 10× denominator in reproductive tissue or a supported embryo-lineage mean.
88 proteins whose normal home is a tissue no panel resolves are
Other:Unscreened tissue unless reproductive evidence qualifies.

## 2. The six filter layers

A protein is kept only if it passes all six. `filter_failed` names the first one it
fails. `max_*` knobs mean **below** the value, `min_*` knobs mean **at least** the
value, and a missing measurement never passes.

| Layer | Question | Knobs | Removed by `strict`, alone |
|---|---|---|---|
| 1 `off_target` | How low must expression be outside tumors and reproductive tissue? | `max_off_target_tpm` | 78 |
| 2 `curation` | Is the protein's recorded normal home in a tissue of this tier or above? | `normal_function` | 8 |
| 3 `host` | Do protected compartments have sufficient coverage and bounded exposure? | `host` (HostContext), or explicit profile/tissue choices | 0 |
| 4 `sets` | How much biological enrichment/replication, and which set? | `min_tumor_fold`, `min_reproductive_fold`, `min_tumor_patients`, `membership` | 196 |
| 5 `papers` | How many papers nominate it or report it in tumors? | `min_papers` | 338 |
| 6 `exclusions` | Which gene families do you rule out by name? | `exclude_prefixes`, `except_labels` | 17 |

**Layer 3 tightens the two sets independently.** Both fold knobs raise the catalog's 10x;
a protein that loses Cancer Antigen membership keeps its Reproductive Antigen
membership, so a CTA with a weaker tumor becomes Reproductive Antigen only and stays in
the list. `membership` then picks the set: `membership="cta"` requires both.

## 3. Each knob alone

Target candidates kept with only that knob changed. Categories are the ones after
layer 3 has tightened Cancer Antigens.

| Knob | Value | Kept | CTA | Cancer Antigen only | Reproductive Antigen only |
|---|---|---|---|---|---|
| `max_off_target_tpm` | 5.0 | 468 | 306 | 53 | 109 |
| `max_off_target_tpm` | 2.0 | 372 | 254 | 37 | 81 |
| `max_off_target_tpm` | 1.0 | 253 | 168 | 19 | 66 |
| `max_off_target_tpm` | 0.5 | 149 | 89 | 13 | 47 |
| `max_off_target_tpm` | 0.1 | 32 | 18 | 4 | 10 |
| `normal_function` | no_vital | 519 | 324 | 68 | 127 |
| `normal_function` | no_protected | 515 | 321 | 68 | 126 |
| `normal_function` | no_manageable | 506 | 320 | 67 | 119 |
| `min_tumor_fold` | 20.0 | 490 | 245 | 36 | 209 |
| `min_tumor_fold` | 30.0 | 476 | 201 | 22 | 253 |
| `min_tumor_fold` | 50.0 | 468 | 151 | 14 | 303 |
| `min_tumor_fold` | 100.0 | 462 | 111 | 8 | 343 |
| `min_reproductive_fold` | 20.0 | 491 | 287 | 109 | 95 |
| `min_reproductive_fold` | 30.0 | 477 | 262 | 134 | 81 |
| `min_reproductive_fold` | 50.0 | 457 | 214 | 182 | 61 |
| `min_reproductive_fold` | 100.0 | 439 | 174 | 222 | 43 |
| `min_tumor_patients` | 1 | 506 | 286 | 52 | 168 |
| `min_tumor_patients` | 2 | 486 | 212 | 32 | 242 |
| `min_tumor_patients` | 3 | 481 | 182 | 27 | 272 |
| `min_tumor_patients` | 5 | 475 | 140 | 21 | 314 |
| `membership` | cancer_antigen | 396 | 327 | 69 | 0 |
| `membership` | reproductive_antigen | 454 | 327 | 0 | 127 |
| `membership` | cta | 327 | 327 | 0 | 0 |
| `min_papers` | 1 | 385 | 274 | 30 | 81 |
| `min_papers` | 2 | 271 | 211 | 14 | 46 |
| `min_papers` | 3 | 185 | 148 | 3 | 34 |
| `min_papers` | 4 | 69 | 65 | 1 | 3 |

## 4. Filter modes and host toxicity profiles

`all` → `loose` → `standard` → `strict`: each mode keeps a subset of the one before.
The reference proteins PRAME, NY-ESO-1, MAGEA4 pass through the same thresholds as every other protein.

| Filter mode | Kept | CTA | Cancer Antigen only | Reproductive Antigen only |
|---|---|---|---|---|
| `all` | 523 | 327 | 69 | 127 |
| `loose` | 385 | 274 | 30 | 81 |
| `standard` | 247 | 198 | 13 | 36 |
| `strict` | 115 | 115 | 0 | 0 |

### `loose`

- Sets: keep proteins that are then in either antigen set.
- Papers: at least 1 distinct paper(s) nominating it or reporting it in tumors.

Target candidates **523** → off_target **523** → curation **523** → host **523** → sets **523** → papers **385** → exclusions **385**

### `standard`

- Off target: the highest tissue value outside tumors and reproductive tissue, Disposable and any protected host compartments included, is below 5 TPM.
- Curation: no recorded normal home in a Vital tissue.
- Sets: keep proteins that are then in either antigen set.
- Papers: at least 2 distinct paper(s) nominating it or reporting it in tumors.

Target candidates **523** → off_target **468** → curation **465** → host **465** → sets **465** → papers **247** → exclusions **247**

### `strict`

- Off target: the highest tissue value outside tumors and reproductive tissue, Disposable and any protected host compartments included, is below 4 TPM.
- Curation: no recorded normal home in a Vital or Protected tissue.
- Sets: keep proteins that are then in cta.
- Papers: at least 3 distinct paper(s) nominating it or reporting it in tumors.
- Exclusions: drop labels starting with MAGE, except MAGEA4.

Target candidates **523** → off_target **445** → curation **439** → host **439** → sets **286** → papers **126** → exclusions **115**

### Host anatomy and tissue-tolerance profiles

Host protection is independent of biological membership. Presets assume anatomy and
protect no reproductive tissue by default. Explicit `protect`, `accept_loss` and `absent`
choices have different meanings.
Missing coverage and unresolved protein evidence fail the screen. Native TPM and nTPM
are screened separately; protected TPM also enters every off-target limit and cancer fold.
[Precise contract](../../host-selection.md). [CSV](host-profiles.csv).

| Profile | Protected compartments | CTA | CA only | RA only | CA | RA |
|---|---|---|---|---|---|
| all | none | 327 | 69 | 127 | 523 |
| female-patient | none | 327 | 69 | 127 | 523 |
| female-patient (protect ovary) | ovary | 277 | 66 | 112 | 455 |
| female-patient (lactating) | lactating_breast | 320 | 66 | 123 | 509 |
| female-patient (pregnant) | embryo, placenta, lactating_breast | 0 | 0 | 0 | 0 |
| male-patient | none | 327 | 69 | 127 | 523 |
| male-patient (protect testis) | testis | 54 | 29 | 24 | 107 |
| male-patient (protect adult tissues) | testis, epididymis, prostate, seminal vesicle, resting_breast | 40 | 23 | 11 | 74 |

## 5. Where the familiar CTAs land

46 well-known names. In the catalog, 27 are Target
candidates, 18 are Unsafe and
1 Other. A mode column shows ✓ or the first layer
that removes the protein. **11 pass every mode:** MAGEA4, NY-ESO-1, SSX1, GAGE1, GAGE2A, PAGE1, AKAP4, ADAM2, TPTE, CT45A1, NUTM1.
Flags and the nTPM column are reported for reading; no knob acts on them.

| Name | Catalog | Off-target TPM | in | Off-target nTPM | Tumor fold | CTA papers | Flags | `loose` | `standard` | `strict` |
|---|---|---|---|---|---|---|---|---|---|---|
| MAGEA1 | CTA | 0.2 | spleen | 0.0 | 346 | 7 | literature conflict | ✓ | ✓ | exclusions |
| MAGEA2 | CTA | 0.6 | white matter | 0.8 | 234 | 4 | mapping risk | ✓ | ✓ | exclusions |
| MAGEA4 | CTA | 3.0 | lymph node | 0.6 | 1084 | 8 | literature conflict | ✓ | ✓ | ✓ |
| MAGEA6 | CTA | 0.7 | adipose tissue | 0.2 | 539 | 6 | mapping risk | ✓ | ✓ | exclusions |
| MAGEA10 | CTA | 2.2 | white matter | 2.4 | 232 | 6 | mapping risk | ✓ | ✓ | exclusions |
| MAGEB1 | CTA | 0.5 | lung | 0.0 | 93 | 6 | mapping risk | ✓ | ✓ | exclusions |
| MAGEB4 | CTA | 0.1 | lymph node | 0.0 | 29 | 1 |  | ✓ | papers | papers |
| MAGEC1 | CTA | 1.0 | spleen | 0.6 | 304 | 5 |  | ✓ | ✓ | exclusions |
| MAGEC2 | CTA | 0.2 | pituitary gland | 0.1 | 682 | 6 |  | ✓ | ✓ | exclusions |
| NY-ESO-1 | CTA | 0.8 | spleen | 0.5 | 944 | 5 | mapping risk | ✓ | ✓ | ✓ |
| PRAME | CTA | 10.9 | adrenal gland | 20.3 | 338 | 6 | dismissed rare-cell, literature conflict | ✓ | off_target | off_target |
| SSX1 | CTA | 2.7 | thyroid gland | 1.3 | 101 | 5 | mapping risk, literature conflict | ✓ | ✓ | ✓ |
| SSX2 | CTA | 9.4 | thyroid gland | 2.8 | 36 | 6 | mapping risk, literature conflict | ✓ | off_target | off_target |
| SSX4 | CTA | 0.5 | heart muscle | 0.2 | 413 | 2 | mapping risk | ✓ | ✓ | papers |
| GAGE1 | CTA | 0.7 | retina | 0.5 | 1182 | 6 | mapping risk | ✓ | ✓ | ✓ |
| GAGE2A | CTA | 1.1 | retina | 0.8 | 1448 | 6 | mapping risk | ✓ | ✓ | ✓ |
| PAGE1 | CTA | 0.5 | blood | 0.0 | 1324 | 4 |  | ✓ | ✓ | ✓ |
| XAGE1A/B | CTA | 5.7 | lung | 1.6 | 709 | 2 |  | ✓ | off_target | off_target |
| SYCP1 | CTA | 4.7 | skin | 1.8 | 63 | 4 |  | ✓ | ✓ | off_target |
| AKAP4 | CTA | 0.2 | adrenal gland | 0.0 | 16 | 3 |  | ✓ | ✓ | ✓ |
| ODF1 | CTA | 1.6 | basal ganglia | 3.4 | 26 | 2 |  | ✓ | ✓ | papers |
| ADAM2 | CTA | 0.6 | midbrain | 0.5 | 100 | 4 |  | ✓ | ✓ | ✓ |
| TPTE | CTA | 2.0 | white matter | 2.5 | 61 | 3 | mapping risk | ✓ | ✓ | ✓ |
| CT45A1 | CTA | 3.7 | lung | 0.8 | 107 | 3 | mapping risk | ✓ | ✓ | ✓ |
| CT83 | CTA | 8.6 | salivary gland | 2.6 | 22 | 6 | literature conflict | ✓ | off_target | off_target |
| CSAG1 | CTA | 2.6 | cerebral cortex | 1.9 | 346 | 1 | mapping risk, literature conflict | ✓ | papers | papers |
| NUTM1 | CTA | 0.9 | cerebellum | 1.1 | 56 | 3 |  | ✓ | ✓ | ✓ |

Not candidates in the catalog:

| Name | Catalog | Off-target TPM | Off-target nTPM | Unsafe reasons |
|---|---|---|---|---|
| MAGEA3 | Unsafe |  | 53.1 | HPA tissue clusters |
| MAGEA12 | Unsafe |  | 0.5 | Established rare normal-cell RNA |
| MAGEB2 | Unsafe |  | 118.4 | HPA tissue clusters |
| CTAG2 | Other | 10.7 | 24.9 |  |
| LAGE3 | Unsafe |  | 136.4 | HPA normal protein | Normal bulk | HPA tissue clusters | HPA immune | Monaco immune | DICE immune | Adult cornea |
| PAGE4 | Unsafe |  | 23.3 | Normal bulk |
| TEX14 | Unsafe |  | 28.4 | HPA tissue clusters |
| TDRD1 | Unsafe |  | 10.8 | Established rare normal-cell RNA |
| HORMAD1 | Unsafe |  | 14.8 | HPA tissue clusters | DICE immune | Adult cornea | Established rare normal-cell RNA |
| MAEL | Unsafe |  | 57.0 | Normal bulk | HPA tissue clusters |
| PIWIL2 | Unsafe |  | 16.0 | Normal bulk | HPA tissue clusters | Adult cornea | Established rare normal-cell RNA |
| DDX43 | Unsafe |  | 6.6 | HPA normal protein | HPA tissue clusters | Monaco immune | DICE immune | Established rare normal-cell RNA |
| CTCFL | Unsafe |  | 15.0 | Normal bulk |
| SPA17 | Unsafe |  | 345.8 | HPA normal protein | Normal bulk | HPA tissue clusters | Adult cornea | Established rare normal-cell RNA |
| SPAG9 | Unsafe |  | 480.1 | HPA normal protein | Normal bulk | HPA tissue clusters | HPA immune | Monaco immune | DICE immune | Adult cornea | Established rare normal-cell RNA |
| PBK | Unsafe |  | 157.7 | Normal bulk | HPA tissue clusters | Adult cornea | Established rare normal-cell RNA |
| CEP55 | Unsafe |  | 106.0 | HPA normal protein | Normal bulk | HPA tissue clusters | DICE immune | Adult cornea | Established rare normal-cell RNA |
| PLAC1 | Unsafe |  | 2.4 | Established rare normal-cell RNA |
| LY6K | Unsafe |  | 106.2 | Normal bulk | HPA tissue clusters | Adult cornea | Established rare normal-cell RNA |

## 6. The MAGE family

`strict` drops every label starting with `MAGE` except MAGEA4:
17 of 18 MAGE candidates. **This is a judgment by family name, not
a measurement.** The last column shows what `strict` would decide from each protein's
own data with the name rule switched off: 11 of the 17 it drops
would pass (MAGEA1, MAGEA10, MAGEA11, MAGEA2/B, MAGEA6, MAGEA8, MAGEA9/B, MAGEB1, MAGEB16, MAGEC1, MAGEC2).

| MAGE | Catalog | Off-target TPM | Off-target nTPM | Tumor fold | CTA papers | Mapping risk | `standard` | `strict` on its own data |
|---|---|---|---|---|---|---|---|---|
| MAGEA1 | CTA | 0.2 | 0.0 | 346.0 | 7 |  | ✓ | ✓ |
| MAGEA10 | CTA | 2.2 | 2.4 | 232.0 | 6 | yes | ✓ | ✓ |
| MAGEA11 | CTA | 1.0 | 0.1 | 109.0 | 6 |  | ✓ | ✓ |
| MAGEA2/B | CTA | 0.6 | 0.8 | 234.0 | 4 | yes | ✓ | ✓ |
| MAGEA4 | CTA | 3.0 | 0.6 | 1,084.0 | 8 |  | ✓ | ✓ |
| MAGEA6 | CTA | 0.7 | 0.2 | 539.0 | 6 | yes | ✓ | ✓ |
| MAGEA8 | CTA | 1.0 | 1.0 | 205.0 | 4 |  | ✓ | ✓ |
| MAGEA9/B | CTA | 0.7 | 0.2 | 321.0 | 4 |  | ✓ | ✓ |
| MAGEB1 | CTA | 0.5 | 0.0 | 93.0 | 6 | yes | ✓ | ✓ |
| MAGEB16 | CTA | 0.0 | 0.0 | 90.0 | 4 |  | ✓ | ✓ |
| MAGEB17 | Cancer Antigen only | 7.8 | 2.8 | 19.0 | 0 |  | off_target | off_target |
| MAGEB18 | CTA | 0.0 | 0.0 | 11.0 | 2 |  | ✓ | papers |
| MAGEB3 | Reproductive Antigen only | 0.2 | 0.3 | 8.0 | 3 |  | ✓ | sets |
| MAGEB4 | CTA | 0.1 | 0.0 | 29.0 | 1 |  | papers | papers |
| MAGEB5 | Reproductive Antigen only | 0.2 | 0.0 | 2.0 | 0 |  | papers | sets |
| MAGEB6B | CTA | 0.0 | 0.0 | 15.0 | 0 | yes | papers | papers |
| MAGEC1 | CTA | 1.0 | 0.6 | 304.0 | 5 |  | ✓ | ✓ |
| MAGEC2 | CTA | 0.2 | 0.1 | 682.0 | 6 |  | ✓ | ✓ |

## 7. Disposable tissue and the "below 1 everywhere" list

- **Disposable.** 19 of 523 candidates reach 1 TPM or 1 nTPM in
  appendix, gallbladder or tonsil. The highest is
  NLRP4 at 7.8 TPM.
  [List](disposable-at-least-1.csv).
- **Below 1 everywhere off target.** 230 candidates are below 1 TPM at tissue
  level and below 1 HPA nTPM, single-cell clusters included.
  [List](off-target-below-1.csv). It is a list, not a knob or a catalog category.

## 8. Curation: what the catalog decides and what the filters add

**The catalog and the filters read the same records.** The literature review
(`normal-expression-literature.json`) holds one entry for each of 912
proteins. The tissue-tier registry tiers the tissue. Nothing below lives in a separate
filter file. [Sweep findings and method](normal-function-sweep.md).

The catalog acts on the review's rubric. A gene-specific report in normal adult tissue
tiered Vital or Protected makes the protein Unsafe, and needs a quote retrieved from
the source. The 2026-11-12 candidate sweep adds 12 such proteins:

| Protein | Tissue | Tier | Basis |
|---|---|---|---|
| ADAM30 | brain | Vital | RNA measured in human frontal cortex of non-demented controls by array and probe assay, with neuronal immunostaining. Absolute level is not given; controls are aged post-mortem brains. |
| ASCL4 | skin | Protected | Gene-specific RT-PCR restricted to skin, quantified in adult skin. Read at abstract level; the adult level is one seventh of fetal and its absolute value is not given. |
| BECN2 | brain | Vital | Protein detected in adult human brain lysate by western blot with a peptide antibody. One laboratory; no cell-type resolution. |
| CCN6 | cartilage | Protected | Protein in chondrocytes of normal adult articular cartilage by immunohistochemistry. Read at abstract level; one 2004 study. |
| DEFB106A/B | lung | Vital | PCR on human tissue cDNA panels reports expression in lung, with in situ hybridisation. Read at abstract level; the lung cell type is not named. |
| GALR3 | adrenal gland | Vital | The human cloning paper reports mRNA detectable in adrenal gland and pancreas. Read at abstract level; a 1998 assay whose method is not given there. |
| MRGPRX1 | dorsal root ganglion | Vital | qPCR of whole human dorsal root ganglia from four donors and single-cell RT-PCR of 85 human sensory neurons. Donor ages are not stated; the neurons picked were mostly small. |
| OBP2A | lung | Vital | The cloning paper reports strong expression in lung, nasal structures, salivary and lachrymal glands, and separates it from the 95%-identical OBP2B. Read at abstract level. RNA-seq of four human olfactory samples and a lung sample shows it below 0.1 FPKM, which the rubric does not count against a gene-specific report. |
| SCGB1C1 | nasal mucosa | Protected | Sequence-confirmed RT-PCR finds it in normal inferior turbinate mucosa of people without sinus disease. Gene-level RNA-seq also gives 322 to 3,126 FPKM in all four human olfactory epithelium samples. |
| SRY | brain | Vital | Immunostaining in dopamine neurons of the adult male substantia nigra and absent in female. One laboratory and one antibody; no RNA. |
| SULT6B1 | kidney | Vital | RT-PCR reports expression in human kidney and testis. Read at abstract level; one non-quantitative report, RNA only. |
| TCHHL1 | skin | Protected | Protein in the basal layer of normal human epidermis by immunohistochemistry, with mRNA in normal skin. Read at abstract level; one in-house antibody. |

The catalog's sequence check moves read-through loci to Other. 262 loci in the
catalog have every coding transcript joining two separate genes. 17 are
Other:Sequence; the rest are Unsafe or Other for another reason. Other:Sequence:
CCDC169-SOHLH2, EFNA4-EFNA3, ENSG00000241690, ENSG00000248235, ENSG00000255432, ENSG00000258417, ENSG00000258529, ENSG00000258691, ENSG00000260836, ENSG00000270299, ENSG00000273155, ENSG00000274944, EPPIN-WFDC6, NT5C1B-RDH14, TRIM39-RPP21, TRIM6-TRIM34, ZNF816-ZNF321P.

Layer 2 acts on the recorded normal home of the candidates that remain: a confirmed
finding that falls short of the rubric, a family rule, or the nose and skin reference.

| Recorded tier | Candidates | Tissues | Removed by `normal_function` |
|---|---|---|---|
| Vital | 4 | blood, eye, pancreas, skeletal muscle | `no_vital` and stricter |
| Protected | 4 | blood, lymph node, stomach | `no_protected` and stricter |
| Manageable | 9 | breast, olfactory epithelium, thyroid gland | `no_manageable` |
| Disposable | 2 | hair follicle, skin | never |

Insufficient coverage of a recorded family home blocks CA membership. Members with qualifying reproductive evidence can enter RA. The table counts members that remain candidates.

| Family | Candidates | Tissue | Tier | Used by filters | Sources |
|---|---|---|---|---|---|
| Olfactory receptors | 6 | olfactory epithelium | Manageable | yes | PMID 24800820 |
| Bitter taste receptors | 0 | taste bud | Protected | yes | PMID 18003842 |
| Hair keratins and keratin-associated proteins | 1 | hair follicle | Disposable | yes | PMID 15028290; PMID 14962103; PMID 15144888; PMID 17235325; PMID 15086535; PMID 17301834 |
| Trace-amine receptors | 0 | olfactory epithelium | Manageable | no | UniProt Q9P1P5; UniProt Q96RI9 |
| Sweet and umami taste receptors | 0 | taste bud | Protected | no |  |
| Late cornified envelope proteins | 0 | skin | Protected | no | PMID 19169253; PMID 15854049; PMID 26556599 |

Tissues added to the tissue-tier registry for these findings
(`normal-tissue-tier-updates.yaml`). Olfactory epithelium is Manageable since
2026-11-15, on the oncology precedent ([tissue tiers](../../tissue-tiers.md#how-a-tier-is-chosen)):

| Tissue | Tier | Confidence | References |
|---|---|---|---|
| Olfactory epithelium and olfactory sensory neurons | Manageable | medium | [ctcae5_anosmia](https://www.cismef.org/page/detail/en/NCI_CO_C143230), [amezaga2023](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9916934/), [rt_smell2025](https://www.epistemonikos.org/en/documents/389d1578f48ebdba37920f7c8460c14bd6e3dc1c), [npc_bulb2014](https://pubmed.ncbi.nlm.nih.gov/24746667/), [buttiron2021](https://pubmed.ncbi.nlm.nih.gov/32734392/), [fda_zinc2009](https://webarchive.library.unt.edu/web/20090807210616/http://www.fda.gov/ForConsumers/ConsumerUpdates/ucm166931.htm), [croy2014](https://pubmed.ncbi.nlm.nih.gov/24429163/), [santos2004](https://pubmed.ncbi.nlm.nih.gov/15023839/), [pinto2014](https://pubmed.ncbi.nlm.nih.gov/25271633/), [child2018](https://pubmed.ncbi.nlm.nih.gov/29934351/) |
| Taste buds and taste receptor cells | Protected | medium | [gunn2021](https://pubmed.ncbi.nlm.nih.gov/33545253/), [kumari2017](https://pubmed.ncbi.nlm.nih.gov/29133390/), [kinnamon2012](https://pubmed.ncbi.nlm.nih.gov/21481196/), [stevie_vismodegib](https://academic.oup.com/oncolo/article-abstract/21/10/1218/6401021), [amezaga2023](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9916934/) |
| Tooth, dental pulp and odontoblasts | Manageable | medium | [ng2010](https://pubmed.ncbi.nlm.nih.gov/20158529/), [aminoshariae2024](https://pubmed.ncbi.nlm.nih.gov/38945200/), [gupta2019](https://pubmed.ncbi.nlm.nih.gov/29573048/) |
| Dorsal root and trigeminal ganglion sensory neurons | Vital | high | [sanchosaldana2021](https://pubmed.ncbi.nlm.nih.gov/33491284/), [graus2001](https://pubmed.ncbi.nlm.nih.gov/11353730/), [camdessanche2009](https://pubmed.ncbi.nlm.nih.gov/19506068/) |
| Nasal respiratory mucosa and its glands | Protected | low | [moore2001](https://pubmed.ncbi.nlm.nih.gov/11777241/), [hussain2025](https://pubmed.ncbi.nlm.nih.gov/39126290/), [sumaily2023](https://pubmed.ncbi.nlm.nih.gov/37449389/) |
| Lacrimal gland | Protected | medium | [stevenson2013](https://pubmed.ncbi.nlm.nih.gov/24204119/), [zhang2021](https://pubmed.ncbi.nlm.nih.gov/34805617/), [veernala2022](https://pubmed.ncbi.nlm.nih.gov/35753665/) |
| Bone, osteoblasts and osteocytes | Vital | medium | [tatsumi2007](https://pubmed.ncbi.nlm.nih.gov/17550781/), [weinstein2000](https://pubmed.ncbi.nlm.nih.gov/10946902/) |

Limits: candidate-sweep findings other than the exclusions carry the second reader's paraphrase, not a quote; a candidate-sweep null is one agent's search, not proof of absence; FSCB was not searched; olfactory, taste and hair-keratin family rules rest on family-level sources and were not reviewed gene by gene; tiers for olfactory epithelium, taste buds, nasal mucosa, lacrimal gland, sensory ganglia, bone and tooth come from the 2026-11-13 registry updates, written from abstracts.

## 9. What to keep in mind

- **One tumor sample can admit a protein.** The rule asks for a single sample at 10x,
  as requested; it does not ask for replication. The most frequent top sample is
  `TREEHOUSE_POLYA_25_01:TCGA-2G-AAFL-01`, the best tumor for 22 Cancer
  Antigens. `samples_at_tenfold`, `independent_patients_at_tenfold`,
  `anonymous_samples_at_tenfold` and `cancer_types_at_tenfold` record qualifying support.
  `min_tumor_patients` optionally requires known-patient replication at the selected
  `min_tumor_fold` (10× by default).
- **The paper layer does most of the cutting.** 356 of 523 candidates
  are nominated by any CTA-list paper. `min_papers=3` keeps 185, and
  148 of those are named by all of the three largest lists
  (Chang et al. (2019), Wang et al. (2016), da Silva et al. (2017)).
- **Paper counts are not independent replication.** The nominating sources are
  published catalogs and screens, and their lists overlap. Several lists naming a gene
  shows that several screens picked it, not that separate experiments confirmed it.
  [Sources and recorded scopes](paper-families.csv).
