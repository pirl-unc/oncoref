# Target filters: categories, knobs and where proteins land

**The catalog assigns five categories and they are not changed here. Filters only
remove proteins from the 478 Target candidates, through six layers that
follow the catalog's own shape.** Filter policy `oncoref.target-filters.2026-11-17.5`, applied to
catalog `oncoref.target-selection.2026-11-18.2`. Every number below is computed by
`scripts/build_target_filter_audit.py`. [Rules reference](../../target-filters.md).

These are outputs of screening rules. They are not validated tumor antigens.

## 1. The categories

| Category | Proteins | Rule |
|---|---|---|
| CTA | 203 | Passes the shared screen; cancer evidence and reproductive evidence both qualify |
| Cancer Antigen only | 64 | Passes the shared screen; only cancer evidence qualifies |
| Reproductive Antigen only | 211 | Passes the shared screen; only reproductive evidence qualifies |
| Other | 1,162 | Not excluded, but fails coverage, identity, expression bounds, sequence, the home-tissue rule or both evidence routes |
| Unsafe | 18,253 | An accepted normal-tissue or clinical exclusion |

Target candidates are the first three rows: 478 proteins. A filter never moves
a protein into or out of Other or Unsafe.

**One cancer rule.** A protein is a Cancer Antigen when one tumor sample, from any
QC-passing comparable TPM cohort including germ cell tumors, is at least 10 times its highest
non-reproductive normal value (normal tissue, immune cells and Disposable tissue; values
below 1 TPM count as 1). The 267 Cancer Antigens: 83 at 10 to 20x, 72 at 20 to 50x, 35 at 50 to 100x, 70 at 100 to 1000x, 7 at 1000x or more.
73 have exactly one qualifying sample;
169 have qualifying samples from
at least two known patients. Anonymous samples are reported separately.
RA uses the same 10× denominator in reproductive tissue or a supported embryo-lineage mean.
63 proteins whose normal home is a tissue no panel resolves are
Other:Unscreened tissue unless reproductive evidence qualifies.

## 2. The six filter layers

A protein is kept only if it passes all six. `filter_failed` names the first one it
fails. `max_*` knobs mean **below** the value, `min_*` knobs mean **at least** the
value, and a missing measurement never passes.

| Layer | Question | Knobs | Removed by `strict`, alone |
|---|---|---|---|
| 1 `off_target` | How low must expression be outside tumors and reproductive tissue? | `max_off_target_tpm` | 212 |
| 2 `curation` | Is the protein's recorded normal home in a tissue of this tier or above? | `normal_function` | 7 |
| 3 `host` | Do protected compartments have sufficient coverage and bounded exposure? | `host` (HostContext), or explicit profile/tissue choices | 0 |
| 4 `sets` | How much biological enrichment/replication, and which set? | `min_tumor_fold`, `min_reproductive_fold`, `min_tumor_patients`, `membership` | 275 |
| 5 `papers` | How many papers nominate it or report it in tumors? | `min_papers` | 306 |
| 6 `exclusions` | Which gene families do you rule out by name? | `exclude_prefixes`, `except_labels` | 15 |

**Layer 3 tightens the two sets independently.** Both fold knobs raise the catalog's 10x;
a protein that loses Cancer Antigen membership keeps its Reproductive Antigen
membership, so a CTA with a weaker tumor becomes Reproductive Antigen only and stays in
the list. `membership` then picks the set: `membership="cta"` requires both.

## 3. Each knob alone

Target candidates kept with only that knob changed. Categories are the ones after
layer 3 has tightened Cancer Antigens.

| Knob | Value | Kept | CTA | Cancer Antigen only | Reproductive Antigen only |
|---|---|---|---|---|---|
| `max_off_target_tpm` | 5.0 | 290 | 164 | 39 | 87 |
| `max_off_target_tpm` | 2.0 | 175 | 108 | 21 | 46 |
| `max_off_target_tpm` | 1.0 | 82 | 52 | 11 | 19 |
| `max_off_target_tpm` | 0.5 | 41 | 24 | 9 | 8 |
| `max_off_target_tpm` | 0.1 | 2 | 1 | 1 | 0 |
| `normal_function` | no_vital | 475 | 202 | 63 | 210 |
| `normal_function` | no_protected | 471 | 200 | 63 | 208 |
| `normal_function` | no_manageable | 465 | 197 | 63 | 205 |
| `min_tumor_fold` | 20.0 | 449 | 149 | 35 | 265 |
| `min_tumor_fold` | 30.0 | 437 | 133 | 23 | 281 |
| `min_tumor_fold` | 50.0 | 427 | 99 | 13 | 315 |
| `min_tumor_fold` | 100.0 | 421 | 70 | 7 | 344 |
| `min_reproductive_fold` | 20.0 | 403 | 150 | 117 | 136 |
| `min_reproductive_fold` | 30.0 | 318 | 114 | 153 | 51 |
| `min_reproductive_fold` | 50.0 | 297 | 86 | 181 | 30 |
| `min_reproductive_fold` | 100.0 | 288 | 67 | 200 | 21 |
| `min_tumor_patients` | 1 | 467 | 175 | 53 | 239 |
| `min_tumor_patients` | 2 | 449 | 134 | 35 | 280 |
| `min_tumor_patients` | 3 | 445 | 117 | 31 | 297 |
| `min_tumor_patients` | 5 | 432 | 100 | 18 | 314 |
| `membership` | cancer_antigen | 267 | 203 | 64 | 0 |
| `membership` | reproductive_antigen | 414 | 203 | 0 | 211 |
| `membership` | cta | 203 | 203 | 0 | 0 |
| `min_papers` | 1 | 361 | 163 | 35 | 163 |
| `min_papers` | 2 | 254 | 120 | 20 | 114 |
| `min_papers` | 3 | 172 | 81 | 9 | 82 |
| `min_papers` | 4 | 67 | 47 | 6 | 14 |

## 4. Filter modes and host toxicity profiles

`all` → `loose` → `standard` → `strict`: each mode keeps a subset of the one before.
The reference proteins PRAME, NY-ESO-1, MAGEA4 pass through the same thresholds as every other protein.

| Filter mode | Kept | CTA | Cancer Antigen only | Reproductive Antigen only |
|---|---|---|---|---|
| `all` | 478 | 203 | 64 | 211 |
| `loose` | 361 | 163 | 35 | 163 |
| `standard` | 134 | 96 | 12 | 26 |
| `strict` | 47 | 47 | 0 | 0 |

### `loose`

- Sets: keep proteins that are then in either antigen set.
- Papers: at least 1 distinct paper(s) nominating it or reporting it in tumors.

Target candidates **478** → off_target **478** → curation **478** → host **478** → sets **478** → papers **361** → exclusions **361**

### `standard`

- Off target: the highest tissue value outside tumors and reproductive tissue, Disposable and any protected host compartments included, is below 5 TPM.
- Curation: no recorded normal home in a Vital tissue.
- Sets: keep proteins that are then in either antigen set.
- Papers: at least 2 distinct paper(s) nominating it or reporting it in tumors.

Target candidates **478** → off_target **290** → curation **288** → host **288** → sets **288** → papers **134** → exclusions **134**

### `strict`

- Off target: the highest tissue value outside tumors and reproductive tissue, Disposable and any protected host compartments included, is below 4 TPM.
- Curation: no recorded normal home in a Vital or Protected tissue.
- Sets: keep proteins that are then in cta.
- Papers: at least 3 distinct paper(s) nominating it or reporting it in tumors.
- Exclusions: drop labels starting with MAGE, except MAGEA4.

Target candidates **478** → off_target **266** → curation **262** → host **262** → sets **148** → papers **54** → exclusions **47**

### Host anatomy and tissue-tolerance profiles

Host protection is independent of biological membership. Presets assume anatomy and
protect no reproductive tissue by default. Explicit `protect`, `accept_loss` and `absent`
choices have different meanings.
Missing coverage and unresolved protein evidence fail the screen. Native TPM and nTPM
are screened separately; protected TPM also enters every off-target limit and cancer fold.
[Precise contract](../../host-selection.md). [CSV](host-profiles.csv).

| Profile | Protected compartments | CTA | CA only | RA only | CA | RA |
|---|---|---|---|---|---|
| all | none | 203 | 64 | 211 | 478 |
| female-patient | none | 203 | 64 | 211 | 478 |
| female-patient (protect ovary) | ovary | 174 | 60 | 183 | 417 |
| female-patient (lactating) | lactating_breast | 199 | 63 | 205 | 467 |
| female-patient (pregnant) | embryo, placenta, lactating_breast | 0 | 0 | 0 | 0 |
| male-patient | none | 203 | 64 | 211 | 478 |
| male-patient (protect testis) | testis | 43 | 28 | 27 | 98 |
| male-patient (protect adult tissues) | testis, epididymis, prostate, seminal vesicle, resting_breast | 35 | 21 | 13 | 69 |

## 5. Where the familiar CTAs land

46 well-known names. In the catalog, 26 are Target
candidates, 18 are Unsafe and
2 Other. A mode column shows ✓ or the first layer
that removes the protein. **2 pass every mode:** GAGE1, PAGE1.
Flags and the nTPM column are reported for reading; no knob acts on them.

| Name | Catalog | Off-target TPM | in | Off-target nTPM | Tumor fold | CTA papers | Flags | `loose` | `standard` | `strict` |
|---|---|---|---|---|---|---|---|---|---|---|
| MAGEA1 | CTA | 3.0 | adipose tissue | 0.0 | 115 | 7 | literature conflict | ✓ | ✓ | exclusions |
| MAGEA2 | CTA | 0.6 | adipose tissue | 0.8 | 234 | 4 | mapping risk | ✓ | ✓ | exclusions |
| MAGEA4 | CTA | 12.8 | pituitary gland | 0.6 | 250 | 8 | literature conflict | ✓ | off_target | off_target |
| MAGEA6 | CTA | 4.1 | adipose tissue | 0.2 | 131 | 6 | mapping risk | ✓ | ✓ | off_target |
| MAGEA10 | CTA | 2.2 | white matter | 2.4 | 232 | 6 | mapping risk | ✓ | ✓ | exclusions |
| MAGEB1 | CTA | 2.6 | esophagus | 0.0 | 36 | 6 | mapping risk | ✓ | ✓ | exclusions |
| MAGEB4 | CTA | 1.4 | esophagus | 0.0 | 21 | 1 |  | ✓ | papers | papers |
| MAGEC1 | Cancer Antigen only | 16.7 | spleen | 0.6 | 19 | 5 |  | ✓ | off_target | off_target |
| MAGEC2 | CTA | 8.1 | adipose tissue | 0.1 | 84 | 6 |  | ✓ | off_target | off_target |
| NY-ESO-1 | Cancer Antigen only | 3.2 | blood vessel | 0.5 | 295 | 5 | mapping risk | ✓ | ✓ | sets |
| PRAME | Cancer Antigen only | 61.9 | adipose tissue | 20.3 | 60 | 6 | dismissed rare-cell, literature conflict | ✓ | off_target | off_target |
| SSX1 | Cancer Antigen only | 25.1 | thyroid gland | 1.3 | 11 | 5 | mapping risk, literature conflict | ✓ | off_target | off_target |
| SSX2 | CTA | 10.0 | thyroid gland | 2.8 | 34 | 6 | mapping risk, literature conflict | ✓ | off_target | off_target |
| SSX4 | CTA | 0.6 | heart muscle | 0.2 | 413 | 2 | mapping risk | ✓ | ✓ | papers |
| GAGE1 | CTA | 1.2 | adipose tissue | 0.5 | 1014 | 6 | mapping risk | ✓ | ✓ | ✓ |
| GAGE2A | CTA | 4.0 | thyroid gland | 0.8 | 395 | 6 | mapping risk | ✓ | ✓ | off_target |
| PAGE1 | CTA | 3.5 | peripheral nerve | 0.0 | 375 | 4 |  | ✓ | ✓ | ✓ |
| XAGE1A/B | CTA | 8.8 | liver | 1.6 | 462 | 2 |  | ✓ | off_target | off_target |
| SYCP1 | CTA | 5.3 | bone marrow | 1.8 | 56 | 4 |  | ✓ | off_target | off_target |
| AKAP4 | Reproductive Antigen only | 84.9 | esophagus | 0.0 | 0 | 3 |  | ✓ | off_target | off_target |
| ODF1 | Reproductive Antigen only | 99.3 | esophagus | 3.4 | 0 | 2 |  | ✓ | off_target | off_target |
| ADAM2 | CTA | 7.0 | esophagus | 0.5 | 14 | 4 |  | ✓ | off_target | off_target |
| TPTE | Reproductive Antigen only | 13.0 | esophagus | 2.5 | 9 | 3 | mapping risk | ✓ | off_target | off_target |
| CT45A1 | Cancer Antigen only | 6.3 | lung | 0.8 | 63 | 3 | mapping risk | ✓ | off_target | off_target |
| CSAG1 | CTA | 8.2 | midbrain | 1.9 | 108 | 1 | mapping risk, literature conflict | ✓ | off_target | off_target |
| NUTM1 | CTA | 4.1 | thyroid gland | 1.1 | 14 | 3 |  | ✓ | ✓ | off_target |

Not candidates in the catalog:

| Name | Catalog | Off-target TPM | Off-target nTPM | Unsafe reasons |
|---|---|---|---|---|
| MAGEA3 | Unsafe |  | 53.1 | HPA tissue clusters |
| MAGEA12 | Unsafe |  | 0.5 | Established rare normal-cell RNA |
| MAGEB2 | Unsafe |  | 118.4 | HPA tissue clusters |
| CTAG2 | Other | 194.7 | 24.9 |  |
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
| CT83 | Other | 77.0 | 2.6 |  |
| PLAC1 | Unsafe |  | 2.4 | Established rare normal-cell RNA |
| LY6K | Unsafe |  | 106.2 | Normal bulk | HPA tissue clusters | Adult cornea | Established rare normal-cell RNA |

## 6. The MAGE family

`strict` drops every label starting with `MAGE` except MAGEA4:
15 of 16 MAGE candidates. **This is a judgment by family name, not
a measurement.** The last column shows what `strict` would decide from each protein's
own data with the name rule switched off: 7 of the 15 it drops
would pass (MAGEA1, MAGEA10, MAGEA2/B, MAGEA8, MAGEA9/B, MAGEB1, MAGEB16).

| MAGE | Catalog | Off-target TPM | Off-target nTPM | Tumor fold | CTA papers | Mapping risk | `standard` | `strict` on its own data |
|---|---|---|---|---|---|---|---|---|
| MAGEA1 | CTA | 3.0 | 0.0 | 115.0 | 7 |  | ✓ | ✓ |
| MAGEA10 | CTA | 2.2 | 2.4 | 232.0 | 6 | yes | ✓ | ✓ |
| MAGEA2/B | CTA | 0.6 | 0.8 | 234.0 | 4 | yes | ✓ | ✓ |
| MAGEA4 | CTA | 12.8 | 0.6 | 250.0 | 8 |  | off_target | off_target |
| MAGEA6 | CTA | 4.1 | 0.2 | 131.0 | 6 | yes | ✓ | off_target |
| MAGEA8 | CTA | 1.6 | 1.0 | 125.0 | 4 |  | ✓ | ✓ |
| MAGEA9/B | CTA | 1.8 | 0.2 | 177.0 | 4 |  | ✓ | ✓ |
| MAGEB1 | CTA | 2.6 | 0.0 | 36.0 | 6 | yes | ✓ | ✓ |
| MAGEB16 | CTA | 0.4 | 0.0 | 90.0 | 4 |  | ✓ | ✓ |
| MAGEB18 | CTA | 0.2 | 0.0 | 11.0 | 2 |  | ✓ | papers |
| MAGEB3 | Reproductive Antigen only | 0.7 | 0.3 | 8.0 | 3 |  | ✓ | sets |
| MAGEB4 | CTA | 1.4 | 0.0 | 21.0 | 1 |  | papers | papers |
| MAGEB5 | Reproductive Antigen only | 0.5 | 0.0 | 2.0 | 0 |  | papers | sets |
| MAGEB6B | CTA | 0.4 | 0.0 | 15.0 | 0 | yes | papers | papers |
| MAGEC1 | Cancer Antigen only | 16.7 | 0.6 | 19.0 | 5 |  | off_target | off_target |
| MAGEC2 | CTA | 8.1 | 0.1 | 84.0 | 6 |  | off_target | off_target |

## 7. Disposable tissue and the "below 1 everywhere" list

- **Disposable.** 15 of 478 candidates reach 1 TPM or 1 nTPM in
  appendix, gallbladder or tonsil. The highest is
  NLRP4 at 8.3 TPM.
  [List](disposable-at-least-1.csv).
- **Below 1 everywhere off target.** 80 candidates are below 1 TPM at tissue
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
| Vital | 3 | blood, eye, pancreas | `no_vital` and stricter |
| Protected | 4 | blood, lymph node, stomach | `no_protected` and stricter |
| Manageable | 6 | breast, olfactory epithelium | `no_manageable` |
| Disposable | 1 | skin | never |

Family rules. A member whose family's home is a tissue no panel resolves, and that has
no reproductive evidence, is Other:Unscreened tissue in the catalog; the rows count the members
that remain candidates:

| Family | Candidates | Tissue | Tier | Used by filters | Sources |
|---|---|---|---|---|---|
| Olfactory receptors | 4 | olfactory epithelium | Manageable | yes | PMID 24800820 |
| Bitter taste receptors | 0 | taste bud | Protected | yes | PMID 18003842 |
| Hair keratins and keratin-associated proteins | 0 | hair follicle | Disposable | yes | PMID 15028290; PMID 14962103; PMID 15144888; PMID 17235325; PMID 15086535; PMID 17301834 |
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
  `TREEHOUSE_POLYA_25_01:THR24_2154_S01`, the best tumor for 19 Cancer
  Antigens. `samples_at_tenfold`, `independent_patients_at_tenfold`,
  `anonymous_samples_at_tenfold` and `cancer_types_at_tenfold` record qualifying support.
  `min_tumor_patients` optionally requires known-patient replication at the selected
  `min_tumor_fold` (10× by default).
- **The paper layer does most of the cutting.** 335 of 478 candidates
  are nominated by any CTA-list paper. `min_papers=3` keeps 172, and
  137 of those are named by all of the three largest lists
  (Chang et al. (2019), Wang et al. (2016), da Silva et al. (2017)).
- **Paper counts are not independent replication.** The nominating sources are
  published catalogs and screens, and their lists overlap. Several lists naming a gene
  shows that several screens picked it, not that separate experiments confirmed it.
  [Sources and recorded scopes](paper-families.csv).
