# Target filters: categories, knobs and where proteins land

**The catalog assigns five categories and they are not changed here. Filters only
remove proteins from the 668 Target candidates, through five layers that
follow the catalog's own shape.** Filter policy `oncoref.target-filters.2026-11-15.1`, applied to
catalog `oncoref.target-selection.2026-11-16.1`. Every number below is computed by
`scripts/build_target_filter_audit.py`. [Rules reference](../../target-filters.md).

These are outputs of screening rules. They are not validated tumor antigens.

## 1. The categories

| Category | Proteins | Rule |
|---|---|---|
| CTA | 437 | Passes the shared screen; cancer evidence and reproductive evidence both qualify |
| Cancer Antigen only | 34 | Passes the shared screen; only cancer evidence qualifies |
| Reproductive Antigen only | 197 | Passes the shared screen; only reproductive evidence qualifies |
| Other | 972 | Not excluded, but fails coverage, identity, expression bounds, sequence, the home-tissue rule or both evidence routes |
| Unsafe | 18,253 | An accepted normal-tissue or clinical exclusion |

Target candidates are the first three rows: 668 proteins. A filter never moves
a protein into or out of Other or Unsafe.

**One cancer rule.** A protein is a Cancer Antigen when one tumor sample, from any
comparable TPM cohort except germ cell tumors, is at least 10 times its highest
non-reproductive normal value (normal tissue, immune cells and Disposable tissue; values
below 1 TPM count as 1). The 471 Cancer Antigens: 121 at 10 to 20x, 138 at 20 to 50x, 68 at 50 to 100x, 126 at 100 to 1000x, 18 at 1000x or more.
83 proteins whose normal home is a tissue no panel resolves are
Other:Tissue unless reproductive evidence qualifies.

## 2. The five filter layers

A protein is kept only if it passes all five. `filter_failed` names the first one it
fails. `max_*` knobs mean **below** the value, `min_*` knobs mean **at least** the
value, and a missing measurement never passes.

| Layer | Question | Knobs | Removed by `strict`, alone |
|---|---|---|---|
| 1 `off_target` | How low must expression be outside tumors and reproductive tissue? | `max_off_target_tpm` | 49 |
| 2 `curation` | Is the protein's recorded normal home in a tissue of this tier or above? | `normal_function` | 14 |
| 3 `sets` | How far above normal must the best tumor be, and which set do you want? | `min_tumor_fold`, `membership` | 231 |
| 4 `papers` | How many papers nominate it or report it in tumors? | `min_papers` | 457 |
| 5 `exclusions` | Which gene families do you rule out by name? | `exclude_prefixes`, `except_labels` | 18 |

**Layer 3 tightens Cancer Antigens only.** `min_tumor_fold` raises the catalog's 10x;
a protein that loses Cancer Antigen membership keeps its Reproductive Antigen
membership, so a CTA with a weaker tumor becomes Reproductive Antigen only and stays in
the list. `membership` then picks the set: `membership="cta"` requires both.

## 3. Each knob alone

Target candidates kept with only that knob changed. Categories are the ones after
layer 3 has tightened Cancer Antigens.

| Knob | Value | Kept | CTA | Cancer Antigen only | Reproductive Antigen only |
|---|---|---|---|---|---|
| `max_off_target_tpm` | 5.0 | 640 | 426 | 31 | 183 |
| `max_off_target_tpm` | 2.0 | 504 | 353 | 23 | 128 |
| `max_off_target_tpm` | 1.0 | 357 | 247 | 13 | 97 |
| `max_off_target_tpm` | 0.5 | 243 | 157 | 8 | 78 |
| `max_off_target_tpm` | 0.1 | 86 | 58 | 2 | 26 |
| `normal_function` | no_vital | 664 | 433 | 34 | 197 |
| `normal_function` | no_protected | 654 | 427 | 34 | 193 |
| `normal_function` | no_manageable | 622 | 409 | 33 | 180 |
| `min_tumor_fold` | 20.0 | 654 | 330 | 20 | 304 |
| `min_tumor_fold` | 30.0 | 650 | 270 | 16 | 364 |
| `min_tumor_fold` | 50.0 | 644 | 202 | 10 | 432 |
| `min_tumor_fold` | 100.0 | 640 | 138 | 6 | 496 |
| `membership` | cancer_antigen | 471 | 437 | 34 | 0 |
| `membership` | reproductive_antigen | 634 | 437 | 0 | 197 |
| `membership` | cta | 437 | 437 | 0 | 0 |
| `min_papers` | 1 | 453 | 333 | 6 | 114 |
| `min_papers` | 2 | 307 | 241 | 1 | 65 |
| `min_papers` | 3 | 211 | 165 | 0 | 46 |
| `min_papers` | 4 | 75 | 71 | 0 | 4 |

## 4. Filter modes

`all` → `loose` → `standard` → `strict`: each mode keeps a subset of the one before.
Every mode keeps the positive controls PRAME, NY-ESO-1, MAGEA4; a test fails
if a change to a mode drops one.

| Filter mode | Kept | CTA | Cancer Antigen only | Reproductive Antigen only |
|---|---|---|---|---|
| `all` | 668 | 437 | 34 | 197 |
| `loose` | 453 | 333 | 6 | 114 |
| `standard` | 296 | 235 | 1 | 60 |
| `strict` | 141 | 141 | 0 | 0 |

### `loose`

- Sets: keep proteins that are then in either antigen set.
- Papers: at least 1 distinct paper(s) nominating it or reporting it in tumors.

Target candidates **668** → off_target **668** → curation **668** → sets **668** → papers **453** → exclusions **453**

### `standard`

- Off target: the highest tissue value outside tumors and reproductive tissue, Disposable tissue included, is below 5 TPM.
- Curation: no recorded normal home in a Vital tissue.
- Sets: keep proteins that are then in either antigen set.
- Papers: at least 2 distinct paper(s) nominating it or reporting it in tumors.

Target candidates **668** → off_target **640** → curation **636** → sets **636** → papers **296** → exclusions **296**

### `strict`

- Off target: the highest tissue value outside tumors and reproductive tissue, Disposable tissue included, is below 4 TPM.
- Curation: no recorded normal home in a Vital or Protected tissue.
- Sets: keep proteins that are then in cta.
- Papers: at least 3 distinct paper(s) nominating it or reporting it in tumors.
- Exclusions: drop labels starting with MAGE, except MAGEA4.

Target candidates **668** → off_target **619** → curation **605** → sets **405** → papers **152** → exclusions **141**

## 5. Where the familiar CTAs land

46 well-known names. In the catalog, 28 are Target
candidates, 18 are Unsafe and
0 Other. A mode column shows ✓ or the first layer
that removes the protein. **16 pass every mode:** MAGEA4, NY-ESO-1, CTAG2, PRAME, SSX1, SSX2, GAGE1, GAGE2A, PAGE1, SYCP1, AKAP4, ADAM2, TPTE, CT45A1, CT83, NUTM1.
Flags and the nTPM column are reported for reading; no knob acts on them.

| Name | Catalog | Off-target TPM | in | Off-target nTPM | Tumor fold | CTA papers | Flags | `loose` | `standard` | `strict` |
|---|---|---|---|---|---|---|---|---|---|---|
| MAGEA1 | CTA | 0.0 |  | 0.0 | 346 | 7 | literature conflict | ✓ | ✓ | exclusions |
| MAGEA2 | CTA | 0.6 | white matter | 0.8 | 234 | 4 | mapping risk | ✓ | ✓ | exclusions |
| MAGEA4 | CTA | 0.7 | lymph node | 0.6 | 3202 | 8 | literature conflict | ✓ | ✓ | ✓ |
| MAGEA6 | CTA | 0.2 | adipose tissue | 0.2 | 539 | 6 | mapping risk | ✓ | ✓ | exclusions |
| MAGEA10 | CTA | 2.2 | white matter | 2.4 | 232 | 6 | mapping risk | ✓ | ✓ | exclusions |
| MAGEB1 | CTA | 0.0 |  | 0.0 | 93 | 6 | mapping risk | ✓ | ✓ | exclusions |
| MAGEB4 | CTA | 0.0 |  | 0.0 | 29 | 1 |  | ✓ | papers | papers |
| MAGEC1 | CTA | 0.5 | white matter | 0.6 | 313 | 5 |  | ✓ | ✓ | exclusions |
| MAGEC2 | CTA | 0.1 | amygdala | 0.1 | 682 | 6 |  | ✓ | ✓ | exclusions |
| NY-ESO-1 | CTA | 0.3 | cerebellum | 0.5 | 944 | 5 | mapping risk | ✓ | ✓ | ✓ |
| CTAG2 | CTA | 3.9 | thymus | 24.9 | 155 | 3 | dismissed rare-cell, mapping risk, literature conflict | ✓ | ✓ | ✓ |
| PRAME | CTA | 3.0 | adrenal gland | 20.3 | 1231 | 6 | dismissed rare-cell, literature conflict | ✓ | ✓ | ✓ |
| SSX1 | CTA | 1.2 | cerebellum | 1.3 | 223 | 5 | mapping risk, literature conflict | ✓ | ✓ | ✓ |
| SSX2 | CTA | 2.4 | hypothalamus | 2.8 | 142 | 6 | mapping risk, literature conflict | ✓ | ✓ | ✓ |
| SSX4 | CTA | 0.2 | cerebral cortex | 0.2 | 413 | 2 | mapping risk | ✓ | ✓ | papers |
| GAGE1 | CTA | 0.7 | retina | 0.5 | 1182 | 6 | mapping risk | ✓ | ✓ | ✓ |
| GAGE2A | CTA | 1.1 | retina | 0.8 | 1448 | 6 | mapping risk | ✓ | ✓ | ✓ |
| PAGE1 | CTA | 0.0 |  | 0.0 | 1324 | 4 |  | ✓ | ✓ | ✓ |
| XAGE1A/B | CTA | 1.8 | lung | 1.6 | 2260 | 2 |  | ✓ | ✓ | papers |
| SYCP1 | CTA | 1.9 | blood | 1.8 | 153 | 4 |  | ✓ | ✓ | ✓ |
| AKAP4 | CTA | 0.1 | blood | 0.0 | 16 | 3 |  | ✓ | ✓ | ✓ |
| ODF1 | CTA | 0.6 | basal ganglia | 3.4 | 31 | 2 |  | ✓ | ✓ | papers |
| ADAM2 | CTA | 0.6 | midbrain | 0.5 | 100 | 4 |  | ✓ | ✓ | ✓ |
| TPTE | CTA | 2.0 | white matter | 2.5 | 61 | 3 | mapping risk | ✓ | ✓ | ✓ |
| CT45A1 | CTA | 0.9 | retina | 0.8 | 397 | 3 | mapping risk | ✓ | ✓ | ✓ |
| CT83 | CTA | 1.7 | salivary gland | 2.6 | 114 | 5 | literature conflict | ✓ | ✓ | ✓ |
| CSAG1 | CTA | 1.4 | basal ganglia | 1.9 | 659 | 1 | mapping risk, literature conflict | ✓ | papers | papers |
| NUTM1 | CTA | 0.9 | cerebellum | 1.1 | 56 | 3 |  | ✓ | ✓ | ✓ |

Not candidates in the catalog:

| Name | Catalog | Off-target TPM | Off-target nTPM | Unsafe reasons |
|---|---|---|---|---|
| MAGEA3 | Unsafe |  | 53.1 | HPA tissue clusters |
| MAGEA12 | Unsafe |  | 0.5 | Established rare normal-cell RNA |
| MAGEB2 | Unsafe |  | 118.4 | HPA tissue clusters |
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
18 of 19 MAGE candidates. **This is a judgment by family name, not
a measurement.** The last column shows what `strict` would decide from each protein's
own data with the name rule switched off: 11 of the 18 it drops
would pass (MAGEA1, MAGEA10, MAGEA11, MAGEA2/B, MAGEA6, MAGEA8, MAGEA9/B, MAGEB1, MAGEB16, MAGEC1, MAGEC2).

| MAGE | Catalog | Off-target TPM | Off-target nTPM | Tumor fold | CTA papers | Mapping risk | `standard` | `strict` on its own data |
|---|---|---|---|---|---|---|---|---|
| MAGEA1 | CTA | 0.0 | 0.0 | 346.0 | 7 |  | ✓ | ✓ |
| MAGEA10 | CTA | 2.2 | 2.4 | 232.0 | 6 | yes | ✓ | ✓ |
| MAGEA11 | CTA | 0.1 | 0.1 | 109.0 | 6 |  | ✓ | ✓ |
| MAGEA2/B | CTA | 0.6 | 0.8 | 234.0 | 4 | yes | ✓ | ✓ |
| MAGEA4 | CTA | 0.7 | 0.6 | 3,202.0 | 8 |  | ✓ | ✓ |
| MAGEA6 | CTA | 0.2 | 0.2 | 539.0 | 6 | yes | ✓ | ✓ |
| MAGEA8 | CTA | 1.0 | 1.0 | 205.0 | 4 |  | ✓ | ✓ |
| MAGEA9/B | CTA | 0.2 | 0.2 | 321.0 | 4 |  | ✓ | ✓ |
| MAGEB1 | CTA | 0.0 | 0.0 | 93.0 | 6 | yes | ✓ | ✓ |
| MAGEB10 | Reproductive Antigen only | 1.9 | 2.4 | 4.0 | 2 | yes | ✓ | sets |
| MAGEB16 | CTA | 0.0 | 0.0 | 90.0 | 4 |  | ✓ | ✓ |
| MAGEB17 | CTA | 1.9 | 2.8 | 80.0 | 0 |  | papers | papers |
| MAGEB18 | CTA | 0.0 | 0.0 | 11.0 | 2 |  | ✓ | papers |
| MAGEB3 | Reproductive Antigen only | 0.2 | 0.3 | 8.0 | 3 |  | ✓ | sets |
| MAGEB4 | CTA | 0.0 | 0.0 | 29.0 | 1 |  | papers | papers |
| MAGEB5 | Reproductive Antigen only | 0.0 | 0.0 | 2.0 | 0 |  | papers | sets |
| MAGEB6B | Cancer Antigen only | 0.0 | 0.0 | 15.0 | 0 | yes | papers | sets |
| MAGEC1 | CTA | 0.5 | 0.6 | 313.0 | 5 |  | ✓ | ✓ |
| MAGEC2 | CTA | 0.1 | 0.1 | 682.0 | 6 |  | ✓ | ✓ |

## 7. Disposable tissue and the "below 1 everywhere" list

- **Disposable.** 8 of 668 candidates reach 1 TPM or 1 nTPM in
  appendix, gallbladder or tonsil. The highest is
  LCE3C at 21.4 TPM.
  [List](disposable-at-least-1.csv).
- **Below 1 everywhere off target.** 318 candidates are below 1 TPM at tissue
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
catalog have every coding transcript joining two separate genes. 30 are
Other:Sequence; the rest are Unsafe or Other for another reason. Other:Sequence:
ATF7-NPFF, CCDC169-SOHLH2, EFNA4-EFNA3, ENSG00000198211, ENSG00000241690, ENSG00000248235, ENSG00000249773, ENSG00000249967, ENSG00000255432, ENSG00000255835, ENSG00000256349, ENSG00000256825, ENSG00000258417, ENSG00000258529, ENSG00000258691, ENSG00000260272, ENSG00000260836, ENSG00000266202, ENSG00000267618, ENSG00000270299, ENSG00000273155, ENSG00000274944, ENSG00000279983, EPPIN-WFDC6, NT5C1B-RDH14, PRR5-ARHGAP8, TRIM39-RPP21, TRIM6-TRIM34, URGCP-MRPS24, ZNF816-ZNF321P.

Layer 2 acts on the recorded normal home of the candidates that remain: a confirmed
finding that falls short of the rubric, a family rule, or the nose and skin reference.

| Recorded tier | Candidates | Tissues | Removed by `normal_function` |
|---|---|---|---|
| Vital | 4 | blood, eye, pancreas, skeletal muscle | `no_vital` and stricter |
| Protected | 10 | blood, lymph node, nasal mucosa, skin, stomach, taste bud | `no_protected` and stricter |
| Manageable | 32 | adipose tissue, breast, olfactory epithelium, thyroid gland | `no_manageable` |
| Disposable | 10 | hair follicle, skin | never |

Family rules. A member whose family's home is a tissue no panel resolves, and that has
no reproductive evidence, is Other:Tissue in the catalog; the rows count the members
that remain candidates:

| Family | Candidates | Tissue | Tier | Used by filters | Sources |
|---|---|---|---|---|---|
| Olfactory receptors | 27 | olfactory epithelium | Manageable | yes | PMID 24800820 |
| Bitter taste receptors | 1 | taste bud | Protected | yes | PMID 18003842 |
| Hair keratins and keratin-associated proteins | 9 | hair follicle | Disposable | yes | PMID 15028290; PMID 14962103; PMID 15144888; PMID 17235325; PMID 15086535; PMID 17301834 |
| Trace-amine receptors | 0 | olfactory epithelium | Manageable | no | UniProt Q9P1P5; UniProt Q96RI9 |
| Sweet and umami taste receptors | 0 | taste bud | Protected | no |  |
| Late cornified envelope proteins | 1 | skin | Protected | no | PMID 19169253; PMID 15854049; PMID 26556599 |

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
  Antigens. `tumor_top_sample` and `independent_patients_above_normal` on each row show
  how much rests on one sample.
- **The paper layer does most of the cutting.** 413 of 668 candidates
  are nominated by any CTA-list paper. `min_papers=3` keeps 211, and
  170 of those are named by all of the three largest lists
  (Wang et al. (2016), Chang et al. (2019), da Silva et al. (2017)).
- **Paper counts are not independent replication.** The nominating sources are
  published catalogs and screens, and their lists overlap. Several lists naming a gene
  shows that several screens picked it, not that separate experiments confirmed it.
  [Sources and recorded scopes](paper-families.csv).
