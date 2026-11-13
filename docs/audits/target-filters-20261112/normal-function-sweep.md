# Normal-tissue function the screen may miss: a literature sweep of all 902 candidates

**The largest gap is tissue coverage, not thresholds. 242 of the 902 Target candidates
belong to gene families whose normal site of action is a tissue the expression panels
do not sample: nose, tongue, hair follicle and mucosal surfaces.**

This is a screening aid produced by automated literature search. No category was
changed. Every finding needs a human read of its source before it is used.

## How it was done

37 research agents each took 25 candidates and searched the literature and curated
databases for normal adult non-reproductive expression or function. A second agent
then tried to refute each flagged finding against its cited source. 901 of 902 labels
were covered (FSCB was skipped). 343 had essentially no literature.

| Outcome after the refutation pass | Findings |
|---|---:|
| Confirmed: human, adult, normal tissue, gene-specific assay | 53 (51 proteins) |
| Confirmed with a limit (not adult human, or not gene-specific) | 11 |
| Weakened: mouse only, fetal, family-level or paralog-ambiguous | 184 |
| Refuted | 8 |

[Every finding, verdict, source and reason](normal-function-sweep.csv).

## 1. Whole families sit in tissues the panels do not sample

| Family | Candidates | Normal site | In `standard` | In `strict` |
|---|---:|---|---:|---:|
| Olfactory receptors (OR) | 164 | Olfactory epithelium | 7 | 1 |
| Hair keratins and keratin-associated proteins (KRT3x, KRTAP) | 41 | Hair follicle | 0 | 0 |
| Beta-defensins (DEFB) | 22 | Epididymis; some in airway, gut and oral mucosa, some inducible | 2 | 1 |
| Taste receptors (TAS1R, TAS2R) | 9 | Taste buds; reported in gut and airway | 0 | 0 |
| Trace-amine and vomeronasal receptors (TAAR, VN1R) | 4 | Olfactory epithelium | 0 | 0 |
| Late cornified envelope (LCE) | 2 | Epidermis | 0 | 0 |

Most are Cancer Antigen only. They pass because bulk RNA of these tissues is absent
from the panels, so their normal maximum is near zero and a low tumor value clears it.
The paper layer already removes nearly all of them: `standard` keeps 9 of the 242.
The catalog itself does not.

## 2. Confirmed findings, by where the protein acts

| Where | Proteins |
|---|---|
| Cartilage, bone and muscle | CCN6, LRRC52 |
| Endocrine and metabolic | DEFB118, GALR3, GFRA4, KLF14, LALBA |
| Hair follicle | KRT37, KRT39, KRTAP10-6, KRTAP10-8, KRTAP12-4, KRTAP24-1, KRTAP29-1, KRTAP5-3, KRTAP5-4, KRTAP5-5, KRTAP5-7, KRTAP5-9 |
| Immune and blood | DEFB104A/B, DEFB130A/B, NLRP11, PRKACG, SCGB1D4, ZNRF4 |
| Kidney, gut, lung and other | DEFB106A/B, SPANXN3, SULT6B1 |
| Nose and olfactory | BPIFB3, OBP2A, SCGB1C1 |
| Sensory ganglia and brain | ADAM30, BECN2, GALP, MRGPRX1, SRY |
| Skin and sweat gland | AADACL3, ASCL4, OR51A7, PRAME, TCHHL1 |
| Taste and tongue | BPIFB6, RTP2, TAS2R1, TAS2R13, TAS2R41, TAS2R50, TAS2R60, TAS2R7 |
| Tooth | DSPP, PANX3 |

## 3. Confirmed findings that survive `standard`

| Protein | Where | Finding | Source | In `strict` |
|---|---|---|---|---|
| LRRC52 | skeletal muscle | Human: UniProt curates LRRC52 tissue specificity as 'mainly expressed in testis and skeletal muscle', citing Yan and Aldrich 2012 (PMID 22547800). The abstract of that paper, which I read, names only testis, so the skeletal muscle statement rests on the UniProt annotation of the paper's tissue panel. | PMID 22547800, full text PMC3356614, Fig. 4 and Discussion | yes |
| NLRP11 | blood / immune system | In human macrophages NLRP11 binds ASC and NLRP3 and is required for NLRP3 inflammasome assembly; NLRP11 deletion blocks NLRP3 activation and siRNA knockdown in primary human macrophages reduces IL-1beta release (Nat Immunol 2022). The gene is absent in mice. | PMID 35624206 (PMC9174058; DOI 10.1038/s41590-022-01220-3) | yes |
| PRKACG | blood / bone marrow | Human: UniProt P22612 curates germline PRKACG variants (Ile74Met) as the cause of Bleeding disorder, platelet-type, 19 (MIM 616176), a macrothrombocytopenia with platelet dysfunction, citing Manchev et al., Blood 2014. This implies a functional role in the megakaryocyte/platelet lineage, although UniProt's tissue-specificity not | Manchev VT et al., Blood 2014;124(16):2554-63, PMID 25061177, PMCID PMC4199957 ( | yes |
| ZNRF4 | blood / innate immune system | ZNRF4 is an ER-anchored E3 ubiquitin ligase that degrades RIP2 and negatively regulates NOD2 signalling: Znrf4-knockdown monocytes showed sustained NF-kB activation, and Znrf4-knockdown mice had reduced NOD2-induced tolerance to MDP and cleared Listeria monocytogenes better. The abstract does not state the species of the monocyt | PMID 28656966 (full text PMC5493756) | yes |
| SPANXN3 | several nongametogenic adult tissues (not named in abstract) | Human: SPANX-N family reported as expressed at a low level in several nongametogenic adult tissues, in contrast to SPANX-A/D. Family-level statement from the abstract; member-specific tissues not verified. | PMID 17406683 full text, Table 2 (PMC1831492) | no |
| ADAM30 | brain | Human: Letronne et al. 2016 (EBioMedicine, PMID 27333034) report from brain transcriptomics that ADAM30 expression is about 50% lower in Alzheimer's disease brains and inversely correlates with amyloid load, which implies measurable ADAM30 mRNA in non-diseased adult brain. The functional work (cathepsin D activation, APP lysosom | https://pmc.ncbi.nlm.nih.gov/articles/PMC4972530/ (PMID 27333034) | no |
| PRAME | skin | Human: nuclear PRAME staining in non-neoplastic skin is usually absent or limited to a few scattered melanocytes, most often in actinically damaged skin; the authors also note occasional weak nuclear staining of stromal cells in scars. This fits the catalog's existing dismissed rare-cell hit in eye melanocytes and hints at damag | PMID 34023098 (https://pmc.ncbi.nlm.nih.gov/articles/PMC8152939/) | yes |
| PANX3 | tooth (dental pulp) | Human dental pulp tissue and odontoblast-like cells express PANX3 mRNA and protein (RT-PCR and immunofluorescence), where it mediates ATP release on cold stimulus. | PMID 26263540 | no |

PRAME's finding is scattered nuclear staining in a few melanocytes of normal skin,
mostly sun-damaged. It is recorded here for completeness and changes nothing.

## 4. A sequence gap: read-through proteins

Three candidates are read-through transcripts whose protein is mostly a ubiquitous
parent protein: ATF7-NPFF (identical to ATF7 over 411 of 463 residues),
ENSG00000273155 (the first 118 residues of MRPL30) and ENSG00000267618 (RAD51D and
RFFL). The catalog's containment rule did measure them: 88.6%, 89.4% and 90.8% of
their 9-mers lie in Unsafe proteins, just under the 95% cutoff, so all three are CTAs.

## What was done with it

- **A curation layer** now carries these findings into the filters. Each confirmed
  finding is assigned a tissue compartment and that compartment's tier, with its
  citation kept. `normal_function="no_vital"` (in `standard`) and `"no_protected"`
  (in `strict`) act on it. [Layer and tiers](index.md#8-the-curation-layer).
- **Family rules** extend three findings to whole families: olfactory receptors
  (olfactory epithelium, Protected), bitter taste receptors (taste buds, Protected) and
  hair keratins (hair follicle, Disposable). Trace-amine receptors, sweet and umami
  receptors and late cornified envelope proteins are recorded without acting, because
  the evidence is by similarity or not gene-specific.
- **Read-through loci** are removed by `readthrough="exclude"` in every mode but `all`.
- The catalog's categories are unchanged. Tiers for seven tissues the hazard registry
  had no call for are proposals in `oncoref/data/target-curation.yaml`, with
  references, and are not yet part of the registry.
