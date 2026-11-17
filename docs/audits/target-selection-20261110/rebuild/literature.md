# Literature review of normal expression

**912 proteins reviewed on 2026-11-15; 19 are Unsafe on published evidence.** Scope: Every protein that would be a Target candidate while its rare-cell signal is dismissed, or after release from the pooled 1% detection flag, and every Target candidate of the 2026-11-10 catalog.

## Rubric

| Rule | Requirement |
|---|---|
| R1 | A primary source reports the human gene in a panel organ by a gene-specific assay, and donor-level single-cell data detect it in that organ in at least 1% of one donor population. The rare-cell grade becomes established. |
| R2 | A primary source reports the human gene, RNA or protein, by a gene-specific assay in normal adult tissue that the tissue-tier registry or its updates tier Vital or Protected, and no gene-specific multi-tissue panel or reliable IHC contradicts it. Unsafe. |
| R3 | Human loss of function causes disease of a Vital tissue, and the ortholog is expressed in that adult tissue. Unsafe. |
| none | No publication, restriction to germline, placenta or embryo, conflicting reports, assays shared with paralogs, cultured or diseased samples, other species alone, or database annotation alone. Recorded; no category effect. |
| recorded | A finding that a second reader confirmed as human, adult, normal-tissue and gene-specific but that does not meet R1 to R3: cultured or induced cells, database annotation alone, a family-level source, no retrievable quote, or a Manageable or Disposable tissue. Recorded with its tissue and tier. No category effect; the candidate filters can act on it. |

## Limits

Web literature search by Claude subagents on 2026-11-10, and on 2026-11-11 for the 21 proteins whose signal the whole-body panel dismissed, and for the 3 proteins released when two lone donors in different cell types stopped counting as replication; each claim carries the URL opened and the quote returned. Candidate sweep: on 2026-11-12, 37 Claude research agents searched the literature for normal-tissue function or expression of every Target candidate of policy oncoref.target-selection.2026-11-10.1, 25 proteins each, and 36 second readers each tried to refute one batch's findings. On 2026-11-13 the abstract, open full text or supplementary table behind each proposed exclusion was retrieved from PubMed, PMC or the publisher and the quote recorded; an exclusion without a retrieved quote was not made.

- Many pages were read through a fetch tool that summarizes, so quotes are as that tool relayed them.
- Several papers were read at abstract level because full text was blocked; each evidence item says so.
- PubMed was reached through Europe PMC or E-utilities.
- A null result means nothing was found in this search, not that nothing exists.
- HPA and GTEx bulk values are numeric filters elsewhere and were not counted as literature.
- Candidate-sweep findings other than the exclusions carry the second reader's paraphrase, not a quote.
- A candidate-sweep null is one agent's search, not proof of absence; FSCB was not searched.
- Olfactory, taste and hair-keratin family rules rest on family-level sources and were not reviewed gene by gene.
- Tiers for olfactory epithelium, taste buds, nasal mucosa, lacrimal gland, sensory ganglia, bone and tooth come from the 2026-11-13 registry updates, written from abstracts.

## Decisions

| Protein | Search verdict | Effect | Reason | Sources |
|---|---|---|---|---|
| ADAM30 | vital expression supported | Unsafe (R2, brain) | RNA measured in human frontal cortex of non-demented controls by array and probe assay, with neuronal immunostaining. Absolute level is not given; controls are aged post-mortem brains. | [Letronne et al. 2016, EBioMedicine 9:278: ADAM30 Downregulat](https://pmc.ncbi.nlm.nih.gov/articles/PMC4972530/); [Same paper, immunohistochemistry](https://pmc.ncbi.nlm.nih.gov/articles/PMC4972530/); [https://pmc.ncbi.nlm.nih.gov/articles/PMC4972530/ (PMID 2733](https://pmc.ncbi.nlm.nih.gov/articles/PMC4972530/) |
| AMBN | vital expression supported | Unsafe (R2, small intestine) | Protein and RNA reported in normal adult duodenum and other gut mucosa by gene-specific assays. One pilot report with limited controls. | [Reseland JE, Heyward CA, Samara A. Revisiting ameloblastin; ](https://www.frontiersin.org/journals/cell-and-developmental-biology/articles/10.3389/fcell.2023.1251540/full); [Same paper, esophagus localization](https://www.frontiersin.org/journals/cell-and-developmental-biology/articles/10.3389/fcell.2023.1251540/full); [Tamburstuen MV et al. Ameloblastin expression and putative a](https://pubmed.ncbi.nlm.nih.gov/20854943/) |
| ASCL4 | protected expression supported | Unsafe (R2, skin) | Gene-specific RT-PCR restricted to skin, quantified in adult skin. Read at abstract level; the adult level is one seventh of fetal and its absolute value is not given. | [Jonsson et al. 2004, Genomics 84:859: Hash4, a novel human a](https://pubmed.ncbi.nlm.nih.gov/15475265/); [PMID 15475265 (Jonsson et al. 2004, Genomics 84:859-866); Un](https://pubmed.ncbi.nlm.nih.gov/15475265/) |
| BECN2 | vital expression supported | Unsafe (R2, brain) | Protein detected in adult human brain lysate by western blot with a peptide antibody. One laboratory; no cell-type resolution. | [He et al. 2013, Cell 154:1085: Beclin 2 functions in autopha](https://pmc.ncbi.nlm.nih.gov/articles/PMC4231430/); [PMID 23954414 (PMC4231430, Fig. 1C); UniProt A8MW95](https://pubmed.ncbi.nlm.nih.gov/23954414/); [PMID 23954414 (PMC4231430)](https://pubmed.ncbi.nlm.nih.gov/23954414/) |
| BTNL2 | vital expression supported | Unsafe (R2, small intestine) | Gene-specific real-time PCR in healthy adult intestine; a BTNL2-specific antibody study reports digestive tract. | [Lebrero-Fernandez et al. 2016, Immun Inflamm Dis 4:191 - Alt](https://pmc.ncbi.nlm.nih.gov/articles/PMC4879465/); [Panea et al. 2021, Commun Biol 4:913 - Butyrophilin-like 2 r](https://pmc.ncbi.nlm.nih.gov/articles/PMC8313535/); [UniProt Q9UIR0 (BTNL2_HUMAN), Tissue specificity comment (ci](https://rest.uniprot.org/uniprotkb/search?query=gene_exact:BTNL2+AND+organism_id:9606&fields=accession,cc_tissue_specificity&format=tsv) |
| CCN6 | protected expression supported | Unsafe (R2, cartilage) | Protein in chondrocytes of normal adult articular cartilage by immunohistochemistry. Read at abstract level; one 2004 study. | [Sen et al. 2004, Arthritis Rheum 50:488: WISP3-dependent reg](https://pubmed.ncbi.nlm.nih.gov/14872491/); [PMID 10471507](https://pubmed.ncbi.nlm.nih.gov/10471507/); [PMID 14872491](https://pubmed.ncbi.nlm.nih.gov/14872491/) |
| DEFB106A/B | vital expression supported | Unsafe (R2, lung) | PCR on human tissue cDNA panels reports expression in lung, with in situ hybridisation. Read at abstract level; the lung cell type is not named. | [Kao et al. 2003, Am J Respir Cell Mol Biol 29:71: ORFeome-ba](https://pubmed.ncbi.nlm.nih.gov/12600824/); [PMID 12600824 (Kao et al. 2003, Am J Respir Cell Mol Biol 29](https://pubmed.ncbi.nlm.nih.gov/12600824/); [PMID 24189797 (Xin et al. 2014, Mol Med Rep 9:149, DOI 10.38](https://pubmed.ncbi.nlm.nih.gov/24189797/) |
| GALR3 | vital expression supported | Unsafe (R2, adrenal gland) | The human cloning paper reports mRNA detectable in adrenal gland and pancreas. Read at abstract level; a 1998 assay whose method is not given there. | [Kolakowski et al. 1998, J Neurochem 71:2239: Molecular chara](https://pubmed.ncbi.nlm.nih.gov/9832121/); [PMID 9832121; rat distribution PMID 9722565](https://pubmed.ncbi.nlm.nih.gov/9832121/); [PMID 36693600 (PMC9918719)](https://pubmed.ncbi.nlm.nih.gov/36693600/) |
| GDF9 | vital expression supported | Unsafe (R2, bone marrow) | Northern blot and RT-PCR report RNA in bone marrow and pituitary. One 1998 report read at abstract level; RNA only. | [Fitzpatrick SL et al. Expression of growth differentiation f](https://pubmed.ncbi.nlm.nih.gov/9564873/); [Aaltonen J et al. Human GDF-9 and its novel homolog GDF-9B a](https://pubmed.ncbi.nlm.nih.gov/10443672/); [UniProt O60383 (GDF9_HUMAN), Tissue specificity comment](https://rest.uniprot.org/uniprotkb/search?query=gene_exact:GDF9+AND+organism_id:9606+AND+reviewed:true&fields=accession,gene_names,cc_tissue_specificity) |
| ITIH6 | vital expression supported | Unsafe (R2, kidney) | A gene-specific protein study reports abundant expression in kidney, colon and stomach. Read at abstract level; antibody validation unverified. | [Shi et al. 2025, Int J Biol Macromol 331:148360 - Characteri](https://pubmed.ncbi.nlm.nih.gov/41109368/); [Human Protein Atlas ITIH6 (summary and tissue pages)](https://www.proteinatlas.org/ENSG00000102313-ITIH6/tissue); [UniProt Q6UXX5 (ITIH6_HUMAN)](https://rest.uniprot.org/uniprotkb/search?query=gene_exact:ITIH6+AND+organism_id:9606&fields=accession,cc_tissue_specificity&format=tsv) |
| MAGEA12 | vital expression supported | Unsafe (R1, brain) | Sequenced qRT-PCR products in normal brain (Morgan 2013); our donor-level data detect it in brain. A later in situ study found no signal and proposed another peptide, and concedes MAGE-A12 cannot be excluded. | [Morgan RA et al. Cancer regression and neurological toxicity](https://pmc.ncbi.nlm.nih.gov/articles/PMC3581823); [Morgan RA et al. 2013 (PMC3581823) - NanoString](https://pmc.ncbi.nlm.nih.gov/articles/PMC3581823); [Morgan RA et al. 2013 (PMC3581823) - RNA-seq](https://pmc.ncbi.nlm.nih.gov/articles/PMC3581823) |
| MRGPRX1 | vital expression supported | Unsafe (R2, dorsal root ganglion) | qPCR of whole human dorsal root ganglia from four donors and single-cell RT-PCR of 85 human sensory neurons. Donor ages are not stated; the neurons picked were mostly small. | [Castro et al. 2019, JCI Insight 4:e131712: Activation of pru](https://pmc.ncbi.nlm.nih.gov/articles/PMC6824308/); [PMID 31536477 (PMC6824308), Results section 'Human DRG neuro](https://pubmed.ncbi.nlm.nih.gov/31536477/); [https://rest.uniprot.org/uniprotkb/Q96LB2 ; underlying prima](https://rest.uniprot.org/uniprotkb/Q96LB2) |
| MYH13 | specialized tissue expression | Unsafe (R2, skeletal muscle) | Isoform-specific qRT-PCR and protein gels place it in human extraocular and laryngeal muscle, which are skeletal muscles. | [Expression and identification of 10 sarcomeric MyHC isoforms](https://europepmc.org/article/MED/26970499); [Phylogenetic implications of the superfast myosin in extraoc](https://europepmc.org/article/MED/12110653); [MYH13, a superfast myosin expressed in extraocular, laryngea](https://europepmc.org/article/MED/38160435) |
| OBP2A | vital expression supported | Unsafe (R2, lung) | The cloning paper reports strong expression in lung, nasal structures, salivary and lachrymal glands, and separates it from the 95%-identical OBP2B. Read at abstract level. RNA-seq of four human olfactory samples and a lung sample shows it below 0.1 FPKM, which the rubric does not count against a gene-specific report. | [Lacazette et al. 2000, Hum Mol Genet 9:289: A novel human od](https://pubmed.ncbi.nlm.nih.gov/10607840/); [PMID 10607840 (DOI 10.1093/hmg/9.2.289) for nasal/salivary/l](https://pubmed.ncbi.nlm.nih.gov/10607840/) |
| OTOG | specialized tissue expression | Unsafe (R3, inner ear) | Human loss of function causes deafness (DFNB18B); the mouse ortholog is produced in the adult inner ear and nowhere else tested. | [Cohen-Salmon et al. 1997, PNAS: Otogelin: a glycoprotein spe](https://pmc.ncbi.nlm.nih.gov/articles/PMC25017/); [Cohen-Salmon et al. 1997 (PubMed abstract)](https://pubmed.ncbi.nlm.nih.gov/9405633/); [El-Amraoui et al. 2001, Hear Res: Spatiotemporal expression ](https://research.pasteur.fr/en/b/T4u) |
| SCGB1C1 | protected expression supported | Unsafe (R2, nasal mucosa) | Sequence-confirmed RT-PCR finds it in normal inferior turbinate mucosa of people without sinus disease. Gene-level RNA-seq also gives 322 to 3,126 FPKM in all four human olfactory epithelium samples. | [Olender et al. 2016, BMC Genomics 17:619: The human olfactor](https://static-content.springer.com/esm/art%3A10.1186%2Fs12864-016-2960-3/MediaObjects/12864_2016_2960_MOESM2_ESM.xlsx); [Lu et al. 2011, Respir Res 12:28: Expression and regulation ](https://pmc.ncbi.nlm.nih.gov/articles/PMC3063214/); [PMID 21385388 / PMC3063214 (DOI 10.1186/1465-9921-12-28)](https://pubmed.ncbi.nlm.nih.gov/21385388/) |
| SCGB1D1 | specialized tissue expression | Unsafe (R1, eye) | The protein was purified and sequenced from normal human tears (Lehrer 1998); our donor-level data detect it in acinar cells of the eye. Read at abstract level; the producing cells are inferred from the secretion. | [Lehrer RI et al. Lipophilin, a novel heterodimeric protein o](https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=EXT_ID:9720917%20AND%20SRC:MED&resultType=core&format=json); [Stoeckelhuber M et al. Immunohistochemical analysis of secre](https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=EXT_ID:16395610%20AND%20SRC:MED&resultType=core&format=json); [Zhao C et al. Lipophilins: human peptides homologous to rat ](https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=EXT_ID:10066439%20AND%20SRC:MED&resultType=core&format=json) |
| SLC25A2 | vital expression supported | Unsafe (R2, liver) | Two primary papers that characterize this isoform separately from SLC25A15 report it in liver, pancreas and kidney. Read at abstract level. | [Fiermonte G et al. The mitochondrial ornithine transporter. ](https://pubmed.ncbi.nlm.nih.gov/12807890/); [Camacho JA et al. Cloning and characterization of human ORNT](https://pubmed.ncbi.nlm.nih.gov/12948741/); [UniProt Q9BXI2 (ORNT2_HUMAN), Tissue specificity comment](https://rest.uniprot.org/uniprotkb/search?query=gene_exact:SLC25A2+AND+organism_id:9606+AND+reviewed:true&fields=accession,gene_names,cc_tissue_specificity) |
| SRY | vital expression supported | Unsafe (R2, brain) | Immunostaining in dopamine neurons of the adult male substantia nigra and absent in female. One laboratory and one antibody; no RNA. | [Czech et al. 2012, J Neurochem 122:260: The human testis-det](https://pubmed.ncbi.nlm.nih.gov/22568433/); [PMID 22568433 (PMC3529967), Czech et al., J Neurochem 2012](https://pubmed.ncbi.nlm.nih.gov/22568433/); [PMID 19914267 (PMC2891862), Ely et al., Steroids 2010 (rat r](https://pubmed.ncbi.nlm.nih.gov/19914267/) |
| SULT6B1 | vital expression supported | Unsafe (R2, kidney) | RT-PCR reports expression in human kidney and testis. Read at abstract level; one non-quantitative report, RNA only. | [Takahashi et al. 2009, Chem Biol Interact 181:253: Identific](https://pubmed.ncbi.nlm.nih.gov/19505954/); [PMID 19505954; UniProt Q6IMI4](https://pubmed.ncbi.nlm.nih.gov/19505954/) |
| TCHHL1 | protected expression supported | Unsafe (R2, skin) | Protein in the basal layer of normal human epidermis by immunohistochemistry, with mRNA in normal skin. Read at abstract level; one in-house antibody. | [Yamakoshi et al. 2013, Biochem Biophys Res Commun 432:66: Tr](https://pubmed.ncbi.nlm.nih.gov/23376073/); [PMID 23376073](https://pubmed.ncbi.nlm.nih.gov/23376073/); [PMID 21614011 (title and MeSH only)](https://pubmed.ncbi.nlm.nih.gov/21614011/) |
| AADACL3 | specialized tissue expression | No effect | Atlas annotation alone; no gene-specific assay. | [https://www.proteinatlas.org/ENSG00000188984-AADACL3 (summar](https://www.proteinatlas.org/ENSG00000188984-AADACL3) |
| ACCSL | no literature | No effect | No publication on normal adult expression was found. |  |
| ACTL7A | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| ACTL7B | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| ACTL8 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| ACTL9 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| ACTRT2 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| ADAD1 | no normal expression evidence | No effect | Sources describe germline, placental or embryonic restriction, or did not test adult somatic tissue. | [Schumacher et al. 1995, Biol Reprod: Distribution of Tenr, a](https://pubmed.ncbi.nlm.nih.gov/7543294/); [Snyder et al. 2020, Sci Rep: ADAD1 and ADAD2, testis-specifi](https://pmc.ncbi.nlm.nih.gov/articles/PMC7360552/); [Human Protein Atlas ADAD1](https://www.proteinatlas.org/ENSG00000164113-ADAD1) |
| ADAM18 | no normal expression evidence | No effect | Sources describe germline, placental or embryonic restriction, or did not test adult somatic tissue. | [UniProt Q9Y3Q7 (ADAM18_HUMAN), Tissue specificity comment](https://rest.uniprot.org/uniprotkb/search?query=gene_exact:ADAM18+AND+organism_id:9606+AND+reviewed:true&fields=accession,gene_names,cc_tissue_specificity); [Human Protein Atlas ENSG00000168619-ADAM18](https://www.proteinatlas.org/ENSG00000168619-ADAM18); [Frayne J, Hurd EA, Hall L. Human tMDC III: a sperm protein w](https://pubmed.ncbi.nlm.nih.gov/12200459/) |
| ADAM2 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| ADIG | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [https://pmc.ncbi.nlm.nih.gov/articles/PMC7966854/ (mouse); G](https://pmc.ncbi.nlm.nih.gov/articles/PMC7966854/) |
| AK4P3 | no literature | No effect | No publication on normal adult expression was found. |  |
| AKAP4 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| AMELX | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [https://pmc.ncbi.nlm.nih.gov/articles/PMC4059656 (PMID 24933](https://pmc.ncbi.nlm.nih.gov/articles/PMC4059656) |
| AMELY | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [https://pmc.ncbi.nlm.nih.gov/articles/PMC4059656 (PMID 24933](https://pmc.ncbi.nlm.nih.gov/articles/PMC4059656) |
| ANHX | no literature | No effect | No publication on normal adult expression was found. |  |
| ANKRD60 | no literature | No effect | No publication on normal adult expression was found. |  |
| ANP32D | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| ANTXRL | no literature | No effect | No publication on normal adult expression was found. |  |
| APOL5 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| ARL13A | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| ARMH2 | no normal expression evidence | No effect | Sources describe germline, placental or embryonic restriction, or did not test adult somatic tissue. | [ARMH2 is a cytosolic component of CatSper crucial for sperm ](https://pmc.ncbi.nlm.nih.gov/articles/PMC12638846/); [Human Protein Atlas gene page ARMH2](https://www.proteinatlas.org/ENSG00000260286-ARMH2) |
| ARSH | no literature | No effect | No publication on normal adult expression was found. |  |
| ATF7-NPFF | no literature | No effect | No gene-specific publication found. | [Human Protein Atlas gene record ATF7-NPFF (JSON for ENSG0000](https://www.proteinatlas.org/ENSG00000267281.json); [NCBI Gene 114108587 ATF7-NPFF readthrough (esummary)](https://www.ncbi.nlm.nih.gov/gene/114108587); [UniProt K7ELQ4 (ATF7-NPFF, unreviewed)](https://rest.uniprot.org/uniprotkb/search?query=gene_exact:ATF7-NPFF+AND+organism_id:9606) |
| ATXN3L | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| BHLHE23 | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [PMID 15363390; UniProt Q8NDY6 (by similarity only)](https://pubmed.ncbi.nlm.nih.gov/15363390/) |
| BMP15 | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [PMID 12446622](https://pubmed.ncbi.nlm.nih.gov/12446622/) |
| BOD1L2 | no literature | No effect | No publication on normal adult expression was found. |  |
| BPIFA3 | no literature | No effect | No publication on normal adult expression was found. |  |
| BPIFB3 | specialized tissue expression | No effect | Gene-level RNA-seq: 242 to 1,234 FPKM in all four human olfactory epithelium samples, 0.06 in respiratory epithelium. Donor ages are not stated. The registry tiers olfactory epithelium manageable; only Vital or Protected tissues exclude. | [Olender et al. 2016, BMC Genomics 17:619: The human olfactor](https://static-content.springer.com/esm/art%3A10.1186%2Fs12864-016-2960-3/MediaObjects/12864_2016_2960_MOESM2_ESM.xlsx); [PMID 27515280, Additional files 2 and 3 (12864_2016_2960_MOE](https://pubmed.ncbi.nlm.nih.gov/27515280/) |
| BPIFB6 | specialized tissue expression | No effect | Gene-level RNA-seq: 5 to 224 FPKM in all four human olfactory epithelium samples, 0.07 in respiratory epithelium. Three of the four samples are below 10 FPKM. The registry tiers olfactory epithelium manageable; only Vital or Protected tissues exclude. | [Olender et al. 2016, BMC Genomics 17:619: The human olfactor](https://static-content.springer.com/esm/art%3A10.1186%2Fs12864-016-2960-3/MediaObjects/12864_2016_2960_MOESM2_ESM.xlsx); [PMID 27515280, Additional files 2 and 3 (12864_2016_2960_MOE](https://pubmed.ncbi.nlm.nih.gov/27515280/); [PMID 26962226 (PMC4859712)](https://pubmed.ncbi.nlm.nih.gov/26962226/) |
| BPY2/B/C | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| BTN1A1 | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [PMID 15226505 (Ogg et al., PNAS 2004); PMID 20208008](https://pubmed.ncbi.nlm.nih.gov/15226505/); [PMID 20208008; UniProt Q13410 (by similarity)](https://pubmed.ncbi.nlm.nih.gov/20208008/) |
| C10orf120 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| C10orf62 | no literature | No effect | No publication on normal adult expression was found. |  |
| C12orf50 | no literature | No effect | No gene-specific publication found. | [Human Protein Atlas gene page C12orf50](https://www.proteinatlas.org/ENSG00000165805-C12orf50); [Wikipedia: C12orf50](https://en.wikipedia.org/wiki/C12orf50); [UniProt Q8NA57 (ZC3H11D / C12orf50), REST query](https://rest.uniprot.org/uniprotkb/search?query=(gene_exact:C12orf50%20OR%20gene_exact:DCAF8L2%20OR%20gene_exact:IZUMO2%20OR%20gene_exact:OR4N2%20OR%20gene_exact:PPP3R2%20OR%20gene_exact:RGSL1%20OR%20gene_exact:SPATA3%20OR%20gene_exact:TUBB8)%20AND%20organism_id:9606%20AND%20reviewed:true&fields=accession,gene_names,cc_tissue_specificity,cc_function,cc_developmental_stage&format=tsv) |
| C12orf71 | no literature | No effect | No publication on normal adult expression was found. |  |
| C16orf78 | no literature | No effect | No publication on normal adult expression was found. |  |
| C16orf82 | no literature | No effect | No publication on normal adult expression was found. |  |
| C16orf90 | no literature | No effect | No publication on normal adult expression was found. |  |
| C19orf67 | no literature | No effect | No publication on normal adult expression was found. |  |
| C19orf85 | no literature | No effect | No publication on normal adult expression was found. |  |
| C1orf202 | no literature | No effect | No publication on normal adult expression was found. |  |
| C1orf232 | no literature | No effect | No publication on normal adult expression was found. |  |
| C20orf141 | no literature | No effect | No publication on normal adult expression was found. |  |
| C20orf144 | no literature | No effect | No publication on normal adult expression was found. |  |
| C20orf173 | no literature | No effect | No publication on normal adult expression was found. |  |
| C2orf78 | no literature | No effect | No publication on normal adult expression was found. |  |
| C3orf22 | no literature | No effect | No publication on normal adult expression was found. |  |
| C4orf17 | no literature | No effect | No publication on normal adult expression was found. |  |
| C7orf33 | no literature | No effect | No gene-specific publication found. | [Human Protein Atlas: C7orf33](https://www.proteinatlas.org/ENSG00000170279-C7orf33); [UniProtKB Q8WU49 (C7orf33)](https://rest.uniprot.org/uniprotkb/search?query=gene_exact:C7orf33+AND+organism_id:9606&fields=accession,reviewed,gene_names,cc_tissue_specificity,cc_developmental_stage&format=tsv) |
| C8orf74 | no literature | No effect | No publication on normal adult expression was found. |  |
| CABS1 | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [PMID 38753592, PMC11098423, DOI 10.1371/journal.pone.0301855](https://pubmed.ncbi.nlm.nih.gov/38753592/); [PMID 25632019, PMC4386001, DOI 10.1152/ajpregu.00153.2014 (a](https://pubmed.ncbi.nlm.nih.gov/25632019/); [PMID 28381457, PMC5495925, DOI 10.1152/ajpregu.00393.2016 (a](https://pubmed.ncbi.nlm.nih.gov/28381457/) |
| CAGE1 | no normal expression evidence | No effect | Sources describe germline, placental or embryonic restriction, or did not test adult somatic tissue. | [Park S et al. 2003, Identification and characterization of a](https://pubmed.ncbi.nlm.nih.gov/12531476/); [Alsheimer M et al. 2005, The cancer/testis antigen CAGE-1 is](https://pubmed.ncbi.nlm.nih.gov/15819420/); [UniProt Q8TC20 (CAGE1_HUMAN), Tissue specificity comment](https://rest.uniprot.org/uniprotkb/Q8TC20) |
| CALR3 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| CAPZA3 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| CATSPER4 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| CBLL2 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| CBY2 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| CCDC116 | no literature | No effect | No publication on normal adult expression was found. |  |
| CCDC166 | no literature | No effect | No publication on normal adult expression was found. |  |
| CCDC179 | no literature | No effect | No publication on normal adult expression was found. |  |
| CCDC182 | no literature | No effect | No publication on normal adult expression was found. |  |
| CCDC185 | no literature | No effect | No publication on normal adult expression was found. |  |
| CCDC195 | no literature | No effect | No publication on normal adult expression was found. |  |
| CCDC197 | no literature | No effect | No publication on normal adult expression was found. |  |
| CCDC201 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| CCDC27 | no literature | No effect | No publication on normal adult expression was found. |  |
| CCDC42 | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. |  |
| CCDC63 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| CCDC70 | no literature | No effect | No publication on normal adult expression was found. |  |
| CCDC87 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| CCER1 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| CCT8L2 | no literature | No effect | No publication on normal adult expression was found. |  |
| CDCP2 | conflicting or nonspecific | No effect | Reports conflict, cannot separate the gene from paralogs, or come from cultured, diseased or fetal samples. | [Cortical gene expression: prognostic value for seizure outco](https://pmc.ncbi.nlm.nih.gov/articles/PMC5063688/); [Human Protein Atlas CDCP2 gene and tissue pages](https://www.proteinatlas.org/ENSG00000157211-CDCP2/tissue) |
| CDRT15L2 | no literature | No effect | No publication on normal adult expression was found. |  |
| CDX4 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| CDY1 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| CDY1B | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| CDY2A/B | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| CER1 | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [PMID 19113921](https://pubmed.ncbi.nlm.nih.gov/19113921/) |
| CETN1 | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [PMID 23641067 (PMC3711207)](https://pubmed.ncbi.nlm.nih.gov/23641067/) |
| CHCT1 | no normal expression evidence | No effect | Sources describe germline, placental or embryonic restriction, or did not test adult somatic tissue. | [Dai et al. 2021, Mol Biol Rep: Expression profile analysis o](https://pubmed.ncbi.nlm.nih.gov/33566224/); [Human Protein Atlas CHCT1](https://www.proteinatlas.org/ENSG00000141371-CHCT1) |
| CIB3 | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [PMID 37001993 (PMC10162464)](https://pubmed.ncbi.nlm.nih.gov/37001993/) |
| CIMAP1A | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| CIMIP2A | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| CIMIP3 | no literature | No effect | No gene-specific publication found. | [Human Protein Atlas gene page GUCA1ANB (ENSG00000287363)](https://www.proteinatlas.org/ENSG00000287363-CIMIP3); [UniProtKB X6R8D5 (CIMIP3 / GUCA1ANB)](https://rest.uniprot.org/uniprotkb/search?query=gene_exact:CIMIP3+AND+organism_id:9606) |
| CIMIP4 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| CIROP | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| CLDN34 | no literature | No effect | No publication on normal adult expression was found. |  |
| CPHXL | no normal expression evidence | No effect | Sources describe germline, placental or embryonic restriction, or did not test adult somatic tissue. | [Madissoon E et al. Characterization and target genes of nine](https://pmc.ncbi.nlm.nih.gov/articles/PMC4944136); [UniProt A0A1W2PPM1 (CPHXL), Developmental stage comment](https://rest.uniprot.org/uniprotkb/search?query=gene_exact:CPHXL+AND+organism_id:9606+AND+reviewed:true&fields=accession,gene_names,cc_tissue_specificity,cc_developmental_stage); [Human Protein Atlas ENSG00000283755-CPHXL](https://www.proteinatlas.org/ENSG00000283755-CPHXL) |
| CPSF4L | no literature | No effect | No publication on normal adult expression was found. |  |
| CPXCR1 | no literature | No effect | No publication on normal adult expression was found. |  |
| CRIPTO3 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| CRISP1 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| CSAG1 | specialized tissue expression | No effect | Evidence is from another species, a database annotation or an assay shared with a paralog, or the tissue is reproductive. | [Lin C et al. Cancer/testis antigen CSAGE is concurrently exp](https://pubmed.ncbi.nlm.nih.gov/12039054/); [UniProt Q6PB30 (CSAG1_HUMAN) tissue specificity](https://rest.uniprot.org/uniprotkb/search?query=gene_exact:CSAG1+AND+organism_id:9606+AND+reviewed:true&fields=accession,gene_names,cc_tissue_specificity,cc_developmental_stage&format=tsv); [Liggins AP et al. A panel of cancer-testis genes exhibiting ](https://pubmed.ncbi.nlm.nih.gov/20726502/) |
| CSAG2 | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [PMID 12039054; UniProt Q9Y5P2](https://pubmed.ncbi.nlm.nih.gov/12039054/) |
| CSAG3 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| CST11 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| CST8 | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [PMID 16837735 (J Androl 2006); PMID 15713831 (J Androl 2005)](https://pubmed.ncbi.nlm.nih.gov/16837735/) |
| CST9L | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| CSTL1 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| CT45A1 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| CT45A10 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| CT45A2/8/9 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| CT47A1/10/11/12/2/3/4/5/6/7/8/9 | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [PMID 16382448 (Chen et al., Genes Chromosomes Cancer 2006, d](https://pubmed.ncbi.nlm.nih.gov/16382448/) |
| CT47B1 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| CT47C1 | no literature | No effect | No publication on normal adult expression was found. |  |
| CT55 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| CT83 | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [PMID 36484013 (J Cancer 2022, PMC9723985); PMID 29666402 (Sc](https://pubmed.ncbi.nlm.nih.gov/36484013/) |
| CTAG2 | conflicting or nonspecific | No effect | Reports conflict, cannot separate the gene from paralogs, or come from cultured, diseased or fetal samples. | [Lethe B et al. LAGE-1, a new gene with tumor specificity. In](https://pubmed.ncbi.nlm.nih.gov/9626360/); [de Carvalho F et al. Evaluation of LAGE-1 and NY-ESO-1 expre](https://pmc.ncbi.nlm.nih.gov/articles/PMC3077292/); [UniProt O75638 (CTAG2_HUMAN) tissue specificity](https://www.uniprot.org/uniprotkb/O75638/entry) |
| CTAGE1 | no literature | No effect | No publication on normal adult expression was found. |  |
| CTAGE15 | no literature | No effect | No publication on normal adult expression was found. |  |
| CTAGE6 | no literature | No effect | No publication on normal adult expression was found. |  |
| CTXND2 | no literature | No effect | No publication on normal adult expression was found. |  |
| CXorf49/B | no literature | No effect | No publication on normal adult expression was found. |  |
| CXorf51A/B | no literature | No effect | No publication on normal adult expression was found. |  |
| CXorf58 | no normal expression evidence | No effect | Sources describe germline, placental or embryonic restriction, or did not test adult somatic tissue. | [Human Protein Atlas gene page CXorf58](https://www.proteinatlas.org/ENSG00000165182-CXorf58); [Wikipedia: CXorf58 (secondary source; cites Vandenbrouck et ](https://en.wikipedia.org/wiki/CXorf58); [Functionathon: a manual data mining workflow to generate fun](https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=TITLE_ABS:(OR2F1%20OR%20CXorf58%20OR%20ARMH2%20OR%20C6orf229)&resultType=core&format=json&pageSize=25) |
| CXorf66 | no literature | No effect | No publication on normal adult expression was found. |  |
| CYLC1 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| CYLC2 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| CYP26C1 | conflicting or nonspecific | No effect | Reports conflict, cannot separate the gene from paralogs, or come from cultured, diseased or fetal samples. | [Taimi et al. 2004, J Biol Chem 279:77 - A novel human cytoch](https://pubmed.ncbi.nlm.nih.gov/14532297/); [UniProt Q6V0L0 (CP26C_HUMAN), Tissue specificity (cites PubM](https://rest.uniprot.org/uniprotkb/search?query=gene_exact:CYP26C1+AND+organism_id:9606&fields=accession,cc_tissue_specificity&format=tsv); [Xi & Yang 2008, Gene Expr Patterns 8:438 - Expression of RAL](https://pubmed.ncbi.nlm.nih.gov/18502188/) |
| DAZL | no normal expression evidence | No effect | Sources describe germline, placental or embryonic restriction, or did not test adult somatic tissue. | [Shan Z et al. 1996, A SPGY copy homologous to the mouse gene](https://pubmed.ncbi.nlm.nih.gov/8968755/); [Yen PH et al. 1996, The human autosomal gene DAZLA: testis s](https://pubmed.ncbi.nlm.nih.gov/8968756/); [Saxena R et al. 1996, The DAZ gene cluster on the human Y ch](https://pubmed.ncbi.nlm.nih.gov/8896558/) |
| DBX1 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| DCAF4L2 | no literature | No effect | No publication on normal adult expression was found. |  |
| DCAF8L1 | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [https://www.proteinatlas.org/ENSG00000226372-DCAF8L1/tissue ](https://www.proteinatlas.org/ENSG00000226372-DCAF8L1/tissue) |
| DCAF8L2 | no literature | No effect | No gene-specific publication found. | [Human Protein Atlas gene page DCAF8L2](https://www.proteinatlas.org/ENSG00000189186-DCAF8L2); [CRL4-DCAF8L2 E3 ligase promotes ubiquitination and degradati](https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=EXT_ID:35487060%20OR%20EXT_ID:38484656&resultType=core&format=json) |
| DCST1 | conflicting or nonspecific | No effect | Reports conflict, cannot separate the gene from paralogs, or come from cultured, diseased or fetal samples. | [Human Protein Atlas gene page DCST1 (ENSG00000163357)](https://www.proteinatlas.org/ENSG00000163357-DCST1); [Nair S et al. Global functional profiling of human ubiquitom](https://pmc.ncbi.nlm.nih.gov/articles/PMC5080589); [Noda T et al. Sperm membrane proteins DCST1 and DCST2 are re](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC8989947/fullTextXML) |
| DDI1 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| DEFB104A/B | specialized tissue expression | No effect | Database annotation alone; the primary paper could not be opened. | [UniProt Q8WTQ1 (https://rest.uniprot.org/uniprotkb/Q8WTQ1.tx](https://rest.uniprot.org/uniprotkb/Q8WTQ1.txt); [PMID 14871351 (Premratanachai et al. 2004, Oral Microbiol Im](https://pubmed.ncbi.nlm.nih.gov/14871351/) |
| DEFB107A/B | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [PMID 14871351](https://pubmed.ncbi.nlm.nih.gov/14871351/) |
| DEFB108B | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [UniProt Q8NET1 (https://rest.uniprot.org/uniprotkb/Q8NET1.tx](https://rest.uniprot.org/uniprotkb/Q8NET1.txt); [PMID 14871351](https://pubmed.ncbi.nlm.nih.gov/14871351/) |
| DEFB110 | no literature | No effect | No publication on normal adult expression was found. |  |
| DEFB112 | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [PMID 14871351](https://pubmed.ncbi.nlm.nih.gov/14871351/) |
| DEFB113 | no literature | No effect | No publication on normal adult expression was found. |  |
| DEFB114 | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [PMID 14871351](https://pubmed.ncbi.nlm.nih.gov/14871351/) |
| DEFB115 | no literature | No effect | No publication on normal adult expression was found. |  |
| DEFB116 | no literature | No effect | No publication on normal adult expression was found. |  |
| DEFB118 | specialized tissue expression | No effect | Database annotation alone; the cited abstract does not mention pancreas. | [UniProt Q96PH6 (https://rest.uniprot.org/uniprotkb/Q96PH6.tx](https://rest.uniprot.org/uniprotkb/Q96PH6.txt) |
| DEFB121 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| DEFB123 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| DEFB125 | no literature | No effect | No publication on normal adult expression was found. |  |
| DEFB126 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| DEFB127 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| DEFB128 | no literature | No effect | No publication on normal adult expression was found. |  |
| DEFB129 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| DEFB130A/B | specialized tissue expression | No effect | Macrophages differentiated and stimulated in culture, not tissue. | [Terkawi et al. 2017, Sci Rep 7:41772, PMID 28181499, https:/](https://pmc.ncbi.nlm.nih.gov/articles/PMC5299406) |
| DEFB132 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| DEFB135 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| DEFB136 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| DMRT1 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| DMRTB1 | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [Ottolenghi et al. 2002, Genomics 79:333-43, PMID 11863363 (p](https://pmc.ncbi.nlm.nih.gov/articles/PMC4197572) |
| DMRTC2 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| DNAJB8 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| DPEP2NB | no literature | No effect | No publication on normal adult expression was found. |  |
| DPPA2 | conflicting or nonspecific | No effect | Reports conflict, cannot separate the gene from paralogs, or come from cultured, diseased or fetal samples. | [John T et al. ECSA/DPPA2 is an embryo-cancer antigen that is](https://europepmc.org/article/MED/18519755); [Du J et al. [Molecular cloning and characterization analysis](https://europepmc.org/article/MED/15583978); [UniProtKB Q7Z7J5 (DPPA2)](https://rest.uniprot.org/uniprotkb/search?query=gene_exact:DPPA2+AND+organism_id:9606&fields=accession,reviewed,gene_names,cc_tissue_specificity,cc_developmental_stage&format=tsv) |
| DPPA3 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| DPPA5 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| DPRX | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| DSPP | specialized tissue expression | No effect | Reported in dental pulp, which the registry tiers manageable; only Vital or Protected tissues exclude. | [Xiao et al. 2001, Nat Genet 27:201-4, PMID 11175790, https:/](https://www.nature.com/articles/ng0201_201); [Xiao et al. 2001, Nat Genet 27:201-4, PMID 11175790](https://pubmed.ncbi.nlm.nih.gov/11175790/); [Ogbureke et al. 2016, J Histochem Cytochem 64:623-36, PMID 2](https://pubmed.ncbi.nlm.nih.gov/27666430/) |
| DUSP21 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| DUXA | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| EBLN1 | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [He et al. 2016, Int J Mol Sci, https://pmc.ncbi.nlm.nih.gov/](https://pmc.ncbi.nlm.nih.gov/articles/PMC4848891) |
| EDDM3A | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| EDDM3B | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| EFCAB3 | no literature | No effect | No gene-specific publication found. | [Human Protein Atlas EFCAB3 gene page](https://www.proteinatlas.org/ENSG00000172421-EFCAB3) |
| EFCAB9 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| ENSG00000237378 | no literature | No effect | No publication on normal adult expression was found. |  |
| ENSG00000241690 | no literature | No effect | No publication on normal adult expression was found. |  |
| ENSG00000248235 | no literature | No effect | No publication on normal adult expression was found. |  |
| ENSG00000249773 | no literature | No effect | No gene-specific publication found. | [Human Protein Atlas ENSG00000249773](https://www.proteinatlas.org/ENSG00000249773) |
| ENSG00000255432 | no literature | No effect | No publication on normal adult expression was found. |  |
| ENSG00000255835 | no literature | No effect | No gene-specific publication found. | [Ensembl REST lookup ENSG00000255835](https://rest.ensembl.org/lookup/id/ENSG00000255835?content-type=application/json); [Human Protein Atlas ENSG00000255835](https://www.proteinatlas.org/ENSG00000255835); [UniProt search xref:ensembl-ENSG00000255835](https://rest.uniprot.org/uniprotkb/search?query=xref:ensembl-ENSG00000255835&fields=accession,gene_names,protein_name,reviewed,cc_tissue_specificity&format=tsv) |
| ENSG00000256349 | no literature | No effect | No gene-specific publication found. | [Human Protein Atlas gene record ENSG00000256349 (JSON)](https://www.proteinatlas.org/ENSG00000256349.json) |
| ENSG00000256825 | no literature | No effect | No gene-specific publication found. | [Human Protein Atlas gene page ENSG00000256825](https://www.proteinatlas.org/ENSG00000256825) |
| ENSG00000258417 | no literature | No effect | No publication on normal adult expression was found. |  |
| ENSG00000258691 | no literature | No effect | No publication on normal adult expression was found. |  |
| ENSG00000266202 | no literature | No effect | No gene-specific publication found. | [Human Protein Atlas gene page ENSG00000266202](https://www.proteinatlas.org/ENSG00000266202) |
| ENSG00000267140 | no literature | No effect | No publication on normal adult expression was found. |  |
| ENSG00000267618 | no literature | No effect | No gene-specific publication found. | [Human Protein Atlas ENSG00000267618](https://www.proteinatlas.org/ENSG00000267618); [Liu, Huse & Kannan 2025, Front Genet 16:1672015 (abstract on](https://pubmed.ncbi.nlm.nih.gov/41356219/); [https://rest.ensembl.org/sequence/id/ENSG00000267618?type=pr](https://rest.ensembl.org/sequence/id/ENSG00000267618?type=protein) |
| ENSG00000270299 | no literature | No effect | No publication on normal adult expression was found. |  |
| ENSG00000273155 | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [https://rest.ensembl.org/sequence/id/ENSG00000273155?type=pr](https://rest.ensembl.org/sequence/id/ENSG00000273155?type=protein) |
| ENSG00000274944 | no literature | No effect | No publication on normal adult expression was found. |  |
| ENSG00000279983 | no literature | No effect | No gene-specific publication found. | [Human Protein Atlas gene page ENSG00000279983](https://www.proteinatlas.org/ENSG00000279983) |
| ENSG00000283205 | no literature | No effect | No publication on normal adult expression was found. |  |
| ENSG00000284779 | no literature | No effect | No gene-specific publication found. | [Ensembl REST lookup ENSG00000284779](https://rest.ensembl.org/lookup/id/ENSG00000284779?content-type=application/json;expand=1); [Human Protein Atlas, ENSG00000284779](https://www.proteinatlas.org/ENSG00000284779) |
| ENSG00000284934 | no literature | No effect | No publication on normal adult expression was found. |  |
| ENSG00000285566 | no literature | No effect | No publication on normal adult expression was found. |  |
| ENSG00000286224 | no literature | No effect | No publication on normal adult expression was found. |  |
| ENSG00000287585 | no literature | No effect | No publication on normal adult expression was found. |  |
| ENSG00000289809 | no literature | No effect | No publication on normal adult expression was found. |  |
| EPPIN-WFDC6 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| ERVV-1 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| ERVV-2 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| ESX1 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| ETDA/B | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| EVX1 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| EVX2 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| EXD1 | no normal expression evidence | No effect | Sources describe germline restriction or did not test adult somatic tissue; the only brain statement is HPA database RNA annotation. Human loss of function causes male infertility (testis, not a vital tissue), so R3 does not fire. | [Yang Z et al. PIWI Slicing and EXD1 Drive Biogenesis of Nucl](https://pubmed.ncbi.nlm.nih.gov/26669262/); [Olotu O et al. Intermitochondrial cement (IMC) harbors piRNA](https://pubmed.ncbi.nlm.nih.gov/36624638/); [Hu K et al. Novel homozygous variants in piRNA pathway facto](https://pubmed.ncbi.nlm.nih.gov/40645105/) |
| EXOC1L | no literature | No effect | No publication on normal adult expression was found. |  |
| FAM170A | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| FAM170B | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| FAM187B | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| FAM237A | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [PMID 36853120 (Li HZ et al., FEBS J 2023;290:3461-3479, doi:](https://pubmed.ncbi.nlm.nih.gov/36853120/); [PMID 35742932 (Narimatsu Y et al., Int J Mol Sci 2022;23:648](https://pubmed.ncbi.nlm.nih.gov/35742932/) |
| FAM240B | no literature | No effect | No publication on normal adult expression was found. |  |
| FAM246B | no literature | No effect | No publication on normal adult expression was found. |  |
| FAM24A | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| FAM47A | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| FAM47B | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| FAM47C | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| FAM90A10 | no literature | No effect | No publication on normal adult expression was found. |  |
| FAM90A26 | conflicting or nonspecific | No effect | Reports conflict, cannot separate the gene from paralogs, or come from cultured, diseased or fetal samples. | [Bosch N et al. Characterization and evolution of the novel g](https://europepmc.org/article/MED/17684299); [Zheng C et al. Evaluation of PE5 Entries on the Human Proteo](https://europepmc.org/article/MED/40503744); [Human Protein Atlas: FAM90A26](https://www.proteinatlas.org/ENSG00000229924-FAM90A26) |
| FAM90A27P | no literature | No effect | No publication on normal adult expression was found. |  |
| FAM90A9 | no literature | No effect | No publication on normal adult expression was found. |  |
| FAM9A | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| FBXO47 | conflicting or nonspecific | No effect | Reports conflict, cannot separate the gene from paralogs, or come from cultured, diseased or fetal samples. | [Simon-Kayser B et al. Molecular cloning and characterization](https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=EXT_ID:15723337%20AND%20SRC:MED&resultType=core&format=json); [UniProtKB Q5MNV8 (FBXO47)](https://rest.uniprot.org/uniprotkb/search?query=gene_exact:FBXO47+AND+organism_id:9606); [Human Protein Atlas gene page FBXO47 (ENSG00000204952)](https://www.proteinatlas.org/ENSG00000204952-FBXO47) |
| FKBP6 | conflicting or nonspecific | No effect | Reports conflict, cannot separate the gene from paralogs, or come from cultured, diseased or fetal samples. | [A novel human gene FKBP6 is deleted in Williams syndrome (Me](https://europepmc.org/article/MED/9782077); [Is a genetic defect in Fkbp6 a common cause of azoospermia i](https://pmc.ncbi.nlm.nih.gov/articles/PMC6275806/); [Human Protein Atlas FKBP6 tissue page](https://www.proteinatlas.org/ENSG00000077800-FKBP6/tissue) |
| FMR1NB | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| FNDC8 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| FOXL3 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| FOXR1 | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [https://pmc.ncbi.nlm.nih.gov/articles/PMC8559929/ (Mota et a](https://pmc.ncbi.nlm.nih.gov/articles/PMC8559929/) |
| FOXR2 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| FREY1 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| FSCN3 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| FTHL17 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| FUT5 | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [PMID 7650030; UniProt Q11128 tissue-specificity note](https://pubmed.ncbi.nlm.nih.gov/7650030/) |
| GAB4 | no literature | No effect | No publication on normal adult expression was found. |  |
| GAGE1 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| GAGE10 | no normal expression evidence | No effect | Sources describe germline, placental or embryonic restriction, or did not test adult somatic tissue. | [Gjerstorff et al. 2006, Br J Cancer: Restriction of GAGE pro](https://pmc.ncbi.nlm.nih.gov/articles/PMC2361341/); [Gjerstorff & Ditzel 2008, Tissue Antigens: An overview of th](https://pubmed.ncbi.nlm.nih.gov/18179644/); [Gjerstorff et al. 2013, BMC Cancer: Analysis of GAGE, NY-ESO](https://pubmed.ncbi.nlm.nih.gov/24103781/) |
| GAGE12C/D/E | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| GAGE12H | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| GAGE12J | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| GAGE13 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| GAGE2A | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| GALP | specialized tissue expression | No effect | The assay detects alarin, a splice product sharing only the peptide's first residues with GALP; antibody specificity is not established. | [PMID 23178871 (doi 10.1016/j.exer.2012.11.009)](https://pubmed.ncbi.nlm.nih.gov/23178871/); [PMID 11861493 (macaque ISH); PMID 14962070 (rat ISH)](https://pubmed.ncbi.nlm.nih.gov/11861493/) |
| GARIN3 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| GARIN5B | no literature | No effect | No publication on normal adult expression was found. |  |
| GFRA4 | specialized tissue expression | No effect | Reported in thyroid gland: parafollicular cells, which the registry tiers manageable; only Vital or Protected tissues exclude. | [PMID 11116144 (doi 10.1074/jbc.M008279200)](https://pubmed.ncbi.nlm.nih.gov/11116144/); [PMID 16497798](https://pubmed.ncbi.nlm.nih.gov/16497798/) |
| GFY | specialized tissue expression | No effect | Evidence is from another species, a database annotation or an assay shared with a paralog, or the tissue is reproductive. | [Kaneko-Goto T et al. Goofy coordinates the acuity of olfacto](https://pmc.ncbi.nlm.nih.gov/articles/PMC6619734); [Same paper, abstract](https://pubmed.ncbi.nlm.nih.gov/23926254/); [Human Protein Atlas gene page GFY](https://www.proteinatlas.org/ENSG00000261949-GFY) |
| GJE1 | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [https://pmc.ncbi.nlm.nih.gov/articles/PMC2719720](https://pmc.ncbi.nlm.nih.gov/articles/PMC2719720) |
| GK2 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| GLT6D1 | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [PMID 19897590 (doi 10.1093/hmg/ddp508)](https://pubmed.ncbi.nlm.nih.gov/19897590/) |
| GLYATL3 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| GOLGA6A | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| GOLGA6B | no literature | No effect | No publication on normal adult expression was found. |  |
| GOLGA6C | no literature | No effect | No publication on normal adult expression was found. |  |
| GOLGA6D | no literature | No effect | No publication on normal adult expression was found. |  |
| GOLGA6L1 | no literature | No effect | No publication on normal adult expression was found. |  |
| GOLGA6L2 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| GOLGA6L22 | no literature | No effect | No publication on normal adult expression was found. |  |
| GOLGA6L24 | no literature | No effect | No publication on normal adult expression was found. |  |
| GOLGA6L25 | no literature | No effect | No publication on normal adult expression was found. |  |
| GOLGA6L26 | no literature | No effect | No publication on normal adult expression was found. |  |
| GOLGA6L6 | no literature | No effect | No publication on normal adult expression was found. |  |
| GOLGA8F | no literature | No effect | No publication on normal adult expression was found. |  |
| GOLGA8G | no literature | No effect | No publication on normal adult expression was found. |  |
| GOLGA8J | no literature | No effect | No publication on normal adult expression was found. |  |
| GOLGA8K | no literature | No effect | No publication on normal adult expression was found. |  |
| GOLGA8T | no literature | No effect | No publication on normal adult expression was found. |  |
| GOT1L1 | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [PMID 25287256; doi:10.1007/s00726-014-1847-3 (Got1l1 knockou](https://pubmed.ncbi.nlm.nih.gov/25287256/) |
| GPR152 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| GPX5 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| GSC2 | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [PMID 20921407; doi:10.1073/pnas.1012764107 (abstract read vi](https://pubmed.ncbi.nlm.nih.gov/20921407/) |
| GSTT4 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| GSX2 | no normal expression evidence | No effect | Sources describe germline, placental or embryonic restriction, or did not test adult somatic tissue. | [De Mori R et al. Agenesis of the putamen and globus pallidus](https://pubmed.ncbi.nlm.nih.gov/31412107/); [Lopez-Juarez A et al. Gsx2 controls region-specific activati](https://pubmed.ncbi.nlm.nih.gov/23723414/); [Dunwell TL, Holland PWH. Diversity of human and mouse homeob](https://pmc.ncbi.nlm.nih.gov/articles/PMC5094009/) |
| H1-8 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| H2AB1 | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [PMID 28234895; doi:10.1371/journal.pgen.1006633 (abstract vi](https://pubmed.ncbi.nlm.nih.gov/28234895/) |
| H2AL1Q | no literature | No effect | No publication on normal adult expression was found. |  |
| H2AP | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| H2BC1 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| H2BW1 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| H2BW2 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| H3Y1 | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [PMID 20819935 (PMC2935562)](https://pubmed.ncbi.nlm.nih.gov/20819935/) |
| HDGFL1 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| HES3 | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [PMID 25371201 (PMC4271235), abstract only](https://pubmed.ncbi.nlm.nih.gov/25371201/); [PMID 30054579 (doi:10.1038/s41598-018-29434-2)](https://pubmed.ncbi.nlm.nih.gov/30054579/) |
| HHLA1 | conflicting or nonspecific | No effect | Reports conflict, cannot separate the gene from paralogs, or come from cultured, diseased or fetal samples. | [Intergenic splicing between a HERV-H endogenous retrovirus a](https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=EXT_ID:10329003%20AND%20SRC:MED&resultType=core&format=json); [Human Protein Atlas gene page HHLA1](https://www.proteinatlas.org/ENSG00000132297-HHLA1); [Differential Expression and Co-Variation of HERV-Associated ](https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=(TITLE:HERV-H%20AND%20TITLE:intergenic)%20OR%20(ABSTRACT:HHLA1)&resultType=core&format=json&pageSize=15) |
| HIGD2B | no literature | No effect | No publication on normal adult expression was found. |  |
| HMGB1P1 | no literature | No effect | No publication on normal adult expression was found. |  |
| HMGB4 | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [PMID 27608812 (PMC5036535)](https://pubmed.ncbi.nlm.nih.gov/27608812/) |
| HNRNPCL1 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| HNRNPCL2 | no literature | No effect | No publication on normal adult expression was found. |  |
| HNRNPCL3/4 | no literature | No effect | No publication on normal adult expression was found. |  |
| HOXB1 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| HOXC12 | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [PMID 11960714 (Shang, Pruett, Awgulewitsch, Mech Dev 2002)](https://pubmed.ncbi.nlm.nih.gov/11960714/) |
| HOXD12 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| HSFX3 | no literature | No effect | No publication on normal adult expression was found. |  |
| HSFX4 | no literature | No effect | No publication on normal adult expression was found. |  |
| HSFY1/2 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| HUS1B | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [PMID 11944979 (Hang et al., Genomics 2002), abstract plus Un](https://pubmed.ncbi.nlm.nih.gov/11944979/) |
| IGBP1C | no literature | No effect | No publication on normal adult expression was found. |  |
| IL36B | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [PMID 16646978 (PMC1526623), Fig 4d](https://pubmed.ncbi.nlm.nih.gov/16646978/); [PMID 21881584 (Carrier et al., J Invest Dermatol 2011) for k](https://pubmed.ncbi.nlm.nih.gov/21881584/) |
| IQCA1L | no literature | No effect | No publication on normal adult expression was found. |  |
| IQCF1 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| IQCF2 | no literature | No effect | No publication on normal adult expression was found. |  |
| IQCF5 | no literature | No effect | No publication on normal adult expression was found. |  |
| IQCF6 | no literature | No effect | No publication on normal adult expression was found. |  |
| IRGC | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| IZUMO1R | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [PMID 36972453 (mouse)](https://pubmed.ncbi.nlm.nih.gov/36972453/) |
| IZUMO2 | no normal expression evidence | No effect | Sources describe germline, placental or embryonic restriction, or did not test adult somatic tissue. | [Izumo is part of a multiprotein family whose members form la](https://pmc.ncbi.nlm.nih.gov/articles/PMC3779078); [Human Protein Atlas gene page IZUMO2](https://www.proteinatlas.org/ENSG00000161652-IZUMO2); [Implication of IZUMO2 in the cell-in-cell phenomenon: a pote](https://pmc.ncbi.nlm.nih.gov/articles/PMC10912533/) |
| KCNU1 | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [PMID 31270758 (mouse)](https://pubmed.ncbi.nlm.nih.gov/31270758/) |
| KDM4F | no literature | No effect | No publication on normal adult expression was found. |  |
| KHDC3L | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| KIAA1210 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| KIF4B | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| KLF14 | specialized tissue expression | No effect | Reported in adipose tissue, which the registry tiers manageable; only Vital or Protected tissues exclude. | [PMID 29632379 (functional follow-up); PMID 21572415](https://pubmed.ncbi.nlm.nih.gov/29632379/) |
| KLF18 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| KRT37 | specialized tissue expression | No effect | Reported in hair follicle, which the registry tiers disposable; only Vital or Protected tissues exclude. | [PMID 15086535](https://pubmed.ncbi.nlm.nih.gov/15086535/) |
| KRT39 | specialized tissue expression | No effect | Reported in hair follicle, which the registry tiers disposable; only Vital or Protected tissues exclude. | [PMID 17301834 (via UniProt Q6A163 tissue-specificity annotat](https://pubmed.ncbi.nlm.nih.gov/17301834/) |
| KRTAP1-4 | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [PMID 11279113 (UniProt P0C5Y4)](https://pubmed.ncbi.nlm.nih.gov/11279113/) |
| KRTAP10-6 | specialized tissue expression | No effect | Reported in hair follicle, which the registry tiers disposable; only Vital or Protected tissues exclude. | [PMID 15028290 with PMID 14962103 (UniProt P60371)](https://pubmed.ncbi.nlm.nih.gov/15028290/) |
| KRTAP10-8 | specialized tissue expression | No effect | Reported in hair follicle, which the registry tiers disposable; only Vital or Protected tissues exclude. | [PMID 15028290 with PMID 14962103 (UniProt P60410)](https://pubmed.ncbi.nlm.nih.gov/15028290/) |
| KRTAP12-4 | specialized tissue expression | No effect | Reported in hair follicle, which the registry tiers disposable; only Vital or Protected tissues exclude. | [PMID 15028290 (names KRTAP12.4 explicitly) with PMID 1496210](https://pubmed.ncbi.nlm.nih.gov/15028290/) |
| KRTAP13-3 | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [PMID 12359730 (abstract only; member-level table in full tex](https://pubmed.ncbi.nlm.nih.gov/12359730/) |
| KRTAP13-4 | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [PMID 12359730 (abstract only)](https://pubmed.ncbi.nlm.nih.gov/12359730/) |
| KRTAP19-2 | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [PMID 12359730 (abstract only)](https://pubmed.ncbi.nlm.nih.gov/12359730/) |
| KRTAP19-4 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| KRTAP19-6 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| KRTAP19-7 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| KRTAP19-8 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| KRTAP20-1 | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [PMID 12359730 (full-text table needed)](https://pubmed.ncbi.nlm.nih.gov/12359730/) |
| KRTAP20-2 | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [PMID 12359730 (full-text table needed)](https://pubmed.ncbi.nlm.nih.gov/12359730/) |
| KRTAP20-3 | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [none found; UniProt Q3LI60 has no expression data](https://pubmed.ncbi.nlm.nih.gov/12359730/) |
| KRTAP20-4 | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [none found; UniProt Q3LI62 has no expression data](https://pubmed.ncbi.nlm.nih.gov/12359730/) |
| KRTAP21-1 | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [PMID 12359730 (full-text table needed)](https://pubmed.ncbi.nlm.nih.gov/12359730/) |
| KRTAP21-2 | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [PMID 12359730 (full-text table needed)](https://pubmed.ncbi.nlm.nih.gov/12359730/) |
| KRTAP21-3 | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [none found; UniProt Q3LHN1 has no expression data](https://pubmed.ncbi.nlm.nih.gov/12359730/) |
| KRTAP22-1 | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [PMID 12359730 (full-text table needed)](https://pubmed.ncbi.nlm.nih.gov/12359730/) |
| KRTAP22-2 | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [none found; UniProt Q3LI68 has no expression data](https://pubmed.ncbi.nlm.nih.gov/12359730/) |
| KRTAP23-1 | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [PMID 12359730 (full-text table needed)](https://pubmed.ncbi.nlm.nih.gov/12359730/) |
| KRTAP24-1 | specialized tissue expression | No effect | Reported in hair follicle, which the registry tiers disposable; only Vital or Protected tissues exclude. | [PMID 17235325](https://pubmed.ncbi.nlm.nih.gov/17235325/) |
| KRTAP25-1 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| KRTAP29-1 | specialized tissue expression | No effect | Reported in hair follicle, which the registry tiers disposable; only Vital or Protected tissues exclude. | [https://www.proteinatlas.org/ENSG00000212658](https://www.proteinatlas.org/ENSG00000212658) |
| KRTAP4-1 | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [UniProt Q9BYQ7 (PMID 15955084); PMID 11279113](https://pubmed.ncbi.nlm.nih.gov/15955084/) |
| KRTAP4-16 | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [PMID 11279113 (family-level only)](https://pubmed.ncbi.nlm.nih.gov/11279113/) |
| KRTAP5-10 | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [PMID 15144888; UniProt Q6L8G5](https://pubmed.ncbi.nlm.nih.gov/15144888/) |
| KRTAP5-3 | specialized tissue expression | No effect | Reported in hair follicle, which the registry tiers disposable; only Vital or Protected tissues exclude. | [PMID 15144888; UniProt Q6L8H2](https://pubmed.ncbi.nlm.nih.gov/15144888/) |
| KRTAP5-4 | specialized tissue expression | No effect | Reported in hair follicle, which the registry tiers disposable; only Vital or Protected tissues exclude. | [PMID 15144888; UniProt Q6L8H1](https://pubmed.ncbi.nlm.nih.gov/15144888/) |
| KRTAP5-5 | specialized tissue expression | No effect | Reported in hair follicle, which the registry tiers disposable; only Vital or Protected tissues exclude. | [PMID 15144888; UniProt Q701N2](https://pubmed.ncbi.nlm.nih.gov/15144888/) |
| KRTAP5-6 | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [PMID 15144888; UniProt Q6L8G9](https://pubmed.ncbi.nlm.nih.gov/15144888/) |
| KRTAP5-7 | specialized tissue expression | No effect | Reported in hair follicle, which the registry tiers disposable; only Vital or Protected tissues exclude. | [PMID 15144888; UniProt Q6L8G8](https://pubmed.ncbi.nlm.nih.gov/15144888/) |
| KRTAP5-9 | specialized tissue expression | No effect | Reported in hair follicle, which the registry tiers disposable; only Vital or Protected tissues exclude. | [PMID 15144888; UniProt P26371](https://pubmed.ncbi.nlm.nih.gov/15144888/) |
| KRTAP6-1 | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [PMID 12359730 (full-text table needed)](https://pubmed.ncbi.nlm.nih.gov/12359730/) |
| KRTAP6-2 | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [PMID 12359730 (full-text table needed)](https://pubmed.ncbi.nlm.nih.gov/12359730/) |
| KRTAP6-3 | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [PMID 12359730 (full-text table needed)](https://pubmed.ncbi.nlm.nih.gov/12359730/) |
| KRTAP9-1 | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [PMID 11279113 (family-level only)](https://pubmed.ncbi.nlm.nih.gov/11279113/) |
| KRTAP9-6 | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [PMID 11279113 (family-level only)](https://pubmed.ncbi.nlm.nih.gov/11279113/) |
| KRTAP9-7 | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [PMID 11279113 (family-level only)](https://pubmed.ncbi.nlm.nih.gov/11279113/) |
| LALBA | specialized tissue expression | No effect | Reported in breast, which the registry tiers manageable; only Vital or Protected tissues exclude. | [PMID 27322324; UniProt P00709](https://pubmed.ncbi.nlm.nih.gov/27322324/) |
| LCE3B | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [PMID 19169253; PMID 15854049; PMID 26556599](https://pubmed.ncbi.nlm.nih.gov/19169253/) |
| LCE3C | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [PMID 19169253 (induction, group-level); PMID 15854049 (group](https://pubmed.ncbi.nlm.nih.gov/19169253/) |
| LCN9 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| LELP1 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| LEUTX | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| LGALS16 | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [PMID 34944551 (atlas and microarray survey only)](https://pubmed.ncbi.nlm.nih.gov/34944551/) |
| LIPJ | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| LITAFD | no literature | No effect | No publication on normal adult expression was found. |  |
| LRRC52 | specialized tissue expression | No effect | The abstract names testis only; the skeletal-muscle statement is in a full text that could not be retrieved for a quote. | [PMID 22547800, full text PMC3356614, Fig. 4 and Discussion](https://pubmed.ncbi.nlm.nih.gov/22547800/); [PMID 31451634 (Lingle 2019, PNAS); supporting PMID 31348683 ](https://pubmed.ncbi.nlm.nih.gov/31451634/) |
| LRRIQ4 | no literature | No effect | No gene-specific publication on expression in any tissue; only database RNA annotation. No loss-of-function disease. | [UniProt A6NIV6 (LRRIQ4_HUMAN)](https://rest.uniprot.org/uniprotkb/search?query=gene_exact:LRRIQ4+AND+organism_id:9606+AND+reviewed:true&fields=accession,gene_names,cc_tissue_specificity,cc_disease); [Human Protein Atlas gene record LRRIQ4 (JSON)](https://www.proteinatlas.org/ENSG00000188306.json); [NCBI Gene 344657 LRRIQ4 (summary)](https://www.ncbi.nlm.nih.gov/gene/344657) |
| LUZP4 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| LYPD4 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| LYZL1 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| LYZL2 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| LYZL6 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| MAGEA1 | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [PMID 14734521, full text PMC2211762 (Fig. 2; Materials and M](https://pubmed.ncbi.nlm.nih.gov/14734521/); [PMID 23377668, full text PMC3581823 (Nanostring and Table 6 ](https://pubmed.ncbi.nlm.nih.gov/23377668/) |
| MAGEA10 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| MAGEA11 | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [PMID 15684378, full text PMC548016 (Fig. 12B and Discussion)](https://pubmed.ncbi.nlm.nih.gov/15684378/) |
| MAGEA2/B | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| MAGEA4 | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [PMID 14734521, full text PMC2211762 (Fig. 2)](https://pubmed.ncbi.nlm.nih.gov/14734521/) |
| MAGEA6 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| MAGEA8 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| MAGEA9/B | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| MAGEB1 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| MAGEB10 | no literature | No effect | No publication on normal adult expression was found. |  |
| MAGEB16 | no literature | No effect | No publication on normal adult expression was found. |  |
| MAGEB17 | no literature | No effect | No publication on normal adult expression was found. |  |
| MAGEB18 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| MAGEB3 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| MAGEB4 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| MAGEB5 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| MAGEB6B | no literature | No effect | No publication on normal adult expression was found. |  |
| MAGEC1 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| MAS1L | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| MBD3L2 | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [PMID 12504854 (Jiang et al., Genomics 2002) and PMID 1570160](https://rest.uniprot.org/uniprotkb/Q8NHZ7) |
| MBD3L2B | no literature | No effect | No publication on normal adult expression was found. |  |
| MBD3L3 | no literature | No effect | No publication on normal adult expression was found. |  |
| MBD3L4 | no literature | No effect | No publication on normal adult expression was found. |  |
| MBD3L5 | no literature | No effect | No publication on normal adult expression was found. |  |
| MEIOSIN | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| MINDY4B | no literature | No effect | No gene-specific publication found. | [Human Protein Atlas, MINDY4B gene and tissue pages](https://www.proteinatlas.org/ENSG00000214237-MINDY4B/tissue); [NCBI Gene 646951 MINDY4B summary](https://www.ncbi.nlm.nih.gov/gene/646951) |
| MOS | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| MPC1L | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| MRGPRG | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [https://rest.uniprot.org/uniprotkb/Q86SM5 (function field, b](https://rest.uniprot.org/uniprotkb/Q86SM5) |
| MROH2B | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| MSANTD5 | no literature | No effect | No publication on normal adult expression was found. |  |
| MSGN1 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| NANOGNB | no literature | No effect | No publication on normal adult expression was found. |  |
| NANOS2 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| NEU2 | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [PMID 23068092 (BBRC 2012, title only); original claim PMID 1](https://pubmed.ncbi.nlm.nih.gov/23068092/) |
| NKX2-6 | conflicting or nonspecific | No effect | Reports conflict, cannot separate the gene from paralogs, or come from cultured, diseased or fetal samples. | [Human Protein Atlas NKX2-6](https://www.proteinatlas.org/ENSG00000180053-NKX2-6); [Heathcote et al. 2005, Hum Mol Genet: Common arterial trunk ](https://pubmed.ncbi.nlm.nih.gov/15649947/); [Ta-Shma et al. 2014, J Med Genet: Conotruncal malformations ](https://pubmed.ncbi.nlm.nih.gov/24421281/) |
| NLRP11 | specialized tissue expression | No effect | Cultured, primed macrophages; not an in situ tissue measurement. | [PMID 35624206 (PMC9174058; DOI 10.1038/s41590-022-01220-3)](https://pubmed.ncbi.nlm.nih.gov/35624206/) |
| NLRP13 | no literature | No effect | No publication on normal adult expression was found. |  |
| NLRP4 | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [https://rest.uniprot.org/uniprotkb/Q96MN2.txt; PMID 22388039](https://rest.uniprot.org/uniprotkb/Q96MN2.txt) |
| NLRP8 | no normal expression evidence | No effect | Sources describe germline, placental or embryonic restriction, or did not test adult somatic tissue. | [Zhang P et al. Expression analysis of the NLRP gene family s](https://pmc.ncbi.nlm.nih.gov/articles/PMC2447171); [Ponsuksili S et al. Bovine NALP5, NALP8, and NALP9 genes ...](https://pubmed.ncbi.nlm.nih.gov/16339045/); [Human Protein Atlas ENSG00000179709-NLRP8](https://www.proteinatlas.org/ENSG00000179709-NLRP8) |
| NOTO | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| NUTM1 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| NY-ESO-1 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| OCM2 | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [PMID 42161597 (PMC13238953) and PMID 42123419 (PMC13163798),](https://pubmed.ncbi.nlm.nih.gov/42161597/) |
| OCSTAMP | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [PMID 41369694 (PMC12891103), mouse knockout; https://rest.un](https://rest.uniprot.org/uniprotkb/Q9BR26.txt) |
| ODF1 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| ODF4 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| OOSP1 | no normal expression evidence | No effect | Sources describe germline, placental or embryonic restriction, or did not test adult somatic tissue. | [Yan C et al. Oosp1 encodes a novel mouse oocyte-secreted pro](https://pubmed.ncbi.nlm.nih.gov/11747200/); [Tremblay K et al. Identification and characterization of a n](https://pubmed.ncbi.nlm.nih.gov/16580794/); [Human Protein Atlas gene page OOSP1](https://www.proteinatlas.org/ENSG00000284873-OOSP1) |
| OOSP2 | no literature | No effect | No publication on normal adult expression was found. |  |
| OOSP4B | no literature | No effect | No gene-specific publication found. | [Human Protein Atlas gene page OOSP4B (ENSG00000255393)](https://www.proteinatlas.org/ENSG00000255393-OOSP4B) |
| OR10A3 | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [PMID 24800820 (PMC4011832); gene-level data would be in its ](https://pubmed.ncbi.nlm.nih.gov/24800820/) |
| OR10A7 | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [PMID 24800820 (PMC4011832); supplementary tables not checked](https://pubmed.ncbi.nlm.nih.gov/24800820/) |
| OR10AG1 | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [PMID 24800820 (PMC4011832); supplementary tables not checked](https://pubmed.ncbi.nlm.nih.gov/24800820/) |
| OR10D3 | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [PMID 24800820 (PMC4011832); supplementary tables not checked](https://pubmed.ncbi.nlm.nih.gov/24800820/) |
| OR10G3 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| OR10G6 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| OR10G8 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| OR10G9 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| OR10H2 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| OR10H3 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| OR10H4 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| OR10J1 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| OR10Q1 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| OR10S1 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| OR10T2 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| OR10V1 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| OR10W1 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| OR10X1 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| OR11A1 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| OR11H4 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| OR11L1 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| OR12D3 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| OR13A1 | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [https://www.proteinatlas.org/ENSG00000256574-OR13A1 (JSON: h](https://www.proteinatlas.org/ENSG00000256574-OR13A1) |
| OR13C2 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| OR13C3 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| OR13C4 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| OR13C5 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| OR13C8 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| OR13C9 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| OR13F1 | no literature | No effect | No publication on normal adult expression was found. |  |
| OR13G1 | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [PMID 24800820 (Verbeurgt et al. 2014, PLoS One; PMC4011832),](https://pubmed.ncbi.nlm.nih.gov/24800820/) |
| OR14A16 | no literature | No effect | No publication on normal adult expression was found. |  |
| OR14A2 | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [PMID 38474024 (Bhomia et al. 2024, Int J Mol Sci; PMC1093159](https://pubmed.ncbi.nlm.nih.gov/38474024/) |
| OR14C36 | no literature | No effect | No publication on normal adult expression was found. |  |
| OR14K1 | no literature | No effect | No publication on normal adult expression was found. |  |
| OR1A2 | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [PMID 36503361 (Weidinger et al. 2022, Mol Med; PMC9743598)](https://pubmed.ncbi.nlm.nih.gov/36503361/); [PMID 38965271 (Danzer et al. 2024, Sci Rep; PMC11224223)](https://pubmed.ncbi.nlm.nih.gov/38965271/) |
| OR1B1 | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [PMID 38474024 (Bhomia et al. 2024, Int J Mol Sci; PMC1093159](https://pubmed.ncbi.nlm.nih.gov/38474024/) |
| OR1D5 | no literature | No effect | No publication on normal adult expression was found. |  |
| OR1E2 | no literature | No effect | No publication on normal adult expression was found. |  |
| OR1J1 | no literature | No effect | No publication on normal adult expression was found. |  |
| OR1J4 | no literature | No effect | No publication on normal adult expression was found. |  |
| OR1K1 | no literature | No effect | No publication on normal adult expression was found. |  |
| OR1L6 | no literature | No effect | No publication on normal adult expression was found. |  |
| OR1N1 | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [PMID 41317264 (Lin et al. 2026, Cardiovasc Drugs Ther; PMC13](https://pubmed.ncbi.nlm.nih.gov/41317264/) |
| OR1N2 | no literature | No effect | No publication on normal adult expression was found. |  |
| OR1R1P | no literature | No effect | No publication on normal adult expression was found. |  |
| OR1S1 | no literature | No effect | No publication on normal adult expression was found. |  |
| OR1S2 | no literature | No effect | No publication on normal adult expression was found. |  |
| OR2A14 | no literature | No effect | No publication on normal adult expression was found. |  |
| OR2A2 | no literature | No effect | No publication on normal adult expression was found. |  |
| OR2AG1 | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [PMID 27540365 (Kalbe et al. 2016, Front Physiol; PMC4972829)](https://pubmed.ncbi.nlm.nih.gov/27540365/) |
| OR2B2 | no literature | No effect | No publication on normal adult expression was found. |  |
| OR2B3 | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [PMID 42066673 (Ferdousi et al. 2026, Biomed Pharmacother; do](https://pubmed.ncbi.nlm.nih.gov/42066673/) |
| OR2C1 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| OR2D2 | no literature | No effect | No publication on normal adult expression was found. |  |
| OR2D3 | no literature | No effect | No publication on normal adult expression was found. |  |
| OR2F1 | no literature | No effect | No gene-specific publication found. | [Human Protein Atlas gene page OR2F1](https://www.proteinatlas.org/ENSG00000213215-OR2F1); [Expression Profile of Ectopic Olfactory Receptors Determined](https://pmc.ncbi.nlm.nih.gov/articles/PMC3566163) |
| OR2G2 | no literature | No effect | No publication on normal adult expression was found. |  |
| OR2G3 | no literature | No effect | No publication on normal adult expression was found. |  |
| OR2H1 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| OR2J3 | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [McRae JF et al. Chem Senses 2012;37(7):585-93, PMID 22714804](https://pubmed.ncbi.nlm.nih.gov/22714804/) |
| OR2K2 | no literature | No effect | No publication on normal adult expression was found. |  |
| OR2M7 | no literature | No effect | No publication on normal adult expression was found. |  |
| OR2T27 | no literature | No effect | No publication on normal adult expression was found. |  |
| OR2T29 | no literature | No effect | No publication on normal adult expression was found. |  |
| OR2T3 | no literature | No effect | No publication on normal adult expression was found. |  |
| OR2T33 | no literature | No effect | No publication on normal adult expression was found. |  |
| OR2T34 | no literature | No effect | No publication on normal adult expression was found. |  |
| OR2T4 | no literature | No effect | No publication on normal adult expression was found. |  |
| OR2T5 | no literature | No effect | No publication on normal adult expression was found. |  |
| OR2Y1 | no literature | No effect | No publication on normal adult expression was found. |  |
| OR4A15 | no literature | No effect | No publication on normal adult expression was found. |  |
| OR4A5 | no literature | No effect | No publication on normal adult expression was found. |  |
| OR4B1 | no literature | No effect | No publication on normal adult expression was found. |  |
| OR4C11 | no literature | No effect | No publication on normal adult expression was found. |  |
| OR4C13 | no literature | No effect | No publication on normal adult expression was found. |  |
| OR4C15 | no literature | No effect | No publication on normal adult expression was found. |  |
| OR4C3 | no literature | No effect | No gene-specific publication found. | [Human Protein Atlas OR4C3](https://www.proteinatlas.org/ENSG00000176547-OR4C3); [Khan et al. 2023, Diagnostics 13:1017 (abstract only)](https://pubmed.ncbi.nlm.nih.gov/36980325/) |
| OR4C46 | no literature | No effect | No publication on normal adult expression was found. |  |
| OR4C5 | no literature | No effect | No publication on normal adult expression was found. |  |
| OR4C6 | no literature | No effect | No publication on normal adult expression was found. |  |
| OR4D11 | no literature | No effect | No publication on normal adult expression was found. |  |
| OR4D2 | no literature | No effect | No publication on normal adult expression was found. |  |
| OR4D5 | no literature | No effect | No publication on normal adult expression was found. |  |
| OR4D6 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| OR4E1 | no literature | No effect | No publication on normal adult expression was found. |  |
| OR4E2 | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [PMID 38069004 (Hirata et al., Int J Mol Sci 2023; PMC1070674](https://pubmed.ncbi.nlm.nih.gov/38069004/) |
| OR4F3 | no literature | No effect | No publication on normal adult expression was found. |  |
| OR4F5 | no literature | No effect | No publication on normal adult expression was found. |  |
| OR4K1 | no literature | No effect | No publication on normal adult expression was found. |  |
| OR4L1 | no literature | No effect | No publication on normal adult expression was found. |  |
| OR4M1 | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [PMID 40806011 (PMC12348634) for the human biopsy data; the u](https://pubmed.ncbi.nlm.nih.gov/40806011/); [PMID 23241557 (Zhao et al., J Alzheimers Dis 2013; PMC396832](https://pubmed.ncbi.nlm.nih.gov/23241557/) |
| OR4M2 | no literature | No effect | No publication on normal adult expression was found. |  |
| OR4M2B | no literature | No effect | No publication on normal adult expression was found. |  |
| OR4N2 | no literature | No effect | No gene-specific publication found. | [Human Protein Atlas gene page OR4N2](https://www.proteinatlas.org/ENSG00000176294-OR4N2); [Human Protein Atlas OR4N2 single cell page](https://www.proteinatlas.org/ENSG00000176294-OR4N2/single+cell) |
| OR4Q3 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| OR4X1 | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [PMID 17964544 and PMID 37451613 (abstracts only; full texts ](https://pubmed.ncbi.nlm.nih.gov/17964544/) |
| OR51A7 | specialized tissue expression | No effect | One healthy donor; the probe was not tested against OR51 paralogs. | [PMID 37533582 (Murayama et al., JID Innov 2023; PMC10392076)](https://pubmed.ncbi.nlm.nih.gov/37533582/) |
| OR51B4 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| OR51C1 | no literature | No effect | No publication on normal adult expression was found. |  |
| OR51F1 | no literature | No effect | No publication on normal adult expression was found. |  |
| OR51F2 | no literature | No effect | No publication on normal adult expression was found. |  |
| OR51I1 | no literature | No effect | No publication on normal adult expression was found. |  |
| OR51I2 | no literature | No effect | No publication on normal adult expression was found. |  |
| OR51M1 | no literature | No effect | No publication on normal adult expression was found. |  |
| OR51Q1 | no literature | No effect | No publication on normal adult expression was found. |  |
| OR51T1 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| OR52B4 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| OR52B6 | no literature | No effect | No publication on normal adult expression was found. |  |
| OR52D1 | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [PMID 15647465 (abstract only)](https://pubmed.ncbi.nlm.nih.gov/15647465/) |
| OR52E2 | no literature | No effect | No publication on normal adult expression was found. |  |
| OR52J3 | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [PMID 23778302 / PMC3805820 (full text)](https://pubmed.ncbi.nlm.nih.gov/23778302/) |
| OR52L1 | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [PMID 23656994 (abstract only; full text not checked)](https://pubmed.ncbi.nlm.nih.gov/23656994/); [PMID 34070492 / PMC8198027 (full text, section 2.5)](https://pubmed.ncbi.nlm.nih.gov/34070492/) |
| OR52M1 | no literature | No effect | No publication on normal adult expression was found. |  |
| OR52R1 | no literature | No effect | No publication on normal adult expression was found. |  |
| OR52W1 | no literature | No effect | No publication on normal adult expression was found. |  |
| OR56A3 | no literature | No effect | No publication on normal adult expression was found. |  |
| OR56B2 | no literature | No effect | No publication on normal adult expression was found. |  |
| OR56B4 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| OR5AC2 | no literature | No effect | No publication on normal adult expression was found. |  |
| OR5AP2 | no literature | No effect | No publication on normal adult expression was found. |  |
| OR5AR1 | no literature | No effect | No publication on normal adult expression was found. |  |
| OR5B12 | no literature | No effect | No publication on normal adult expression was found. |  |
| OR5B17 | no literature | No effect | No publication on normal adult expression was found. |  |
| OR5B21 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| OR5BS1 | no literature | No effect | No publication on normal adult expression was found. |  |
| OR5F1 | no literature | No effect | No publication on normal adult expression was found. |  |
| OR5H1 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| OR5H15 | no literature | No effect | No publication on normal adult expression was found. |  |
| OR5H2 | no literature | No effect | No publication on normal adult expression was found. |  |
| OR5H6 | no literature | No effect | No publication on normal adult expression was found. |  |
| OR5J2 | no literature | No effect | No publication on normal adult expression was found. |  |
| OR5K2 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| OR5K3 | no literature | No effect | No publication on normal adult expression was found. |  |
| OR5K4 | no literature | No effect | No publication on normal adult expression was found. |  |
| OR5L1 | no literature | No effect | No publication on normal adult expression was found. |  |
| OR5M1 | no literature | No effect | No publication on normal adult expression was found. |  |
| OR5M11 | no literature | No effect | No publication on normal adult expression was found. |  |
| OR5M8 | no literature | No effect | No publication on normal adult expression was found. |  |
| OR5T3 | no literature | No effect | No publication on normal adult expression was found. |  |
| OR5V1 | no literature | No effect | No gene-specific publication found. | [Ge XY et al. Identification of Or5v1/Olfr110 as an oxylipin ](https://europepmc.org/article/MED/41570820); [Human Protein Atlas: OR5V1](https://www.proteinatlas.org/ENSG00000243729-OR5V1); [UniProtKB Q9UGF6 (OR5V1)](https://rest.uniprot.org/uniprotkb/search?query=gene_exact:OR5V1+AND+organism_id:9606&fields=accession,reviewed,gene_names,cc_tissue_specificity,cc_developmental_stage&format=tsv) |
| OR6C65 | no literature | No effect | No publication on normal adult expression was found. |  |
| OR6C68 | no literature | No effect | No publication on normal adult expression was found. |  |
| OR6C70 | no literature | No effect | No publication on normal adult expression was found. |  |
| OR6F1 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| OR6K2 | no literature | No effect | No publication on normal adult expression was found. |  |
| OR6M1 | no literature | No effect | No publication on normal adult expression was found. |  |
| OR6P1 | no literature | No effect | No publication on normal adult expression was found. |  |
| OR6T1 | no literature | No effect | No publication on normal adult expression was found. |  |
| OR7C2 | no literature | No effect | No publication on normal adult expression was found. |  |
| OR7E24 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| OR7G1 | no literature | No effect | No publication on normal adult expression was found. |  |
| OR7G3 | no literature | No effect | No publication on normal adult expression was found. |  |
| OR8B2 | no literature | No effect | No publication on normal adult expression was found. |  |
| OR8B3 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| OR8B4 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| OR8D2 | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [PMID 11705801](https://pubmed.ncbi.nlm.nih.gov/11705801/) |
| OR8D4 | no literature | No effect | No gene-specific publication found. | [Human Protein Atlas OR8D4 gene page](https://www.proteinatlas.org/ENSG00000181518-OR8D4) |
| OR8G1 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| OR8I2 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| OR8K1 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| OR8K5 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| OR8S1 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| OR8U1 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| OTOP1 | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [PMID 24379350](https://pubmed.ncbi.nlm.nih.gov/24379350/); [PMID 20554841; PMID 12651873](https://pubmed.ncbi.nlm.nih.gov/20554841/); [PMID 31543453 (mouse Otop1 knockout, type III taste cells); ](https://pubmed.ncbi.nlm.nih.gov/31543453/) |
| OTUD6A | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [PMID 38105755](https://pubmed.ncbi.nlm.nih.gov/38105755/) |
| P2RY4 | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [PMID 12644577](https://pubmed.ncbi.nlm.nih.gov/12644577/); [PMID 22437266](https://pubmed.ncbi.nlm.nih.gov/22437266/); [PMID 11690642; PMID 8617367 for the pancreas-restricted coun](https://pubmed.ncbi.nlm.nih.gov/11690642/) |
| P3R3URF | no literature | No effect | No publication on normal adult expression was found. |  |
| PABPC3 | no normal expression evidence | No effect | Sources describe germline, placental or embryonic restriction, or did not test adult somatic tissue. | [Feral C, Guellaen G, Pawlak A. Human testis expresses a spec](https://pmc.ncbi.nlm.nih.gov/articles/PMC37253); [UniProt Q9H361 (PABP3_HUMAN), Tissue specificity comment](https://rest.uniprot.org/uniprotkb/search?query=gene_exact:PABPC3+AND+organism_id:9606+AND+reviewed:true&fields=accession,gene_names,cc_tissue_specificity); [Human Protein Atlas ENSG00000151846-PABPC3](https://www.proteinatlas.org/ENSG00000151846-PABPC3) |
| PABPN1L | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [PMID 18483763 (abstract contradicts the claim)](https://pubmed.ncbi.nlm.nih.gov/18483763/) |
| PAGE1 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| PAGE2 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| PAGE3 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| PANX3 | specialized tissue expression | No effect | Reported in dental pulp, which the registry tiers manageable; only Vital or Protected tissues exclude. | [PMID 26263540](https://pubmed.ncbi.nlm.nih.gov/26263540/); [PMID 22947051 (human skin immunostaining); PMID 36813158 (mo](https://pubmed.ncbi.nlm.nih.gov/22947051/); [PMID 26138248; PMID 33426805](https://pubmed.ncbi.nlm.nih.gov/26138248/) |
| PASD1 | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [PMID 15162151](https://pubmed.ncbi.nlm.nih.gov/15162151/) |
| PATE1 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| PATE3 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| PDHA2 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| PFN3 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| PGAM4 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| PGK2 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| PIWIL3 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| PKDREJ | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| PMIS2 | no literature | No effect | No publication on normal adult expression was found. |  |
| PNMA6E | no literature | No effect | No publication on normal adult expression was found. |  |
| POM121L12 | no literature | No effect | No publication on normal adult expression was found. |  |
| POM121L2 | no literature | No effect | No gene-specific publication found. | [Human Protein Atlas POM121L2 tissue page](https://www.proteinatlas.org/ENSG00000158553-POM121L2/tissue); [Human Protein Atlas gene page POM121L2](https://www.proteinatlas.org/ENSG00000158553-POM121L2); [Aberg KA et al. A comprehensive family-based replication stu](https://pubmed.ncbi.nlm.nih.gov/23894747/) |
| POTEB/2 | no literature | No effect | No publication on normal adult expression was found. |  |
| POTEB3 | no literature | No effect | No publication on normal adult expression was found. |  |
| POTEC | no normal expression evidence | No effect | Sources describe germline, placental or embryonic restriction, or did not test adult somatic tissue. | [Barger CJ et al. Expression of the POTE gene family in human](https://pmc.ncbi.nlm.nih.gov/articles/PMC6244393/); [Bera TK et al. POTE, a highly homologous gene family ... exp](https://pmc.ncbi.nlm.nih.gov/articles/PMC139254/); [Bera TK et al. PNAS 2002 (same paper, paralog cloning)](https://pmc.ncbi.nlm.nih.gov/articles/PMC139254/) |
| POTED | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [PMID 12475935 (https://pmc.ncbi.nlm.nih.gov/articles/PMC1392](https://pmc.ncbi.nlm.nih.gov/articles/PMC139254/) |
| POTEF | conflicting or nonspecific | No effect | Reports conflict, cannot separate the gene from paralogs, or come from cultured, diseased or fetal samples. | [POTE, a highly homologous gene family located on numerous ch](https://pmc.ncbi.nlm.nih.gov/articles/PMC139254); [Evolution and expression of chimeric POTE-actin genes in the](https://pmc.ncbi.nlm.nih.gov/articles/PMC1693842/); [Barger et al., Expression of the POTE gene family in human o](https://www.med.upenn.edu/drapkinlab/assets/user-content/documents/Barger_et_al_Scientific_Reports 2018.pdf) |
| POTEH | no literature | No effect | No publication on normal adult expression was found. |  |
| PPIAL4D | no literature | No effect | No publication on normal adult expression was found. |  |
| PPIAL4E/F | no literature | No effect | No publication on normal adult expression was found. |  |
| PPP1R2B | no normal expression evidence | No effect | Sources describe germline, placental or embryonic restriction, or did not test adult somatic tissue. | [Korrodi-Gregorio et al. 2013, BMC Cell Biol 14:15 - Two dist](https://pmc.ncbi.nlm.nih.gov/articles/PMC3606321/); [UniProt Q6NXS1 (IPP2B_HUMAN), Tissue specificity (cites PubM](https://rest.uniprot.org/uniprotkb/search?query=gene_exact:PPP1R2B+AND+organism_id:9606&fields=accession,cc_tissue_specificity&format=tsv); [Human Protein Atlas PPP1R2B](https://www.proteinatlas.org/ENSG00000231989-PPP1R2B) |
| PPP3R2 | no normal expression evidence | No effect | Sources describe germline, placental or embryonic restriction, or did not test adult somatic tissue. | [Characterization of a human regulatory subunit of protein ph](https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=EXT_ID:15865209%20AND%20SRC:MED&resultType=core&format=json); [UniProt Q96LZ3 (PPP3R2)](https://rest.uniprot.org/uniprotkb/Q96LZ3.txt); [PPP3R2 sperm maturation/capacitation study in mouse (2020, P](https://pmc.ncbi.nlm.nih.gov/articles/PMC7493031/) |
| PPP4R3C | no literature | No effect | No publication on normal adult expression was found. |  |
| PRAME | conflicting or nonspecific | No effect | Conflicting reports, so no rule fires. Brain, kidney and skin are attributed to the 1997 discovery paper only through UniProt's annotation; its full text was not opened, so tissue levels are unverified. A 70-tissue IHC series (Kaczorowski 2022) found no staining in brain, kidney or adrenal, and Boon 2003 found little or no RNA in normal brain. Low adrenal staining (Lezcano 2021, abstract) is contradicted by the same series. No report in normal eye. PRAME-directed T-cell therapy trials report no dose-limiting on-target toxicity in these tissues. Neither negative series excludes a 1-4% neuronal subpopulation. | [UniProt P78395 (PRAME_HUMAN) tissue specificity](https://rest.uniprot.org/uniprotkb/search?query=gene_exact:PRAME+AND+organism_id:9606+AND+reviewed:true&fields=accession,gene_names,cc_tissue_specificity,cc_disease); [Ikeda H et al. Characterization of an antigen that is recogn](https://pubmed.ncbi.nlm.nih.gov/9047241/); [Kaczorowski M et al. PRAME Expression in Cancer. A Systemati](https://pmc.ncbi.nlm.nih.gov/articles/PMC9588667/) |
| PRAMEF1 | no literature | No effect | No publication on normal adult expression was found. |  |
| PRAMEF10 | no literature | No effect | No publication on normal adult expression was found. |  |
| PRAMEF11 | no literature | No effect | No publication on normal adult expression was found. |  |
| PRAMEF13 | no literature | No effect | No publication on normal adult expression was found. |  |
| PRAMEF14 | no literature | No effect | No publication on normal adult expression was found. |  |
| PRAMEF15 | no literature | No effect | No publication on normal adult expression was found. |  |
| PRAMEF17 | no literature | No effect | No publication on normal adult expression was found. |  |
| PRAMEF19 | no literature | No effect | No publication on normal adult expression was found. |  |
| PRAMEF2 | no literature | No effect | No publication on normal adult expression was found. |  |
| PRAMEF20 | no literature | No effect | No publication on normal adult expression was found. |  |
| PRAMEF25/6 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| PRAMEF4 | no literature | No effect | No publication on normal adult expression was found. |  |
| PRAMEF5 | no literature | No effect | No publication on normal adult expression was found. |  |
| PRAMEF6 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| PRDM14 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| PRDM7 | conflicting or nonspecific | No effect | Reports conflict, cannot separate the gene from paralogs, or come from cultured, diseased or fetal samples. | [Fumasoni I et al. Family expansion and gene rearrangements c](https://pmc.ncbi.nlm.nih.gov/articles/PMC2082429); [Fumasoni I et al. 2007 (PMC2082429) - qPCR panel](https://pmc.ncbi.nlm.nih.gov/articles/PMC2082429); [Fumasoni I et al. 2007 (PMC2082429) - abstract](https://pmc.ncbi.nlm.nih.gov/articles/PMC2082429) |
| PRDM9 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| PRKACG | specialized tissue expression | No effect | R3 needs expression in the affected tissue; megakaryocyte expression could not be checked. One family with one missense variant. | [Manchev VT et al., Blood 2014;124(16):2554-63, PMID 25061177](https://pubmed.ncbi.nlm.nih.gov/25061177/) |
| PRM1 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| PRM2 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| PROKR1 | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [Lin DC et al., J Biol Chem 2002;277(22):19276-80, PMID 11886](https://pubmed.ncbi.nlm.nih.gov/11886876/) |
| PRR20A/B/C/D/E | no literature | No effect | No publication on normal adult expression was found. |  |
| PRR20G | no literature | No effect | No publication on normal adult expression was found. |  |
| PRR23A | no literature | No effect | No publication on normal adult expression was found. |  |
| PRR23B | no literature | No effect | No publication on normal adult expression was found. |  |
| PRR23C | no literature | No effect | No publication on normal adult expression was found. |  |
| PRR23D1/2 | no literature | No effect | No publication on normal adult expression was found. |  |
| PRR23E | no literature | No effect | No publication on normal adult expression was found. |  |
| PRR30 | no literature | No effect | No publication on normal adult expression was found. |  |
| PRR5-ARHGAP8 | conflicting or nonspecific | No effect | Reports conflict, cannot separate the gene from paralogs, or come from cultured, diseased or fetal samples. | [Expression of acidosis-dependent genes in human cancer nests](https://europepmc.org/article/MED/25279216); [Bipolar disorder with binge eating behavior: a genome-wide a](https://europepmc.org/article/PMC/PMC5804024); [Human Protein Atlas PRR5-ARHGAP8 gene page](https://www.proteinatlas.org/ENSG00000248405-PRR5-ARHGAP8) |
| PRSS37 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| PRSS38 | no literature | No effect | No publication on normal adult expression was found. |  |
| PRSS54 | no literature | No effect | No publication on normal adult expression was found. |  |
| PRSS58 | no literature | No effect | No publication on normal adult expression was found. |  |
| PSG11 | conflicting or nonspecific | No effect | Reports conflict, cannot separate the gene from paralogs, or come from cultured, diseased or fetal samples. | [Houston et al. 2016, mAbs: Pregnancy-specific glycoprotein e](https://pmc.ncbi.nlm.nih.gov/articles/PMC4966852); [Arnold et al. 1999, Am J Reprod Immunol: Pregnancy-specific ](https://pubmed.ncbi.nlm.nih.gov/10326620/); [Endoh et al. 2009, Biogerontology: Coordinate expression of ](https://pubmed.ncbi.nlm.nih.gov/18792801/) |
| PSG3 | conflicting or nonspecific | No effect | Reports conflict, cannot separate the gene from paralogs, or come from cultured, diseased or fetal samples. | [Houston A et al. Pregnancy-specific glycoprotein expression ](https://pmc.ncbi.nlm.nih.gov/articles/PMC4966852); [Human Protein Atlas ENSG00000221826-PSG3](https://www.proteinatlas.org/ENSG00000221826-PSG3); [UniProt Q16557 (PSG3_HUMAN)](https://rest.uniprot.org/uniprotkb/search?query=gene_exact:PSG3+AND+organism_id:9606+AND+reviewed:true&fields=accession,gene_names,cc_tissue_specificity,cc_developmental_stage) |
| PSG7 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| PSG8 | conflicting or nonspecific | No effect | Reports conflict, cannot separate the gene from paralogs, or come from cultured, diseased or fetal samples. | [Houston A et al. Pregnancy-specific glycoprotein expression ](https://pmc.ncbi.nlm.nih.gov/articles/PMC4966852); [Same paper, antibody cross-reactivity](https://pmc.ncbi.nlm.nih.gov/articles/PMC4966852); [Same paper, expression level](https://pmc.ncbi.nlm.nih.gov/articles/PMC4966852) |
| PSMA8 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| PTTG1IP2 | no literature | No effect | No publication on normal adult expression was found. |  |
| RAB6D | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| RAD21L1 | no normal expression evidence | No effect | Sources describe germline, placental or embryonic restriction, or did not test adult somatic tissue. | [A new meiosis-specific cohesin complex implicated in the coh](https://pmc.ncbi.nlm.nih.gov/articles/PMC3059921); [RAD21L is a novel kleisin subunit of the cohesin complex (co](https://pmc.ncbi.nlm.nih.gov/articles/PMC3184609); [Europe PMC abstract record, Genes (Basel) 2024, PMID 3933675](https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=TITLE_ABS:(TRIM60%20OR%20RNF129%20OR%20RAD21L1%20OR%20SPAG11A)&resultType=core&format=json&pageSize=40) |
| RBMXL2 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| RBMXL3 | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [PMID 33719563, DOI 10.1152/ajplung.00022.2020 (Europe PMC ab](https://rest.uniprot.org/uniprotkb/Q8N7X1.txt) |
| RBMY1A1 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| RBMY1B | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| RBMY1D | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| RBMY1E | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| RBMY1F/J | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| RFPL3 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| RFPL4AL1 | no literature | No effect | No publication on normal adult expression was found. |  |
| RFPL4B | no literature | No effect | No publication on normal adult expression was found. |  |
| RGPD3 | conflicting or nonspecific | No effect | Reports conflict, cannot separate the gene from paralogs, or come from cultured, diseased or fetal samples. | [Ciccarelli et al. 2005, Genome Res 15:343 - Complex genomic ](https://pmc.ncbi.nlm.nih.gov/articles/PMC551560); [Human Protein Atlas RGPD3 (summary and tissue pages)](https://www.proteinatlas.org/ENSG00000153165-RGPD3/tissue); [UniProt A6NKT7 (RGPD3_HUMAN)](https://rest.uniprot.org/uniprotkb/search?query=gene_exact:RGPD3+AND+organism_id:9606&fields=accession,gene_names,cc_tissue_specificity&format=tsv) |
| RGPD4 | conflicting or nonspecific | No effect | Reports conflict, cannot separate the gene from paralogs, or come from cultured, diseased or fetal samples. | [Ciccarelli FD et al. 2005, Complex genomic rearrangements le](https://pmc.ncbi.nlm.nih.gov/articles/PMC551560/); [Ciccarelli FD et al. 2005 (same paper)](https://pmc.ncbi.nlm.nih.gov/articles/PMC551560/); [Human Protein Atlas, RGPD4 single cell page](https://www.proteinatlas.org/ENSG00000196862-RGPD4/single+cell) |
| RGSL1 | no literature | No effect | No gene-specific publication found. | [Human Protein Atlas gene page RGSL1](https://www.proteinatlas.org/ENSG00000121446-RGSL1); [Exploring Sex Differences in Type 2 Diabetes via a Male-Domi](https://pmc.ncbi.nlm.nih.gov/articles/PMC12602021/) |
| RHOXF2B | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| RIPPLY3 | conflicting or nonspecific | No effect | Reports conflict, cannot separate the gene from paralogs, or come from cultured, diseased or fetal samples. | [Shibuya K et al. Isolation of two novel genes, DSCR5 and DSC](https://europepmc.org/article/MED/10814524); [UniProtKB P57055 (RIPPLY3)](https://rest.uniprot.org/uniprotkb/search?query=gene_exact:RIPPLY3+AND+organism_id:9606&fields=accession,reviewed,gene_names,cc_tissue_specificity,cc_developmental_stage&format=tsv); [Human Protein Atlas: RIPPLY3](https://www.proteinatlas.org/ENSG00000183145-RIPPLY3) |
| RLN1 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| RNASE11 | no literature | No effect | No publication on normal adult expression was found. |  |
| RNASE12 | no literature | No effect | No publication on normal adult expression was found. |  |
| RNASE13 | no literature | No effect | No publication on normal adult expression was found. |  |
| RNASE9 | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [PMID 18992174 / PMC2669477 (the RT-PCR claim); PMID 19137000](https://pubmed.ncbi.nlm.nih.gov/18992174/) |
| RNF151 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| RPS4Y2 | no literature | No effect | No publication on normal adult expression was found. |  |
| RSPH6A | no normal expression evidence | No effect | Sources describe germline, placental or embryonic restriction, or did not test adult somatic tissue. | [RSPH6A is required for sperm flagellum formation and male fe](https://ir.library.osaka-u.ac.jp/repo/ouka/all/78569/JCellSci_131_19_jcs.221648.pdf); [A mammalian radial spokehead-like gene, RSHL1, at the myoton](https://europepmc.org/article/MED/11237735); [Human Protein Atlas RSPH6A gene page](https://www.proteinatlas.org/ENSG00000104941-RSPH6A) |
| RTP2 | specialized tissue expression | No effect | The abstract reports the RTP and REEP families in taste papillae; the gene-specific statement is a database annotation. | [PMID 16720576 (DOI 10.1074/jbc.M513637200), abstract only; U](https://pubmed.ncbi.nlm.nih.gov/16720576/); [PMID 28262096 / PMC5362263 (DOI 10.7554/eLife.21895); UniPro](https://pubmed.ncbi.nlm.nih.gov/28262096/) |
| SAGE1 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| SATL1 | conflicting or nonspecific | No effect | Reports conflict, cannot separate the gene from paralogs, or come from cultured, diseased or fetal samples. | [Alayoubi et al. 2024, Sci Rep: Loss-of-function variant in s](https://pmc.ncbi.nlm.nih.gov/articles/PMC10923806/); [Human Protein Atlas SATL1](https://www.proteinatlas.org/ENSG00000184788-SATL1) |
| SCGB1D4 | specialized tissue expression | No effect | Induced expression; the resting pattern conflicts between sources. | [PMID 15034037 (DOI 10.4049/jimmunol.172.7.4245), abstract; c](https://pubmed.ncbi.nlm.nih.gov/15034037/) |
| SCP2D1 | no literature | No effect | No publication on normal adult expression was found. |  |
| SCYGR2 | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [UniProt A0A286YFB4; HPA ENSG00000284643 (skin enriched, 1.6 ](https://rest.uniprot.org/uniprotkb/A0A286YFB4) |
| SCYGR3 | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [UniProt A0A286YF60; HPA ENSG00000284704 (skin enriched)](https://rest.uniprot.org/uniprotkb/A0A286YF60) |
| SCYGR5 | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [UniProt A0A286YF46; HPA ENSG00000284667 (not detected)](https://rest.uniprot.org/uniprotkb/A0A286YF46) |
| SCYGR6 | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [UniProt A0A286YF77; HPA ENSG00000284725 (skin enriched)](https://rest.uniprot.org/uniprotkb/A0A286YF77) |
| SCYGR7 | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [UniProt A0A286YF01; HPA ENSG00000284718 (skin enriched)](https://rest.uniprot.org/uniprotkb/A0A286YF01) |
| SCYGR8 | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [UniProt A0A286YFG1; HPA ENSG00000284635 (skin enriched, 3.1 ](https://rest.uniprot.org/uniprotkb/A0A286YFG1) |
| SEBOX | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| SEMG2 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| SIGLECL1 | no literature | No effect | No publication on normal adult expression was found. |  |
| SLC25A52 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| SLC35G3 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| SLC35G4 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| SLC35G5 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| SLC35G6 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| SLC36A3 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| SMCP | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| SMIM21 | no literature | No effect | No publication on normal adult expression was found. |  |
| SMIM28 | no literature | No effect | No publication on normal adult expression was found. |  |
| SMIM42 | no literature | No effect | No publication on normal adult expression was found. |  |
| SMIM9 | no literature | No effect | No gene-specific publication found. | [Oyama Y et al. CRISPR/Cas9-mediated genome editing reveals 1](https://pmc.ncbi.nlm.nih.gov/articles/PMC9226692); [Human Protein Atlas gene page SMIM9](https://www.proteinatlas.org/ENSG00000203870-SMIM9) |
| SP7 | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [PMID 20579626 (Lapunzina 2010, Am J Hum Genet); UniProt Q8TD](https://pubmed.ncbi.nlm.nih.gov/20579626/); [PMID 20615976 (Zhou 2010, PNAS)](https://pubmed.ncbi.nlm.nih.gov/20615976/); [PMID 26599722 (He 2016, Int Endod J), a review of mouse prim](https://pubmed.ncbi.nlm.nih.gov/26599722/) |
| SPACA1 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| SPADH | no literature | No effect | No publication on normal adult expression was found. |  |
| SPAG11A | specialized tissue expression | No effect | Evidence is from another species, a database annotation or an assay shared with a paralog, or the tissue is reproductive. | [Human Protein Atlas SPAG11A gene and tissue pages](https://www.proteinatlas.org/ENSG00000178287-SPAG11A/tissue); [Yamaguchi & Ouchi, Antimicrobial peptide defensin: identific](https://pmc.ncbi.nlm.nih.gov/articles/PMC3406309); [Europe PMC abstract record, JBRA Assist Reprod 2024, PMID 38](https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=TITLE_ABS:(TRIM60%20OR%20RNF129%20OR%20RAD21L1%20OR%20SPAG11A)&resultType=core&format=json&pageSize=40) |
| SPAG11B | specialized tissue expression | No effect | Evidence is from another species, a database annotation or an assay shared with a paralog, or the tissue is reproductive. | [Human Protein Atlas SPAG11B tissue / primary data](https://www.proteinatlas.org/ENSG00000164871-SPAG11B/tissue/primary+data); [Osterhoff et al. 1994, Biol Reprod 50:516 - novel human sper](https://pubmed.ncbi.nlm.nih.gov/8167223/); [Hamil et al. 2000, Endocrinology 141:1245 - HE2beta and HE2g](https://pubmed.ncbi.nlm.nih.gov/10698202/) |
| SPAM1 | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [PMID 17446714 (Grigorieva 2007, mouse); PMID 8575780 (Jones ](https://pubmed.ncbi.nlm.nih.gov/17446714/) |
| SPANXA1/2 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| SPANXB1 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| SPANXC | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| SPANXD | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| SPANXN1 | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [PMID 17406683 full text, Table 2 (PMC1831492)](https://pubmed.ncbi.nlm.nih.gov/17406683/) |
| SPANXN2 | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [PMID 17406683 full text, Table 2 and Figure 6b (PMC1831492)](https://pubmed.ncbi.nlm.nih.gov/17406683/) |
| SPANXN3 | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [PMID 17406683 full text, Table 2 (PMC1831492)](https://pubmed.ncbi.nlm.nih.gov/17406683/) |
| SPANXN4 | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [PMID 17406683 full text, Table 2 (PMC1831492)](https://pubmed.ncbi.nlm.nih.gov/17406683/) |
| SPANXN5 | no literature | No effect | No publication on normal adult expression was found. |  |
| SPATA19 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| SPATA3 | no normal expression evidence | No effect | Sources describe germline, placental or embryonic restriction, or did not test adult somatic tissue. | [Roles of human SPATA3 in cell proliferation and expression p](https://www.spandidos-publications.com/10.3892/mmr.2025.13620); [Human Protein Atlas gene page SPATA3](https://www.proteinatlas.org/ENSG00000173699-SPATA3) |
| SPATA31A1 | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [PMID 28264649 (PMC5338094), Bekpen et al., BMC Genomics 2017](https://pubmed.ncbi.nlm.nih.gov/28264649/) |
| SPATA31A3 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| SPATA31A5 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| SPATA31A6 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| SPATA31A7 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| SPATA31D1 | no literature | No effect | No publication on normal adult expression was found. |  |
| SPATA31D3 | no literature | No effect | No publication on normal adult expression was found. |  |
| SPATA31D4 | no literature | No effect | No publication on normal adult expression was found. |  |
| SPATA31F1 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| SPATA31F3 | no literature | No effect | No publication on normal adult expression was found. |  |
| SPEM1 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| SPEM2 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| SPEM3 | no literature | No effect | No publication on normal adult expression was found. |  |
| SPMAP1 | no literature | No effect | No publication on normal adult expression was found. |  |
| SPMIP10 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| SPMIP9 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| SPO11 | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [PMID 11807770 (Klein et al., Eur J Immunol 2002) for the kno](https://pubmed.ncbi.nlm.nih.gov/11807770/) |
| SPPL2C | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| SRMS | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [PMID 37523817 (Kikuchi et al., BBRC 2023); PMID 7935409 (Koh](https://pubmed.ncbi.nlm.nih.gov/37523817/) |
| SSMEM1 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| SSU72L1 | no literature | No effect | No publication on normal adult expression was found. |  |
| SSU72L2 | no literature | No effect | No publication on normal adult expression was found. |  |
| SSU72L3 | no literature | No effect | No publication on normal adult expression was found. |  |
| SSU72L4 | no literature | No effect | No publication on normal adult expression was found. |  |
| SSU72L5 | no literature | No effect | No publication on normal adult expression was found. |  |
| SSU72L6 | no literature | No effect | No publication on normal adult expression was found. |  |
| SSX1 | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [PMID 15781633](https://pubmed.ncbi.nlm.nih.gov/15781633/) |
| SSX2/B | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [UniProt Q16385 / Q16384, citing PMID 7539744](https://pubmed.ncbi.nlm.nih.gov/7539744/) |
| SSX3 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| SSX4/B | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| SSX5 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| STARD6 | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [PMID 17885612](https://pubmed.ncbi.nlm.nih.gov/17885612/) |
| STPG3 | no literature | No effect | No publication on normal adult expression was found. |  |
| SUN5 | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [UniProt Q8TC36; PMID 12621555](https://pubmed.ncbi.nlm.nih.gov/12621555/) |
| SUPT20HL1 | no literature | No effect | No publication on normal adult expression was found. |  |
| SYCP1 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| TAAR2 | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [PMID 23315425](https://pubmed.ncbi.nlm.nih.gov/23315425/); [UniProt Q9P1P5](https://rest.uniprot.org/uniprotkb/search?query=organism_id:9606+AND+reviewed:true+AND+gene_exact:TAAR2) |
| TAAR9 | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [UniProt Q96RI9](https://rest.uniprot.org/uniprotkb/search?query=organism_id:9606+AND+reviewed:true+AND+gene_exact:TAAR9) |
| TAS1R2 | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [PMID 19352508](https://pubmed.ncbi.nlm.nih.gov/19352508/); [PMID 17724330; UniProt Q8TE23 for taste function](https://pubmed.ncbi.nlm.nih.gov/17724330/) |
| TAS2R1 | specialized tissue expression | No effect | Evidence is from another species, a database annotation or an assay shared with a paralog, or the tissue is reproductive. | [UniProtKB Q9NYW7 (TAS2R1)](https://rest.uniprot.org/uniprotkb/search?query=gene_exact:TAS2R1+AND+organism_id:9606&fields=accession,reviewed,gene_names,cc_tissue_specificity,cc_developmental_stage&format=tsv); [Upadhyaya JD et al. Dextromethorphan mediated bitter taste r](https://pmc.ncbi.nlm.nih.gov/articles/PMC4207743); [Jaggupilli A et al. Analysis of the expression of human bitt](https://europepmc.org/article/MED/28012014) |
| TAS2R13 | specialized tissue expression | No effect | Family-level source for all 25 bitter receptors; not opened per gene. | [PMID 18003842; UniProt Q9NYV9](https://pubmed.ncbi.nlm.nih.gov/18003842/); [PMID 20972434 full text, Table 1 (PMC3066567)](https://pubmed.ncbi.nlm.nih.gov/20972434/) |
| TAS2R41 | specialized tissue expression | No effect | Family-level source for all 25 bitter receptors; not opened per gene. | [PMID 18003842; UniProt P59536](https://pubmed.ncbi.nlm.nih.gov/18003842/) |
| TAS2R50 | specialized tissue expression | No effect | Family-level source for all 25 bitter receptors; not opened per gene. | [PMID 18003842; UniProt P59544](https://pubmed.ncbi.nlm.nih.gov/18003842/); [PMID 33461297](https://pubmed.ncbi.nlm.nih.gov/33461297/) |
| TAS2R60 | specialized tissue expression | No effect | Family-level source for all 25 bitter receptors; not opened per gene. | [PMID 18003842; UniProt P59551](https://pubmed.ncbi.nlm.nih.gov/18003842/) |
| TAS2R7 | specialized tissue expression | No effect | Family-level source for all 25 bitter receptors; not opened per gene. | [PMID 18003842; UniProt Q9NYW3](https://pubmed.ncbi.nlm.nih.gov/18003842/); [PMID 31632299 (PMC6783802)](https://pubmed.ncbi.nlm.nih.gov/31632299/) |
| TAS2R8 | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [PMID 18003842 (Behrens et al., J Neurosci 2007; PMCID PMC667](https://pubmed.ncbi.nlm.nih.gov/18003842/) |
| TAS2R9 | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [PMID 18003842 (Behrens et al., J Neurosci 2007) for expressi](https://pubmed.ncbi.nlm.nih.gov/18003842/) |
| TBC1D21 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| TBC1D28 | no literature | No effect | No publication on normal adult expression was found. |  |
| TBX22 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| TBXT | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| TCP11 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| TCP11X1 | no literature | No effect | No publication on normal adult expression was found. |  |
| TCP11X2 | no literature | No effect | No publication on normal adult expression was found. |  |
| TDRD15 | conflicting or nonspecific | No effect | Reports conflict, cannot separate the gene from paralogs, or come from cultured, diseased or fetal samples. | [Tdrd15 is dispensable for male fertility and spermatogenesis](https://europepmc.org/article/PMC/PMC12548254); [Human Protein Atlas TDRD15 tissue page](https://www.proteinatlas.org/ENSG00000218819-TDRD15/tissue); [Human Protein Atlas TDRD15 JSON record](https://www.proteinatlas.org/ENSG00000218819.json) |
| TEDDM1 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| TENT5D | no normal expression evidence | No effect | Sources describe germline, placental or embryonic restriction, or did not test adult somatic tissue. | [Bettoni et al. 2009, Genomics: Identification of FAM46D as a](https://pubmed.ncbi.nlm.nih.gov/19540335/); [UniProt Q8NEK8 (TENT5D_HUMAN), tissue specificity field via ](https://rest.uniprot.org/uniprotkb/search?query=gene_exact:TENT5D+AND+organism_id:9606&fields=accession,gene_names,cc_tissue_specificity&format=tsv); [Cong et al. 2022, Cell Discov: Deficiency of X-linked TENT5D](https://pmc.ncbi.nlm.nih.gov/articles/PMC8901658/) |
| TEX13A | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| TEX13B | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| TEX13C | no literature | No effect | No publication on normal adult expression was found. |  |
| TEX13D | no literature | No effect | No publication on normal adult expression was found. |  |
| TEX19 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| TEX28 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| TEX38 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| TEX44 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| TEX46 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| TEX50 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| TEX51 | no literature | No effect | No publication on normal adult expression was found. |  |
| TEX55 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| TFDP3 | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [PMID 12097419 (Wang et al., J Immunol 2002)](https://pubmed.ncbi.nlm.nih.gov/12097419/); [PMID 20559320 (PMCID PMC3131880)](https://pubmed.ncbi.nlm.nih.gov/20559320/) |
| TGIF2LX | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| TGIF2LY | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| TGM4 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| TKTL2 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| TLE7 | no literature | No effect | No publication on normal adult expression was found. |  |
| TMCO2 | no literature | No effect | No publication on normal adult expression was found. |  |
| TMEM202 | no literature | No effect | No publication on normal adult expression was found. |  |
| TMEM225 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| TMEM239 | no literature | No effect | No publication on normal adult expression was found. |  |
| TMEM262 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| TMEM269 | no literature | No effect | No gene-specific publication found. | [Human Protein Atlas ENSG00000274386-TMEM269](https://www.proteinatlas.org/ENSG00000274386-TMEM269); [Wikipedia: TMEM269](https://en.wikipedia.org/wiki/TMEM269) |
| TMEM270 | no literature | No effect | No publication on normal adult expression was found. |  |
| TMEM31 | no literature | No effect | No publication on normal adult expression was found. |  |
| TMEM89 | no literature | No effect | No publication on normal adult expression was found. |  |
| TMEM95 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| TNP2 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| TOMM20L | no literature | No effect | No publication on normal adult expression was found. |  |
| TPD52L3 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| TPRX1 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| TPRX2 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| TPTE | no normal expression evidence | No effect | Sources describe germline, placental or embryonic restriction, or did not test adult somatic tissue. | [Chen H et al. A testis-specific gene, TPTE, encodes a putati](https://pubmed.ncbi.nlm.nih.gov/10598804/); [Dong X-Y et al. Identification of two novel CT antigens and ](https://pmc.ncbi.nlm.nih.gov/articles/PMC2394243); [Tapparel C et al. The TPTE gene family: cellular expression,](https://pubmed.ncbi.nlm.nih.gov/14659893/) |
| TRIM42 | no literature | No effect | No publication on normal adult expression was found. |  |
| TRIM43 | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [PMID 30420784 (PMC6294671), Full et al., Nat Microbiol 2019;](https://pubmed.ncbi.nlm.nih.gov/30420784/) |
| TRIM43B | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [PMID 30420784 (PMC6294671), methods RNA-seq section; UniProt](https://pubmed.ncbi.nlm.nih.gov/30420784/) |
| TRIM48 | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [PMID 40609779 (Wu et al., Cell Signal 2025, doi 10.1016/j.ce](https://pubmed.ncbi.nlm.nih.gov/40609779/) |
| TRIM49 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| TRIM49B | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| TRIM49C | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| TRIM51G | no literature | No effect | No publication on normal adult expression was found. |  |
| TRIM60 | no normal expression evidence | No effect | Sources describe germline, placental or embryonic restriction, or did not test adult somatic tissue. | [Human Protein Atlas TRIM60 gene and tissue pages](https://www.proteinatlas.org/ENSG00000176979-TRIM60/tissue); [TRIM60 negatively regulates TAK1 signalosome (Cell Mol Immun](https://pmc.ncbi.nlm.nih.gov/articles/PMC8322076/); [PMID 33184450 (Gu et al., Cell Mol Immunol 2021;18:1981, PMC](https://pubmed.ncbi.nlm.nih.gov/33184450/) |
| TRIM64C | no literature | No effect | No publication on normal adult expression was found. |  |
| TRIM77 | no literature | No effect | No gene-specific publication found. | [Human Protein Atlas TRIM77 (summary and tissue pages)](https://www.proteinatlas.org/ENSG00000214414-TRIM77/tissue); [Peng et al. 2021, Front Oncol 11:626190 (abstract only)](https://pubmed.ncbi.nlm.nih.gov/34993129/) |
| TRIML1 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| TSGA10IP | no literature | No effect | No gene-specific publication found. | [Human Protein Atlas gene record TSGA10IP (JSON)](https://www.proteinatlas.org/ENSG00000175513.json); [Levine TP. Structural bioinformatics predicts that the Retin](https://pmc.ncbi.nlm.nih.gov/articles/PMC7551519/); [Wang M et al. Casein kinase-2 inhibition promotes retinal ga](https://pmc.ncbi.nlm.nih.gov/articles/PMC10749609/) |
| TSGA13 | no literature | No effect | No publication on normal adult expression was found. |  |
| TSPY1 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| TSPY10/3 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| TSPY2 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| TSPY4 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| TSPY8 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| TSPY9 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| TSSK1B | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| TSSK2 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| TTLL13 | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [PMID 40593606 (Sun et al., Nat Commun 2025;16:5591, PMC12219](https://pubmed.ncbi.nlm.nih.gov/40593606/) |
| TTLL2 | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [UniProt A4Q9E4 citing PMID 17499049 (van Dijk et al., Mol Ce](https://pubmed.ncbi.nlm.nih.gov/17499049/) |
| TTLL8 | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [UniProt A4Q9F1 tissue specificity citing PMID 23897886 (epen](https://pubmed.ncbi.nlm.nih.gov/23897886/) |
| TUBB8 | no normal expression evidence | No effect | Sources describe germline, placental or embryonic restriction, or did not test adult somatic tissue. | [Mutations in TUBB8 and Human Oocyte Meiotic Arrest (Feng et ](https://pmc.ncbi.nlm.nih.gov/articles/PMC4767273); [UniProt Q3ZCM7 (TUBB8)](https://rest.uniprot.org/uniprotkb/Q3ZCM7.txt); [Human Protein Atlas gene page TUBB8 (summary and tissue page](https://www.proteinatlas.org/ENSG00000261456-TUBB8/tissue) |
| TXNDC8 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| UBQLN3 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| UBTFL1 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| URGCP-MRPS24 | no literature | No effect | No gene-specific publication found. | [Human Protein Atlas: URGCP-MRPS24](https://www.proteinatlas.org/ENSG00000270617-URGCP-MRPS24) |
| USP17L1 | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [PMID 14699124 (DUB-3/USP17L2); UniProt Q6R6M4 induction anno](https://pubmed.ncbi.nlm.nih.gov/14699124/) |
| USP17L3 | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [PMID 14699124; PMID 20388806 (both USP17L2/family-level)](https://pubmed.ncbi.nlm.nih.gov/14699124/) |
| USP17L4 | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [PMID 14699124 (USP17L2/family-level)](https://pubmed.ncbi.nlm.nih.gov/14699124/) |
| USP17L7 | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [PMID 14699124 (USP17L2/family-level)](https://pubmed.ncbi.nlm.nih.gov/14699124/) |
| USP26 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| VCX | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| VCX2 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| VCX3A | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [PMID 18076704 and PMID 18194880 (contradicting follow-ups); ](https://pubmed.ncbi.nlm.nih.gov/18076704/) |
| VCX3B | no normal expression evidence | No effect | Sources describe germline, placental or embryonic restriction, or did not test adult somatic tissue. | [Lahn BT, Page DC 2000, A human sex-chromosomal gene family e](https://pubmed.ncbi.nlm.nih.gov/10607842/); [UniProt Q9H321 (VCX3B_HUMAN), Tissue specificity comment](https://rest.uniprot.org/uniprotkb/Q9H321); [NCBI Gene 425054 VCX3B RefSeq summary](https://www.ncbi.nlm.nih.gov/gene/425054) |
| VCY/1B | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| VHLL | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [PMID 14757845](https://pubmed.ncbi.nlm.nih.gov/14757845/) |
| VN1R2 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| VN1R4 | no literature | No effect | No publication on normal adult expression was found. |  |
| VRTN | conflicting or nonspecific | No effect | Reports conflict, cannot separate the gene from paralogs, or come from cultured, diseased or fetal samples. | [VRTN is Required for the Development of Thoracic Vertebrae i](https://europepmc.org/article/PMC/PMC6001657); [Identification of marker genes to monitor residual iPSCs in ](https://europepmc.org/article/MED/36319564); [Human Protein Atlas VRTN tissue page](https://www.proteinatlas.org/ENSG00000133980-VRTN/tissue) |
| WEE2 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| WFDC10A | no literature | No effect | No publication on normal adult expression was found. |  |
| WFDC13 | specialized tissue expression | No effect | Evidence is from another species, a database annotation or an assay shared with a paralog, or the tissue is reproductive. | [Human Protein Atlas gene page WFDC13 (ENSG00000168634)](https://www.proteinatlas.org/ENSG00000168634-WFDC13); [Jain A et al. Large-scale CRISPR/Cas9 deletions within the W](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC11665854/fullTextXML); [Clauss A et al. A locus on human chromosome 20 contains seve](https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=EXT_ID:12206714%20AND%20SRC:MED&resultType=core&format=json) |
| WFDC8 | specialized tissue expression | No effect | Evidence is from another species, a database annotation or an assay shared with a paralog, or the tissue is reproductive. | [Human Protein Atlas WFDC8 tissue page (antibodies HPA042710,](https://www.proteinatlas.org/ENSG00000158901-WFDC8/tissue); [UniProt Q8IUA0 (WFDC8_HUMAN) tissue specificity](https://www.uniprot.org/uniprotkb/Q8IUA0/entry); [Clauss A, Lilja H, Lundwall A. A locus on human chromosome 2](https://pubmed.ncbi.nlm.nih.gov/12206714/) |
| WFDC9 | no literature | No effect | No publication on normal adult expression was found. |  |
| WNT8B | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| XAGE1A/B | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| XAGE3 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| XAGE5 | no literature | No effect | No publication on normal adult expression was found. |  |
| XKR3 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| XKR5 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| ZACN | conflicting or nonspecific | No effect | No finding survived a second reading as human, adult, normal-tissue and gene-specific. | [PMID 16083862 (primary); PMID 33962015 for the statement tha](https://pubmed.ncbi.nlm.nih.gov/16083862/); [PMID 12381728](https://pubmed.ncbi.nlm.nih.gov/12381728/) |
| ZAR1L | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| ZC3H11B | conflicting or nonspecific | No effect | Reports conflict, cannot separate the gene from paralogs, or come from cultured, diseased or fetal samples. | [Fan et al. 2012, PLoS Genet: Genetic Variants on Chromosome ](https://pmc.ncbi.nlm.nih.gov/articles/PMC3369958/); [Fan et al. 2012, PLoS Genet (same paper, methods)](https://pmc.ncbi.nlm.nih.gov/articles/PMC3369958/); [Human Protein Atlas ZC3H11B](https://www.proteinatlas.org/ENSG00000215817-ZC3H11B) |
| ZC3H11C | no literature | No effect | No publication on normal adult expression was found. |  |
| ZCCHC13 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| ZIM3 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| ZNF280A | conflicting or nonspecific | No effect | Reports conflict, cannot separate the gene from paralogs, or come from cultured, diseased or fetal samples. | [Wang X et al. 2019, ZNF280A Promotes Proliferation and Tumor](https://pmc.ncbi.nlm.nih.gov/articles/PMC6389780/); [Liu H et al. 2021, ZNF280A promotes lung adenocarcinoma deve](https://pmc.ncbi.nlm.nih.gov/articles/PMC7791122/); [Zheng W et al. 2025, ZNF280A promotes malignant melanoma dev](https://pmc.ncbi.nlm.nih.gov/articles/PMC12008104/) |
| ZNF479 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| ZNF492 | conflicting or nonspecific | No effect | Reports conflict, cannot separate the gene from paralogs, or come from cultured, diseased or fetal samples. | [PMID 16714237, Curr Eye Res 2006 (doi 10.1080/02713680600678](https://pubmed.ncbi.nlm.nih.gov/16714237/); [Human Protein Atlas ENSG00000229676-ZNF492](https://www.proteinatlas.org/ENSG00000229676-ZNF492); [Kim YJ et al. ZNF492 and GPR149 methylation patterns as prog](https://pubmed.ncbi.nlm.nih.gov/31115548/) |
| ZNF560 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| ZNF679 | no literature | No effect | No publication on normal adult expression was found. |  |
| ZNF705A | no literature | No effect | No publication on normal adult expression was found. |  |
| ZNF705B | no literature | No effect | No publication on normal adult expression was found. |  |
| ZNF705D | no literature | No effect | No publication on normal adult expression was found. |  |
| ZNF705G | no literature | No effect | No gene-specific publication found. | [Human Protein Atlas gene page ZNF705G (ENSG00000215372)](https://www.proteinatlas.org/ENSG00000215372-ZNF705G) |
| ZNF716 | no literature | No effect | No publication on normal adult expression was found. |  |
| ZNF722 | no literature | No effect | No publication on normal adult expression was found. |  |
| ZNF723 | no literature | No effect | No publication on normal adult expression was found. |  |
| ZNF729 | no literature | No effect | No publication on normal adult expression was found. |  |
| ZNF735 | no literature | No effect | No publication on normal adult expression was found. |  |
| ZNF816-ZNF321P | no literature | No effect | No publication on normal adult expression was found. |  |
| ZNF878 | no literature | No effect | No publication on normal adult expression was found. |  |
| ZNF99 | no literature | No effect | No gene-specific publication found. | [Human Protein Atlas gene page ZNF99](https://www.proteinatlas.org/ENSG00000213973-ZNF99) |
| ZNRF4 | specialized tissue expression | No effect | Cultured human monocytes and macrophages; expression is induced. | [PMID 28656966 (full text PMC5493756)](https://pubmed.ncbi.nlm.nih.gov/28656966/) |
| ZP4 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |
| ZSWIM2 | no normal expression evidence | No effect | No report of expression or function in normal adult non-reproductive tissue. |  |

Summaries, every evidence item and its quote are in `oncoref/data/normal-expression-literature.json`.
