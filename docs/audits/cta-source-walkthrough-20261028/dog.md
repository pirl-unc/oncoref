# Dog CTA evidence funnel

**14 native nominations, all unassigned/unassessed.** There is no completed canine quantitative selection or genome-wide denominator. None is automatically called Neither. Human RNA ceilings, HPA grades and clinical evidence do not transfer by orthology.

| citation | mapped_entities | mapped_gene_loci | source_url |
| --- | --- | --- | --- |
| [Hindriks et al. (2025)](https://doi.org/10.1016/j.vetimm.2025.110946) | 10 | 9 | [Source](https://doi.org/10.1016/j.vetimm.2025.110946) |
| [Nemec et al. (2019)](https://doi.org/10.1111/vco.12475) | 5 | 5 | [Source](https://doi.org/10.1111/vco.12475) |

Nemec et al. contributes 5 nominations; Hindriks et al. contributes 10; NT5C1B is shared, giving 14. Thirteen have unique gene mappings; one MAGE family assay remains ambiguous. A unique gene mapping does not establish identical protein isoforms.

![Dog funnel](figures/dog-funnel.png)

| stage | status | note |
| --- | --- | --- |
| Curated native nominations | available | 14 native nominations; not a genome-wide screen |
| Reference and assay identity | partial | 13 unique canine mappings and one ambiguous MAGE-family assay |
| Normal and reproductive/developmental RNA | not_evaluated | Bgee TPM and BarkBase TMM-CPM require native reference reconciliation |
| Tumor RNA | not_evaluated | GSE304066 count/FPKM and source-reference conflicts unresolved |
| Normal protein | not_evaluated | Study IHC is limited; Pfizer pooled iBAQ is not cell-resolved IHC |
| Quantitative selection | not_evaluated | No native species selection policy; no biological category assignment |

## Sources and gaps

| name | evidence_category | native_unit | status | source_url |
| --- | --- | --- | --- | --- |
| [bgee-srp009687](https://www.ncbi.nlm.nih.gov/sra/?term=SRP009687) | normal_bulk_rna | TPM | native_matrix | [Source](https://www.ncbi.nlm.nih.gov/sra/?term=SRP009687) |
| [bgee-srp009687](https://www.ncbi.nlm.nih.gov/sra/?term=SRP009687) | reproductive_expression | TPM | native_matrix | [Source](https://www.ncbi.nlm.nih.gov/sra/?term=SRP009687) |
| [bgee-srp181218](https://www.ncbi.nlm.nih.gov/sra/?term=SRP181218) | normal_bulk_rna | TPM | native_matrix | [Source](https://www.ncbi.nlm.nih.gov/sra/?term=SRP181218) |
| [bgee-srp181218](https://www.ncbi.nlm.nih.gov/sra/?term=SRP181218) | reproductive_expression | TPM | native_matrix | [Source](https://www.ncbi.nlm.nih.gov/sra/?term=SRP181218) |
| [barkbase](https://pubmed.ncbi.nlm.nih.gov/31181663/) | normal_bulk_rna | TMM-CPM | native_matrix | [Source](https://pubmed.ncbi.nlm.nih.gov/31181663/) |
| [barkbase](https://pubmed.ncbi.nlm.nih.gov/31181663/) | reproductive_expression | TMM-CPM | native_matrix | [Source](https://pubmed.ncbi.nlm.nih.gov/31181663/) |
| [barkbase](https://pubmed.ncbi.nlm.nih.gov/31181663/) | developmental_expression | TMM-CPM | native_matrix | [Source](https://pubmed.ncbi.nlm.nih.gov/31181663/) |
| [pfizer-dog-protein](https://doi.org/10.1038/s41597-026-07512-9) | normal_protein | iBAQ | native_matrix | [Source](https://doi.org/10.1038/s41597-026-07512-9) |
| [doga](https://doi.org/10.1038/s41467-024-52798-1) | normal_bulk_rna | unresolved | metadata_only | [Source](https://doi.org/10.1038/s41467-024-52798-1) |
| [doga](https://doi.org/10.1038/s41467-024-52798-1) | reproductive_expression | unresolved | metadata_only | [Source](https://doi.org/10.1038/s41467-024-52798-1) |
| [doga](https://doi.org/10.1038/s41467-024-52798-1) | developmental_expression | unresolved | metadata_only | [Source](https://doi.org/10.1038/s41467-024-52798-1) |
| [pfizer-dog-rna](https://doi.org/10.1038/s41597-026-07512-9) | normal_bulk_rna | unresolved | metadata_only | [Source](https://doi.org/10.1038/s41597-026-07512-9) |
| [gse304066-osteosarcoma](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE304066) | tumor_expression | unresolved | blocked_provenance | [Source](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE304066) |
| [hindriks2025](https://research-portal.uu.nl/ws/portalfiles/portal/264124617/1-s2.0-S0165242725000662-main.pdf) | nomination | source-native observations; no absolute expression unit | native_observations | [Source](https://research-portal.uu.nl/ws/portalfiles/portal/264124617/1-s2.0-S0165242725000662-main.pdf) |
| [hindriks2025](https://research-portal.uu.nl/ws/portalfiles/portal/264124617/1-s2.0-S0165242725000662-main.pdf) | normal_bulk_rna | source-native observations; no absolute expression unit | native_observations | [Source](https://research-portal.uu.nl/ws/portalfiles/portal/264124617/1-s2.0-S0165242725000662-main.pdf) |
| [hindriks2025](https://research-portal.uu.nl/ws/portalfiles/portal/264124617/1-s2.0-S0165242725000662-main.pdf) | normal_protein | source-native observations; no absolute expression unit | native_observations | [Source](https://research-portal.uu.nl/ws/portalfiles/portal/264124617/1-s2.0-S0165242725000662-main.pdf) |
| [hindriks2025](https://research-portal.uu.nl/ws/portalfiles/portal/264124617/1-s2.0-S0165242725000662-main.pdf) | reproductive_expression | source-native observations; no absolute expression unit | native_observations | [Source](https://research-portal.uu.nl/ws/portalfiles/portal/264124617/1-s2.0-S0165242725000662-main.pdf) |
| [hindriks2025](https://research-portal.uu.nl/ws/portalfiles/portal/264124617/1-s2.0-S0165242725000662-main.pdf) | tumor_expression | source-native observations; no absolute expression unit | native_observations | [Source](https://research-portal.uu.nl/ws/portalfiles/portal/264124617/1-s2.0-S0165242725000662-main.pdf) |
| [nemec2019](https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=EXT_ID:30854786%20AND%20SRC:MED&resultType=core&format=json) | nomination | source-native observations; no absolute expression unit | native_observations | [Source](https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=EXT_ID:30854786%20AND%20SRC:MED&resultType=core&format=json) |
| [nemec2019](https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=EXT_ID:30854786%20AND%20SRC:MED&resultType=core&format=json) | normal_bulk_rna | source-native observations; no absolute expression unit | native_observations | [Source](https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=EXT_ID:30854786%20AND%20SRC:MED&resultType=core&format=json) |
| [nemec2019](https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=EXT_ID:30854786%20AND%20SRC:MED&resultType=core&format=json) | tumor_expression | source-native observations; no absolute expression unit | native_observations | [Source](https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=EXT_ID:30854786%20AND%20SRC:MED&resultType=core&format=json) |

Bgee and BarkBase include 46 overlapping libraries, not independent replication. BarkBase has 37,149 features × 150 TMM-CPM columns; Pfizer MS has 9,929 reported groups × 38 tissue pools and missing matrix cells. Pooled iBAQ is not cell-specific IHC or antigen presentation. Tumor GSE304066 count/FPKM and reference provenance remain unresolved. These quantities must retain their native scales.

[RNA/reference/lineage integration issue #574](https://github.com/pirl-unc/oncoref/issues/574) and [Pfizer MS integration issue #585](https://github.com/pirl-unc/oncoref/issues/585) track the main work. Systematic canine IHC and normal-cell atlases remain coverage gaps. [Pairwise nomination overlaps](dog-overlaps.csv.gz).

The same API exposes sources, memberships, overlaps, criteria, funnel and explanations for human, dog and mouse; unavailable evaluations have null counts and null biological categories.
