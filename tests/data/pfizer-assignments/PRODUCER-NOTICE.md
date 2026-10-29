# Pfizer dog Proteome Discoverer source slices

`native.json` contains exact selected rows/columns retrieved on 2026-10-26 from
the public MSV000099744 / PXD070290 deposit's
`search/search/IQP_1269_Dog_P1toP4.pdResult` (86,023,450,624 bytes;
MDTM 20251028233329), using bounded explicit-FTPS range reads through SQLite.
The database's full-file checksum is unavailable. The fixture records the
source identity and schema fingerprint, not a full-file integrity claim.

Tables and columns are named as deposited. The samples include the first three
target PSMs/proteins/peptide groups; three quantified and two shared peptide
groups from the first 100 row IDs; the corresponding selected distribution
metadata; and all 97 WorkflowInputFiles rows (96 raw files and one processing
MSF). Numeric BLOBs retain their original bytes in base64 with a byte length
and SHA-256. Other metadata rows are retained exactly.

These factual slices are not a full extraction or a consumer release. The
fixture's candidate examples do not establish unique protein attribution,
individual-dog expression or DLA presentation. The full searched FASTA contains
more exact sequence candidates for peptide group 8 than the vendor's reported
protein list; a complete importer must retain those distinct facts.

Attribution: Oziolor et al., Joel D. Federspiel and Pfizer. MassIVE dataset terms:
CC0-1.0. Source article: https://doi.org/10.1038/s41597-026-07512-9.
Deposit: https://massive.ucsd.edu/ProteoSAFe/dataset.jsp?task=062b109e9adf461ba619369dd0d85aa9.

The numeric value/presence storage convention is independently documented in
https://github.com/BenBruyneel/proteinDiscover/blob/main/R/transformRaw.R.
No implementation from that project is included here; native type/label
metadata and exact source bytes govern the decoding tests.

`definitions.json` retains the complete selected DataTypes (30), DataTypesColumns
(577), EnumDataTypes (20) and EnumDataTypeValues (141) metadata records, plus one
native Acetyl modification definition and its PSM link. These factual records
come from the same deposited database and bounded read workflow. They test the
distinct confidence, ambiguity, quantification and master-protein enum domains,
source modification specificity, and unaltered native position/record IDs.

`memberships.json` adds exact examples from the first three deposited
protein-group/protein and peptide-group/protein link rows, their six target
proteins, three protein groups and peptide groups 1, 2, 3 and 8. Provenance names
the original sample SQL/file hashes and the staging source checkpoints. Local
test SQLite rowids are fixture IDs, not claims about the vendor link rowids.
The fixture is not a complete set of native links. Its seven candidate peptide
placements were independently checked by exhaustive substring search across the
complete pinned 20,649-record FASTA. `membership-search-slice.fasta` contains the
ten complete source entries needed by these examples; its hash is in the JSON.
Missing endpoints, sequence mismatches, repeated I/L placements and direct
quantification links in the tests are explicitly synthetic cases.

`acquisition-missingness.json` contains three exact native quantification
spectra, their 22 reporter links and 22 reporter records. Spectrum 8546281 has
no reporter-link rows; 8546311 has eleven links despite null raw and normalized
abundance vectors. Both carry native `QuanInfo=1` (`No Quan Values`). Spectrum
8546263 retains its measured vectors. The fixture preserves original source
rowids, typed slice digests, pinned source identity and the extraction checkpoints
at capture; the reporter table was still being extracted beyond these records.
This is a bounded regression slice, not a claim of full-source completion.
Tests normalize its exact scalar/vector fields and retain the original missing
lookup key/reason separately from records containing null vectors.

`parent-status.json` retains four additional native peptide records (160079,
516803, 697816 and 993593), all 13 of their protein links, and the corresponding
13 protein records. The complete source-table checkpoints and original rowids
are retained separately from the local test rowids. Each displayed parent list
omits one linked protein with `ExcludedBy=3`; its declared peptide count is zero.
Both the displayed values and the complete native links must survive. A resolved
link is not an accepted identification. `parent-status-search-slice.fasta` holds
the 15 complete searched-reference entries for these peptides. All 15 candidate
placements were independently enumerated over the full pinned FASTA, including
the excluded native parents and candidates absent from the native parent links.

`search-slice.fasta` contains twelve complete records from the checksum-pinned
20,649-record searched FASTA. They are the candidate accessions for the eight
native peptide examples. `search-slice.json` records the full-reference audit,
every matching offset (zero-based, end-exclusive), and the slice's own hash.
The full-reference candidate set was independently compared against exhaustive
substring searches. This offline slice preserves that regression, but is not
itself the complete searched reference. Isoleucine/leucine-equivalent matches
are distinguished from exact residue matches; neither changes the original
vendor-reported memberships or supplies evidence of individual protein detection.

`psm-parent-sets.json` retains source PSMs 1890, 15032 and 16602, their three
peptide groups, four proteins and all 11 links connecting these selected records.
PSM and group parent sets remain independent, including when their native
sequence IDs agree. It preserves selected source fields, original row IDs and
completed table checkpoints. `psm-parent-sets-search-slice.fasta` contains all
four searched entries matching any of the five literal PSM/group sequences;
all ten placements were independently enumerated against the complete 20,649
entry deposited FASTA. Fixture row IDs are local to the regression database.

The FASTA entries in these slices are from the UniProt Consortium, distributed
under CC BY 4.0 (https://creativecommons.org/licenses/by/4.0/), with attribution
to UniProt (https://www.uniprot.org/help/license). They retain the original
headers and sequences; the change here is selection of entries for tests.
Their source is the deposited `Dog_UP000002254_9615_noisoform.fasta`, SHA-256
`d8488bf283de653efcff2e5de6d413b116e1ba4a2eb89298bd31205e72fef682`.
The deposit's CC0 terms do not replace the searched FASTA's UniProt license.

`study-slice.xml` retains the complete FileSets, QuanMethods and Samples elements
from the deposited 20210614_zoomapdog_nodbsearch.pdStudy. This includes both cyno
and dog, to check that the dog reconciliation selects only its four file sets.
The full source is 3,528,792 bytes, SHA-256
`fd0a494cd0a73d120045e092668f7ce02cc23ec57722cc48d309e57ca6394bb8`.
The XML slice is 60,419 bytes, SHA-256
`ad15f7d631559d48d85c4e3d41302c363ba38cc8fbebe50a6ec7a647e34ab2a4`.
`summary-samples.json` contains all 38 distinct seven-field sample records from
the already pinned protein summary; it contains no abundance rows.

`dog-panel.json` records literal dog channel/tissue facts and cell locations
from Supplemental Table S2 (TableS2_tissues_proteomics!M1:Q13). The original
316,225-byte workbook has SHA-256
`a95b7a2514679ba6f7283573059f27e6182b5620f6d0beead601e47daf85e740` and MD5
`e354f58534b5c41c7819285605a2ad97` (the latter independently supplied by the
article's Europe PMC full-text XML). All 61 retained cells were independently
checked with openpyxl 3.1.5 against the untouched original workbook.

Supplement source: Oziolor et al., DOI 10.1038/s41597-026-07512-9,
41597_2026_7512_MOESM15_ESM.xlsx, retrieved via
https://www.ebi.ac.uk/europepmc/webservices/rest/PMC13482810/supplementaryFiles?includeInlineImage=false.
The publication/supplement license is CC BY-NC-ND 4.0
(https://creativecommons.org/licenses/by-nc-nd/4.0/), distinct from the CC0
deposited study/data inputs above. The complete workbook and its formatting are
not redistributed here. Tests use these factual cells in a synthetic minimal
OOXML wrapper to check reading without recalculation or modification.
