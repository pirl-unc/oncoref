# Species and data groups

Species is the first data boundary. Select it explicitly, then inspect conceptual
collections within it. All species expose the same groups: `references`,
`expression`, `cancer`, `antigens`, `normalization`, `governance`, and `legacy`.
A known species does not imply that every group has data or complete coverage.

```python
import oncoref

dog = oncoref.for_species("dog")  # also 9615, "canine", or scientific name
oncoref.available_species()
dog.references.datasets()
dog.expression.sources(context="normal")
dog.antigens.coverage()  # unavailable until canine antigen facts are registered
human = oncoref.for_species("human")
human.expression.sources()  # existing human registry
```

```bash
oncoref species
oncoref species dog
oncoref data list --species dog --group expression
oncoref expression-sources --species dog --context normal
```

The existing flat APIs, cache paths and data bundle retain their historical human
behavior. Their default is not a claim of species independence. The new species
collections never substitute HPA, human cancer facts or human CTA membership for
missing dog/mouse facts. Unknown species fail explicitly. `catalog.inventory()`
defaults to human and adds `species`, `taxon_id`, and `group` to each row;
`catalog.inventory(species="dog")` selects only dog facts. Availability describes
registered/local data, not complete normal-tissue or antigen coverage.

## Expression contracts

`ExpressionSource.expression_contract` returns versioned, immutable metadata:

- Species, source ID/version, normal/tumor/cell-line/mixed/unknown biological scope.
- The shared `ReferenceIdentity` (taxon, assembly, annotation, namespaces, hashes),
  or an explicit absent reference with `unresolved`/`legacy_human` status.
- Measured feature level (`gene`, `transcript`, `promoter`, `protein`,
  `protein_group`, or `unknown`), native unit, scale and proxy status.
- Source asset SHA-256s, provenance notes and ordered transformations, including
  their unit changes, assumptions and input/output feature counts.

A verified reference is an explicit source declaration supported by provenance;
syntax/digest checks alone cannot establish the correct biological reference.
Historical human rows use `legacy_human`; missing historical reference/feature
information is not reconstructed by guessing. They remain readable by existing
APIs. Nonhuman sources must declare their contract explicitly.

`ExpressionContract.as_dict()` / `from_dict()` provide strict JSON round trips,
version/digest verification and unknown-field rejection. The `contract_id`
fingerprints all metadata, so new derived artifacts must include it in their cache
identity. Existing human bundle fingerprints and replay are unchanged.

`validate_expression_join(left, right)` requires verified, matching species,
references, namespaces and measured feature levels. A directed,
checksum-supported `ReferenceCrosswalk` can authorize an annotation join between
different references within one species. It preserves one-to-many alternatives
and unresolved IDs. It does not authorize abundance allocation. With
`combine_values=True`, matching references and complete identical measurement
semantics are mandatory; joins do not infer missing values as zero.

`contract.require_feature_level("transcript")` rejects gene or promoter totals.
Exact protein identity also cannot establish individual isoform abundance.
Native counts, CPM, ERCC-normalized counts, FPKM, TPM, pseudo-TPM and iBAQ retain
their unit labels; transformations never turn a proxy into a measured quantity.

## Staged canine support

Version 1.9.0 supplies species collections, source/measurement contracts and join
validation. Canine source manifests, offline hash verification and native-unit
readers follow in separate versioned changes. No canine CTA admission or target
selection follows from registering references or expression sources.

## Pinned canine references (1.9.1)

`dog.references.registered_references()` exposes RefSeq **GCF_011100685.1 /
annotation 106** and Ensembl **CanFam3.1 / release 98 (GCA_000002285.2)**. The
former pins genome, GTF, CDS, RNA, complete annotated protein FASTA and feature
table; the latter pins GTF and peptide FASTA. Different assemblies/annotations
remain different reference identities. Registration does not claim a complete
validated translation subset or establish an expression-reference crosswalk.

```python
reference, paths = dog.references.ensure_reference(
    "refseq106", offline=True,
    local_overrides={"9615:refseq106-genome": "/existing/genome.fna.gz",
                     # supply each reference asset or reuse its managed cache
                    },
)
# One asset, online if absent; an explicit transfer budget is required for large files:
path = dog.references.ensure("ensembl98-gtf", max_bytes=25_000_000)
```

Local overrides use **`taxon_id:asset_id`** keys and are verified in place, without
copying or overwriting the original. `ONCOREF_ASSET_OVERRIDES` can point to a JSON
mapping of those keys to paths. `ONCOREF_ASSET_CACHE` selects the datacache root;
the default uses datacache's `oncoref` application directory, partitioned by taxon,
group and content digest. Existing human cache locations are unchanged.

```bash
oncoref data metadata refseq106-gtf --species dog
oncoref data status references --species dog --asset-overrides /path/to/overrides.json
oncoref data fetch ensembl98-gtf --species dog --offline --asset-overrides /path/to/overrides.json
```

`status`, `path`, and offline replay verify the pinned size and SHA-256. Inventory
listing only checks local presence and reports `integrity: not_checked`. Online
fetches use datacache's atomic publication and bounded HTTP resume, with a default
100 MiB total transfer budget (`--max-bytes` to change it). Reference fetches check
the total missing bytes before starting. An invalid override fails even online;
managed-cache corruption can be repaired online but fails offline. The wheel
contains small manifests only, with upstream URLs, byte pins and source terms.

## Native canine expression (1.9.2)

```python
dataset = dog.expression.read(
    "bgee-srp181218", offline=True,
    local_overrides={"9615:bgee-srp181218-matrix": "/existing/bgee-srp181218.tsv.gz"},
)
dataset.values       # source features × libraries; missing remains NaN
dataset.features     # original IDs/aliases; no human canonicalization
dataset.samples      # tissues, donor scope, pooling, runs, original facts/evidence
dataset.contract     # gene TPM; quantification reference explicitly unresolved
dataset.tissue_summary()  # native values + measured-library/known-donor coverage
path = dataset.save("new-artifact-directory")
from oncoref.species_expression import ExpressionDataset
replayed = ExpressionDataset.load(path)  # offline file/digest/contract verification
```

| Source ID | Native feature space | Unit | Source libraries |
| --- | --- | --- | --- |
| `bgee-srp009687` | Ensembl source gene IDs | TPM | 13 |
| `bgee-srp181218` | Ensembl source gene IDs | TPM | 75 |
| `barkbase` | All 37,149 StringTie gene rows, retaining Ensembl/symbol annotations | TMM-CPM | 150 |

Bgee's sparse gene/library rows remain missing where absent. BarkBase retains all
source features, including unannotated StringTie IDs and ambiguous annotation
strings; no length proxy or feature filter is applied. Exact quantification
references remain unresolved for these sources, so cross-reference joins and
cross-source abundance combination remain blocked by their contracts.

Sample metadata preserves source fields, BioSample/run provenance and unknown
age/sex/health/treatment/pooling facts. Five Bgee SRP181218 colon labels are corrected
to heart using BioSample `source_name`, with both labels, a correction reason and
source hashes retained. Donors are scoped to their study. BarkBase adult donor
labels follow the source columns and the publication's five-adult design;
embryonic donor identities and matrix-column-to-run mappings remain unresolved.
Pooled observations never contribute multiple independent donors. A zero count
of **known** donors means unresolved coverage, not proof that no animals exist.

`dog.expression.lineage_overlaps()` exposes the 46 shared runs between the global
Bgee inventory and BarkBase's RNA-seq inventory. This is inventory-level lineage;
it does not assert overlap of the two imported Bgee matrices or resolve BarkBase
column-to-run mappings. `combine_expression()` rejects known reused runs,
unresolved run lineage, unresolved reference identities and incompatible units.

Saved artifacts include all feature/sample facts, the expression contract and
checksummed Parquet/JSON files. New destinations are published atomically; existing
results are not overwritten. The artifact identity changes with its metadata or
file contents. `scripts/build_canine_source_manifests.py` replays the intake from
an acquired source snapshot; no downstream package is imported at runtime.
