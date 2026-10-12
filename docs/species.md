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
