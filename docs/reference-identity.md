# OpenVax reference and protein identity contract, version 1

These primitives live in OncoRef's reference layer and use the Python standard
library. Their contract is shared across consumers: no species, application,
assay, target panel, tissue threshold or RNA allocation policy is implicit.

## Three distinct identifiers

| Identifier | Meaning | Form |
| --- | --- | --- |
| `reference_id` | Exact declared taxon, assembly, annotation, namespaces and source assets | `openvax:reference:v1:sha256:<hex>` |
| `protein_sequence_id` | Normalized exact amino-acid sequence across references | `openvax:protein-sequence:v1:sha256:<hex>` |
| `group_id` | That sequence within that reference | `openvax:protein-group:v1:sha256:<hex>` |

These are opaque, case-sensitive strings. Persist the complete identifier, including
kind, version and algorithm. Identical proteins from different references share a
sequence ID, but have different group IDs. Different reference IDs can describe
biologically equivalent material; equivalence requires separate evidence.

## Reference schema

`ReferenceIdentity.as_dict()` serializes these exact fields:

- `schema`: `"openvax.reference"`; `schema_version`: integer `1`.
- `taxon_id`: a positive integer NCBI taxonomy ID; booleans/strings are rejected.
- `assembly_accession`: the exact provider accession, including its version when
  available. Custom assemblies require an explicit stable identifier and asset hashes.
- `annotation_source`: explicit naming authority, e.g. `"NCBI RefSeq"` or `"Ensembl"`.
- `annotation_release`: the provider's annotation release identifier, as text.
- `source_version`: the upstream reference snapshot identifier. This is not the
  version of OncoRef, a consumer, a local run or a normalization pipeline.
- `identifier_namespaces`: entity kinds mapped to their identifier authorities.
  Protein grouping requires `gene`, `transcript` and `protein` entries; other
  reference consumers may declare additional or different kinds.
- `content_provenance`: stable logical asset roles mapped to lowercase SHA-256 of
  exact source bytes. Hexadecimal input is case-insensitive. Role names use
  `[A-Za-z0-9][A-Za-z0-9_.-]*`; local paths and retrieval URLs belong in external
  provenance, not the identity. Fix role names when registering a reference.
- `reference_id`: the derived identifier.

All fields except the derived ID participate in identity. The complete declared
asset inventory is significant; adding an asset or renaming a role changes the ID.
Include reference-defining assets, not unrelated reports or pipeline outputs.
Compression differences change byte digests. No sequence-equivalence or annotation
crosswalk is inferred from matching names, symbols or orthology. Hash validation
checks syntax; acquisition adapters must verify the actual file bytes.

The object copies and freezes its mappings. `from_dict()` requires the exact v1
field set, verifies the supplied ID, and rejects unsupported versions. Extra
provenance belongs beside the identity record. Labels retain case and Unicode
spelling without normalization; empty labels, edge whitespace and ASCII controls
are rejected. There are no implicit namespace aliases.

## Normative digest encoding

For the reference digest, remove only `reference_id` from the serialized record.
Encode the remaining object as ASCII JSON with keys sorted, no insignificant
whitespace, ASCII escaping, and no NaN/Infinity. The Python spelling is:

```python
json.dumps(payload, sort_keys=True, separators=(",", ":"),
           ensure_ascii=True, allow_nan=False).encode("ascii")
```

SHA-256 these bytes and append the lowercase hexadecimal digest to the typed ID
prefix. Reference and group digest payloads contain strings, dictionaries and integers;
arbitrary floating-point provenance is outside those payloads.

Protein v1 normalization first validates the input against ASCII upper/lowercase
AA20 letters, the six ASCII whitespace characters (`space`, `tab`, `LF`, `CR`,
`VT`, `FF`), and `*`. Remove those whitespace characters, remove at most one final
`*`, then uppercase. Require a nonempty sequence of `ACDEFGHIKLMNPQRSTVWY`.
Hash those ASCII residues directly. I/L remain different. U/O and ambiguity codes
are explicitly unsupported by this version. Non-ASCII characters are rejected
before conversion: `Mß`, `Mı`, ligatures and Unicode whitespace cannot become valid
residue strings. Completeness is a separate source assertion (`complete=True`).

The group digest uses the same JSON encoding with exactly four keys:
`schema="openvax.protein-group"`, `schema_version=1`, `reference_id`,
`protein_sequence_id`. Occurrence order and provenance do not change a group ID.
Occurrences are sorted by `occurrence_id`; groups are sorted by sequence ID.

[Interchange vectors](https://github.com/pirl-unc/oncoref/blob/main/tests/fixtures/reference-identity-v1.json)
freeze human, dog and mouse examples with synthetic asset digests. Consumers can
use them without real genomes or an OncoRef data download. Schema/normalization
changes require a new identity version; package releases do not change these IDs.

## Explicit legacy import

The unreleased first draft of #573 used untyped digests and a smaller reference
schema. Those keys, and consumer-specific existing keys, are not v1 aliases.
Keep their original records. The generic adapter requires all of:

1. A fully declared target `ReferenceIdentity`.
2. The original source-reference record and its opaque source ID. The caller
   verifies that ID using the source producer's algorithm; it is never guessed.
3. A field map from target reference fields to source fields, and explicit values
   for missing declarations. Every target field must be supplied once. Taxon,
   assembly, annotation release and content hashes must be mapped from the source.
4. Checksum-pinned migration evidence, such as the reviewed source manifest.
5. The name of the source reference-ID field on occurrences. Every row must match.

For example, a consumer whose schema uses `species`, `assembly`, `release`,
`upstream_version`, `files`, and occurrence field `old_reference` can use:

```python
from oncoref import adapt_reference_occurrences, protein_sequence_groups

migration = adapt_reference_occurrences(
    source_occurrences,
    reference=reference,
    source_reference=source_manifest,
    source_reference_id=verified_legacy_id,
    reference_id_field="old_reference",
    field_map={
        "taxon_id": "species", "assembly_accession": "assembly",
        "annotation_release": "release", "source_version": "upstream_version",
        "content_provenance": "files",
    },
    added_fields={
        "annotation_source": reference.annotation_source,
        "identifier_namespaces": dict(reference.identifier_namespaces),
    },
    evidence={"reviewed_source_manifest": verified_manifest_sha256},
)
groups = protein_sequence_groups(migration["occurrences"], reference=reference)
```

The mapped source declaration must reproduce the target identity exactly. Conflicting
per-row declarations are rejected. The adapter performs a schema migration only;
it cannot convert species, lift coordinates, translate gene IDs, or reconcile
reference assemblies. Each output preserves its original `source_record` and adds
`reference_id` and `reference_migration_id`. Existing fields are retained. Reserved
migration-field collisions fail rather than discarding provenance.

Persist the returned `openvax.reference-migration` v1 audit envelope alongside
its grouped occurrences. Its typed migration ID hashes all envelope fields except
`migration_id` and `occurrences` using the same JSON encoding. The envelope records
both references, the original ID field, field mappings, added declarations and
evidence hashes. Consumers own their source adapters and biological validation;
OncoRef's schema does not import or depend on those consumers.
