# Verified CTA sequence records

The current sequence bundle binds source proteins and every occurrence to the
selection, translation and consumer-role policies. It preserves original sequence
bytes, ambiguous residues, incomplete records and explicit coverage gaps. The
[bundle audit](audits/sequence-bundle-20261104/index.md) gives its exact identity,
counts and comparison with the historical v1.20.0 forbidden assets.

Acquisition is explicit and bounded by the packaged asset hash and size. Reading
the descriptor or loading an existing bundle performs no network access.

```python
from oncoref import sequence_bundle

pin = sequence_bundle.descriptor()
path = sequence_bundle.ensure(download=True)  # explicit one-asset acquisition
bundle = sequence_bundle.load(path)          # verifies current package dependencies
records = bundle.records()                  # original bytes and all occurrences
segments = bundle.segments(role="exclude", min_length=8)
coverage = bundle.coverage
```

`ONCOREF_SEQUENCE_CACHE_DIR` or `ensure(cache_dir=...)` selects an isolated cache.
A missing or corrupt cache raises an error; an explicit `download=True` permits
repair using the pinned asset. Mixed policy/reference bundles fail validation.
For historical use, supply a separately trusted `expected_id` and
`require_current=False`; an identity read from an untrusted bundle is not a pin.

Each segment has zero-based, half-open coordinates in its original sequence.
Segments never bridge an `X`, internal stop, selenocysteine or separate record.
The current AA20 contract requires review for unsupported residues; that is not
a claim of absent biological production. Original sequences remain available.

`groups()` groups exact original sequence strings and retains all records and
occurrences. Their combined role follows **exclude > review > allow**. A consumer
matching peptides must likewise aggregate every matching segment occurrence;
an allowed CTA occurrence never cancels another occurrence's exclusion or review.
Exact sequence grouping is not a validated approximate similarity rule.

An Unsafe gene conservatively excludes its source isoforms, while positive
allowance requires the qualified selected representative. Other isoforms require
review. Gene-level evidence does not demonstrate isoform-specific production.
Missing translations, source genes without decisions, the empty legacy outside-inventory list and
ncORF nominations without supplied peptide bytes remain explicit in `coverage`.

The optional `oncoref[quantitative]` extra supports `bundle.as_orfdata()`. It loads
the preserved native ORFdata proteome and verifies its original content identity.
The returned raw proteome does not apply CTA policy roles; retain the bundle when
performing exclusion/review decisions. The current source is the full Ensembl 115
peptide inventory, not a complete empirical ncORF sequence database.

Historical Unsafe plus additional-forbidden FASTAs retain their original policy
and reference identities. The downloadable transition table compares exact full
sequence hashes; an old sequence absent from the new source requires coverage
review rather than implicit allowance. Neither product establishes patient-HLA
presentation, receptor recognition or clinical safety.
