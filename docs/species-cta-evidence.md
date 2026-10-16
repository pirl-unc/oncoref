# Species CTA candidates and evidence

Use the same methods on `oncoref.human`, `oncoref.canine`, and `oncoref.mouse`:

```python
from oncoref import human, canine, mouse

human.cta_candidates()       # existing human curation, wrapped without losing fields
canine.cta_candidates()      # canine study nominations and assay follow-ups
mouse.cta_candidates()       # () -- no curated mouse evidence yet

canine.cta_evidence()        # all native observations, sources and mapping details
canine.cta_evidence(candidate_id="CTCFL")
mouse.cta_snapshot().coverage  # {"status": "unavailable", "limitations": [...]}
```

Each query returns a tuple of detached JSON records. Candidate and evidence
records have the same fields across species. The same methods are also available
under each species' existing `antigens` collection. There are no membership or
ambiguity switches on these queries: all nominations and counterevidence remain
visible. These are evidence APIs, not patient admission or vaccine selection.

For an annotation-bound set of **candidate** gene IDs, choose the reference:

```python
canine.cta_gene_ids(reference="refseq106")
human.cta_gene_ids(reference="ensembl115")
```

This set contains only complete, unique **gene** mappings. Ambiguous alternatives
remain in `cta_candidates()` under `mapping.targets` and `gene_ids`; no preferred
gene is chosen. Assay follow-ups can have unique gene mappings even when their
CTA interpretation remains unresolved. Neither mapping uniqueness nor inclusion
in this set establishes normal-tissue restriction.

`cta_candidates()` and `cta_evidence()` use the snapshot's sole native reference
for display; the returned records include its full identity. Their optional
`reference=` argument also accepts a `ReferenceIdentity`. The ID query requires
an explicit reference. Unknown annotations, changed content identities and
wrong-species references raise `ValueError`. Empty mouse evidence does not permit
a reference-bound ID query. Inspect coverage before querying its IDs.

The original flat `oncoref.cta_gene_ids()` keeps its historical human default
panel and options. The new `human.cta_gene_ids(reference=...)` enumerates the
broader candidate set; it does not apply the historical default-panel decision.
The original `oncoref.cta_evidence()` DataFrame is unchanged.

## One evidence flow, different sources

Each species uses the same flow: source records become common facts, facts remain
linked to their original observations, and candidates report their evidence and
reference coverage. Source adapters handle the historical human HPA summary and
native study observations. Validation, queries, coverage reporting and archive
replay share one implementation. A new species with native observations uses the
same implementation without a new consumer policy.

```python
for species in (human, canine, mouse):
    for record in species.cta_evidence():
        for fact in record["facts"]:
            print(record["candidate_ids"], fact["tissue_role"], fact["tissue"])
            print(fact["measurement"], fact["observation"])

candidate = canine.cta_candidates()[0]
candidate["evidence_coverage"]["missing"]
```

`facts` have the same fields for every source. Each fact identifies its parent
evidence record and, for a summary, its original `source_field`. The parent
retains source/version/checksum, study lineage, native assay, all candidate IDs,
sample/donor counts and technical replicates. A family or multi-gene observation
remains one observation linked to all of its candidates.

`statistic` distinguishes `native_observation`, `tissue_summary`,
`tissue_maximum`, `somatic_maximum` and `tissue_detection_summary`. Human RNA
values retain nTPM; canine PCR retains percent of testis-normalized expression;
IHC remains categorical or a count of positive/examined samples. The common
representation does **not** make these values numerically interchangeable.
Human organ maxima can overlap the overall somatic maximum. The brain group
includes retina and spinal cord; it is not a measurement of one brain tissue.

Human IHC tissue lists supply qualitative detection summaries with protein-group
scope. Missing IHC and false derived protein flags never become measured negative
observations. Historical RNA summary zeros are retained as summary values: the
historical builder could fill absent tissues with zero, so these values do not
establish measured absence. Native censored bounds and no-amplification results
remain distinguishable from both zero and unreported measurements.

`curation` holds explicitly named, versioned historical policy records separately
from observations. The human `oncoref-human-cta-core` policy preserves its original
decision, decision field, rule and source rationale. Its included genes reproduce
the existing flat `oncoref.cta_gene_ids()` exactly. Shared consumers can inspect
the same fact types for any species while owning their own thresholds and
overrides; historical membership does not override normal-tissue counterevidence.

`evidence_coverage` groups fact IDs by tissue role and reports specific missing
capabilities: observations, a complete tissue panel, independent sample/donor
totals, quantification-reference identity, sequence reference assets, or
transcript/protein pairs. No missing tissue list or independent total is inferred
from summary columns, assay counts or technical replicates. `unreported_tissues`
is `null` because no complete tissue-panel denominator is declared. Per-observation
counts remain available on the evidence record.

A unique gene mapping is still separate from measured tissue restriction and
complete sequence coverage. Human mappings currently lack a full sequence-asset
identity and exhaustive transcript/protein mappings; those gaps are explicit.
Canine mappings retain their pinned reference assets and sequence pairs. A
measurement verified against another reference reports
`matching_quantification_reference` as missing. Matching assembly labels or
symbols never establish reference equivalence. The human sequence identity bridge
and initial mouse curation remain follow-up work in #599 and #600.

## Reading a record

- `nomination_kind` records provenance: `published_nomination`, `assay_followup`,
  `expression_derived`, or `legacy_human`. No expression-derived canine panel is
  supplied in this initial snapshot. Orthology is recorded as a nomination
  derivation, never as evidence of canine tissue restriction.
- `evidence_status` describes limited, conflicting or unresolved evidence;
  `derivation` and source locators explain why the candidate was included.
- `mapping` contains the taxon/annotation-scoped gene namespace, all supported
  target IDs, associated transcript/protein accession pairs, completeness and
  mapping method. These sequence pairs describe the annotation, not measured
  isoform abundance. Family assays retain family scope.
- Evidence includes the original assay, tissue and biological context, sample
  identifiers, sample/donor counts, technical replicates, study lineage, source
  version/checksum and an existing `ExpressionContract`.
- `observation.status` distinguishes reported values, upper bounds, no
  amplification, qualitative findings and unreported observations. An upper bound
  is not an exact value. Missing values/counts are `null`, never filled with zero.
  Unlisted tissues have no evidence, not evidence of absence.
- `source_fields` preserves a heterogeneous source record when needed. The human
  adapter retains every field of the original CTA summary here, including native
  RNA nTPM values, categorical protein calls, source labels and curation decisions.
  Its expression contract remains `legacy_human`; donor counts and upstream
  quantification-reference equivalence are not reconstructed.

The human `ensembl115` reference identifies the bundled primary gene-space
projection of Ensembl 115 / GRCh38.p14 (GCA_000001405.29), hashed as logical JSON
records. It is distinct from a full GTF/sequence-backed reference identity. The
canine `refseq106` mapping reuses the existing registered reference identity,
including its source asset checksums. Historical canine assay references remain
`unresolved`; mapping a target to annotation 106 does not retroactively establish
the assay's annotation or authorize an expression join.

## Frozen offline replay

```python
snapshot = canine.cta_snapshot()
snapshot.write("canine-cta.json")
receipt = snapshot.snapshot_id

# Later, with no network or reference assets, possibly a different package version:
replay = canine.cta_snapshot(path="canine-cta.json", expected_snapshot_id=receipt)
assert replay.snapshot_id == receipt
replay.candidates()
replay.evidence(candidate_id="CTCFL")
replay.gene_ids(reference="refseq106")
```

Human and mouse snapshots use the same export/replay contract. The complete
snapshot is immutable and content-addressed; exported metadata, evidence and
mapping changes invalidate its SHA-256 identity. A caller-saved expected identity
pins the selected version. An internal checksum detects changes, not authorship.
No source PDFs, expression matrices, GTFs or network access are needed to replay.
New default snapshots use `oncoref.cta-evidence`, version 2, with frozen `facts`
and `curation` sections. Unknown fields and inconsistent species, references,
contracts, links, source values or counts are rejected. Replay uses the archived
facts, without running adapters or reading current data.

Version 1 archives remain readable with their original identity and contents.
Their query results report missing normalized observations rather than silently
reinterpreting the archive. To explicitly derive a new v2 archive using only an
old archive's records:

```python
from oncoref.cta_evidence_adapters import normalize_cta_snapshot

upgraded = normalize_cta_snapshot(replay)
upgraded.write("cta-v2.json")
```

The canine builder below continues to reproduce the pinned v1 source snapshot.
Default package queries adapt it to v2; explicit file replay preserves whichever
version was saved.

## Initial canine curation

Version `2026-10-13.1` includes 12 published candidates plus MAGEA13P and MAGEB10
assay follow-ups. It contains 209 observations and reports partial coverage.

[Nemec et al. 2019](https://pubmed.ncbi.nlm.nih.gov/30854786/) contributes five
nominations and their abstract-level counterevidence. Only the indexed abstract
was imported: quantitative units, sample counts and DLA peptide sequences remain
unavailable. The AKR1E2/SPECC1/TPX2 normal-tissue finding is retained as one
multi-gene statement; it cannot establish which individual gene was elevated in
each named tissue. This limitation is in both the source receipt and evidence.

[Hindriks et al. 2025](https://doi.org/10.1016/j.vetimm.2025.110946) contributes
Table 1 target identities, Table 2 RT-PCR observations, Table 3 tumor IHC counts,
and normal/reproductive tissue findings from sections 3.1–3.4. PCR values retain
percent of testis-normalized expression, reference-gene normalization and
technical triplicates. Kidney and placental expression are retained. IHC remains
a protein-group assay. MAGE family mapping lists the two reported subtype
follow-ups as incomplete alternatives; it does not claim to enumerate or resolve
all primer/antibody targets. Supplementary material is not imported.

The paper's PCR threshold descriptions disagree between section 3.1 and Table 2;
some table signs also disagree with their printed numeric values. Both signs and
values are preserved, with the conflict documented. No positivity cutoff is
recomputed. Repeated assays and observations share their study lineage, and
unknown donor overlap remains unknown. Reproductive tissue labels do not imply
an exclusion, safety waiver or transfer of human CTA rules.

### Rebuilding the canine snapshot

The deterministic builder requires the exact source bytes identified in the
snapshot. It validates all three inputs before generating any output and works
with a fresh destination directory:

```bash
python scripts/build_canine_cta_evidence.py \
  --feature-table /data/GCF_011100685.1_UU_Cfam_GSD_1.0_feature_table.txt.gz \
  --hindriks-pdf /data/hindriks2025.pdf \
  --nemec-record /data/nemec2019.json \
  --output /tmp/fresh-cta/cta-evidence.json
```

The PDF is the [Utrecht publisher copy](https://research-portal.uu.nl/ws/portalfiles/portal/264124617/1-s2.0-S0165242725000662-main.pdf).
The Nemec receipt identifies the Europe PMC JSON response containing the primary
article abstract. Source URLs, versions and full byte checksums are recorded in
`snapshot.as_dict()["sources"]`; RefSeq input checksums come from the existing
species reference registry. If an upstream response changes, reacquiring a URL
does not authorize silently updating its checksum. Raw study files are not
redistributed in the wheel. The packaged snapshot supports offline replay
independently of those upstream files.
