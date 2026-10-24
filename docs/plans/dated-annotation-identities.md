# OncoRef #609: pinned numbered and dated annotation identities

## Contract and scope

Allow an explicitly selected positive numbered release or exact ASCII YYYY_MM
annotation date, retaining int/string type and value in GeneIdentity.as_dict().
Reject booleans, zero/negative releases, numeric strings, noncanonical/invalid
months and year zero. Do not coerce dates or interpret "latest".

Use genome.assembly when the provider exposes it; an explicitly unknown assembly
must fail rather than falling back to a convenient cache label. Providers without
that attribute retain reference_name compatibility. Require explicit homo_sapiens
and resolved GRCh37 or GRCh38. Assembly names/patches and accessions remain the
upstream resolver's responsibility; do not implement another mapping catalogue.

Keep GeneIdentity's fields and identity contract version unchanged, except widen
annotation_release's annotation to int | str. Do not add optional fields to old
records or change reference-table hashes, CTA membership, alias resolution,
source identities, data versions, cache paths or cached files.

## Implementation and tests

1. Change only the annotation boundary, release typing and relevant docs/version.
2. Test real PyEnsembl metadata selections for numbered GRCh37/GRCh38, dated
   GRCh37 and GRCh38, explicit accession and custom display label. Perform a
   public resolve_gene_identity lookup against a real tiny GTF-backed
   EnsemblRelease (temporary mirror), not a mocked annotation resolver.
3. Reject invalid releases, missing/nonhuman species and unknown assemblies before
   lookup. Keep older annotation-provider compatibility for reference_name-only
   objects. Test a present-but-unknown assembly doesn't use the fallback.
4. Capture a numbered GeneIdentity from unmodified OncoRef 1.19.0 as a compatibility
   fixture. Load/serialize it unchanged and compare current numbered results and
   reference hashes, excluding the expected package-version field change.
5. Test exact dated JSON provenance round-trips. Snapshot temporary installed
   files before and after lookup/reopen; preserve all stable file metadata and
   bytes. No writes to user caches or the shared environment.

## Shipping

Use an isolated feature checkout. Baseline tests run against unchanged main first.
Bump the next available patch version after checking concurrent main and PyPI.
Run repository lint and full tests, docs, pip check, package checks, and CI on
Python 3.10–3.14. Merge and deploy.sh from clean main, verify both public artifacts.
Then update existing Vaxrank PR #665 to this published OncoRef minimum, remove its
upstream mocks, and add the unmocked dated admission test before releasing it.
