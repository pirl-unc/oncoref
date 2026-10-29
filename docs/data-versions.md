# Data versions and reproducibility

Every managed source has a stable name and a concrete release. `version=None` and `version="latest"` select the newest release registered by the installed OncoRef build. An explicit version selects that release only. OncoRef does not replace a missing historical release with newer data or discover upstream updates during a load.

```python
from oncoref import catalog, gtex

catalog.versions("gtex_tissue_rna")              # ('v11', 'v10', 'v8')
catalog.resolve_version("gtex_tissue_rna")       # 'v11'
catalog.ensure("gtex_tissue_rna")                # newest registered release
catalog.ensure("gtex_tissue_rna", version="v10") # exact older release
catalog.path("gtex_tissue_rna", version="v8")    # inspect cache; no download
catalog.provenance("gtex_tissue_rna", "v11", verify_content=True)

matrix = gtex.tissue_rna(version="v11")
matrix.attrs["version"]                         # 'v11'
matrix.attrs["reference"]                       # 'GENCODE v47'
```

The same selector applies to cancer data and native species assets:

```python
catalog.ensure("per-sample:LUAD", version="5.22.13")
catalog.ensure("cancer-reference-expression-percentiles", version="5.23.25")
catalog.versions("barkbase-matrix", species="dog")
catalog.ensure("barkbase-matrix", species="dog", version="latest")
catalog.provenance("cta-selection-human")
catalog.ensure("cta-evidence", species="dog", version="oncoref-1.20.0")
```

The CLI exposes the same choices:

```sh
oncoref data versions gtex_tissue_rna
oncoref data path gtex_tissue_rna --version v10
oncoref data provenance gtex_tissue_rna --version v11
oncoref data versions barkbase-matrix --species dog
```

| Data kind | What `version` identifies | Historical access |
|---|---|---|
| Reference download, such as GTEx, HPA or an atlas | Upstream release, or accession plus retrieval date when the provider has no release number | Registered versions have separate paths and provenance; new snapshots have independent content hashes |
| Cancer per-sample matrix | OncoRef matrix build, separate from the study's accession and source release | Exact `source-vX.Y.Z` release URL; absent assets fail explicitly |
| Expression bundle | OncoRef data build | Exact `vX.Y.Z` release and matching checksum manifest; incompatible old bundle layouts fail rather than substituting current artifacts |
| Curated human table shipped in the wheel | `oncoref-X.Y.Z` package release | Historical data are extracted from the exact checksum-verified PyPI wheel without installing or executing old code |
| Native dog/mouse asset | Its species manifest's `source_version` | Multiple releases can share a stable asset ID; the newest registered row is the default, older rows remain selectable; reference registrations retain their exact content hashes |
| Packaged native CTA evidence | `oncoref-X.Y.Z`, with the curated snapshot date recorded separately as `source_version` | The same verified-wheel access as human tables, plus the existing reader for separately exported CTA evidence snapshots |

`versions()` lists registered releases and, for bundles/matrices, discovered local historical releases. It is not an inventory of all releases on the internet. Explicit published bundle, matrix and package pins can be requested even before they are cached. An upstream source with only one registered version exposes one version; OncoRef does not invent historical availability. Species with no registered source return an unavailable/unknown-source error and never fall back to human data.

## Provenance has separate identities

The artifact version, genome/annotation release and processing pipeline are different fields. GTEx v11 uses GENCODE v47 and RNA-SeQC 2.4.3; v10 uses GENCODE v39 and RNA-SeQC 2.4.2. Values remain native TPM. [GTEx pipeline release table](https://github.com/broadinstitute/gtex-pipeline/blob/master/rnaseq/README.md).

The common provenance fields identify the source, concrete version, species, local path, URL where available, content hash and verification state. `content_recorded` means the local bytes were hashed; it does not claim that they match an independent publication pin. `verified` means comparison with an expected checksum succeeded. An unavailable checksum remains explicit. Version-separated caches prevent an older selection from reading a current file.

HPA consensus and single-cell summary loaders now default to v25.1; both retain v23. HPA's legacy `normal_tissue.tsv` IHC export remains at v23, the newest registered release of that particular export. Its version is not silently promoted to the newer website release. The current RNA files are preserved in an [immutable HPA 25.1 snapshot release](https://github.com/pirl-unc/oncoref/releases/tag/reference-hpa-v25.1-20261028), with original URLs and archive/content hashes. [HPA release history](https://www.proteinatlas.org/about/releases).

## Saved analyses pin their inputs

Resolve a default once, save the concrete version and checksum, and pass the concrete version when reproducing the analysis. Do not save the word `latest` as a scientific input identity.

The [current CTA analysis](audits/cta-source-walkthrough-20261028/index.md) pins GTEx v11, HPA v23 normal panels, and the Tyser and Petropoulos snapshots. It does not pool v8, v10 and v11 as independent assays. Public loader defaults can advance without changing these saved choices. HPA nTPM and direct GTEx donor-median TPM remain separate evidence streams; a donor median is not a bound on every donor. The GTEx reader retains all 68 native contexts; the current whole-organ adapter uses 52 and explicitly reports the remaining microdissected/cultured contexts.

## Registering another release

Append a concrete version under the existing source name in `oncoref/data/reference-releases.json`, recording its URL, content and archive hashes, byte count, taxon, annotation, units and processing metadata where known. Set `default_version` to the newest validated release and retain the older records. The loader merges these releases with legacy source manifests, so existing source names and pins continue to work. An unversioned provider URL requires a dated snapshot and content pin; archive it when possible so an upstream update cannot make the old version unavailable.

For native species assets, append a release with the same ID and a distinct `source_version`; keep the manifest in release order. Existing assembly/reference registrations resolve their own hashes, even after a newer asset is added. For a new cancer data build, update its build registry and publish the matching release manifest. The catalog resolves these backend details through the same `version=` API.
