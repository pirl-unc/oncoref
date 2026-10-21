# Licensed under the Apache License, Version 2.0.
"""Species is the top-level data boundary; reference identity is a finer boundary."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class Species:
    taxon_id: int
    name: str
    scientific_name: str
    aliases: tuple[str, ...] = ()


_SPECIES = (
    Species(9606, "human", "Homo sapiens", ("homo_sapiens",)),
    Species(
        9615,
        "dog",
        "Canis lupus familiaris",
        ("canine", "canis_lupus_familiaris", "canis familiaris"),
    ),
    Species(10090, "mouse", "Mus musculus", ("mus_musculus",)),
)
DATA_GROUPS = (
    "references",
    "expression",
    "cancer",
    "antigens",
    "normalization",
    "governance",
    "legacy",
)


def available_species() -> tuple[Species, ...]:
    """Known species, including those with incomplete or unavailable data."""
    return _SPECIES


def resolve_species(value: Species | str | int) -> Species:
    """Resolve an explicit taxon/name; unknown species never fall back to human."""
    for species in _SPECIES:
        if isinstance(value, Species) and value == species:
            return species
        if type(value) is int and value == species.taxon_id:
            return species
        if type(value) is str and value.strip().lower() in {
            species.name,
            species.scientific_name.lower(),
            str(species.taxon_id),
            *species.aliases,
        }:
            return species
    raise ValueError(f"Unknown species {value!r}; choose human, dog, mouse or their taxonomy IDs")


def _human_data_group(row):
    category = row["category"]
    if category == "gene-id":
        return "references"
    if category in {"expression", "hpa", "protein-calibration"}:
        return "expression"
    if category == "cta":
        return "antigens"
    if category == "normalization":
        return "normalization"
    if category == "legacy-compat":
        return "legacy"
    if category == "governance":
        return "governance"
    return "cancer"


@dataclass(frozen=True)
class SpeciesDataGroup:
    species: Species
    name: str

    def __post_init__(self):
        object.__setattr__(self, "species", resolve_species(self.species))
        if self.name not in DATA_GROUPS:
            raise ValueError(f"Unknown data group {self.name!r}; choose from {DATA_GROUPS}")

    def datasets(self) -> tuple[dict, ...]:
        """Inspect this species/group without downloading or substituting human data."""
        from . import catalog

        return tuple(
            row for row in catalog.inventory(species=self.species) if row["group"] == self.name
        )

    def sources(self, *, context: str | None = None):
        """Expression sources in native measurement space; optional normal/tumor scope."""
        if self.name != "expression":
            raise ValueError("sources() is available on the expression group")
        from .expression_registry import expression_sources

        return expression_sources(species=self.species, context=context)

    def read(self, source_id, **options):
        if self.name != "expression":
            raise ValueError("read() requires the expression group")
        from .species_expression import read_species_expression

        return read_species_expression(source_id, species=self.species, **options)

    def lineage_overlaps(self):
        if self.name != "expression":
            raise ValueError("lineage_overlaps() requires the expression group")
        from .species_expression import expression_lineage_overlaps

        return expression_lineage_overlaps(self.species)

    def cta_snapshot(self, *, path=None, expected_snapshot_id=None):
        """Frozen species CTA evidence; no downloads or implicit human fallback."""
        if self.name != "antigens":
            raise ValueError("cta_snapshot() requires the antigens group")
        from .species_antigens import cta_snapshot

        return cta_snapshot(self.species, path=path, expected_snapshot_id=expected_snapshot_id)

    def cta_candidates(self, *, reference=None):
        """All CTA nominations, including counterevidence and unresolved mappings."""
        return self.cta_snapshot().candidates(reference=reference)

    def cta_evidence(self, *, reference=None, candidate_id=None):
        """Native study observations, optionally for one returned candidate ID."""
        return self.cta_snapshot().evidence(reference=reference, candidate_id=candidate_id)

    def cta_tissue_annotations(self, *, candidate_id=None):
        """Current cited tissue scopes alongside the frozen CTA evidence."""
        return self.cta_snapshot().tissue_annotations(candidate_id=candidate_id)

    def normal_tissue_scopes(self):
        """Anatomical scope unions applicable to this species, independent of data coverage."""
        from .normal_tissues import normal_tissue_scopes

        return normal_tissue_scopes(species=self.species)

    def resolve_normal_tissue(self, tissue=None, **options):
        """Map a normal source context under this explicit species."""
        from .normal_tissues import resolve_normal_tissue

        return resolve_normal_tissue(tissue, species=self.species, **options)

    def normal_tissue_inventory(self):
        """All canonical concepts, with species applicability and five-way categories."""
        from .normal_tissues import normal_tissue_inventory

        return normal_tissue_inventory(species=self.species)

    def normal_tissue_mapping_inventory(self, *, source_id=None, source_version=None):
        """Bundled native mappings; an empty result never implies absent expression."""
        from .normal_tissues import normal_tissue_mapping_inventory

        return normal_tissue_mapping_inventory(
            species=self.species, source_id=source_id, source_version=source_version
        )

    def cta_gene_ids(self, *, reference):
        """Unambiguously mapped candidate IDs for an explicit annotation; not admission."""
        return self.cta_snapshot().gene_ids(reference=reference)

    def assets(self):
        from .species_assets import species_assets

        return species_assets(self.species, group=self.name)

    def registered_references(self):
        if self.name != "references":
            raise ValueError("registered_references() requires the references group")
        from .species_assets import registered_references

        return registered_references(self.species)

    def ensure_reference(self, reference_id, **options):
        if self.name != "references":
            raise ValueError("ensure_reference() requires the references group")
        from .species_assets import ensure_reference

        return ensure_reference(reference_id, species=self.species, **options)

    def ensure(self, asset_id, **options):
        from .species_assets import ensure_asset, species_asset

        asset = species_asset(asset_id, species=self.species)
        if asset.group != self.name:
            raise ValueError(f"{asset_id} belongs to {asset.group}, not {self.name}")
        return ensure_asset(asset, **options)

    def coverage(self) -> dict:
        """Inventory availability is not a claim of complete biological coverage."""
        rows = self.datasets()
        return {
            "taxon_id": self.species.taxon_id,
            "species": self.species.name,
            "group": self.name,
            "registered_datasets": len(rows),
            "available_datasets": sum(bool(row["available"]) for row in rows),
            "coverage": "partial" if rows else "unavailable",
        }


@dataclass(frozen=True)
class SpeciesData:
    species: Species

    def __post_init__(self):
        object.__setattr__(self, "species", resolve_species(self.species))

    def group(self, name: str) -> SpeciesDataGroup:
        return SpeciesDataGroup(self.species, name)

    @property
    def groups(self) -> tuple[SpeciesDataGroup, ...]:
        return tuple(self.group(name) for name in DATA_GROUPS)

    @property
    def references(self):
        return self.group("references")

    @property
    def expression(self):
        return self.group("expression")

    @property
    def cancer(self):
        return self.group("cancer")

    @property
    def antigens(self):
        return self.group("antigens")

    @property
    def normalization(self):
        return self.group("normalization")

    def inventory(self):
        from .catalog import inventory

        return inventory(species=self.species)


def for_species(species: Species | str | int) -> SpeciesData:
    """Enter a species explicitly, then select references/expression/cancer/etc."""
    return SpeciesData(resolve_species(species))
