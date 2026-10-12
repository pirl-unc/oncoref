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
