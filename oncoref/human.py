# Licensed under the Apache License, Version 2.0.
"""Human datasets and conceptual data groups.

Dataset modules are shared with their top-level imports, for example
``human.human_protein_atlas`` is ``oncoref.human_protein_atlas``.
"""

from . import human_protein_atlas
from .species import for_species as _for_species

_data = _for_species(9606)
species = _data.species
references = _data.references
expression = _data.expression
cancer = _data.cancer
antigens = _data.antigens
normalization = _data.normalization
governance = _data.group("governance")
legacy = _data.group("legacy")
groups = (references, expression, cancer, antigens, normalization, governance, legacy)
inventory = _data.inventory

# Named dataset modules exposed here; inventory() describes the complete catalog.
datasets = (human_protein_atlas,)

__all__ = [
    "antigens",
    "cancer",
    "datasets",
    "expression",
    "governance",
    "groups",
    "human_protein_atlas",
    "inventory",
    "legacy",
    "normalization",
    "references",
    "species",
]
