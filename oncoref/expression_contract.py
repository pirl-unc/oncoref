# Licensed under the Apache License, Version 2.0.
"""Immutable, serializable source and measurement contracts for expression facts.

Reference equivalence, measurement comparability and species are separate checks.
Unknown historical provenance remains unknown; a digest does not repair it.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from types import MappingProxyType

from .reference_identity import ReferenceIdentity, _identity, _label, _string_map

FEATURE_LEVELS = {"gene", "transcript", "promoter", "protein", "protein_group", "unknown"}
REFERENCE_STATUSES = {"verified", "unresolved", "legacy_human"}
EXPRESSION_CONTEXTS = {"normal", "tumor", "cell_line", "mixed", "unknown"}


def _fields(record, expected, kind):
    if not isinstance(record, Mapping) or set(record) != set(expected):
        raise ValueError(f"{kind} requires exact schema fields: {sorted(expected)}")


def _optional_label(value, name):
    if value is not None:
        _label(value, name)


@dataclass(frozen=True)
class MeasurementTransform:
    operation: str
    input_unit: str
    output_unit: str
    assumptions: tuple[str, ...]
    input_features: int | None = None
    output_features: int | None = None

    def __post_init__(self):
        for name in ("operation", "input_unit", "output_unit"):
            _label(getattr(self, name), name)
        if not isinstance(self.assumptions, (tuple, list)) or not self.assumptions:
            raise ValueError("Transform assumptions must be explicitly recorded")
        for item in self.assumptions:
            _label(item, "assumption")
        object.__setattr__(self, "assumptions", tuple(self.assumptions))
        for name in ("input_features", "output_features"):
            value = getattr(self, name)
            if value is not None and (type(value) is not int or value < 0):
                raise ValueError(f"{name} must be a nonnegative integer or unknown")
        if (self.input_features is None) != (self.output_features is None):
            raise ValueError("Feature coverage requires both input and output counts")

    def as_dict(self):
        return {
            "operation": self.operation,
            "input_unit": self.input_unit,
            "output_unit": self.output_unit,
            "assumptions": list(self.assumptions),
            "input_features": self.input_features,
            "output_features": self.output_features,
        }

    @classmethod
    def from_dict(cls, record):
        _fields(record, cls.__dataclass_fields__, "Measurement transform")
        return cls(**record)


@dataclass(frozen=True)
class ExpressionMeasurement:
    feature_level: str
    unit: str | None
    scale: str | None = None
    is_proxy: bool | None = None
    transformations: tuple[MeasurementTransform, ...] = ()

    def __post_init__(self):
        if self.feature_level not in FEATURE_LEVELS:
            raise ValueError(f"Unsupported measured feature level {self.feature_level!r}")
        _optional_label(self.unit, "unit")
        _optional_label(self.scale, "scale")
        if self.is_proxy is not None and type(self.is_proxy) is not bool:
            raise ValueError("is_proxy must be boolean or unknown")
        if not isinstance(self.transformations, (tuple, list)) or not all(
            isinstance(t, MeasurementTransform) for t in self.transformations
        ):
            raise ValueError("transformations must contain MeasurementTransform records")
        transforms = tuple(self.transformations)
        for left, right in zip(transforms, transforms[1:]):
            if left.output_unit != right.input_unit:
                raise ValueError("Transformation units do not form a consistent chain")
        if transforms and transforms[-1].output_unit != self.unit:
            raise ValueError("Final transformation output must match the measurement unit")
        object.__setattr__(self, "transformations", transforms)

    def as_dict(self):
        return {
            "feature_level": self.feature_level,
            "unit": self.unit,
            "scale": self.scale,
            "is_proxy": self.is_proxy,
            "transformations": [t.as_dict() for t in self.transformations],
        }

    @classmethod
    def from_dict(cls, record):
        _fields(record, cls.__dataclass_fields__, "Expression measurement")
        return cls(
            **{
                **record,
                "transformations": tuple(
                    MeasurementTransform.from_dict(t) for t in record["transformations"]
                ),
            }
        )


@dataclass(frozen=True)
class ExpressionContract:
    taxon_id: int
    source_id: str
    source_version: str | None
    context: str
    reference: ReferenceIdentity | None
    reference_status: str
    feature_namespace: str | None
    measurement: ExpressionMeasurement
    content_provenance: Mapping[str, str]
    provenance_notes: str

    def __post_init__(self):
        if type(self.taxon_id) is not int or self.taxon_id < 1:
            raise ValueError("Positive integer taxonomy identifier required")
        _label(self.source_id, "source_id")
        _optional_label(self.source_version, "source_version")
        _optional_label(self.feature_namespace, "feature_namespace")
        _label(self.provenance_notes, "provenance_notes")
        if self.context not in EXPRESSION_CONTEXTS:
            raise ValueError(f"Unknown expression context {self.context!r}")
        if self.reference_status not in REFERENCE_STATUSES:
            raise ValueError("Unknown reference provenance status")
        if not isinstance(self.measurement, ExpressionMeasurement):
            raise ValueError("ExpressionMeasurement required")
        if self.reference is not None:
            if (
                not isinstance(self.reference, ReferenceIdentity)
                or self.reference.taxon_id != self.taxon_id
            ):
                raise ValueError("Expression/reference taxon mismatch")
            if self.reference_status != "verified":
                raise ValueError("An unresolved reference must remain absent")
        elif self.reference_status == "verified":
            raise ValueError("Verified provenance requires a complete ReferenceIdentity")
        if self.reference_status == "legacy_human" and self.taxon_id != 9606:
            raise ValueError("Legacy human compatibility cannot be used for another species")
        if not isinstance(self.content_provenance, Mapping):
            raise ValueError("content_provenance must be an asset digest mapping")
        hashes = (
            _string_map(self.content_provenance, "content_provenance", digests=True)
            if self.content_provenance
            else {}
        )
        if self.reference_status == "verified" and (
            not hashes or self.source_version is None or self.feature_namespace is None
        ):
            raise ValueError(
                "Verified expression requires pinned source bytes, version and namespace"
            )
        object.__setattr__(self, "content_provenance", MappingProxyType(hashes))

    def _payload(self):
        return {
            "schema": "oncoref.expression-contract",
            "schema_version": 1,
            "taxon_id": self.taxon_id,
            "source_id": self.source_id,
            "source_version": self.source_version,
            "context": self.context,
            "reference": self.reference.as_dict() if self.reference else None,
            "reference_status": self.reference_status,
            "feature_namespace": self.feature_namespace,
            "measurement": self.measurement.as_dict(),
            "content_provenance": dict(self.content_provenance),
            "provenance_notes": self.provenance_notes,
        }

    @property
    def contract_id(self):
        """Include all provenance/measurement semantics in derived-artifact cache keys."""
        return _identity("expression-contract", 1, self._payload())

    def as_dict(self):
        return {**self._payload(), "contract_id": self.contract_id}

    @classmethod
    def from_dict(cls, record):
        _fields(
            record,
            {*cls.__dataclass_fields__, "schema", "schema_version", "contract_id"},
            "Expression contract",
        )
        if (
            record["schema"] != "oncoref.expression-contract"
            or type(record["schema_version"]) is not int
            or record["schema_version"] != 1
        ):
            raise ValueError("Unsupported expression contract schema/version")
        values = {key: record[key] for key in cls.__dataclass_fields__}
        values["measurement"] = ExpressionMeasurement.from_dict(values["measurement"])
        if values["reference"] is not None:
            values["reference"] = ReferenceIdentity.from_dict(values["reference"])
        result = cls(**values)
        if record["contract_id"] != result.contract_id:
            raise ValueError("Expression contract digest does not match its contents")
        return result

    def require_reference(self):
        if self.reference_status != "verified":
            raise ValueError(
                f"{self.source_id}: reference provenance is {self.reference_status}; an explicit verified reference is required"
            )
        return self.reference

    def require_feature_level(self, feature_level):
        if feature_level not in FEATURE_LEVELS or feature_level == "unknown":
            raise ValueError("A known measured feature level is required")
        if self.measurement.feature_level != feature_level:
            raise ValueError(
                f"Measured {self.measurement.feature_level} values cannot establish {feature_level} abundance"
            )


@dataclass(frozen=True)
class ReferenceCrosswalk:
    """Directed annotation join evidence, never permission to sum/allocate abundance.

    Pairs retain all alternatives. Unmapped IDs remain unresolved; symbols and
    orthology cannot implicitly construct this record. Cross-species joins are
    intentionally unsupported by expression APIs.
    """

    source: ReferenceIdentity
    target: ReferenceIdentity
    feature_level: str
    source_namespace: str
    target_namespace: str
    version: str
    evidence: Mapping[str, str]
    pairs: tuple[tuple[str, str], ...]

    def __post_init__(self):
        if not isinstance(self.source, ReferenceIdentity) or not isinstance(
            self.target, ReferenceIdentity
        ):
            raise ValueError("Crosswalk requires two reference identities")
        if self.source.taxon_id != self.target.taxon_id:
            raise ValueError("Cross-species expression crosswalks are unsupported")
        if self.feature_level not in FEATURE_LEVELS - {"unknown"}:
            raise ValueError("Crosswalk requires a known feature level")
        for name in ("source_namespace", "target_namespace", "version"):
            _label(getattr(self, name), name)
        object.__setattr__(
            self, "evidence", MappingProxyType(_string_map(self.evidence, "evidence", digests=True))
        )
        if not isinstance(self.pairs, (tuple, list)):
            raise ValueError("Crosswalk pairs must be a sequence")
        pairs = []
        for pair in self.pairs:
            if not isinstance(pair, (tuple, list)) or len(pair) != 2:
                raise ValueError("Each crosswalk pair needs source and target feature IDs")
            for value in pair:
                _label(value, "crosswalk feature ID")
            pairs.append(tuple(pair))
        object.__setattr__(self, "pairs", tuple(sorted(set(pairs))))

    def as_dict(self):
        payload = {
            "schema": "oncoref.reference-crosswalk",
            "schema_version": 1,
            "source": self.source.as_dict(),
            "target": self.target.as_dict(),
            "feature_level": self.feature_level,
            "source_namespace": self.source_namespace,
            "target_namespace": self.target_namespace,
            "version": self.version,
            "evidence": dict(self.evidence),
            "pairs": [list(p) for p in self.pairs],
        }
        return {**payload, "crosswalk_id": _identity("reference-crosswalk", 1, payload)}

    @classmethod
    def from_dict(cls, record):
        _fields(
            record,
            {*cls.__dataclass_fields__, "schema", "schema_version", "crosswalk_id"},
            "Reference crosswalk",
        )
        if (
            record["schema"] != "oncoref.reference-crosswalk"
            or type(record["schema_version"]) is not int
            or record["schema_version"] != 1
        ):
            raise ValueError("Unsupported crosswalk schema/version")
        values = {key: record[key] for key in cls.__dataclass_fields__}
        for key in ("source", "target"):
            values[key] = ReferenceIdentity.from_dict(values[key])
        result = cls(**values)
        if result.as_dict()["crosswalk_id"] != record["crosswalk_id"]:
            raise ValueError("Crosswalk digest does not match its contents")
        return result

    def targets(self, feature_id):
        """All mapping alternatives; an empty tuple means unresolved."""
        return tuple(target for source, target in self.pairs if source == feature_id)


def validate_expression_join(left, right, *, crosswalk=None, combine_values=False):
    """Guard a feature join or derived abundance summary before touching the values.

    A crosswalk permits annotation joins only. Combining values additionally
    requires matching complete measurement semantics and the same exact reference.
    Legacy/unresolved sources can be inspected but cannot pass this guard.
    """
    if left.taxon_id != right.taxon_id:
        raise ValueError(f"Expression taxon mismatch: {left.taxon_id} != {right.taxon_id}")
    a, b = left.require_reference(), right.require_reference()
    if (
        left.measurement.feature_level == "unknown"
        or left.measurement.feature_level != right.measurement.feature_level
    ):
        raise ValueError("Expression measured feature levels must be known and equal")
    matches = a == b and left.feature_namespace == right.feature_namespace
    if not matches:
        if crosswalk is None or combine_values:
            raise ValueError(
                "Expression reference/namespace mismatch; explicit crosswalk required for annotation joins; abundance allocation is unsupported"
            )
        if not isinstance(crosswalk, ReferenceCrosswalk) or (
            crosswalk.source != a
            or crosswalk.target != b
            or crosswalk.feature_level != left.measurement.feature_level
            or crosswalk.source_namespace != left.feature_namespace
            or crosswalk.target_namespace != right.feature_namespace
        ):
            raise ValueError("Crosswalk does not support this directed expression join")
    if combine_values:
        measurement = left.measurement
        if (
            measurement != right.measurement
            or measurement.unit is None
            or measurement.scale is None
            or measurement.is_proxy is None
        ):
            raise ValueError("Combining abundance requires identical, known measurement semantics")
