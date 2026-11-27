"""Verified, offline sequence records for the explicitly versioned CTA policy.

Acquisition is opt-in. Original records, every occurrence and interrupted valid
segments remain available; exact sequence grouping never discards an exclusion.
"""

from __future__ import annotations

import gzip
import hashlib
import io
import json
import re
import tempfile
import zipfile
from collections import defaultdict
from copy import deepcopy
from dataclasses import dataclass
from importlib.resources import files
from pathlib import Path

from . import data_file_renames

SCHEMA = "oncoref.sequence-bundle"
ROLE_POLICY_ID = "oncoref.target-selection-sequence-roles.2026-11-21.1"
NCORF_TE_CONTRACT_ID = "oncoref.ncorf-te-screening.2026-11-22.1"
PAYLOADS = frozenset(
    {"records.jsonl.gz", "coverage.json", "decisions.csv.gz", "orfdata-proteome.json.gz"}
)
ROLE_ORDER = {"allow": 0, "review": 1, "exclude": 2}


def target_selection_dispositions():
    """Assign representative-protein sequence roles from ``target_selection.results()``.

    ``Unsafe`` exclusions take precedence across occurrences. Qualified representatives
    receive ``shared_sequence == "allow"`` only when the protein passes all six
    ``target_selection.SHARED_CHECKS`` and its recorded normal home is screened
    (``CTA``, ``Cancer Antigen only``, ``Reproductive Antigen only``, or eligible
    ``Other:Expression``). ``Other:Review``, ``Other:Sequence``, ``Other:Unscreened tissue``,
    and ``Other:Expression`` failing expression bounds receive ``"review"``.
    ``direct_target_candidate`` marks confirmed ``CTA`` proteins, while
    ``cancer_antigen_candidate`` and ``reproductive_antigen_candidate`` preserve the
    broader set memberships explicitly.
    """
    import pandas as pd

    from . import target_selection

    frame = target_selection.results()
    shared_ok = frame[list(target_selection.SHARED_CHECKS)].astype(bool).all(axis=1)
    excluded = frame.category.eq("Unsafe")
    unscreened_home = frame.classification.eq("Other:Unscreened tissue")
    allowed = ~excluded & shared_ok & ~unscreened_home

    out = frame[
        [
            "entity_id",
            "protein_group_label",
            "member_gene_ids",
            "representative_sequence_sha256s",
            "category",
            "subtype",
            "classification",
            "is_cta",
            "is_cancer_antigen",
            "is_reproductive_antigen",
            "is_target_candidate",
            "cta_tier",
            "decision_reason",
            "clinical_clearance_scope",
            "clinical_cleared_organs",
            "clinical_cleared_epitopes",
            "clinical_hla_alleles",
            "clinical_cleared_transcripts",
            "clinical_unverified_isoforms",
            "clinical_trial_source_ids",
            "closed_homolog_family_id",
            "co_silenced_coding_sibling_ids",
            "closed_self_paralog_ids",
            "exact_unique_9mers",
            "safe_2mm_9mers",
            "unsafe_near_self_1mm_pct",
            "unsafe_shared_tryptic_peptides",
        ]
    ].copy()
    out["gene_ids"] = out.member_gene_ids.map(lambda raw: tuple(json.loads(raw)))
    out["representative_sequence_sha256s"] = out.representative_sequence_sha256s.map(
        lambda raw: tuple(json.loads(raw)) if isinstance(raw, str) and raw else ()
    )
    out["shared_sequence"] = pd.Series("review", index=out.index, dtype="string")
    out.loc[allowed, "shared_sequence"] = "allow"
    out.loc[excluded, "shared_sequence"] = "exclude"
    out["direct_target_candidate"] = out.is_cta.astype(bool) & allowed
    out["cancer_antigen_candidate"] = out.is_cancer_antigen.astype(bool) & allowed
    out["reproductive_antigen_candidate"] = out.is_reproductive_antigen.astype(bool) & allowed
    out["role_policy_id"] = ROLE_POLICY_ID
    out["selection_policy_id"] = target_selection.POLICY_ID
    return out


def digest(content: bytes) -> str:
    return hashlib.sha256(content).hexdigest()


def content_id(manifest: dict) -> str:
    body = {k: v for k, v in manifest.items() if k != "bundle_id"}
    return "oncoref:sequence-bundle:v1:sha256:" + digest(
        json.dumps(body, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()
    )


def descriptor() -> dict:
    """Packaged current asset pin; reading it performs no download or cache write."""
    return json.loads(files("oncoref").joinpath("data/sequence-bundle.json").read_text())


def _decode_records(payload):
    with gzip.GzipFile(fileobj=io.BytesIO(payload)) as handle:
        for line in handle:
            if line.strip():
                yield json.loads(line)


@dataclass(frozen=True)
class SequenceBundle:
    """Immutable verified bytes; public objects are detached on every access."""

    _manifest: bytes
    _records: bytes
    _coverage: bytes
    _orfdata: bytes
    _project_current_roles: bool = False

    @property
    def identity(self):
        return json.loads(self._manifest)["bundle_id"]

    @property
    def manifest(self):
        return json.loads(self._manifest)

    @property
    def coverage(self):
        return json.loads(self._coverage)

    def records(self):
        """Yield full originals with qualification, roles and all source occurrences."""
        if not self._project_current_roles:
            yield from _decode_records(self._records)
            return
        disp = target_selection_dispositions()
        gene_roles: dict[str, tuple[str, bool, bool, bool, str]] = {}
        for r in disp.itertuples(index=False):
            info = (
                str(r.shared_sequence),
                bool(r.direct_target_candidate),
                bool(r.cancer_antigen_candidate),
                bool(r.reproductive_antigen_candidate),
                str(r.entity_id),
            )
            for gid in r.gene_ids:
                gene_roles[str(gid)] = info
        for record in _decode_records(self._records):
            for occ in record["occurrences"]:
                gid = str(occ.get("gene_id", ""))
                if gid in gene_roles:
                    s_role, is_cta, is_ca, is_ra, eid = gene_roles[gid]
                    occ["shared_sequence"] = s_role
                    occ["direct_target_candidate"] = is_cta
                    occ["cancer_antigen_candidate"] = is_ca
                    occ["reproductive_antigen_candidate"] = is_ra
                    occ["decision_ids"] = [eid]
            record["shared_sequence"] = max(
                (o["shared_sequence"] for o in record["occurrences"]), key=ROLE_ORDER.get
            )
            allow_rec = record["shared_sequence"] == "allow"
            record["direct_target_candidate"] = allow_rec and any(
                bool(o.get("direct_target_candidate")) for o in record["occurrences"]
            )
            record["cancer_antigen_candidate"] = allow_rec and any(
                bool(o.get("cancer_antigen_candidate")) for o in record["occurrences"]
            )
            record["reproductive_antigen_candidate"] = allow_rec and any(
                bool(o.get("reproductive_antigen_candidate")) for o in record["occurrences"]
            )
            yield record

    def segments(self, *, role=None, min_length=1):
        """Yield coordinate-bearing runs; never concatenate across an interruption."""
        if role is not None and role not in ROLE_ORDER:
            raise ValueError("Unknown sequence role")
        if type(min_length) is not int or min_length < 1:
            raise ValueError("min_length must be a positive integer")
        for record in self.records():
            if role is not None and record["shared_sequence"] != role:
                continue
            for segment in record["valid_segments"]:
                if segment["end"] - segment["start"] >= min_length:
                    yield {
                        **segment,
                        "record_id": record["record_id"],
                        "shared_sequence": record["shared_sequence"],
                        "occurrences": deepcopy(record["occurrences"]),
                        "reference_id": record["reference_id"],
                    }

    def groups(self):
        """Exact original-sequence groups; retain records/occurrences and role precedence.

        Full originals containing ambiguity are grouped only by their exact raw
        bytes; use segments for matching. Group identity does not prove production.
        """
        grouped = defaultdict(list)
        for record in self.records():
            grouped[record["sequence_sha256"]].append(record)
        for key, records in sorted(grouped.items()):
            yield {
                "sequence_sha256": key,
                "sequence": records[0]["sequence"],
                "shared_sequence": max((r["shared_sequence"] for r in records), key=ROLE_ORDER.get),
                "records": records,
                "record_count": len(records),
                "occurrence_count": sum(len(r["occurrences"]) for r in records),
            }

    def as_orfdata(self):
        """Load the preserved native ORFdata proteome with its original identity.

        Requires the optional quantitative extra. Policy roles remain in this
        bundle; the returned raw Proteome does not itself apply CTA exclusions.
        """
        try:
            from orfdata import Proteome
        except ImportError as error:
            raise ImportError("Install oncoref[quantitative] to use the ORFdata adapter") from error
        with tempfile.TemporaryDirectory(prefix="oncoref-orfdata-") as folder:
            path = Path(folder) / "proteome.json.gz"
            path.write_bytes(self._orfdata)
            proteome = Proteome.load(path)
        if proteome.identity != self.manifest["orfdata_identity"]:
            raise ValueError("Native ORFdata content identity mismatch")
        return proteome

    def to_sequences_table(self, *, roles: tuple[str, ...] | None = None):
        """Expand bundle records and occurrences into an ``antigen_sets``/``PeptideIndex`` table.

        Preserves every transcript/gene occurrence without requiring ``pyensembl``.
        When ``roles`` is supplied (e.g. ``('exclude', 'review')``), keeps every
        occurrence of any record whose record-level ``shared_sequence`` or
        occurrence-level ``shared_sequence`` belongs to ``roles`` so ``exclude > review > allow``
        precedence is never weakened.
        """
        import pandas as pd

        allowed_roles = None
        if roles is not None:
            allowed_roles = tuple(dict.fromkeys(roles))
            if not allowed_roles or any(r not in ROLE_ORDER for r in allowed_roles):
                raise ValueError(
                    f"Unknown sequence role in {roles!r}; choose from {tuple(ROLE_ORDER)}"
                )
            allowed_set = set(allowed_roles)

        manifest = self.manifest
        coverage = self.coverage
        rows = []
        for record in self.records():
            rec_role = record["shared_sequence"]
            seq = record["sequence"]
            seq_hash = record["sequence_sha256"]
            rid = record["record_id"]
            qual = record.get("qualification") or {}
            for occ in record["occurrences"]:
                occ_role = occ.get("shared_sequence", rec_role)
                if (
                    allowed_roles is not None
                    and rec_role not in allowed_set
                    and occ_role not in allowed_set
                ):
                    continue
                raw_tid = str(occ.get("transcript_id") or qual.get("transcript_id") or rid)
                rows.append(
                    {
                        "gene_id": str(occ.get("gene_id", "")),
                        "transcript_id": raw_tid.split(".", 1)[0],
                        "transcript_accession": raw_tid,
                        "transcript_biotype": str(
                            occ.get("transcript_biotype")
                            or qual.get("transcript_biotype")
                            or "protein_coding"
                        ),
                        "record_id": str(rid),
                        "occurrence_id": str(occ["occurrence_id"]),
                        "sequence": seq,
                        "sequence_sha256": seq_hash,
                        "shared_sequence": str(occ_role),
                        "record_shared_sequence": str(rec_role),
                        "direct_target_candidate": bool(
                            occ.get(
                                "direct_target_candidate",
                                record.get("direct_target_candidate", False),
                            )
                        ),
                        "cancer_antigen_candidate": bool(
                            occ.get(
                                "cancer_antigen_candidate",
                                record.get("cancer_antigen_candidate", False),
                            )
                        ),
                        "reproductive_antigen_candidate": bool(
                            occ.get(
                                "reproductive_antigen_candidate",
                                record.get("reproductive_antigen_candidate", False),
                            )
                        ),
                    }
                )
        columns = [
            "gene_id",
            "transcript_id",
            "transcript_accession",
            "transcript_biotype",
            "record_id",
            "occurrence_id",
            "sequence",
            "sequence_sha256",
            "shared_sequence",
            "record_shared_sequence",
            "direct_target_candidate",
            "cancer_antigen_candidate",
            "reproductive_antigen_candidate",
        ]
        out = pd.DataFrame(rows, columns=columns)
        ref = manifest["reference"]
        from . import target_selection

        effective_policy_id = (
            target_selection.POLICY_ID
            if self._project_current_roles
            else manifest.get("selection_policy_id", "")
        )
        out.attrs.update(
            reference=ref["reference_id"],
            reference_id=ref["reference_id"],
            taxon_id=ref.get("taxon_id", 9606),
            bundle_id=self.identity,
            selection_policy_id=effective_policy_id,
            role_policy_id=manifest.get("role_policy_id", ROLE_POLICY_ID),
            sequence_source="sequence_bundle",
            roles=list(allowed_roles) if allowed_roles is not None else list(ROLE_ORDER),
            missing_translations=list(coverage.get("missing_translations", [])),
        )
        return out


def load(path=None, *, expected_id=None, require_current=True) -> SequenceBundle:
    """Verify an unpacked bundle without network access; default to the current pin.

    For a named historical bundle, explicitly pass its trusted expected_id and
    require_current=False. A manifest from the bundle itself is not a trusted pin.
    Missing coverage cannot be represented by an empty record list.
    """
    pin = descriptor() if require_current else None
    if not require_current and not expected_id:
        raise ValueError("Historical loading requires an explicit trusted bundle identity")
    if pin and expected_id and expected_id != pin["bundle_id"]:
        raise ValueError("Expected bundle is not the current package policy")
    expected_id = pin["bundle_id"] if pin else expected_id
    folder = Path(path) if path is not None else ensure(download=False)
    manifest_bytes = (folder / "manifest.json").read_bytes()
    manifest = json.loads(manifest_bytes)
    if (manifest.get("schema"), manifest.get("schema_version")) != (SCHEMA, 1):
        raise ValueError("Unsupported sequence bundle schema")
    if manifest.get("bundle_id") != expected_id or content_id(manifest) != expected_id:
        raise ValueError("Sequence bundle identity mismatch")
    if set(manifest["payloads"]) != PAYLOADS:
        raise ValueError("Missing or unexpected sequence bundle payload")
    if pin:
        baseline_deps: dict[str, str] = {}
        delta_gz = files("oncoref").joinpath("data/baseline-2026-11-20-deltas.json.gz")
        delta_res = files("oncoref").joinpath("data/baseline-2026-11-20-deltas.json")
        if delta_gz.is_file():
            delta_doc = json.loads(gzip.decompress(delta_gz.read_bytes()).decode("utf-8"))
            if manifest.get("selection_policy_id") == delta_doc.get("baseline_policy_id"):
                baseline_deps = dict(delta_doc.get("baseline_package_dependencies_sha256", {}))
        elif delta_res.is_file():
            delta_doc = json.loads(delta_res.read_text(encoding="utf-8"))
            if manifest.get("selection_policy_id") == delta_doc.get("baseline_policy_id"):
                baseline_deps = dict(delta_doc.get("baseline_package_dependencies_sha256", {}))
        for name, expected in manifest["package_dependencies"].items():
            if Path(name).name != name:
                raise ValueError("Invalid package dependency name")
            live = files("oncoref").joinpath(f"data/{name}")
            if live.is_file() and digest(live.read_bytes()) == expected:
                continue
            if baseline_deps.get(name) == expected:
                continue
            if digest(data_file_renames.as_released(name)) != expected:
                raise ValueError(f"Mixed package/bundle policy: {name}")
    payloads = {}
    for name, entry in manifest["payloads"].items():
        payload = (folder / name).read_bytes()
        if len(payload) != entry["bytes"] or digest(payload) != entry["sha256"]:
            raise ValueError(f"Sequence bundle payload checksum mismatch: {name}")
        payloads[name] = payload
    coverage = json.loads(payloads["coverage.json"])
    if coverage.get("reference_id") != manifest["reference"]["reference_id"]:
        raise ValueError("Coverage reference mismatch")
    ids, occurrences = set(), set()
    for record in _decode_records(payloads["records.jsonl.gz"]):
        rid = record["record_id"]
        if rid in ids or record["reference_id"] != coverage["reference_id"]:
            raise ValueError("Repeated record or mixed reference")
        ids.add(rid)
        if not record["occurrences"] or record["shared_sequence"] not in ROLE_ORDER:
            raise ValueError("Record lacks an occurrence or explicit role")
        sequence = record["sequence"]
        if digest(sequence.encode()) != record["sequence_sha256"]:
            raise ValueError("Original sequence checksum mismatch")
        previous = 0
        for segment in record["valid_segments"]:
            start, end = segment["start"], segment["end"]
            if not 0 <= previous <= start < end <= len(sequence):
                raise ValueError("Invalid segment coordinates")
            value = sequence[start:end]
            if (
                value != segment["sequence"]
                or set(value) - set("ACDEFGHIKLMNPQRSTVWY")
                or digest(value.encode()) != segment["sha256"]
            ):
                raise ValueError("Invalid segment sequence")
            previous = end
        expected_segments = [
            (m.start(), m.end()) for m in re.finditer(r"[ACDEFGHIKLMNPQRSTVWY]+", sequence)
        ]
        if [(s["start"], s["end"]) for s in record["valid_segments"]] != expected_segments:
            raise ValueError("Valid segments do not preserve every complete AA20 run")
        for occurrence in record["occurrences"]:
            oid = occurrence["occurrence_id"]
            if oid in occurrences:
                raise ValueError("Duplicate occurrence identity")
            occurrences.add(oid)
            if occurrence.get("shared_sequence") not in ROLE_ORDER:
                raise ValueError("Occurrence lacks an explicit sequence role")
        role = max((o["shared_sequence"] for o in record["occurrences"]), key=ROLE_ORDER.get)
        if record["shared_sequence"] != role:
            raise ValueError("Record lost occurrence role precedence")
    if len(ids) != manifest["record_count"] or len(occurrences) != manifest["occurrence_count"]:
        raise ValueError("Sequence bundle count mismatch")
    from . import target_selection

    project_live = bool(pin and manifest.get("selection_policy_id") != target_selection.POLICY_ID)
    return SequenceBundle(
        manifest_bytes,
        payloads["records.jsonl.gz"],
        payloads["coverage.json"],
        payloads["orfdata-proteome.json.gz"],
        project_live,
    )


def cache_root(cache_dir=None) -> Path:
    """Base sequence-bundle cache root without creating directories."""
    from . import _downloads

    return _downloads.resolve_cache_subdir(
        "sequence-bundles",
        override=cache_dir,
        envkeys=("ONCOREF_SEQUENCE_CACHE_DIR",),
    )


def is_cached(*, cache_dir=None) -> bool:
    """Return True when the current pinned sequence bundle is staged locally (#789).

    Performs zero network calls and creates no directories on disk.
    """
    pin = descriptor()
    destination = cache_root(cache_dir) / pin["bundle_id"].rsplit(":", 1)[-1]
    manifest_path = destination / "manifest.json"
    if not manifest_path.is_file():
        return False
    try:
        manifest = json.loads(manifest_path.read_bytes())
    except (ValueError, OSError):
        return False
    if (manifest.get("schema"), manifest.get("schema_version")) != (SCHEMA, 1):
        return False
    if (
        manifest.get("bundle_id") != pin["bundle_id"]
        or content_id(manifest) != pin["bundle_id"]
        or set(manifest.get("payloads", ())) != PAYLOADS
    ):
        return False
    try:
        return all(
            (destination / name).is_file() and (destination / name).stat().st_size == entry["bytes"]
            for name, entry in manifest["payloads"].items()
        )
    except (KeyError, OSError, TypeError):
        return False


def ensure(*, cache_dir=None, download=False) -> Path:
    """Acquire only the one pinned ZIP, bounded by its declared byte limits.

    No source discovery, gene downloads or automatic repairs. Corrupt caches fail;
    download=True explicitly permits replacing them with a verified pinned asset.
    """
    from filelock import FileLock

    from . import _downloads

    pin = descriptor()
    root = cache_root(cache_dir)
    destination = root / pin["bundle_id"].rsplit(":", 1)[-1]
    if destination.exists():
        try:
            load(destination)
            return destination
        except (ValueError, OSError):
            if not download:
                raise
    if not download:
        raise FileNotFoundError(
            "Pinned sequence bundle is not cached; use ensure(download=True) explicitly"
        )
    if _downloads.is_offline():
        raise FileNotFoundError(
            "Offline mode is enabled; refusing to download pinned sequence bundle"
        )
    root.mkdir(parents=True, exist_ok=True)
    with FileLock(str(root / (destination.name + ".lock"))):
        import shutil

        import requests

        with tempfile.TemporaryDirectory(prefix=".sequence-", dir=root) as temporary:
            temp = Path(temporary)
            archive = temp / "asset.zip"
            size = 0
            with requests.get(pin["url"], stream=True, timeout=(20, 120)) as response:
                response.raise_for_status()
                with archive.open("wb") as output:
                    for block in response.iter_content(1024 * 1024):
                        size += len(block)
                        if size > pin["archive_bytes"]:
                            raise ValueError("Sequence download exceeds declared bound")
                        output.write(block)
            if (
                size != pin["archive_bytes"]
                or digest(archive.read_bytes()) != pin["archive_sha256"]
            ):
                raise ValueError("Sequence archive checksum mismatch")
            unpacked = temp / "bundle"
            unpacked.mkdir()
            with zipfile.ZipFile(archive) as zipped:
                names = zipped.namelist()
                if set(names) != PAYLOADS | {"manifest.json"} or len(names) != len(set(names)):
                    raise ValueError("Unexpected sequence archive members")
                if sum(i.file_size for i in zipped.infolist()) > pin["unpacked_bytes"]:
                    raise ValueError("Sequence archive exceeds expansion bound")
                for name in names:
                    (unpacked / name).write_bytes(zipped.read(name))
            load(unpacked)
            if destination.exists():
                shutil.rmtree(destination)
            unpacked.rename(destination)
    return destination


def sequences_table(
    bundle: SequenceBundle | None = None,
    *,
    path=None,
    cache_dir=None,
    download: bool = False,
    roles: tuple[str, ...] | None = None,
):
    """Load or convert a verified ``SequenceBundle`` into an occurrence DataFrame."""
    loaded = (
        bundle
        if bundle is not None
        else load(path)
        if path is not None
        else load(ensure(cache_dir=cache_dir, download=download))
    )
    return loaded.to_sequences_table(roles=roles)


def ncorf_te_screening_contract() -> dict:
    """Return the versioned screening contract for non-canonical ORFs and TE/HERV antigens (#736).

    Keeps canonical Ensembl coding proteins (`target_selection`) cleanly separated from
    uORF/dORF/lncRNA-sORF and transposable-element (HERV/LINE-1/LTR) nominations, which
    require empirical translation proof and multicopy-aware normal-atlas screening.
    """
    from . import peptide_subtraction, target_selection

    return {
        "contract_id": NCORF_TE_CONTRACT_ID,
        "catalog_policy_id": target_selection.POLICY_ID,
        "near_self_contract_id": peptide_subtraction.NEAR_SELF_CONTRACT_ID,
        "scope": [
            "uORF",
            "dORF",
            "lncRNA_sORF",
            "intronic_cryptic_ORF",
            "transposable_element_HERV_LINE_LTR",
        ],
        "mandatory_gates": {
            "translation_evidence": (
                "Empirical Ribo-seq triplet periodicity and/or FDR-controlled HLA "
                "immunopeptidomics MS/MS confirmation; in-silico 3-frame translation "
                "alone never qualifies."
            ),
            "multicopy_mappability": (
                "Identical k-mer or locus-subfamily multimapper accounting across all "
                "dispersed TE copies in both tumor and normal RNA-seq; single-copy "
                "unique-mapper filtering on repetitive elements is rejected."
            ),
            "normal_atlas_and_mtec_screen": (
                "Donor-resolved healthy adult bulk and rare-cell single-cell screening "
                "under identical multimapper rules, retaining thymic mTEC promiscuous "
                "transcription caveats separately from peripheral somatic toxicity."
            ),
            "self_and_polymorphic_subtraction": (
                "Exact and bounded near-self subtraction (PeptideIndex + population "
                "CodingVariant augmentation) against the full Unsafe|Other proteome."
            ),
        },
        "canonical_catalog_separation": (
            "ncORF and TE candidates never enter canonical Ensembl protein-group "
            "counts silently; they must carry explicit ORF/repeat coordinates and "
            "pass all four mandatory_gates."
        ),
    }
