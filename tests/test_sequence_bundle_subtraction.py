# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.

from __future__ import annotations

import builtins
import gzip
import json
from pathlib import Path

import pandas as pd
import pytest

from oncoref import antigen_sets, cli, host_context, sequence_bundle, target_selection
from oncoref.peptide_subtraction import MatchPolicy, PeptideIndex, near_self_contract


def _build_synthetic_bundle(folder: Path) -> sequence_bundle.SequenceBundle:
    folder.mkdir(parents=True, exist_ok=True)
    raw_records = [
        (
            "ENSP_UNSAFE.1",
            "ACDEFGHIKLMNPQRST",
            "exclude",
            "ENSG_UNSAFE",
            "ENST_UNSAFE.2",
            False,
        ),
        (
            "ENSP_OTHER.1",
            "VWYAACDEFGHIKWWWW",
            "review",
            "ENSG_OTHER",
            "ENST_OTHER.1",
            False,
        ),
        (
            "ENSP_CTA.1",
            "YYYYYLMNPQRSTYYYY",
            "allow",
            "ENSG_CTA",
            "ENST_CTA.3",
            True,
        ),
    ]
    records = []
    for rid, seq, role, gid, tid, is_cta in raw_records:
        records.append(
            {
                "record_id": rid,
                "reference_id": "synthetic:ref:115",
                "sequence": seq,
                "sequence_sha256": sequence_bundle.digest(seq.encode()),
                "valid_segments": [
                    {
                        "start": 0,
                        "end": len(seq),
                        "sequence": seq,
                        "sha256": sequence_bundle.digest(seq.encode()),
                    }
                ],
                "shared_sequence": role,
                "direct_target_candidate": is_cta,
                "cancer_antigen_candidate": is_cta,
                "reproductive_antigen_candidate": is_cta,
                "occurrences": [
                    {
                        "occurrence_id": f"{rid}:0",
                        "gene_id": gid,
                        "transcript_id": tid,
                        "transcript_biotype": "protein_coding",
                        "shared_sequence": role,
                        "direct_target_candidate": is_cta,
                        "cancer_antigen_candidate": is_cta,
                        "reproductive_antigen_candidate": is_cta,
                    }
                ],
            }
        )
    files = {
        "records.jsonl.gz": gzip.compress(
            ("\n".join(json.dumps(r) for r in records) + "\n").encode(), mtime=0
        ),
        "coverage.json": json.dumps(
            {"reference_id": "synthetic:ref:115", "missing_translations": []}
        ).encode(),
        "decisions.csv.gz": gzip.compress(b"gene,role\n", mtime=0),
        "orfdata-proteome.json.gz": gzip.compress(b"{}", mtime=0),
    }
    manifest = {
        "schema": sequence_bundle.SCHEMA,
        "schema_version": 1,
        "reference": {"reference_id": "synthetic:ref:115", "taxon_id": 9606},
        "selection_policy_id": target_selection.POLICY_ID,
        "role_policy_id": sequence_bundle.ROLE_POLICY_ID,
        "record_count": len(records),
        "occurrence_count": len(records),
        "package_dependencies": {},
        "payloads": {
            n: {"bytes": len(v), "sha256": sequence_bundle.digest(v)} for n, v in files.items()
        },
    }
    manifest["bundle_id"] = sequence_bundle.content_id(manifest)
    for name, value in files.items():
        (folder / name).write_bytes(value)
    (folder / "manifest.json").write_text(json.dumps(manifest))
    return sequence_bundle.load(folder, expected_id=manifest["bundle_id"], require_current=False)


def test_sequence_bundle_to_sequences_table_and_role_index(tmp_path):
    syn = _build_synthetic_bundle(tmp_path / "syn")
    full_table = syn.to_sequences_table()
    assert list(full_table["gene_id"]) == ["ENSG_UNSAFE", "ENSG_OTHER", "ENSG_CTA"]
    assert list(full_table["transcript_id"]) == ["ENST_UNSAFE", "ENST_OTHER", "ENST_CTA"]
    assert full_table.attrs["bundle_id"] == syn.identity
    assert full_table.attrs["role_policy_id"] == sequence_bundle.ROLE_POLICY_ID

    exclude_only = PeptideIndex.from_sequence_bundle(syn, roles=("exclude",), lengths=[9])
    assert (
        exclude_only.matches("LMNPQRSTY").empty is False
        or not exclude_only.matches("ACDEFGHIK").empty
    )
    assert not exclude_only.matches("ACDEFGHIK").empty
    assert exclude_only.matches("VWYAACDEF").empty

    conservative = PeptideIndex.from_sequence_bundle(
        syn, roles=("exclude", "review"), lengths=[9], policy=MatchPolicy(max_distance=1)
    )
    assert not conservative.matches("VWYAACDEG").empty
    assert conservative.subtract("YYYYYWWWWZZZZZ", lengths=[9]).peptide.tolist()


def test_near_self_contract_and_cross_species_guardrails(tmp_path):
    contract = near_self_contract()
    assert contract["catalog_policy_id"] == target_selection.POLICY_ID
    assert contract["sequence_role_policy_id"] == sequence_bundle.ROLE_POLICY_ID
    assert contract["default_background"] == "Unsafe|Other"

    syn = _build_synthetic_bundle(tmp_path / "syn_species")
    human_table = syn.to_sequences_table()
    with pytest.raises(ValueError, match="taxon mismatch"):
        antigen_sets.forbidden_background(
            "Unsafe|Other", species="canine", sequences_table=human_table
        )
    with pytest.raises(ValueError, match="taxon 9606"):
        antigen_sets.sequences("Other", species="canine", sequence_source="sequence_bundle")


def test_cached_bundle_bridges_antigen_sets_and_peptide_index_without_pyensembl(monkeypatch):
    try:
        sequence_bundle.ensure(download=False)
    except FileNotFoundError:
        pytest.skip("Current sequence bundle is not staged in local cache")

    real_import = builtins.__import__

    def _block_pyensembl(name, *args, **kwargs):
        if name == "pyensembl" or name.startswith("pyensembl."):
            raise ModuleNotFoundError("pyensembl intentionally blocked in offline bundle test")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", _block_pyensembl)

    cta_seqs = antigen_sets.sequences("CTA", sequence_source="sequence_bundle")
    assert not cta_seqs.empty
    assert cta_seqs.attrs["sequence_source"] == "sequence_bundle"
    assert cta_seqs.attrs["missing_gene_ids"] == []

    # Regression (#640): Unsafe vs Other stay distinct; Other is never promoted to Unsafe,
    # and Unscreened/Review entities are explicitly disclosed in attrs.
    unsafe_bg = antigen_sets.forbidden_background("Unsafe", sequence_source="sequence_bundle")
    both_bg = antigen_sets.forbidden_background("Unsafe|Other", sequence_source="sequence_bundle")
    assert len(both_bg) > len(unsafe_bg)
    assert unsafe_bg.attrs["omitted_categories"] == ["Other"]
    assert both_bg.attrs["omitted_categories"] == []
    assert both_bg.attrs["cta_exemption_applied"] is False
    assert len(both_bg.attrs["unmeasured_home_entity_ids"]) == 89
    assert len(both_bg.attrs["review_hold_entity_ids"]) > 700

    # Protecting testis moves testis-expressed CTAs (e.g. MAGEA4 / CTAG1B) into the background!
    male_bg = antigen_sets.forbidden_background(
        "Unsafe",
        host=host_context.HostContext(protect=("testis",)),
        sequence_source="sequence_bundle",
    )
    assert len(male_bg.attrs["host_excluded_entity_ids"]) > 0
    assert (male_bg["background_reason"].str.startswith("Host:")).any()

    # CLI subtraction works with --sequence-source sequence_bundle without pyensembl.
    sample_cta_seq = cta_seqs.iloc[0]["sequence"][:18]
    index = PeptideIndex.from_sequence_bundle(background="Unsafe", lengths=[9])
    sub = index.subtract(sample_cta_seq, lengths=[9])
    assert isinstance(sub, pd.DataFrame)
    assert (
        cli.main(
            [
                "antigens",
                "subtract",
                "WWWWWWWWW",
                "--background",
                "Unsafe",
                "--sequence-source",
                "sequence_bundle",
                "--min-length",
                "9",
                "--max-length",
                "9",
                "--format",
                "json",
            ]
        )
        == 0
    )
