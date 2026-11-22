# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0

"""Tests for the controlled/local paired-FASTQ Salmon builder (#665)."""

from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from oncoref import expression_builders, expression_registry

SCRIPTS_DIR = Path(__file__).resolve().parents[1] / "scripts"


def _load_script(stem: str):
    path = SCRIPTS_DIR / f"{stem}.py"
    spec = importlib.util.spec_from_file_location(f"_script_{stem}", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _write_fastq_pair(
    directory: Path,
    stem: str,
    payload_seed: str,
) -> tuple[Path, Path, tuple[str, str], tuple[str, str]]:
    directory.mkdir(parents=True, exist_ok=True)
    r1 = directory / f"{stem}_R1.fastq.gz"
    r2 = directory / f"{stem}_R2.fastq.gz"
    b1 = f"@{stem}/1\nACGTACGT\n+\nIIIIIIII\n#{payload_seed}\n".encode()
    b2 = f"@{stem}/2\nTGCATGCA\n+\nIIIIIIII\n#{payload_seed}\n".encode()
    r1.write_bytes(b1)
    r2.write_bytes(b2)
    sha = (hashlib.sha256(b1).hexdigest(), hashlib.sha256(b2).hexdigest())
    md5 = (hashlib.md5(b1).hexdigest(), hashlib.md5(b2).hexdigest())
    return r1, r2, sha, md5


def _write_quant_dir(
    quant_root: Path,
    sample_ids: list[str],
    transcript_to_gene: dict[str, str],
) -> dict[str, Path]:
    paths: dict[str, Path] = {}
    for index, sample_id in enumerate(sample_ids):
        qdir = quant_root / sample_id
        qdir.mkdir(parents=True, exist_ok=True)
        qpath = qdir / "quant.sf"
        tp53 = 650_000.0 - index * 25_000.0
        pd.DataFrame(
            {
                "Name": list(transcript_to_gene),
                "Length": [1000, 800, 2000],
                "EffectiveLength": [800, 600, 1800],
                "TPM": [tp53 / 2.0, tp53 / 2.0, 1_000_000.0 - tp53],
                "NumReads": [120, 80, 100],
            }
        ).to_csv(qpath, sep="\t", index=False)
        aux = qdir / "aux_info"
        aux.mkdir(exist_ok=True)
        (aux / "meta_info.json").write_text(
            json.dumps(
                {
                    "salmon_version": "1.10.3",
                    "library_types": ["ISR"],
                    "percent_mapped": 94.5 - index,
                    "num_processed": 2000,
                    "num_mapped": 1890 - index * 20,
                }
            )
        )
        paths[sample_id] = qpath
    return paths


def test_unc_nutm1_registry_entry_keeps_legacy_artifact_withheld_while_binding_salmon_builder():
    src = expression_registry.expression_source("unc-nutm1")
    assert src is not None
    assert src.builder == "scripts/build_local_fastq_salmon_source.py"
    assert src.builder_args == ("--source-id", "unc-nutm1")
    assert src.external_build_exemption is None
    assert src.access_level == "controlled"
    assert src.acquisition_status == "local_manifest_required"
    # Legacy v5.22.10 artifact stays explicitly withheld until rebuilt from controlled reads (#604).
    assert src.source_scale_class == "invalid_legacy_tpm"
    assert src.linear_tpm_comparable is False

    entries = expression_builders.local_fastq_salmon_source_entries()
    assert any(entry["id"] == "unc-nutm1" for entry in entries)


def test_unc_nutm1_loads_multilane_yaml_and_csv_manifests_and_strips_private_paths(tmp_path):
    private_dir = tmp_path / "sensitive_clinical_mount" / "patient_reads"
    c1_l1 = _write_fastq_pair(private_dir, "nutm_case1_L001", "c1l1")
    c1_l2 = _write_fastq_pair(private_dir, "nutm_case1_L002", "c1l2")
    c2_l1 = _write_fastq_pair(private_dir, "nutm_case2_L001", "c2l1")
    c3_l1 = _write_fastq_pair(private_dir, "nutm_case3_L001", "c3l1")
    ctrl_l1 = _write_fastq_pair(private_dir, "nutm_ctrl1_L001", "ctrll1")

    manifest_csv = private_dir / "manifest.csv"
    pd.DataFrame(
        [
            {
                "sample_id": "UNC_NUTM_01",
                "case_id": "CASE_01",
                "sample_title": "Primary NUTM tumor 1",
                "role": "tumor",
                "cancer_code": "NUTM",
                "lane_id": "L001",
                "read_1": c1_l1[0].name,
                "read_2": c1_l1[1].name,
                "read_1_sha256": c1_l1[2][0],
                "read_2_sha256": c1_l1[2][1],
            },
            {
                "sample_id": "UNC_NUTM_01",
                "case_id": "CASE_01",
                "sample_title": "Primary NUTM tumor 1",
                "role": "tumor",
                "cancer_code": "NUTM",
                "lane_id": "L002",
                "read_1": c1_l2[0].name,
                "read_2": c1_l2[1].name,
                "read_1_sha256": c1_l2[2][0],
                "read_2_sha256": c1_l2[2][1],
            },
            {
                "sample_id": "UNC_NUTM_02",
                "case_id": "CASE_02",
                "sample_title": "Primary NUTM tumor 2",
                "role": "tumor",
                "cancer_code": "NUTM",
                "lane_id": "L001",
                "read_1": str(c2_l1[0]),
                "read_2": str(c2_l1[1]),
                "read_1_sha256": c2_l1[2][0],
                "read_2_sha256": c2_l1[2][1],
            },
            {
                "sample_id": "UNC_NUTM_03",
                "case_id": "CASE_03",
                "sample_title": "Primary NUTM tumor 3",
                "role": "tumor",
                "cancer_code": "NUTM",
                "lane_id": "L001",
                "read_1": str(c3_l1[0]),
                "read_2": str(c3_l1[1]),
                "read_1_sha256": c3_l1[2][0],
                "read_2_sha256": c3_l1[2][1],
            },
            {
                "sample_id": "UNC_NUTM_01_ADJ",
                "case_id": "CASE_01",
                "sample_title": "Matched normal 1",
                "role": "matched_normal",
                "cancer_code": "",
                "lane_id": "L001",
                "read_1": str(ctrl_l1[0]),
                "read_2": str(ctrl_l1[1]),
                "read_1_md5": ctrl_l1[3][0],
                "read_2_md5": ctrl_l1[3][1],
            },
        ]
    ).to_csv(manifest_csv, index=False)

    source = expression_builders.local_fastq_salmon_source_from_registry(
        "unc-nutm1",
        manifest=manifest_csv,
    )
    assert source.source_scale_class == "salmon_gene_tpm"
    assert source.linear_tpm_comparable is True
    assert source.tpm_proxy is False
    assert len(source.samples) == 4
    assert len(source.samples[0].lanes) == 2

    transcript_to_gene = {
        "ENST_TP53_A": "ENSG00000141510.17",
        "ENST_TP53_B": "ENSG00000141510.18",
        "ENST_NUTM1": "ENSG00000184507.14",
    }
    quant_paths = _write_quant_dir(
        tmp_path / "quant",
        [s.sample_id for s in source.samples],
        transcript_to_gene,
    )

    result = expression_builders.build_local_fastq_salmon_source_matrices(
        source,
        cache_dir=tmp_path / "cache",
        output_dir=tmp_path / "derived",
        quant_paths=quant_paths,
        transcript_to_gene=transcript_to_gene,
        verify_reads=True,
    )

    matrix = result.matrices["NUTM"]
    assert list(matrix.columns) == [
        "Ensembl_Gene_ID",
        "Symbol",
        "UNC_NUTM_01",
        "UNC_NUTM_02",
        "UNC_NUTM_03",
    ]
    assert np.allclose(matrix.iloc[:, 2:].sum(axis=0), 1_000_000.0)
    assert set(result.sample_qc["source_type"]) == {"local-fastq-salmon"}
    assert set(result.sample_qc["source_scale_class"]) == {"salmon_gene_tpm"}
    assert result.sample_qc["linear_tpm_comparable"].all()

    manifest_out = result.sidecar_paths["run_manifest"].read_text()
    audit_out = result.sidecar_paths["salmon_quantification_audit"].read_text()
    assert "sensitive_clinical_mount" not in manifest_out
    assert "patient_reads" not in manifest_out
    assert "fastq.gz" not in manifest_out
    assert "sensitive_clinical_mount" not in audit_out

    manifest_df = pd.read_csv(result.sidecar_paths["run_manifest"])
    assert manifest_df["included"].tolist() == [True, True, True, False]
    assert manifest_df["n_lanes"].tolist() == [2, 1, 1, 1]
    assert manifest_df.loc[3, "exclusion_reason"] == (
        "matched_normal_excluded_from_tumor_reference"
    )


def test_quantify_local_fastq_sample_passes_multilane_mates_and_invalidates_on_digest_change(
    tmp_path, monkeypatch
):
    reads_dir = tmp_path / "reads"
    l1 = _write_fastq_pair(reads_dir, "s1_L1", "seed1")
    l2 = _write_fastq_pair(reads_dir, "s1_L2", "seed2")
    s2 = _write_fastq_pair(reads_dir, "s2_L1", "seed3")
    s3 = _write_fastq_pair(reads_dir, "s3_L1", "seed4")

    manifest_yaml = tmp_path / "manifest.yaml"
    manifest_yaml.write_text(
        json.dumps(
            {
                "samples": [
                    {
                        "sample_id": "NUTM_S1",
                        "case_id": "C1",
                        "sample_title": "Tumor 1",
                        "role": "tumor",
                        "cancer_code": "NUTM",
                        "lanes": [
                            {
                                "lane_id": "L001",
                                "read_1": str(l1[0]),
                                "read_2": str(l1[1]),
                                "sha256": list(l1[2]),
                            },
                            {
                                "lane_id": "L002",
                                "read_1": str(l2[0]),
                                "read_2": str(l2[1]),
                                "sha256": list(l2[2]),
                            },
                        ],
                    },
                    {
                        "sample_id": "NUTM_S2",
                        "case_id": "C2",
                        "sample_title": "Tumor 2",
                        "role": "tumor",
                        "cancer_code": "NUTM",
                        "lanes": [
                            {
                                "lane_id": "L001",
                                "read_1": str(s2[0]),
                                "read_2": str(s2[1]),
                                "sha256": list(s2[2]),
                            }
                        ],
                    },
                    {
                        "sample_id": "NUTM_S3",
                        "case_id": "C3",
                        "sample_title": "Tumor 3",
                        "role": "tumor",
                        "cancer_code": "NUTM",
                        "lanes": [
                            {
                                "lane_id": "L001",
                                "read_1": str(s3[0]),
                                "read_2": str(s3[1]),
                                "sha256": list(s3[2]),
                            }
                        ],
                    },
                ]
            }
        )
    )
    source = expression_builders.local_fastq_salmon_source_from_registry(
        "unc-nutm1",
        manifest=manifest_yaml,
    )

    calls: list[list[str]] = []
    monkeypatch.setattr(
        expression_builders, "_salmon_executable", lambda cmd: "/usr/local/bin/salmon"
    )
    monkeypatch.setattr(expression_builders, "_salmon_version", lambda exe: "1.10.3")

    def _fake_run(cmd, check):
        assert check is True
        calls.append(list(cmd))
        out_idx = cmd.index("--output") + 1
        out_dir = Path(cmd[out_idx])
        out_dir.mkdir(parents=True, exist_ok=True)
        (out_dir / "quant.sf").write_text("Name\tLength\tEffectiveLength\tTPM\tNumReads\n")

    monkeypatch.setattr(expression_builders.subprocess, "run", _fake_run)

    qpath = expression_builders.quantify_local_fastq_salmon_sample(
        source,
        source.samples[0],
        cache_dir=tmp_path / "cache",
        index_dir=tmp_path / "index",
    )
    assert qpath.exists()
    assert len(calls) == 1
    cmd = calls[0]
    m1_idx = cmd.index("--mates1")
    m2_idx = cmd.index("--mates2")
    threads_idx = cmd.index("--threads")
    assert cmd[m1_idx + 1 : m2_idx] == [str(l1[0]), str(l2[0])]
    assert cmd[m2_idx + 1 : threads_idx] == [str(l1[1]), str(l2[1])]

    # Cache hit does not invoke Salmon again, and cached input JSON omits file paths.
    inputs_json = (qpath.parent / "oncoref_quantification_inputs.json").read_text()
    assert str(reads_dir) not in inputs_json
    expression_builders.quantify_local_fastq_salmon_sample(
        source,
        source.samples[0],
        cache_dir=tmp_path / "cache",
        index_dir=tmp_path / "index",
    )
    assert len(calls) == 1

    # Salmon version upgrade invalidates the cache and re-quantifies.
    monkeypatch.setattr(expression_builders, "_salmon_version", lambda exe: "1.11.0")
    expression_builders.quantify_local_fastq_salmon_sample(
        source,
        source.samples[0],
        cache_dir=tmp_path / "cache",
        index_dir=tmp_path / "index",
    )
    assert len(calls) == 2

    # Corrupting a FASTQ byte on disk fails closed before reusing the cache.
    l1[0].write_bytes(b"corrupted fastq content")
    with pytest.raises(RuntimeError, match="sha256 mismatch"):
        expression_builders.quantify_local_fastq_salmon_sample(
            source,
            source.samples[0],
            cache_dir=tmp_path / "cache",
            index_dir=tmp_path / "index",
        )


def test_local_fastq_manifest_rejects_wrong_sample_count_and_reused_reads(tmp_path):
    r1 = _write_fastq_pair(tmp_path, "s1", "a")
    r2 = _write_fastq_pair(tmp_path, "s2", "b")

    # Only 2 NUTM tumor samples when registry expects 3 -> fail closed.
    two_sample_manifest = [
        {
            "sample_id": "S1",
            "role": "tumor",
            "cancer_code": "NUTM",
            "read_1": str(r1[0]),
            "read_2": str(r1[1]),
            "read_1_sha256": r1[2][0],
            "read_2_sha256": r1[2][1],
        },
        {
            "sample_id": "S2",
            "role": "tumor",
            "cancer_code": "NUTM",
            "read_1": str(r2[0]),
            "read_2": str(r2[1]),
            "read_1_sha256": r2[2][0],
            "read_2_sha256": r2[2][1],
        },
    ]
    with pytest.raises(ValueError, match="expected 3"):
        expression_builders.local_fastq_salmon_source_from_registry(
            "unc-nutm1",
            manifest=two_sample_manifest,
        )

    # Reusing the same physical FASTQ across samples -> fail closed.
    reused_manifest = [
        *two_sample_manifest,
        {
            "sample_id": "S3",
            "role": "tumor",
            "cancer_code": "NUTM",
            "read_1": str(r2[0]),
            "read_2": str(r2[1]),
            "read_1_sha256": r2[2][0],
            "read_2_sha256": r2[2][1],
        },
    ]
    with pytest.raises(ValueError, match="reuses read file"):
        expression_builders.local_fastq_salmon_source_from_registry(
            "unc-nutm1",
            manifest=reused_manifest,
        )


def test_build_local_fastq_salmon_script_cli(tmp_path, monkeypatch, capsys):
    mod = _load_script("build_local_fastq_salmon_source")
    r1 = _write_fastq_pair(tmp_path, "s1", "a")
    r2 = _write_fastq_pair(tmp_path, "s2", "b")
    r3 = _write_fastq_pair(tmp_path, "s3", "c")
    manifest = tmp_path / "manifest.json"
    manifest.write_text(
        json.dumps(
            [
                {
                    "sample_id": f"S{idx}",
                    "role": "tumor",
                    "cancer_code": "NUTM",
                    "read_1": str(pair[0]),
                    "read_2": str(pair[1]),
                    "read_1_sha256": pair[2][0],
                    "read_2_sha256": pair[2][1],
                }
                for idx, pair in enumerate((r1, r2, r3), start=1)
            ]
        )
    )

    with pytest.raises(SystemExit) as err:
        mod.main(
            [
                "--source-id",
                "unc-nutm1",
                "--manifest",
                str(manifest),
                "--quant-dir",
                str(tmp_path / "q"),
                "--force-quant",
            ]
        )
    assert err.value.code == 2
    assert "--quant-dir cannot be combined" in capsys.readouterr().err

    def _fake_build(source_obj, *, cache_dir, quant_paths=None, **_kwargs):
        assert source_obj.source_id == "unc-nutm1"
        assert quant_paths == {
            f"S{idx}": tmp_path / "q" / f"S{idx}" / "quant.sf" for idx in (1, 2, 3)
        }
        matrix = pd.DataFrame(
            {
                "Ensembl_Gene_ID": ["ENSG00000184507"],
                "Symbol": ["NUTM1"],
                "S1": [1e6],
                "S2": [1e6],
                "S3": [1e6],
            }
        )
        return expression_builders.SourceMatrixBuildResult(
            source=source_obj,
            matrices={"NUTM": matrix},
            matrix_paths={"NUTM": tmp_path / "NUTM_per_sample_tpm.parquet"},
            summary_rows=pd.DataFrame(
                columns=list(expression_builders.REFERENCE_EXPRESSION_COLUMNS)
            ),
            mapping_audit=pd.DataFrame(),
            parse_diagnostics=pd.DataFrame(),
            sample_qc=pd.DataFrame(),
            sidecar_paths={"run_manifest": tmp_path / "run_manifest.csv"},
        )

    monkeypatch.setattr(mod, "build_local_fastq_salmon_source_matrices", _fake_build)
    assert (
        mod.main(
            [
                "--source-id",
                "unc-nutm1",
                "--manifest",
                str(manifest),
                "--cache-dir",
                str(tmp_path / "cache"),
                "--quant-dir",
                str(tmp_path / "q"),
            ]
        )
        == 0
    )
    out = capsys.readouterr().out
    assert '"source_id": "unc-nutm1"' in out
    assert '"NUTM": 3' in out
