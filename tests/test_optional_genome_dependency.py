# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0

import subprocess
import sys
from pathlib import Path

_PROJECT = Path(__file__).resolve().parents[1]


def test_base_layer_without_pyensembl():
    code = r"""
import builtins

real_import = builtins.__import__

def without_pyensembl(name, *args, **kwargs):
    if name == "pyensembl" or name.startswith("pyensembl."):
        raise ModuleNotFoundError("blocked optional dependency", name="pyensembl")
    return real_import(name, *args, **kwargs)

builtins.__import__ = without_pyensembl

import oncoref

assert oncoref.canonical_gene_id("TP53") == "ENSG00000141510"
assert oncoref.canonical_gene_id("NOT_A_REAL_GENE") is None
assert not oncoref.cancer_tmb_df().empty

try:
    oncoref.genomes()
except ImportError as exc:
    assert "oncoref[genome]" in str(exc)
else:
    raise AssertionError("genome API should require the optional extra")
"""
    result = subprocess.run(
        [sys.executable, "-c", code],
        cwd=_PROJECT,
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode == 0, result.stderr + result.stdout
