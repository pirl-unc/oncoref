"""Run identity compatibility tests against an extracted source distribution."""

import os
import subprocess
import sys
import tarfile
from pathlib import Path
from tempfile import TemporaryDirectory

archive_path = Path(sys.argv[1]).resolve()
with TemporaryDirectory(prefix="oncoref-609-sdist-") as directory:
    with tarfile.open(archive_path) as archive:
        archive.extractall(directory, filter="data")
    root = next(Path(directory).iterdir())
    checkout = Path(__file__).resolve().parents[1]
    tracked = (
        subprocess.check_output(["git", "ls-files", "-z", "tests/"], cwd=checkout)
        .decode()
        .split("\0")
    )
    for name in filter(None, tracked):
        assert (root / name).read_bytes() == (checkout / name).read_bytes(), name
    env = dict(os.environ)
    env.pop("PYTHONPATH", None)
    env["PYENSEMBL_CACHE_DIR"] = str(Path(directory) / "cache")
    runner = """import pathlib,sys
sys.path.insert(0,sys.argv[1])
import oncoref
assert pathlib.Path(oncoref.__file__).resolve().is_relative_to(pathlib.Path(sys.argv[1]).resolve())
import pytest
raise SystemExit(pytest.main(['-q','--import-mode=importlib','tests/test_gene_identity.py','tests/test_dated_gene_identity.py']))
"""
    subprocess.run([sys.executable, "-I", "-c", runner, str(root)], cwd=root, env=env, check=True)
