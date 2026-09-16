# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Where figures land: ``figures/run_<YYYYMMDD-HHMMSS>/``.

One convention for every figure oncoref writes, whether it came from a single
``oncoref plot`` invocation or the whole-batch driver. A run never overwrites an
older one, so a figure in a paper draft can always be traced back to the run that
produced it, and ``figures/latest`` points at the most recent run.

    figures/
      run_20260915-104233/
        cta-curation/...
        expression-provenance/...
      latest -> run_20260915-104233

``resolve()`` is the single entry point. Callers that pass an explicit path get
exactly that path (no timestamping, no surprises); callers that pass nothing get
a fresh run directory.
"""

from __future__ import annotations

import os
from datetime import datetime
from pathlib import Path

#: Env var overriding the ``figures/`` root (useful in tests and CI).
FIGURES_DIR_ENV_VAR = "ONCOREF_FIGURES_DIR"

#: Directory name for the symlink to the newest run.
LATEST = "latest"

_RUN_PREFIX = "run_"
_STAMP_FORMAT = "%Y%m%d-%H%M%S"


def figures_root(base=None) -> Path:
    """The ``figures/`` root: explicit ``base``, the env override, or ``./figures``."""
    if base is not None:
        return Path(base)
    env = os.environ.get(FIGURES_DIR_ENV_VAR)
    return Path(env) if env else Path("figures")


def new_run_dir(base=None, *, stamp=None, link_latest=True) -> Path:
    """Create and return a fresh ``figures/run_<stamp>/`` directory.

    A second run inside the same second gets a ``-2``, ``-3``, ... suffix rather
    than silently sharing a directory with the first.
    """
    root = figures_root(base)
    stamp = stamp or datetime.now().strftime(_STAMP_FORMAT)
    run = root / f"{_RUN_PREFIX}{stamp}"
    n = 2
    while run.exists():
        run = root / f"{_RUN_PREFIX}{stamp}-{n}"
        n += 1
    run.mkdir(parents=True)
    if link_latest:
        _link_latest(root, run)
    return run


def _link_latest(root: Path, run: Path) -> None:
    """Point ``figures/latest`` at ``run``, tolerating filesystems without symlinks."""
    link = root / LATEST
    try:
        if link.is_symlink() or link.exists():
            link.unlink()
        link.symlink_to(run.name)
    except OSError:
        # Windows without developer mode, or a filesystem that forbids symlinks.
        pass


def resolve(out=None, *, name=None, suffix=".png", base=None, stamp=None) -> Path:
    """Resolve a figure destination.

    ``out`` given
        Used verbatim — the caller decides. Parent directories are created.
    ``out`` omitted
        A fresh run directory. With ``name``, the path becomes
        ``figures/run_<stamp>/<name><suffix>`` for a single figure, or
        ``figures/run_<stamp>/<name>/`` when ``suffix`` is empty (a figure family
        that writes several files).
    """
    if out is not None:
        path = Path(out)
        parent = path.parent if path.suffix else path
        parent.mkdir(parents=True, exist_ok=True)
        return path
    run = new_run_dir(base, stamp=stamp)
    if name is None:
        return run
    if suffix:
        return run / f"{name}{suffix}"
    target = run / name
    target.mkdir(parents=True, exist_ok=True)
    return target
