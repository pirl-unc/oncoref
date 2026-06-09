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

"""On-disk cache location for the downloadable data bundle.

Resolution order for the cache root:
  1. ``CANCERDATA_BUNDLED_DATA`` env var (this package's override);
  2. ``PIRLYGENES_BUNDLED_DATA`` env var (back-compat: the bundle was historically
     fetched by pirlygenes, and existing installs already have it cached there);
  3. ``~/.cache/pirlygenes/bundled_data`` (the historical default; preserved so
     no one is forced to re-download ~340 MB after the data moved here).

The per-version bundle lives under ``<root>/v<DATA_VERSION>/``. The fetch/extract
machinery that populates it is added with the expression bundle in a later
milestone; this module just resolves the path so the CLI and accessors agree on
where the bundle lives.
"""

from __future__ import annotations

import os
from pathlib import Path

from .version import DATA_VERSION

#: Environment variable that overrides the cache root for this package.
CACHE_DIR_ENV_VAR = "CANCERDATA_BUNDLED_DATA"
#: Back-compat env var honored when this package's own override is unset.
LEGACY_CACHE_DIR_ENV_VAR = "PIRLYGENES_BUNDLED_DATA"
#: Historical default cache root (kept so existing caches are reused as-is).
DEFAULT_CACHE_ROOT = Path.home() / ".cache" / "pirlygenes" / "bundled_data"


def cache_root() -> Path:
    """Parent of all version-pinned bundle cache dirs (honors env overrides)."""
    override = os.environ.get(CACHE_DIR_ENV_VAR) or os.environ.get(LEGACY_CACHE_DIR_ENV_VAR)
    return Path(override) if override else DEFAULT_CACHE_ROOT


def bundle_cache_dir() -> Path:
    """Directory the ``v<DATA_VERSION>`` data bundle lives in for this install."""
    return cache_root() / f"v{DATA_VERSION}"
