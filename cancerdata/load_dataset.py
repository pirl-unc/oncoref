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

"""Loader for the small bundled reference CSVs.

The cancer-type registry, cohort vocabulary, TMB, and incidence/mortality
tables are tiny and ship inside the wheel (``cancerdata/data/*.csv``). This
module resolves a dataset name to its CSV and returns a cached, defensively
copied DataFrame.

The heavy per-cohort expression bundle (lazy-downloaded from a GitHub Release)
is added in a later milestone; this loader is forward-compatible — callers ask
for a dataset by name and get a DataFrame regardless of where it lives.
"""

from pathlib import Path

import pandas as pd

_DATA_DIR = Path(__file__).parent / "data"
_CACHED_DATAFRAMES: dict[str, pd.DataFrame] = {}


def _clear_cache() -> None:
    """Drop cached frames. Test hook for swapping fixture CSVs via a
    monkey-patched :func:`get_data`."""
    _CACHED_DATAFRAMES.clear()


def get_data(name, _dataframes_dict=None, *, copy=True):
    """Load a bundled dataset as a DataFrame.

    By default returns a defensive ``.copy()`` so callers that mutate in place
    can't corrupt the shared cache. Pass ``copy=False`` for read-only callers
    that slice/copy before mutating.
    """
    stem = name.removesuffix(".csv")
    candidates = [name, name.lower(), f"{stem}.csv", f"{stem.lower()}.csv"]

    if _dataframes_dict is not None:
        for candidate in candidates:
            if candidate in _dataframes_dict:
                df = _dataframes_dict[candidate]
                return df.copy() if copy else df
        raise ValueError(f"Dataset {name} not found")

    if stem not in _CACHED_DATAFRAMES:
        path = _DATA_DIR / f"{stem}.csv"
        if not path.exists():
            raise ValueError(f"Dataset {name} not found at {path}")
        _CACHED_DATAFRAMES[stem] = pd.read_csv(str(path), low_memory=False)
    cached = _CACHED_DATAFRAMES[stem]
    return cached.copy() if copy else cached
