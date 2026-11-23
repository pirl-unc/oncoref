"""Compatibility shim for frozen legacy module ``oncoref.cta_membership`` (#756).

Implementation lives in :mod:`oncoref._legacy.cta_membership`; active target
selection lives in :mod:`oncoref.target_selection` and :mod:`oncoref.antigen_sets`.
"""

from __future__ import annotations

import importlib as _importlib
import sys as _sys

_sys.modules[__name__] = _importlib.import_module("oncoref._legacy.cta_membership")
