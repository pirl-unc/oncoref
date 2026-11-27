# Licensed under the Apache License, Version 2.0.
"""Dog species namespace alias for :mod:`oncoref.canine` (#791).

Exposes ``oncoref.dog`` (and ``import oncoref.dog``) with full symmetry alongside
``oncoref.human``, ``oncoref.mouse``, and ``oncoref.canine``.
"""

import sys as _sys

from . import canine as _canine
from .canine import *  # noqa: F403
from .canine import __all__  # noqa: F401

_sys.modules[__name__] = _canine
