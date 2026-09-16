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

"""One publication style for every oncoref figure.

Two entry points, applied at the two places a figure can pick up styling:

``apply()``
    Global rcParams — opaque white figure/axes/savefig faces (so a PNG dropped on
    a white manuscript page has no grey plate and no transparent halo), hairline
    spines, and a single serif-free type scale.

``finalize(fig)``
    Per-figure cleanup at save time: drop the in-figure titles, strip the top and
    right spines, and thin the tick marks. Titles, captions and callouts belong in
    the manuscript's caption text, not burned into the raster, so the default
    ``minimal`` mode removes them rather than asking every plot function to stop
    setting them (they stay useful for interactive/notebook use, where the figure
    travels without a caption).

Both are idempotent and safe to call repeatedly. ``set_minimal(False)`` restores
the chatty labelling for interactive work.
"""

from __future__ import annotations

#: Whether :func:`finalize` strips in-figure titles and other redundant chrome.
_MINIMAL = True

_APPLIED = False

#: Categorical palette (colour-blind safe), shared by the curation and provenance
#: figure families so a multi-panel display reads as one system.
PALETTE = (
    "#2a7f4f",  # green   - kept / present
    "#b0b0b0",  # grey    - dropped / absent
    "#f0c419",  # amber   - weak / caveated
    "#3b6ea5",  # blue    - secondary series
    "#c0392b",  # red     - thresholds
    "#7b5aa6",  # purple  - tertiary series
)

KEPT = PALETTE[0]
DROP = PALETTE[1]
WEAK = PALETTE[2]
ACCENT = PALETTE[3]
THRESHOLD = PALETTE[4]

_RC = {
    # Opaque white everywhere, including the saved raster.
    "figure.facecolor": "white",
    "figure.edgecolor": "white",
    "axes.facecolor": "white",
    "savefig.facecolor": "white",
    "savefig.edgecolor": "white",
    "savefig.transparent": False,
    "savefig.dpi": 300,
    "figure.dpi": 150,
    # Hairline frame; finalize() removes the top/right pair per figure.
    "axes.linewidth": 0.8,
    "axes.edgecolor": "#333333",
    "axes.grid": False,
    "axes.titlesize": 11,
    "axes.labelsize": 10,
    "axes.labelcolor": "#222222",
    "xtick.color": "#333333",
    "ytick.color": "#333333",
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "xtick.direction": "out",
    "ytick.direction": "out",
    "xtick.major.width": 0.8,
    "ytick.major.width": 0.8,
    "xtick.major.size": 3.0,
    "ytick.major.size": 3.0,
    "font.size": 10,
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "legend.frameon": False,
    "legend.fontsize": 9,
    "legend.handlelength": 1.2,
    "legend.borderaxespad": 0.2,
    "lines.linewidth": 1.5,
    "patch.linewidth": 0.0,
    # Keep text as text in vector exports rather than outlines.
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
}


def set_minimal(minimal: bool) -> None:
    """Turn the label-stripping half of the style on or off.

    ``apply()`` still governs colours and type; this only controls whether
    :func:`finalize` removes titles and redundant chrome.
    """
    global _MINIMAL
    _MINIMAL = bool(minimal)


def minimal() -> bool:
    """Whether :func:`finalize` is currently stripping in-figure labels."""
    return _MINIMAL


def apply(force: bool = False) -> None:
    """Install the publication rcParams (once per process unless ``force``)."""
    global _APPLIED
    if _APPLIED and not force:
        return
    import matplotlib

    matplotlib.use("Agg", force=False)
    matplotlib.rcParams.update(_RC)
    _APPLIED = True


def finalize(fig, *, keep_titles: bool | None = None):
    """Apply the per-figure half of the style and return ``fig``.

    Strips the top/right spines on every axes, forces opaque white faces (a figure
    built before :func:`apply` ran still lands on white), and — in minimal mode —
    removes axes titles and the suptitle.
    """
    keep = (not _MINIMAL) if keep_titles is None else keep_titles
    fig.patch.set_facecolor("white")
    fig.patch.set_alpha(1.0)
    for ax in fig.get_axes():
        ax.set_facecolor("white")
        for side in ("top", "right"):
            if side in ax.spines:
                ax.spines[side].set_visible(False)
        for side in ("left", "bottom"):
            if side in ax.spines:
                ax.spines[side].set_linewidth(0.8)
                ax.spines[side].set_color("#333333")
        if not keep:
            ax.set_title("")
    # ``_suptitle`` is the only handle matplotlib exposes for clearing one.
    if not keep and getattr(fig, "_suptitle", None) is not None:
        fig.suptitle("")
    return fig


def save(fig, path, *, dpi=300, keep_titles: bool | None = None):
    """Finalize and write ``fig`` to ``path`` on an opaque white ground."""
    finalize(fig, keep_titles=keep_titles)
    fig.savefig(path, dpi=dpi, bbox_inches="tight", facecolor="white", transparent=False)
    return fig


#: Largest dimension, in inches, any figure may reach. A figure sized by row count
#: grows without limit otherwise: 85 cancer types at 0.42 in/row is a 36-inch page,
#: which is not a figure but a table, and at 300 dpi it is a 10,000px raster that no
#: layout can place. matplotlib itself hard-fails past 65,536px a side.
MAX_FIGURE_INCHES = 20.0

#: Never shrink a tick label below this, however dense a clamped figure gets.
MIN_TICK_FONTSIZE = 4.0


def stack_size(n, *, per_item, floor, cap=MAX_FIGURE_INCHES):
    """Inches for ``n`` stacked rows (or columns), clamped to a printable maximum.

    Returns ``(inches, density)`` where ``density <= 1`` is how far the per-item
    spacing had to compress to fit ``cap``. Callers scale their tick-label font by
    it, so a clamped figure gets denser rather than overprinting itself.

    A clamped figure is a signal, not a solution: when ``density`` is small the
    honest fix is to filter the rows (``min_regimens``, ``n=``, ``only_multi``),
    because no font size makes 85 rows on one page a readable figure.
    """
    want = max(float(floor), float(per_item) * int(n))
    inches = min(want, float(cap))
    density = (inches / want) if want > 0 else 1.0
    return inches, density


def tick_fontsize(density, base=8.0):
    """Tick-label size for a figure clamped to ``density`` (see :func:`stack_size`)."""
    return max(MIN_TICK_FONTSIZE, base * density)
