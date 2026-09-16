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

#: Active preset name. ``"print"`` is the dense manuscript figure; ``"slide"`` is the
#: same data cut down to what reads from the back of a room.
_PRESET = "print"

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
    # An outside legend overflows a plain tight_layout, so always crop to the real
    # ink bounds when saving rather than to the nominal figure rectangle.
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.08,
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
    """Install the figure rcParams for the active preset (once per process unless
    ``force``). ``use()`` forces a re-apply when the preset changes."""
    global _APPLIED
    if _APPLIED and not force:
        return
    import matplotlib

    matplotlib.use("Agg", force=False)
    matplotlib.rcParams.update(_RC)
    matplotlib.rcParams.update(_preset_rc())
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
    params = PRESETS[_PRESET]
    if cap is MAX_FIGURE_INCHES:  # caller took the default: honour the preset's page
        cap = params["max_inches"]
    per_item = row_height(float(per_item))
    want = max(float(floor), per_item * int(n))
    inches = min(want, float(cap))
    density = (inches / want) if want > 0 else 1.0
    return inches, density


def tick_fontsize(density, base=8.0):
    """Tick-label size for a figure clamped to ``density`` (see :func:`stack_size`)."""
    return max(MIN_TICK_FONTSIZE, base * density)


# ---------------------------------------------------------------- presets ----
#
# A slide is not a small manuscript page. A figure in a paper is read at ~30cm by
# one person who can stop and study it; a figure on a slide is read at 5m by a
# room that gets it for thirty seconds. Bigger type alone does not bridge that —
# the figure also has to carry fewer things. So a preset sets three knobs
# together: type scale, page geometry, and how many rows/points survive.
#
# Presets never change what the data says. A slide figure shows the top slice and
# says so on the axis; it does not re-rank, re-scale, or quietly drop outliers.

#: Per-preset knobs. ``font_scale`` multiplies every text size (and, through it, row
#: spacing and marker area); ``max_inches`` caps any figure dimension; ``max_items``
#: is the row/point budget a plot trims to (``None`` = no trimming).
PRESETS = {
    "print": {"font_scale": 1.0, "max_inches": 20.0, "max_items": None},
    # 1.7 is the scale pirlygenes settled on for projected figures.
    "slide": {"font_scale": 1.7, "max_inches": 12.0, "max_items": 12},
}


def font_scale() -> float:
    """The active text-size multiplier."""
    return PRESETS[_PRESET]["font_scale"]


def fs(pt) -> float:
    """Scale one explicit point size for the active preset.

    Pure multiplication, deliberately: the plot functions already encode a size
    hierarchy (a 4pt label inside a dense grid, a 9pt row label on a bar chart)
    and that hierarchy *is* information about how crowded each panel is. Lifting
    everything to a common floor instead would flatten it and turn the dense
    panels into a pile of overlapping text.
    """
    return round(font_scale() * float(pt), 1)


def _preset_rc() -> dict:
    """rcParams derived from the active preset's font scale."""
    s = font_scale()
    if s == 1.0:
        return {}
    return {
        "font.size": fs(10),
        "axes.titlesize": fs(12),
        "axes.labelsize": fs(11),
        "xtick.labelsize": fs(9),
        "ytick.labelsize": fs(9),
        "legend.fontsize": fs(9),
        "legend.title_fontsize": fs(9),
        "figure.titlesize": fs(13),
        # Bigger type needs proportionally more room around the axes.
        "axes.labelpad": 4.0 * s,
        "axes.titlepad": 6.0 * s,
        "axes.linewidth": 0.8 * min(s, 1.8),
        "xtick.major.width": 0.8 * min(s, 1.8),
        "ytick.major.width": 0.8 * min(s, 1.8),
        "xtick.major.size": 3.0 * min(s, 1.8),
        "ytick.major.size": 3.0 * min(s, 1.8),
        "lines.linewidth": 1.5 * min(s, 2.0),
    }


#: Slide figures cap here rather than at MAX_FIGURE_INCHES.
SLIDE_ASPECT = (12.0, 6.75)


def preset() -> str:
    """The active preset name."""
    return _PRESET


def use(name: str) -> None:
    """Select a preset (``"print"`` or ``"slide"``) and re-apply the style."""
    global _PRESET
    if name not in PRESETS:
        raise ValueError(f"unknown preset {name!r}; expected one of {sorted(PRESETS)}")
    _PRESET = name
    apply(force=True)


def _params() -> dict:
    return PRESETS[_PRESET]


def max_items():
    """Row/point budget for the active preset, or ``None`` when unlimited.

    A plot with more rows than this trims to the top slice and says so in its axis
    label, rather than shrinking type until the figure is unreadable on a wall.
    """
    return _params()["max_items"]


def trim(items, budget=None):
    """Return ``(kept, n_dropped)`` for a ranked sequence under the item budget.

    ``items`` must already be in the order that matters (best first, or the order
    the figure will draw). Trimming keeps the head, so a chart sorted by magnitude
    keeps its largest rows.
    """
    budget = max_items() if budget is None else budget
    items = list(items)
    if budget is None or len(items) <= budget:
        return items, 0
    return items[:budget], len(items) - budget


def marker_size() -> int:
    """Scatter marker *area* for the active preset.

    Area scales with the square of the text scale so a marker keeps the same
    visual weight next to its label rather than shrinking beside it.
    """
    return int(70 * font_scale() ** 2)


def row_height(base: float = 0.30) -> float:
    """Per-row inches for a horizontal bar chart.

    Grows with the text scale, but sub-linearly: taller labels need more room,
    while a strictly proportional growth would make a slide bar chart taller than
    the slide. Mirrors pirlygenes' row_height.
    """
    return base * (1.0 + 0.72 * (font_scale() - 1.0))


def figure_size(width, height, *, wide=False):
    """Clamp a ``(width, height)`` to the active preset's page.

    ``wide`` asks for the full 16:9 slide rather than a clamp, which is what a
    scatter wants: rendering at the slide's own aspect means a point size here
    lands as roughly that point size on a full-bleed slide.
    """
    cap = _params()["max_inches"]
    if _PRESET == "slide":
        if wide:
            return SLIDE_ASPECT
        return (min(width, SLIDE_ASPECT[0]), min(height, SLIDE_ASPECT[1]))
    return (min(width, cap), min(height, cap))


def label(text) -> str:
    """Display form of a machine label for the active preset.

    Slide labels lose their underscores: ``non_hodgkin_lymphoma`` is a column name,
    ``non hodgkin lymphoma`` is a phrase a room can read. Print keeps the literal
    identifier, which is what a methods section and a data table should agree on.
    """
    text = str(text)
    return text.replace("_", " ") if _PRESET == "slide" else text
