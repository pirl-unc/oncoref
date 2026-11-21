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

"""Consistent publication exports for package and analysis figures."""

from __future__ import annotations

from pathlib import Path

import matplotlib


def queue_scatter_labels(texts, ax) -> None:
    """Place labels only after titles, legends, scales and layout are final."""
    groups = getattr(ax.figure, "_pirly_label_groups", [])
    groups.append((ax, list(texts)))
    ax.figure._pirly_label_groups = groups


def _place_scatter_labels(fig) -> None:
    """Deterministic nearest-free label placement in display coordinates.

    Reserve the full painted marker extents, including highlight/outer rings.
    Each label stays inside the axes and clears other labels and fixed notes.
    Point coordinates and marker areas are never moved to make labels fit.
    """
    import numpy as np
    from matplotlib.collections import PathCollection
    from matplotlib.legend import Legend
    from matplotlib.patches import FancyArrowPatch
    from matplotlib.transforms import Bbox

    groups = getattr(fig, "_pirly_label_groups", [])
    if not groups:
        return
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    pad = renderer.points_to_pixels(1.5)
    for ax, texts in groups:
        texts = [t for t in texts if t.get_visible() and t.get_text()]
        bounds = ax.get_window_extent().padded(-pad)
        obstacles = []
        # A scatter can carry multiple paths, sizes and border widths. Compute
        # each marker in pixels rather than treating it as a zero-size point.
        for collection in ax.collections:
            if not isinstance(collection, PathCollection) or not collection.get_visible():
                continue
            offsets = collection.get_offset_transform().transform(collection.get_offsets())
            paths, transforms = collection.get_paths(), collection.get_transforms()
            widths = collection.get_linewidths()
            for i, center in enumerate(offsets):
                if not np.isfinite(center).all() or not paths or not len(transforms):
                    continue
                vertices = paths[i % len(paths)].vertices
                matrix = transforms[i % len(transforms)]
                vertices = vertices @ matrix[:2, :2].T + matrix[:2, 2] + center
                border = (
                    renderer.points_to_pixels(widths[i % len(widths)] / 2) if len(widths) else 0
                )
                box = Bbox.from_extents(*vertices.min(axis=0), *vertices.max(axis=0))
                obstacles.append(box.padded(border + pad))
        for fixed in ax.texts:
            if fixed not in texts and fixed.get_visible() and fixed.get_text():
                obstacles.append(fixed.get_window_extent(renderer).padded(pad))
        for legend in ax.get_children():
            if isinstance(legend, Legend) and legend.get_visible():
                obstacles.append(legend.get_window_extent(renderer).padded(pad))
        anchors = {t: t.get_transform().transform(t.get_position()) for t in texts}
        # Larger highlighted labels get space first; ties retain input order.
        ordered = sorted(texts, key=lambda t: -t.get_fontsize())
        placed = []
        for text in ordered:
            target = anchors[text]
            if not np.isfinite(target).all():
                continue
            size = text.get_window_extent(renderer)
            width, height = size.width + 2 * pad, size.height + 2 * pad
            step = renderer.points_to_pixels(3)
            dx = np.arange(-bounds.width, bounds.width + step, step)
            dy = np.arange(-bounds.height, bounds.height + step, step)
            xx, yy = np.meshgrid(dx, dy)
            centers = np.column_stack((xx.ravel(), yy.ravel())) + target
            boxes = np.column_stack(
                (
                    centers[:, 0] - width / 2,
                    centers[:, 1] - height / 2,
                    centers[:, 0] + width / 2,
                    centers[:, 1] + height / 2,
                )
            )
            inside = (
                (boxes[:, 0] >= bounds.x0)
                & (boxes[:, 1] >= bounds.y0)
                & (boxes[:, 2] <= bounds.x1)
                & (boxes[:, 3] <= bounds.y1)
            )
            centers, boxes = centers[inside], boxes[inside]
            delta = centers - target
            # Prefer near anchors; a small tie-break favours above/right.
            distance = np.sum(delta**2, axis=1) + 0.01 * (-delta[:, 0] - delta[:, 1])
            order = np.argsort(distance, kind="stable")
            reserved = np.array([b.extents for b in obstacles + placed])
            chosen = None
            for start in range(0, len(order), 256):
                candidates = order[start : start + 256]
                b = boxes[candidates]
                if len(reserved):
                    clashes = (
                        (b[:, None, 0] < reserved[None, :, 2])
                        & (b[:, None, 2] > reserved[None, :, 0])
                        & (b[:, None, 1] < reserved[None, :, 3])
                        & (b[:, None, 3] > reserved[None, :, 1])
                    ).any(axis=1)
                else:
                    clashes = np.zeros(len(b), dtype=bool)
                free = np.flatnonzero(~clashes)
                if len(free):
                    chosen = candidates[free[0]]
                    break
            if chosen is None:
                raise ValueError(
                    f"No collision-free space for {text.get_text()!r}; enlarge the figure"
                )
            center = centers[chosen]
            placed.append(Bbox.from_extents(*boxes[chosen]))
            text.set_ha("center")
            text.set_va("center")
            text.set_position(text.get_transform().inverted().transform(center))
            text.set_zorder(max(8, text.get_zorder()))
            if np.linalg.norm(center - target) > renderer.points_to_pixels(7):
                # Leader lines use data coordinates so PNG and PDF retain the
                # same geometry even though they have different renderers/DPI.
                inv = ax.transData.inverted()
                arrow = FancyArrowPatch(
                    inv.transform(target),
                    inv.transform(center),
                    arrowstyle="-",
                    color="0.6",
                    linewidth=0.45,
                    shrinkA=3,
                    shrinkB=2,
                    patchB=text,
                    transform=ax.transData,
                    zorder=2,
                )
                ax.add_patch(arrow)
    fig._pirly_label_groups = []


def save_figure(fig, path, *, dpi=300, **kwargs) -> None:
    """Save a lossless PNG (at least 300 dpi) and a matching vector PDF.

    Text and paths remain vector objects in PDF; raster artists such as imshow
    retain the requested resolution. Never build the PDF from the PNG.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    options = {"bbox_inches": "tight", "facecolor": "white", "transparent": False}
    options.update(kwargs)
    _place_scatter_labels(fig)
    with matplotlib.rc_context({"pdf.fonttype": 42, "ps.fonttype": 42}):
        fig.savefig(path, dpi=max(300, dpi), **options)
        if path.suffix.lower() == ".png":
            fig.savefig(path.with_suffix(".pdf"), dpi=max(300, dpi), **options)
