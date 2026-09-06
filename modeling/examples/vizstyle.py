"""
vizstyle.py -- one visual system for every figure the modeling code produces.

Colour roles rather than raw hex, so the example scripts read as one set.  The
values are the validated default data-viz palette:

* Categorical slots 1-3 (blue / orange / aqua) are used in fixed order, never
  cycled.  Those three clear the all-pairs colour-vision-deficiency and
  normal-vision separation floors in light mode, which is the gate that matters
  for scatter- and locus-style plots where any two series can end up adjacent.
  Aqua sits at 2.74:1 against the surface, below the 3:1 bar, so every figure
  that uses it also carries a visible direct label -- identity is never left to
  colour alone.
* The blue ordinal ramp encodes deployed length, which is an ordered magnitude.
  Its light end (step 250) clears 2:1 against the surface.

Figures are written as PNG: these are thin-line plots with small text, and JPEG
ringing around 2 px strokes is clearly visible at the sizes used here.
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

# -- roles ------------------------------------------------------------------

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_2 = "#52514e"
GRID = "#e5e3de"
SPINE = "#b9b7b0"

#: Categorical hues, assigned in this order and never cycled.
SERIES = ("#2a78d6", "#eb6834", "#1baf7a")

#: Single-hue ordinal ramp, light -> dark, for deployed length.
RAMP = ("#86b6ef", "#5598e7", "#2a78d6", "#1c5cab", "#104281")

#: Reference / annotation marks that must not read as a data series.
NEUTRAL = "#8a8880"


def ramp(n: int):
    """``n`` steps of the ordinal ramp, keeping the largest available gaps."""
    if n <= 1:
        return [RAMP[2]]
    if n > len(RAMP):
        raise ValueError(
            f"{n} ordered levels exceeds the {len(RAMP)}-step ramp; bin the "
            "levels or facet instead of inventing steps"
        )
    idx = [round(i * (len(RAMP) - 1) / (n - 1)) for i in range(n)]
    return [RAMP[i] for i in idx]


def apply() -> None:
    """Install the shared rcParams.  Call once, before creating figures."""
    plt.rcParams.update(
        {
            "figure.facecolor": SURFACE,
            "savefig.facecolor": SURFACE,
            "axes.facecolor": SURFACE,
            "axes.edgecolor": SPINE,
            "axes.labelcolor": INK_2,
            "axes.titlecolor": INK,
            "axes.titlesize": 11,
            "axes.titleweight": "semibold",
            "axes.titlelocation": "left",
            "axes.titlepad": 10,
            "axes.labelsize": 9,
            "axes.linewidth": 0.8,
            "axes.grid": True,
            "axes.axisbelow": True,
            "grid.color": GRID,
            "grid.linewidth": 0.8,
            "xtick.color": INK_2,
            "ytick.color": INK_2,
            "xtick.labelsize": 8.5,
            "ytick.labelsize": 8.5,
            "xtick.direction": "out",
            "ytick.direction": "out",
            "lines.linewidth": 2.0,
            "lines.solid_capstyle": "round",
            "legend.frameon": False,
            "legend.fontsize": 9,
            "legend.labelcolor": INK_2,
            "font.size": 9.5,
            "figure.dpi": 150,
            "savefig.dpi": 150,
            "savefig.bbox": "tight",
        }
    )


def tidy(ax) -> None:
    """Drop the top/right spines so the grid stays recessive."""
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)


def tidy3d(ax) -> None:
    """Make a 3D axes' panes recede instead of boxing the data in."""
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        axis.pane.set_facecolor(SURFACE)
        axis.pane.set_edgecolor(GRID)
        axis.pane.set_alpha(1.0)
        axis._axinfo["grid"]["color"] = GRID
        axis._axinfo["grid"]["linewidth"] = 0.6
    ax.tick_params(colors=INK_2, labelsize=8)


def equal_aspect_3d(ax, pts) -> None:
    """Give a 3D axes a true 1:1:1 aspect around ``pts`` (an (M, 3) array).

    Without this a backbone's shape is distorted by the axis ranges and the
    plotted curvature is not the modelled curvature.
    """
    lo = pts.min(axis=0)
    hi = pts.max(axis=0)
    mid = 0.5 * (lo + hi)
    half = 0.5 * max((hi - lo).max(), 1e-9) * 1.1
    ax.set_xlim(mid[0] - half, mid[0] + half)
    ax.set_ylim(mid[1] - half, mid[1] + half)
    ax.set_zlim(mid[2] - half, mid[2] + half)
    try:
        ax.set_box_aspect((1.0, 1.0, 1.0))
    except AttributeError:  # pragma: no cover - matplotlib < 3.3
        pass


def caption(fig, text: str) -> None:
    """One line of provenance under the figure."""
    fig.text(0.0, -0.02, text, ha="left", va="top", fontsize=8, color=INK_2)
