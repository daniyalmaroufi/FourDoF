#!/usr/bin/env python3
"""
plot_workspace.py -- tip loci swept by the inner tube's rotation, per deployment.

    python3 examples/plot_workspace.py
    python3 examples/plot_workspace.py --deployments 10 20 30 --steps 72

For each deployed length the inner tube is rotated a full turn while the outer
tube is held, and the tip is traced.  Each locus is one deployment, so the
series are an ordered magnitude and take the single-hue ordinal ramp rather
than categorical colours.

The sweep is solved by continuation (``CTSDR.solve_path``), marching each
configuration from its neighbour's converged torsions.  That is not just a
speed-up: once the tubes are wound up the torsional BVP has more than one
solution branch, and a cold start can land on a branch the hardware would not
follow.  The printed windup column is the model's headline result -- the tip
does *not* rotate by the amount the ITR motor was commanded.
"""

from __future__ import annotations

import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import vizstyle as vs  # noqa: E402
from ctr import CTSDR, Joints  # noqa: E402
from ctr.robot import INNER, OUTER  # noqa: E402

import matplotlib.pyplot as plt  # noqa: E402
from mpl_toolkits.mplot3d import Axes3D  # noqa: E402,F401

FIGURES = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "figures")


def sweep(robot, deployment_mm, angles_deg):
    """Trace the tip while ITR turns, holding OTT = ITT = ``deployment_mm``."""
    path = [Joints(ott=deployment_mm, itt=deployment_mm, otr=0.0, itr=a) for a in angles_deg]
    sols = robot.solve_path(path)
    tips = np.array([s.tip_position for s in sols]) * 1e3
    windup = np.array([
        np.degrees(s.alphas[INNER] - s.alphas[OUTER])
        - np.degrees(s.theta[-1, INNER] - s.theta[-1, OUTER])
        for s in sols
    ])
    return tips, windup, sols


def draw(results, angles, path):
    vs.apply()
    colours = vs.ramp(len(results))

    fig = plt.figure(figsize=(13.0, 4.6))
    ax3d = fig.add_subplot(1, 3, 1, projection="3d")
    ax_xy = fig.add_subplot(1, 3, 2)
    ax_w = fig.add_subplot(1, 3, 3)

    all_pts = np.vstack([tips for _, tips, _ in results])
    for (d, tips, windup), colour in zip(results, colours):
        label = f"{d:g} mm"
        ax3d.plot(tips[:, 0], tips[:, 1], tips[:, 2], color=colour, lw=2.0, label=label)
        ax_xy.plot(tips[:, 0], tips[:, 1], color=colour, lw=2.0, label=label)
        ax_w.plot(angles, windup, color=colour, lw=2.0, label=label)
        # Direct label at the widest point of each locus, so the ordering is
        # readable without tracing the legend back to a colour.
        k = int(np.argmax(tips[:, 0]))
        ax_xy.annotate(label, xy=(tips[k, 0], tips[k, 1]), xytext=(6, 0),
                       textcoords="offset points", color=vs.INK_2, fontsize=8.5,
                       va="center")
        # Label at the peak, not the right edge: after a snap the curves all
        # collapse toward zero there and the labels would pile up.
        pk = int(np.argmax(windup))
        ax_w.annotate(label, xy=(angles[pk], windup[pk]), xytext=(-6, 5),
                      textcoords="offset points", color=vs.INK_2, fontsize=8.5,
                      ha="right")
        # Mark the fold where the wound branch ceases to exist and the robot
        # snaps to the surviving equilibrium.
        jump = int(np.argmin(np.diff(windup)))
        if windup[jump] - windup[jump + 1] > 20.0:
            ax_w.plot(angles[jump], windup[jump], "o", ms=6, color=colour,
                      mec=vs.SURFACE, mew=1.4, zorder=5)

    ax3d.set_title("tip locus over a full ITR turn")
    ax3d.set_xlabel("x (mm)")
    ax3d.set_ylabel("y (mm)")
    ax3d.set_zlabel("z (mm)")
    vs.tidy3d(ax3d)
    vs.equal_aspect_3d(ax3d, all_pts)
    ax3d.view_init(elev=24, azim=-60)

    ax_xy.set_title("top view (x-y)")
    ax_xy.set_xlabel("x (mm)")
    ax_xy.set_ylabel("y (mm)")
    ax_xy.set_aspect("equal", adjustable="datalim")
    vs.tidy(ax_xy)

    ax_w.set_title("rotation lost to elastic windup")
    ax_w.set_xlabel("commanded ITR (deg)")
    ax_w.set_ylabel("commanded - achieved relative roll at the tip (deg)")
    ax_w.axhline(0.0, color=vs.NEUTRAL, lw=1.0, ls=(0, (4, 3)), zorder=1)
    ax_w.set_xticks(np.arange(0, 361, 90))
    vs.tidy(ax_w)

    handles, labels = ax_xy.get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper right", ncol=len(labels),
               bbox_to_anchor=(0.99, 1.05), title="deployed length",
               title_fontsize=9)
    fig.suptitle("CT-SDR workspace and torsional windup   (OTT = ITT, OTR = 0)",
                 x=0.005, ha="left", fontsize=12, color=vs.INK, fontweight="semibold")
    vs.caption(
        fig,
        "Geometrically exact Cosserat model (Rucker et al. 2010), solved by continuation over the "
        "ITR sweep.  Dots mark snap-through: the wound branch folds away and the robot jumps to the\n"
        "surviving equilibrium.  Loci are open, not closed, because a full turn leaves the tubes "
        "wound - one revolution of ITR does not return the tip to where it started.",
    )
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(path)
    plt.close(fig)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--deployments", type=float, nargs="+", default=[10, 20, 30, 40, 50],
                    help="OTT = ITT values to sweep (mm); at most 5")
    ap.add_argument("--steps", type=int, default=49, help="ITR samples over one turn")
    ap.add_argument("--config", default=None)
    ap.add_argument("--out", default=None)
    args = ap.parse_args(argv)

    robot = CTSDR.from_yaml(args.config) if args.config else CTSDR.from_yaml()
    angles = np.linspace(0.0, 360.0, args.steps)

    results = []
    print(f"{'deploy':>8}  {'tip x range':>13}  {'tip y range':>13}  "
          f"{'max windup':>11}  {'max |resid|':>11}")
    for d in args.deployments:
        robot.check(Joints(ott=d, itt=d))
        tips, windup, sols = sweep(robot, d, angles)
        results.append((d, tips, windup))
        print(f"{d:8.1f}  {tips[:, 0].ptp():13.3f}  {tips[:, 1].ptp():13.3f}  "
              f"{np.abs(windup).max():10.2f}d  "
              f"{max(s.residual_norm for s in sols):11.2e}")

    os.makedirs(FIGURES, exist_ok=True)
    out = args.out or os.path.join(FIGURES, "workspace_itr_sweep.png")
    draw(results, angles, out)
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
