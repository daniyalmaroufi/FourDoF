#!/usr/bin/env python3
"""
plot_drilling_load.py -- how far the tip moves when the drill bit pushes back.

    python3 examples/plot_drilling_load.py
    python3 examples/plot_drilling_load.py --deployments 20 35 --max-force 5

This is the part of the model that a constant-curvature kinematic model cannot
give you.  During drilling the bit applies a reaction wrench at the tip, and
the resulting tip deviation is the error between where the robot was commanded
to drill and where it actually drills.  The model carries that load as the
internal force/moment states ``n(s)``, ``m(s)`` with the tip boundary condition
``n(L) = F_tip``, ``m(L) = M_tip``.

Two loading directions are shown separately, never on shared axes: a lateral
cutting force perpendicular to the tip tangent, and an axial thrust along it.
They differ by more than an order of magnitude, which is the point -- a
concentric-tube drill is stiff along its axis and soft across it.
"""

from __future__ import annotations

import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import vizstyle as vs  # noqa: E402
from ctr import CTSDR, Joints  # noqa: E402

import matplotlib.pyplot as plt  # noqa: E402

FIGURES = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "figures")


def load_sweep(robot, joints, forces, axial):
    """Tip deviation from the unloaded shape over a range of tip forces [mm].

    The wrench is defined in the *unloaded* tip frame.  Re-solving the frame
    against the deflected shape changes the answer by well under the plotted
    resolution at these loads, so one pass is enough; the second pass is what
    you would add if the deflections ever grew comparable to the tube radius.
    """
    free = robot.solve(joints)
    dev, guess = [], free.u_z0
    for F in forces:
        if axial:
            load = robot.drilling_load(F, solution=free)
        else:
            load = robot.drilling_load(0.0, lateral=(F, 0.0), solution=free)
        sol = robot.solve(joints, load=load, guess=guess)
        guess = sol.u_z0
        dev.append(np.linalg.norm(sol.tip_position - free.tip_position) * 1e3)
    return np.array(dev), free


def deviation_profile(robot, joints, force, free):
    """|loaded - free| at every arc length [mm], for a lateral tip force.

    Plotted instead of an overlay of the two backbones: at these loads the
    deviation is a few tenths of a millimetre against a 35 mm arc, so the two
    shapes are indistinguishable when drawn to true scale, and drawing them to
    any other scale would misrepresent the deflection.
    """
    load = robot.drilling_load(0.0, lateral=(force, 0.0), solution=free)
    loaded = robot.solve(joints, load=load, guess=free.u_z0)
    return free.s * 1e3, np.linalg.norm(loaded.p - free.p, axis=1) * 1e3


def draw(deployments, forces, lateral, axial, shapes, max_force, path):
    vs.apply()
    colours = vs.ramp(len(deployments))

    fig, axes = plt.subplots(1, 3, figsize=(13.0, 4.3))
    ax_lat, ax_ax, ax_shape = axes

    for d, dev, colour in zip(deployments, lateral, colours):
        ax_lat.plot(forces, dev, color=colour, lw=2.0, label=f"{d:g} mm")
        ax_lat.annotate(f"{d:g} mm", xy=(forces[-1], dev[-1]), xytext=(-4, 6),
                        textcoords="offset points", ha="right",
                        color=vs.INK_2, fontsize=8.5)
    for d, dev, colour in zip(deployments, axial, colours):
        ax_ax.plot(forces, dev, color=colour, lw=2.0, label=f"{d:g} mm")
        ax_ax.annotate(f"{d:g} mm", xy=(forces[-1], dev[-1]), xytext=(-4, 6),
                       textcoords="offset points", ha="right",
                       color=vs.INK_2, fontsize=8.5)

    ax_lat.set_title("lateral cutting force")
    ax_ax.set_title("axial thrust")
    for ax in (ax_lat, ax_ax):
        ax.set_xlabel("tip force (N)")
        ax.set_ylabel("tip deviation from the free shape (mm)")
        ax.set_xlim(0, max_force)
        ax.set_ylim(bottom=0)
        vs.tidy(ax)

    for (d, (s_mm, dev)), colour in zip(shapes, colours):
        ax_shape.plot(s_mm, dev, color=colour, lw=2.0, label=f"{d:g} mm")
        # Below the end point: the shorter runs finish underneath the longer
        # ones, so labels placed above would land on another curve.
        ax_shape.annotate(f"{d:g} mm", xy=(s_mm[-1], dev[-1]), xytext=(-3, -13),
                          textcoords="offset points", ha="right",
                          color=vs.INK_2, fontsize=8.5)
    ax_shape.set_title(f"deviation along the arc ({max_force:g} N lateral)")
    ax_shape.set_xlabel("arc length from the guide exit (mm)")
    ax_shape.set_ylabel("deviation from the free shape (mm)")
    ax_shape.set_xlim(left=0)
    ax_shape.set_ylim(bottom=0)
    vs.tidy(ax_shape)

    handles, labels = ax_lat.get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper right", ncol=len(labels),
               bbox_to_anchor=(0.99, 1.05), title="deployed length",
               title_fontsize=9)
    fig.suptitle("CT-SDR tip deviation under a drilling reaction wrench",
                 x=0.005, ha="left", fontsize=12, color=vs.INK, fontweight="semibold")
    vs.caption(fig, "Geometrically exact Cosserat model (Rucker et al. 2010); wrench applied "
                    "at the tip in the unloaded tip frame.")
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(path)
    plt.close(fig)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--deployments", type=float, nargs="+", default=[17.5, 26.0, 35.0],
                    help="OTT = ITT values to compare (mm); at most 5")
    ap.add_argument("--max-force", type=float, default=3.0, help="largest tip force (N)")
    ap.add_argument("--steps", type=int, default=9, help="force samples")
    ap.add_argument("--itr", type=float, default=0.0, help="inner tube rotation (deg)")
    ap.add_argument("--config", default=None)
    ap.add_argument("--out", default=None)
    args = ap.parse_args(argv)

    robot = CTSDR.from_yaml(args.config) if args.config else CTSDR.from_yaml()
    forces = np.linspace(0.0, args.max_force, args.steps)

    lateral, axial, frees = [], [], []
    print(f"{'deploy':>8}  {'lateral dev @ max':>18}  {'axial dev @ max':>16}  "
          f"{'stiffness ratio':>15}")
    for d in args.deployments:
        j = Joints(ott=d, itt=d, itr=args.itr)
        robot.check(j)
        lat, free = load_sweep(robot, j, forces, axial=False)
        ax_, _ = load_sweep(robot, j, forces, axial=True)
        lateral.append(lat)
        axial.append(ax_)
        frees.append((j, free))
        ratio = lat[-1] / ax_[-1] if ax_[-1] > 0 else float("inf")
        print(f"{d:8.1f}  {lat[-1]:15.4f} mm  {ax_[-1]:13.4f} mm  {ratio:15.1f}x")

    profiles = [
        (d, deviation_profile(robot, j, args.max_force, free))
        for d, (j, free) in zip(args.deployments, frees)
    ]

    os.makedirs(FIGURES, exist_ok=True)
    out = args.out or os.path.join(FIGURES, "drilling_load.png")
    draw(args.deployments, forces, lateral, axial, profiles, args.max_force, out)
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
