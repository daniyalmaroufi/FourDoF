#!/usr/bin/env python3
"""
plot_shape.py -- solve and draw the CT-SDR backbone for one joint configuration.

    python3 examples/plot_shape.py --ott 35 --itt 35 --itr 180
    python3 examples/plot_shape.py --ott 35 --itt 70 --itr 180 --thrust 2.0
    python3 examples/plot_shape.py --ott 20 --itt 20 --csv backbone.csv

Writes ``figures/shape_<joints>.png``: the 3D backbone plus the two projections
the NDI tip-tracking analysis already reports (top x-y, side y-z), so a modelled
shape can be laid straight over a measured one.

The backbone is drawn in two sections because they are physically different
beams: where both tubes overlap the bending stiffness is EI_outer + EI_inner,
and past the outer tube's tip only the inner tube carries load.
"""

from __future__ import annotations

import argparse
import csv
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import vizstyle as vs  # noqa: E402
from ctr import CTSDR, ExternalLoad, Joints  # noqa: E402
from ctr.robot import INNER, OUTER  # noqa: E402

import matplotlib.pyplot as plt  # noqa: E402
from mpl_toolkits.mplot3d import Axes3D  # noqa: E402,F401

FIGURES = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "figures")

BOTH = "outer + inner (EI added)"
INNER_ONLY = "inner tube only"


def split_sections(sol):
    """Split the backbone at the outer tube's tip into (both, inner-only)."""
    overlap = ~np.isnan(sol.theta[:, OUTER])
    if overlap.all():
        return [(BOTH, sol.p, vs.SERIES[0])]
    cut = int(np.argmax(~overlap))
    # Repeat the boundary sample so the two polylines meet without a gap.
    lo = slice(0, max(cut, 1))
    hi = slice(max(cut - 1, 0), None)
    out = []
    if cut > 0:
        out.append((BOTH, sol.p[lo], vs.SERIES[0]))
    out.append((INNER_ONLY, sol.p[hi], vs.SERIES[1]))
    return out


def draw(sol, joints, load_label, path):
    vs.apply()
    mm = sol.p * 1e3
    sections = split_sections(sol)

    fig = plt.figure(figsize=(12.5, 4.4))
    ax3d = fig.add_subplot(1, 3, 1, projection="3d")
    ax_xy = fig.add_subplot(1, 3, 2)
    ax_yz = fig.add_subplot(1, 3, 3)

    for label, p, colour in sections:
        q = p * 1e3
        ax3d.plot(q[:, 0], q[:, 1], q[:, 2], color=colour, lw=2.0, label=label)
        ax_xy.plot(q[:, 0], q[:, 1], color=colour, lw=2.0, label=label)
        ax_yz.plot(q[:, 1], q[:, 2], color=colour, lw=2.0, label=label)

    # Base and tip: the two points anyone reads off this figure.
    for ax, cols in ((ax_xy, (0, 1)), (ax_yz, (1, 2))):
        ax.plot(mm[0, cols[0]], mm[0, cols[1]], "o", ms=5, color=vs.NEUTRAL, zorder=3)
        ax.plot(mm[-1, cols[0]], mm[-1, cols[1]], "o", ms=7,
                color=sections[-1][2], mec=vs.SURFACE, mew=1.5, zorder=4)
    ax3d.scatter(*mm[-1], s=45, color=sections[-1][2], edgecolors=vs.SURFACE,
                 linewidths=1.2, depthshade=False, zorder=5)

    # Tip coordinates go on the side view, not the 3D panel: mpl draws 3D text
    # in screen space and it lands on the z-axis tick labels.
    tip = mm[-1]
    ax_yz.annotate(f"tip ({tip[0]:.1f}, {tip[1]:.1f}, {tip[2]:.1f}) mm",
                   xy=(tip[1], tip[2]), xytext=(10, -2), textcoords="offset points",
                   color=vs.INK, fontsize=8.5, va="center")
    # Direct labels on the sections, so identity never rests on colour alone.
    for label, p, _ in sections:
        q = p * 1e3
        k = len(q) // 2
        ax_yz.annotate(label, xy=(q[k, 1], q[k, 2]), xytext=(8, -10),
                       textcoords="offset points", color=vs.INK_2, fontsize=8.5)

    ax3d.set_title("3D backbone")
    ax3d.set_xlabel("x (mm)")
    ax3d.set_ylabel("y (mm)")
    ax3d.set_zlabel("z (mm)")
    vs.tidy3d(ax3d)
    vs.equal_aspect_3d(ax3d, mm)
    ax3d.view_init(elev=22, azim=-58)

    ax_xy.set_title("top view (x-y)")
    ax_xy.set_xlabel("x (mm)")
    ax_xy.set_ylabel("y (mm)")
    ax_yz.set_title("side view (y-z), z along the guide axis")
    ax_yz.set_xlabel("y (mm)")
    ax_yz.set_ylabel("z (mm)")
    for ax in (ax_xy, ax_yz):
        ax.set_aspect("equal", adjustable="datalim")
        vs.tidy(ax)

    # One section needs no legend box -- the direct label already names it.
    if len(sections) > 1:
        handles, labels = ax_xy.get_legend_handles_labels()
        fig.legend(handles, labels, loc="upper right", ncol=len(labels),
                   bbox_to_anchor=(0.99, 1.04))
    fig.suptitle(f"CT-SDR backbone   {joints}", x=0.005, ha="left",
                 fontsize=12, color=vs.INK, fontweight="semibold")
    vs.caption(fig, f"Geometrically exact Cosserat model (Rucker et al. 2010).  {load_label}.  "
                    f"Shooting residual {sol.residual_norm:.1e}.")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(path)
    plt.close(fig)


def report(robot, sol, joints):
    mm = 1e3
    tip = sol.tip_position * mm
    print(robot.summary())
    print()
    print(f"joints            {joints}")
    print(f"deployed length   outer {robot.deployed_lengths(joints)[OUTER] * mm:.2f} mm, "
          f"inner {robot.deployed_lengths(joints)[INNER] * mm:.2f} mm")
    print(f"backbone length   {sol.arc_length * mm:.2f} mm")
    print(f"tip position      x {tip[0]:8.3f}   y {tip[1]:8.3f}   z {tip[2]:8.3f}  mm")
    print(f"tip tangent       {np.round(sol.tip_tangent, 4)}")
    print(f"base torsion      outer {sol.u_z0[OUTER]:+.4f}   inner {sol.u_z0[INNER]:+.4f}  rad/m")

    # Relative roll is only defined where both tubes are present, which is not
    # the tip once the inner tube has been advanced past the outer one.
    overlap = np.flatnonzero(~np.isnan(sol.theta[:, OUTER]))
    cmd = np.degrees(sol.alphas[INNER] - sol.alphas[OUTER])
    last = overlap[-1]
    at_end = np.degrees(sol.theta[last, INNER] - sol.theta[last, OUTER])
    at_base = np.degrees(sol.theta[0, INNER] - sol.theta[0, OUTER])
    where = "the tip" if last == sol.s.size - 1 else f"s={sol.s[last] * mm:.1f} mm (outer tip)"
    print(f"relative roll     commanded {cmd:8.2f} deg -> {at_base:8.2f} at s=0 "
          f"-> {at_end:8.2f} at {where}")
    print(f"                  transmission windup {cmd - at_base:.2f} deg, "
          f"deployed windup {at_base - at_end:.2f} deg")
    print(f"shooting residual {sol.residual_norm:.2e}  ({sol.n_iter} integrations, "
          f"{sol.restarts} restarts)")

    print("\nmaterial utilisation (linear-elastic model -- see README on nitinol):")
    for row in robot.material_limits(sol):
        if row is None:
            continue
        flag = "" if row["shear_ok"] else "   <-- past the linear-elastic range"
        print(f"  {row['tube']:<6s} peak bending strain {row['bending_strain'] * 100:5.2f} %   "
              f"torque {row['torque']:+7.4f} N.m   "
              f"peak shear {row['shear_stress'] / 1e6:6.1f} MPa{flag}")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ott", type=float, default=35.0, help="outer tube translation (mm)")
    ap.add_argument("--itt", type=float, default=35.0, help="inner tube translation (mm)")
    ap.add_argument("--otr", type=float, default=0.0, help="outer tube rotation (deg)")
    ap.add_argument("--itr", type=float, default=180.0, help="inner tube rotation (deg)")
    ap.add_argument("--thrust", type=float, default=0.0,
                    help="drilling thrust along the tip tangent (N)")
    ap.add_argument("--lateral", type=float, nargs=2, default=(0.0, 0.0),
                    metavar=("FX", "FY"), help="transverse tip force in the tip frame (N)")
    ap.add_argument("--torque", type=float, default=0.0,
                    help="drill reaction torque about the tip tangent (N.m)")
    ap.add_argument("--gravity", action="store_true", help="include tube self-weight")
    ap.add_argument("--config", default=None, help="path to ct_sdr.yaml")
    ap.add_argument("--csv", default=None, help="also write the backbone samples here")
    ap.add_argument("--out", default=None, help="output figure path")
    args = ap.parse_args(argv)

    robot = CTSDR.from_yaml(args.config) if args.config else CTSDR.from_yaml()
    joints = Joints(ott=args.ott, itt=args.itt, otr=args.otr, itr=args.itr)
    robot.check(joints)

    free = robot.solve(joints)
    loaded = bool(args.thrust or args.torque or any(args.lateral) or args.gravity)
    if loaded:
        drill = robot.drilling_load(args.thrust, args.lateral, args.torque, solution=free)
        if args.gravity:
            load = robot.gravity_load(joints, tip_force=drill.tip_force,
                                      tip_moment=drill.tip_moment)
        else:
            load = drill
        sol = robot.solve(joints, load=load, guess=free.u_z0)
        bits = []
        if args.thrust:
            bits.append(f"{args.thrust:g} N thrust")
        if any(args.lateral):
            bits.append(f"{tuple(args.lateral)} N lateral")
        if args.torque:
            bits.append(f"{args.torque:g} N.m drill torque")
        if args.gravity:
            bits.append("self-weight")
        label = "Loaded: " + ", ".join(bits)
        shift = np.linalg.norm(sol.tip_position - free.tip_position) * 1e3
    else:
        sol, label, shift = free, "Unloaded (free space)", 0.0

    report(robot, sol, joints)
    if loaded:
        print(f"load-induced tip shift {shift:.4f} mm vs the unloaded shape")

    os.makedirs(FIGURES, exist_ok=True)
    name = (f"shape_ott{args.ott:g}_itt{args.itt:g}_otr{args.otr:g}_itr{args.itr:g}"
            f"{'_loaded' if loaded else ''}.png")
    out = args.out or os.path.join(FIGURES, name)
    draw(sol, joints, label, out)
    print(f"\nwrote {out}")

    if args.csv:
        with open(args.csv, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["s_mm", "x_mm", "y_mm", "z_mm",
                        "theta_outer_deg", "theta_inner_deg",
                        "uz_outer_1pm", "uz_inner_1pm"])
            for k in range(sol.s.size):
                w.writerow([f"{sol.s[k] * 1e3:.5f}"]
                           + [f"{v:.6f}" for v in sol.p[k] * 1e3]
                           + [f"{np.degrees(sol.theta[k, i]):.6f}" for i in (OUTER, INNER)]
                           + [f"{sol.u_z[k, i]:.6f}" for i in (OUTER, INNER)])
        print(f"wrote {args.csv}")


if __name__ == "__main__":
    main()
