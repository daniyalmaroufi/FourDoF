#!/usr/bin/env python3
"""
friction_study.py -- calibrate Coulomb friction on set2, then test it elsewhere.

`set2` is the one experiment where friction is the **only** loss mechanism the
model has.  Only the inner tube is deployed, so there is no second tube to
couple to: the frictionless model predicts a perfect 360 deg of tip rotation
and has no way at all to produce the 8.2 deg the tracker actually recorded.
Everything that shortfall contains is friction against the guide bore.

That makes it a calibration, exactly like `set2`'s role for the tube geometry:

1. **Identify** the guide contact force from `set2`'s 8.2 deg loss.
2. **Predict** with it on `set4a` / `set4b`, where the tubes also rub on each
   other, and compare against those sets' measured rotations.

The prediction is the interesting part, and it is not one-signed.  Held at a
fixed part-wound pose, friction is a pure torque sink between the actuator and
the deployed section and costs delivered rotation.  But across a whole 0 -> 180
deg sweep the deployed section winds *and then unwinds*, and friction resists
the unwinding too -- holding the assembly more wound than it would otherwise
settle.  On these two sets the second effect wins for `set4a` (37 -> 71 deg
against 60 measured) and the first wins for `set4b` (55 -> 39 against 76).

    python3 validation/friction_study.py
    python3 validation/friction_study.py --steps 25
"""

from __future__ import annotations

import argparse
import os
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_MODELING = os.path.dirname(_HERE)
sys.path.insert(0, _MODELING)
sys.path.insert(0, os.path.join(_MODELING, "examples"))
sys.path.insert(0, _HERE)

import ndi_experiments as X  # noqa: E402
import vizstyle as vs  # noqa: E402
from ctr import FrictionModel, Joints  # noqa: E402
from ctr.friction import MU_NITINOL_NITINOL, MU_STEEL_NITINOL  # noqa: E402
from ctr.robot import INNER  # noqa: E402
from run_validation import FITTED, OUT_DIR, build_robot  # noqa: E402

import matplotlib.pyplot as plt  # noqa: E402

#: Length of guide the retracted tube is gripped over, for the a-priori normal
#: force estimate.  Not measured -- which is exactly why it gets calibrated.
GUIDE_GRIP_M = 0.15

#: set2: inner tube alone, ITR 360 deg, measured 351.8 deg (8.2 deg lost).
SET2_LOSS_DEG = 8.2

CASES = {
    "set4a": dict(ott=35.0, itt=35.0, cmd=180.0, measured=60.0),
    "set4b": dict(ott=17.5, itt=17.5, cmd=180.0, measured=76.2),
}


def make_friction(robot, scale, direction=1.0, mu_tt=MU_NITINOL_NITINOL):
    f = FrictionModel(mu_tube_tube=mu_tt, mu_guide=MU_STEEL_NITINOL,
                      direction=direction, force_scale=scale)
    f.guide_normal_force = FrictionModel.estimate_guide_force(
        robot.tubes[INNER], GUIDE_GRIP_M)
    return f


def set2_loss(robot, scale):
    """Rotation lost at the tip over set2's revolution [deg].

    A single deployed tube reaches its own free end with zero torsional moment,
    so the whole of this is the guide contact -- no continuation sweep needed,
    the end pose carries it.
    """
    j = Joints(ott=0.0, itt=35.0, itr=360.0)
    base = robot.solve(j)
    if scale <= 0.0:
        return 0.0
    sol = robot.solve(j, friction=make_friction(robot, scale))
    return float(abs(np.degrees(sol.theta[0, INNER] - base.theta[0, INNER])))


def calibrate(robot, target=SET2_LOSS_DEG):
    """Contact-force multiplier that reproduces set2's measured loss."""
    from scipy.optimize import brentq
    lo, hi = 0.05, 60.0
    if set2_loss(robot, hi) < target:
        return hi, False
    return brentq(lambda x: set2_loss(robot, x) - target, lo, hi, xtol=1e-3), True


def sweep(robot, case, scale, steps=19, reverse=False):
    """Tip path over a case's ITR step, by continuation."""
    angles = np.linspace(0.0, case["cmd"], steps)
    direction = 1.0
    if reverse:
        angles, direction = angles[::-1], -1.0
    guess, tips = None, []
    for a in angles:
        j = Joints(ott=case["ott"], itt=case["itt"], itr=a)
        f = None if scale <= 0.0 else make_friction(robot, scale, direction)
        sol = robot.solve(j, friction=f, guess=guess, strict=False)
        guess = sol.u_z0
        tips.append(sol.position_at(robot.deployed_lengths(j)[INNER]))
    tips = np.asarray(tips) * 1e3
    return tips[::-1] if reverse else tips


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--steps", type=int, default=19)
    ap.add_argument("--out-dir", default=OUT_DIR)
    args = ap.parse_args(argv)

    robot = build_robot(**FITTED)
    os.makedirs(args.out_dir, exist_ok=True)

    w0 = FrictionModel.estimate_guide_force(robot.tubes[INNER], GUIDE_GRIP_M)
    print(f"dry coefficients: nitinol-nitinol {MU_NITINOL_NITINOL}, "
          f"steel-nitinol {MU_STEEL_NITINOL}")
    print(f"a-priori guide contact force {w0:.1f} N/m "
          f"over a {GUIDE_GRIP_M * 1e3:.0f} mm grip\n")

    print("1. calibrate on set2 (inner tube alone -- friction is the only loss)")
    scales = np.array([0.0, 1.0, 2.0, 4.0, 8.0, 16.0])
    losses = [set2_loss(robot, s) for s in scales]
    for s, v in zip(scales, losses):
        print(f"   contact x{s:4.1f}  ({w0 * s:7.1f} N/m)  ->  {v:5.2f} deg lost")
    scale, ok = calibrate(robot)
    print(f"   measured {SET2_LOSS_DEG} deg  ->  contact x{scale:.2f} "
          f"({w0 * scale:.0f} N/m){'' if ok else '  [not bracketed]'}\n")

    print("2. predict with it on the 180 deg inner-tube rotations")
    print(f"   {'set':7s}{'measured':>10s}{'frictionless':>14s}{'calibrated':>12s}")
    preds = {}
    for key, case in CASES.items():
        free = X.arc_metrics(sweep(robot, case, 0.0, args.steps))["swept"]
        rubbed = X.arc_metrics(sweep(robot, case, scale, args.steps))["swept"]
        preds[key] = (free, rubbed)
        print(f"   {key:7s}{case['measured']:9.1f}d{free:13.1f}d{rubbed:11.1f}d")
    print("   -> friction is large here but not consistently signed:\n         it closes most of set4a's gap and opens set4b's.\n")

    print("3. hysteresis at the calibrated contact force")
    hyst = {}
    for key, case in CASES.items():
        fwd = sweep(robot, case, scale, args.steps)
        rev = sweep(robot, case, scale, args.steps, reverse=True)
        gap = np.linalg.norm(fwd - rev, axis=1)
        hyst[key] = (fwd, rev, gap)
        print(f"   {key}: mean {gap.mean():.3f} mm, max {gap.max():.3f} mm")

    path = os.path.join(args.out_dir, "friction_study.png")
    draw(scales, losses, w0, scale, preds, hyst, path)
    print(f"\nwrote {path}")
    return 0


def draw(scales, losses, w0, scale, preds, hyst, path):
    vs.apply()
    fig, axes = plt.subplots(1, 3, figsize=(13.4, 4.4))
    ax_a, ax_b, ax_c = axes

    ax_a.plot(scales * w0, losses, color=vs.SERIES[0], lw=2.0, marker="o", ms=5,
              mec=vs.SURFACE, mew=1.2, label="model")
    ax_a.axhline(SET2_LOSS_DEG, color=vs.INK_2, lw=1.2, ls=(0, (4, 3)), zorder=2)
    ax_a.annotate(f"measured {SET2_LOSS_DEG}°", xy=(scales[-1] * w0, SET2_LOSS_DEG),
                  xytext=(-4, 5), textcoords="offset points", ha="right",
                  color=vs.INK_2, fontsize=8.5)
    ax_a.plot([scale * w0], [SET2_LOSS_DEG], "o", ms=9, color=vs.SERIES[1],
              mec=vs.SURFACE, mew=1.5, zorder=5)
    ax_a.annotate(f"calibrated\n{scale * w0:.0f} N/m", xy=(scale * w0, SET2_LOSS_DEG),
                  xytext=(8, -18), textcoords="offset points",
                  color=vs.SERIES[1], fontsize=8.5)
    ax_a.set_xlabel("guide contact force (N/m)")
    ax_a.set_ylabel("rotation lost over 360° (deg)")
    ax_a.set_title("1. calibrate on set2 (no tube-tube coupling)")
    vs.tidy(ax_a)

    keys = list(preds)
    xs = np.arange(len(keys))
    w = 0.26
    meas = [CASES[k]["measured"] / CASES[k]["cmd"] * 100 for k in keys]
    free = [preds[k][0] / CASES[k]["cmd"] * 100 for k in keys]
    rub = [preds[k][1] / CASES[k]["cmd"] * 100 for k in keys]
    for off, vals, colour, lab in ((-w, meas, vs.SERIES[0], "measured"),
                                   (0.0, free, vs.NEUTRAL, "model, frictionless"),
                                   (w, rub, vs.SERIES[1], "model, calibrated friction")):
        ax_b.bar(xs + off, vals, width=w - 0.03, color=colour, label=lab, zorder=3)
    ax_b.set_xticks(xs)
    ax_b.set_xticklabels([f"{k}\nITR 180°" for k in keys], fontsize=8.5)
    ax_b.set_ylabel("rotation reaching the tip (% of commanded)")
    ax_b.set_title("2. large, but inconsistent in sign")
    ax_b.set_ylim(0, 62)
    ax_b.legend(loc="upper right", fontsize=8)
    vs.tidy(ax_b)

    key = "set4a"
    fwd, rev, gap = hyst[key]
    c, u_ax, v_ax, _, _ = X.ndi.fit_plane(np.vstack([fwd, rev]))
    for p, colour, lab in ((fwd, vs.SERIES[0], "winding up (0 → 180°)"),
                           (rev, vs.SERIES[1], "unwinding (180° → 0)")):
        u, v = X.ndi.project_plane(p, c, u_ax, v_ax)
        ax_c.plot(u, v, color=colour, lw=2.0, label=lab)
        ax_c.plot(u[0], v[0], "o", ms=5, color=colour, mec=vs.SURFACE, mew=1.2)
    ax_c.set_title(f"3. {key} hysteresis: max {gap.max():.2f} mm")
    ax_c.set_xlabel("in-plane u (mm)")
    ax_c.set_ylabel("in-plane v (mm)")
    ax_c.set_aspect("equal", adjustable="datalim")
    ax_c.legend(loc="best", fontsize=8)
    vs.tidy(ax_c)

    fig.suptitle("Coulomb friction: calibrated on set2, tested on the 180° rotations",
                 x=0.005, ha="left", fontsize=12, color=vs.INK, fontweight="semibold")
    vs.caption(fig, "Friction resists both the winding and the unwinding leg of a sweep, so "
                    "its net sign is configuration-dependent -- it closes most of set4a's gap "
                    "and opens set4b's.\nThe hysteresis loop is the prediction worth testing: "
                    "sweep ITR out and back, and the return path should not retrace the "
                    "outbound one.")
    fig.tight_layout(rect=(0, 0, 1, 0.90))
    fig.savefig(path)
    plt.close(fig)


if __name__ == "__main__":
    sys.exit(main())
