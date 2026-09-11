#!/usr/bin/env python3
"""
calibrate_roll.py -- is the tubes' pre-curvature aligned at OTR = ITR = 0?

The rotation joints have **no absolute reference**.  `dxl_control_4dof_cli.py`'s
``home <joint>`` zeroes whatever position the tube is currently at, so
"OTR = ITR = 0" means "homed here", not "the two pre-curvature planes coincide".
Nothing in the recordings forces the two to agree, and if the robot was re-homed
between sessions the offset can differ from one experiment set to the next.

This script treats the initial relative roll

    psi0 = roll_offsets[INNER] - roll_offsets[OUTER]

as one free parameter **per experiment set** and identifies it by scanning it
against that set's measured tip trajectories.

Why it is not obvious a priori that this matters
------------------------------------------------
A misalignment applied at the actuator clamp does *not* arrive at the deployed
section: it is absorbed by the same torsional windup that eats the commanded
rotations (VALIDATION_REPORT.md section 3).  At OTT = ITT = 35 mm a 90 deg clamp
offset delivers only ~14 deg to the deployed overlap, moving the common
curvature from 19.31 to 19.13 /m -- against 14.7 /m if the tubes were
torsionally rigid.  So the advances are nearly blind to psi0.

The *rotation* segments are not: psi0 shifts the whole sweep in phase, and the
tip path is strongly sensitive to it.  That is where the identification comes
from, and it is why the scan is run on whole trials rather than on the advances
alone.

    python3 validation/calibrate_roll.py
    python3 validation/calibrate_roll.py --coarse 30 --n-per-step 15
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
import run_validation as V  # noqa: E402
from run_validation import FITTED, OUT_DIR, build_robot, compare_trial  # noqa: E402

import matplotlib.pyplot as plt  # noqa: E402


def set_score(robot, es, trials, psi_deg, n_per_step):
    """Frame-invariant mismatch for one set at a given psi0 [deg].

    **Not** the registration residual.  A rigid Kabsch transform absorbs a
    great deal of shape error -- that is exactly why geometry is identified on
    turn angles in ``run_validation`` -- and scoring psi0 on it is actively
    misleading: at psi0 = -180 deg `set4a`'s registration RMS *improves* from
    2.28 to 1.57 mm while its traced arc radius goes from 51.7 mm to 132.4 mm
    against a measured 62.5 mm.  The transform simply hides the wrong shape.

    So the score here is built only from quantities that need no registration:
    the accumulated turn of each advance, and the swept angle of each rotation.
    Both are in degrees, so they combine without an arbitrary weight.
    Registration RMS is still returned alongside, as a diagnostic.
    """
    robot.roll_offsets = np.radians([0.0, psi_deg])
    p, slices, _ = X.model_trajectory(robot, es, n_per_step)

    resid, rms_vals = [], []
    for label, _ in es.trials:
        tr = trials[label]
        if not tr.segments:
            continue
        rms_vals.append(compare_trial(tr, es, p, slices)["rms"])
        for k, step in enumerate(es.steps):
            if k >= len(tr.segments):
                break
            if step.kind == "T":
                a = V.turn_profile(tr.segment(k))
                b = V.turn_profile(p[slices[k]])
                ok = np.isfinite(a) & np.isfinite(b)
                resid.extend((a[ok] - b[ok]).tolist())
            else:
                resid.append(X.arc_metrics(tr.segment(k))["swept"]
                             - X.arc_metrics(p[slices[k]])["swept"])
    score = float(np.sqrt(np.mean(np.square(resid)))) if resid else np.nan
    return score, (float(np.mean(rms_vals)) if rms_vals else np.nan)


def scan(robot, es, trials, coarse_step, n_per_step):
    """Coarse scan over a full turn, then a finer pass around the best point."""
    grid = np.arange(-180.0, 180.0, coarse_step)
    out = [set_score(robot, es, trials, a, n_per_step) for a in grid]
    score = np.array([o[0] for o in out])
    rms = np.array([o[1] for o in out])
    best = grid[int(np.nanargmin(score))]

    fine = np.arange(best - coarse_step, best + coarse_step + 1e-9, coarse_step / 3.0)
    out_f = [set_score(robot, es, trials, a, n_per_step) for a in fine]
    score_f = np.array([o[0] for o in out_f])
    j = int(np.nanargmin(score_f))
    return grid, score, rms, float(fine[j]), float(score_f[j]), float(out_f[j][1])


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--coarse", type=float, default=30.0, help="scan step (deg)")
    ap.add_argument("--n-per-step", type=int, default=15)
    ap.add_argument("--out-dir", default=OUT_DIR)
    args = ap.parse_args(argv)

    os.makedirs(args.out_dir, exist_ok=True)
    print("loading NDI trials ...")
    trials = {}
    for es in X.EXPERIMENTS:
        for label, fn in es.trials:
            trials[label] = X.load_trial(label, fn, es.key)

    robot = build_robot(**FITTED)
    print(f"\nscanning psi0 in {args.coarse:.0f} deg steps over a full turn, "
          f"then refining\n")
    print("  score = frame-invariant mismatch (deg); RMS = registration residual (mm)")
    print(f"  {'set':7s}{'score@0':>10s}{'best psi0':>11s}{'score':>9s}"
          f"{'gain':>8s}{'span':>8s}{'RMS@0':>9s}{'RMS@best':>10s}")

    results = {}
    for es in X.EXPERIMENTS:
        grid, score, rms, best, best_score, best_rms = scan(
            robot, es, trials, args.coarse, args.n_per_step)
        zero, zero_rms = set_score(robot, es, trials, 0.0, args.n_per_step)
        span = float(np.nanmax(score) - np.nanmin(score))
        results[es.key] = dict(grid=grid, score=score, rms=rms, best=best,
                               best_score=best_score, best_rms=best_rms,
                               zero=zero, zero_rms=zero_rms, span=span)
        print(f"  {es.key:7s}{zero:9.2f}°{best:10.0f}°{best_score:8.2f}°"
              f"{zero - best_score:7.2f}°{span:7.2f}°{zero_rms:8.2f}mm"
              f"{best_rms:9.2f}mm", flush=True)

    print("\nVerdict per set:")
    for key, r in results.items():
        if r["span"] < 1.0:
            note = "psi0 is not identifiable here (score is flat)"
        elif abs(r["best"]) < 30.0:
            note = "consistent with the tubes being ALIGNED at zero"
        else:
            note = f"prefers psi0 = {r['best']:.0f} deg"
        print(f"  {key:7s} span {r['span']:5.2f} deg  ->  {note}")

    path = os.path.join(args.out_dir, "roll_calibration.png")
    draw(results, path)
    print(f"\nwrote {path}")

    print("\nsuggested ROLL_OFFSETS for run_validation.py:")
    print("ROLL_OFFSETS = {")
    for key, r in results.items():
        print(f'    "{key}": {r["best"]:.0f}.0,')
    print("}")
    return 0


def draw(results, path):
    vs.apply()
    keys = list(results)
    ncol = 3
    nrow = int(np.ceil(len(keys) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(12.6, 3.6 * nrow), squeeze=False)

    for ax, key in zip(axes.ravel(), keys):
        r = results[key]
        ax.plot(r["grid"], r["score"], color=vs.SERIES[0], lw=2.0, marker="o", ms=3.5,
                mec=vs.SURFACE, mew=0.8)
        ax.axvline(0.0, color=vs.NEUTRAL, lw=1.2, ls=(0, (4, 3)), zorder=1)
        ax.plot([r["best"]], [r["best_score"]], "o", ms=9, color=vs.SERIES[1],
                mec=vs.SURFACE, mew=1.5, zorder=5)
        ax.annotate(f"{r['best']:.0f}°", xy=(r["best"], r["best_score"]),
                    xytext=(6, 8), textcoords="offset points",
                    color=vs.SERIES[1], fontsize=9)
        ax.annotate("aligned", xy=(0.0, r["zero"]), xytext=(6, 6),
                    textcoords="offset points", color=vs.INK_2, fontsize=8.5)
        ax.set_title(f"{key}   (span {r['span']:.1f}°)")
        ax.set_xlabel("initial relative roll ψ₀ (deg)")
        ax.set_ylabel("frame-invariant mismatch (deg)")
        ax.set_xticks([-180, -90, 0, 90, 180])
        vs.tidy(ax)
    for ax in axes.ravel()[len(keys):]:
        ax.axis("off")

    fig.suptitle("Is the tubes' pre-curvature aligned when the rotation joints read zero?",
                 x=0.005, ha="left", fontsize=12, color=vs.INK, fontweight="semibold")
    vs.caption(fig, "Scored on turn angles and swept angles only — quantities a rigid "
                    "registration cannot flatter.  A flat curve means ψ₀ is simply not "
                    "identifiable from that set:\ntorsional windup absorbs a clamp-side "
                    "misalignment before it reaches the deployed section.  Dashed line = the "
                    "aligned assumption; dot = best fit.")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(path)
    plt.close(fig)


if __name__ == "__main__":
    sys.exit(main())
