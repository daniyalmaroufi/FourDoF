#!/usr/bin/env python3
"""
windup_vs_superposition.py -- why does a 180 deg ITR command move the tip 60 deg?

``ICRA2027/EVALUATION_REPORT.md`` reports that inner-tube rotation reaches the
tip almost perfectly with the outer tube retracted (97.7 % at OTT = 0) but only
33-42 % once the outer tube is advanced, and attributes the loss to *curvature
superposition* rather than torsional windup -- ruling windup out on the grounds
that the 360 deg trials bound it at 2-8 deg.

The Cosserat model contains both mechanisms, so it can separate them instead of
arguing from bounds.  Running it twice settles it:

* **full model** -- torsionally compliant tubes, as built.
* **torsionally rigid** -- the same tubes with ``GJ -> infinity``, so every
  degree commanded at the actuator arrives at the deployed section and *only*
  curvature superposition can move the answer.

If the rigid model already reproduces the shortfall, superposition is the
cause.  If it does not, the gap between the two runs is what windup contributes.

    python3 validation/windup_vs_superposition.py
"""

from __future__ import annotations

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
from ctr import Tube  # noqa: E402
from ctr.model import CosseratModel  # noqa: E402
from ctr.robot import INNER, OUTER  # noqa: E402
from run_validation import FITTED, OUT_DIR, build_robot  # noqa: E402

import matplotlib.pyplot as plt  # noqa: E402

#: Poisson ratio driving G = E / (2(1+nu)) to ~10^7 x its real value.
RIGID_NU = -0.999999


def make_rigid(robot):
    """Same robot with torsionally rigid tubes."""
    stiff = [
        Tube(t.name, t.outer_diameter, t.inner_diameter, t.length, t.curved_length,
             t.curvature, t.youngs_modulus, RIGID_NU, t.density, t.tip_straight_length)
        for t in robot.tubes
    ]
    rigid = build_robot(**FITTED)
    rigid.tubes = stiff
    rigid.model = CosseratModel(stiff)
    return rigid


def sweep(robot, es, n=41):
    """Swept tip angle over the set's rotation step, plus the roll actually
    delivered to the deployed section."""
    k = [i for i, s in enumerate(es.steps) if s.kind == "R"][0]
    # Steps after the rotation are irrelevant here and are the expensive ones.
    p, slices, sols = X.model_trajectory(robot, es, n, steps=range(k + 1))
    seg = p[slices[k]]
    swept = X.arc_metrics(seg)["swept"]

    end = sols[slices[k].stop - 1]
    if end is None:
        return swept, np.nan
    both = ~np.isnan(end.theta[:, OUTER])
    if not np.any(both):
        return swept, np.nan
    last = np.flatnonzero(both)[-1]
    delivered = np.degrees(end.theta[last, INNER] - end.theta[last, OUTER])
    return swept, delivered


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    print("loading NDI trials ...")
    trials = {}
    for es in X.EXPERIMENTS:
        for label, fn in es.trials:
            trials[label] = X.load_trial(label, fn, es.key)

    full = build_robot(**FITTED)
    rigid = make_rigid(full)

    rows = []
    for es in X.EXPERIMENTS:
        rot = [s for s in es.steps if s.kind == "R"]
        if not rot:
            continue
        k = [i for i, s in enumerate(es.steps) if s.kind == "R"][0]
        cmd = max(abs(rot[0].otr), abs(rot[0].itr))

        meas = []
        for label, _ in es.trials:
            tr = trials[label]
            if len(tr.segments) > k:
                meas.append(X.arc_metrics(tr.segment(k))["swept"])
        sw_full, rel_full = sweep(full, es)
        sw_rigid, rel_rigid = sweep(rigid, es)
        rows.append(dict(key=es.key, moves=rot[0].moves, cmd=cmd,
                         meas=float(np.mean(meas)), meas_sd=float(np.std(meas)),
                         full=sw_full, rigid=sw_rigid,
                         rel_full=rel_full, rel_rigid=rel_rigid))

    print(f"\n{'set':7s}{'rotates':9s}{'cmd':>6s}{'measured':>10s}"
          f"{'model':>8s}{'rigid':>8s}   {'roll at deployed section':>26s}")
    print(f"{'':7s}{'':9s}{'':>6s}{'swept deg':>10s}{'swept':>8s}{'swept':>8s}"
          f"   {'full':>12s}{'rigid':>13s}")
    for r in rows:
        print(f"{r['key']:7s}{r['moves']:9s}{r['cmd']:6.0f}"
              f"{r['meas']:8.1f}±{r['meas_sd']:.1f}{r['full']:8.1f}{r['rigid']:8.1f}"
              f"   {r['rel_full']:12.1f}{r['rel_rigid']:13.1f}")

    print("\nattribution of the shortfall (percentage points of the command):")
    for r in rows:
        loss_total = 100.0 * (r["cmd"] - r["meas"]) / r["cmd"]
        loss_super = 100.0 * (r["cmd"] - r["rigid"]) / r["cmd"]
        loss_wind = 100.0 * (r["rigid"] - r["full"]) / r["cmd"]
        print(f"  {r['key']:7s} measured loss {loss_total:6.1f}  |  "
              f"superposition {loss_super:6.1f}  windup {loss_wind:6.1f}  "
              f"unexplained {loss_total - loss_super - loss_wind:6.1f}")

    fig_attribution(rows, os.path.join(OUT_DIR, "windup_vs_superposition.png"))
    print(f"\nwrote {OUT_DIR}/windup_vs_superposition.png")


def fig_attribution(rows, path):
    vs.apply()
    fig, axes = plt.subplots(1, 2, figsize=(11.8, 4.4))
    ax_a, ax_b = axes

    keys = [f"{r['key']}\n{r['moves']} {r['cmd']:.0f}°" for r in rows]
    xs = np.arange(len(rows))
    super_loss = np.array([100.0 * (r["cmd"] - r["rigid"]) / r["cmd"] for r in rows])
    wind_loss = np.array([100.0 * (r["rigid"] - r["full"]) / r["cmd"] for r in rows])
    meas_loss = np.array([100.0 * (r["cmd"] - r["meas"]) / r["cmd"] for r in rows])

    ax_a.bar(xs, super_loss, width=0.56, color=vs.SERIES[1],
             label="curvature superposition", zorder=3)
    ax_a.bar(xs, wind_loss, width=0.56, bottom=super_loss, color=vs.SERIES[2],
             label="torsional windup", zorder=3)
    ax_a.plot(xs, meas_loss, "o", ms=9, color=vs.SERIES[0], mec=vs.SURFACE, mew=1.5,
              zorder=5, label="measured total loss")
    ax_a.set_xticks(xs)
    ax_a.set_xticklabels(keys, fontsize=8)
    ax_a.set_ylabel("rotation lost at the tip (% of commanded)")
    ax_a.set_title("what the model attributes the loss to")
    ax_a.legend(loc="upper left", fontsize=8.5)
    vs.tidy(ax_a)

    ott = {"set2": 0.0, "set3a": 35.0, "set3b": 35.0, "set4a": 35.0, "set4b": 17.5}
    inner = [r for r in rows if "ITR" in r["moves"] and "OTR" not in r["moves"]]
    xo = [ott[r["key"]] for r in inner]
    order = np.argsort(xo)
    xo = np.array(xo)[order]
    for vals, colour, lab, style in (
            ([r["meas"] / r["cmd"] * 100 for r in inner], vs.SERIES[0], "measured", "-"),
            ([r["rigid"] / r["cmd"] * 100 for r in inner], vs.SERIES[1],
             "model, torsionally rigid", (0, (4, 3))),
            ([r["full"] / r["cmd"] * 100 for r in inner], vs.SERIES[2],
             "model, full (with windup)", "-")):
        v = np.array(vals)[order]
        ax_b.plot(xo, v, color=colour, lw=2.0, ls=style, marker="o", ms=6,
                  mec=vs.SURFACE, mew=1.2, label=lab)
    ax_b.set_xlabel("outer tube extension at the time of rotation (mm)")
    ax_b.set_ylabel("inner-tube rotation reaching the tip (%)")
    ax_b.set_title("inner-tube rotation vs outer-tube overlap")
    ax_b.set_ylim(0, 115)
    ax_b.legend(loc="lower left", fontsize=8.5)
    vs.tidy(ax_b)

    fig.suptitle("Why a 180° inner-tube command moves the tip only ~60°",
                 x=0.005, ha="left", fontsize=12, color=vs.INK, fontweight="semibold")
    vs.caption(fig, "The torsionally rigid run isolates curvature superposition: every "
                    "commanded degree reaches the deployed section, so any remaining "
                    "shortfall is geometric.\nThe gap to the full model is what elastic "
                    "windup adds.")
    fig.tight_layout(rect=(0, 0, 1, 0.91))
    fig.savefig(path)
    plt.close(fig)


if __name__ == "__main__":
    main()
