#!/usr/bin/env python3
"""
run_validation.py -- score the Cosserat model against the NDI free-space trials.

    python3 validation/run_validation.py --fit        # identify geometry, then compare
    python3 validation/run_validation.py              # compare with the shipped config
    python3 validation/run_validation.py --sets set4a

Writes, into ``modeling/figures/validation/``:

    <trial>.png            one per experiment: measured vs modelled, registered
    curvature_profile.png  what the advances say about the tubes' real shape
    rotation_delivery.png  commanded vs delivered rotation, model and measurement
    summary_error.png      per-set registration residual, as-specified vs fitted
    metrics.csv            every number in the report

Method
------
Model and measurement live in unrelated frames -- the NDI tracker was never
registered to the robot's guide -- so shapes are compared after a rigid
(rotation + translation, no scaling) Kabsch alignment of the whole trial.  What
survives that alignment is real shape disagreement.  Alongside it we report
frame-invariant scalars (path length, fitted arc radius, swept angle) that need
no registration at all and so cannot be flattered by it.

Correspondence is by fraction of each commanded step completed, not global arc
length: the steps are separated by operator pauses, and one global parameter
would smear a disagreement in one step across the others.

Geometry identification uses **translation segments only**.  The rotation
segments are then a genuine prediction rather than a fit, which is what makes
them able to settle whether the 180 deg rotation shortfall is torsional windup
or curvature superposition.
"""

from __future__ import annotations

import argparse
import csv
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
from ctr import CTSDR, Tube  # noqa: E402
from ctr.model import CosseratModel  # noqa: E402
from ctr.robot import INNER, OUTER  # noqa: E402

import matplotlib.pyplot as plt  # noqa: E402
from mpl_toolkits.mplot3d import Axes3D  # noqa: E402,F401

OUT_DIR = os.path.join(_MODELING, "figures", "validation")

#: As shipped in config/ct_sdr.yaml: the MEASURED free-tube radius, 57 mm on
#: both tubes (2026-09-08), no lead-in, tips flush at zero.
AS_SPECIFIED = dict(R_outer=57.0, S_outer=0.0, R_inner=57.0,
                    home_outer=0.0, home_inner=0.0)

#: Identified from the data.  Two of these are well determined and two are not
#: -- see modeling/VALIDATION_REPORT.md before using them as a calibration.
#:
#: * ``R_inner`` 39.5 mm and ``home_inner`` 2.90 mm are pinned tightly and
#:   independently by `set2`, where the inner tube acts alone: together they
#:   reproduce that set's advance path length (34.39 vs 34.40 mm measured), its
#:   chord (33.33 vs 33.40) and its rotation circle radius (16.83 vs 16.83).
#: * ``R_outer`` and ``S_outer`` are the best fit to the both-tube advances but
#:   are *not* well determined: those profiles taper -- the assembly starts as
#:   curved as the inner tube alone and progressively straightens -- and no
#:   constant-curvature pair of tubes reproduces that shape at all (4.78 deg
#:   RMS against 5.25 deg as specified).
#: * ``home_outer`` is only weakly identifiable: once the outer tube covers the
#:   inner tube's tip, advancing it further changes nothing the tracker can
#:   see.  It is pinned to zero rather than left free because `set2` is run at
#:   OTT = 0 and is the calibration case for the two parameters above -- a
#:   non-zero value would leave a stub of outer tube stiffening the first few
#:   millimetres there and pull that set's rotation circle from 16.83 mm
#:   (measured, and matched exactly) down to 15.06 mm.
#:
#: Updated 2026-09-08: ``R_outer`` now takes the **measured** 57 mm rather than
#: the 53.63 mm fit, since that fit was the badly-determined one and a direct
#: measurement beats it.  Costs almost nothing -- 3.77 vs 3.74 mm mean RMS --
#: and replaces a fitted number with a measured one.  ``R_inner`` stays at the
#: data-identified 39.5 mm, which conflicts with the same measurement; §10 of
#: VALIDATION_REPORT.md lays out the conflict and the experiment that settles it.
FITTED = dict(R_outer=57.0, S_outer=0.0, R_inner=39.51,
              home_outer=0.0, home_inner=2.90)


# --- robot construction -----------------------------------------------------

def build_robot(R_outer: float, S_outer: float, R_inner: float,
                home_outer: float = 0.0, home_inner: float = 0.0,
                curved_mm: float = 78.54) -> CTSDR:
    """CT-SDR with the two tubes' shape and home offsets overridden.

    ``home_outer`` / ``home_inner`` are the lengths already protruding from
    the guide when the translational joints read zero [mm].  The shipped
    config leaves them ``null`` (meaning zero, tips flush with the guide
    exit); the NDI rotation circles say otherwise.
    """
    base = CTSDR.from_yaml()
    o, i = base.tubes
    base.tubes = [
        Tube(o.name, o.outer_diameter, o.inner_diameter, o.length,
             curved_mm * 1e-3, 1.0 / (R_outer * 1e-3), o.youngs_modulus,
             o.poisson_ratio, o.density, S_outer * 1e-3),
        Tube(i.name, i.outer_diameter, i.inner_diameter, i.length,
             curved_mm * 1e-3, 1.0 / (R_inner * 1e-3), i.youngs_modulus,
             i.poisson_ratio, i.density, 0.0),
    ]
    base.model = CosseratModel(base.tubes)
    base.base_offsets = np.array([-o.length + home_outer * 1e-3,
                                  -i.length + home_inner * 1e-3])
    return base


def rotation_circle_radius(robot, joints) -> float:
    """Radius the tip sweeps when the whole robot is rolled about the guide axis.

    For `set2` (inner tube alone) and `set3b` (both tubes co-rotating) the
    rotation is a rigid roll of the whole shape about +z, so the swept circle's
    radius is just the tip's distance from that axis -- one boundary value
    problem instead of a whole 360 deg sweep.  That is what makes the home
    offsets cheap to identify: the advances fix the curvatures and are blind to
    the offsets, while these radii depend on total deployed length and so fix
    the offsets.
    """
    sol = robot.solve(joints)
    tip = sol.position_at(robot.deployed_lengths(joints)[INNER])
    return float(np.hypot(tip[0], tip[1]) * 1e3)


# --- geometry identification ------------------------------------------------

#: Advance distances at which the tip path's accumulated turn is compared [mm].
PROFILE_SIGMA = np.array([12.0, 16.0, 20.0, 25.0, 30.0])


def turn_profile(p: np.ndarray, sigma=PROFILE_SIGMA) -> np.ndarray:
    """Accumulated turn of a tip path, in degrees, at each advance in ``sigma``.

    Measured by fitting a circle to the sub-arc ``[0, sigma]`` and taking
    ``sigma / R``, which uses every sample in the sub-arc rather than
    differentiating a noisy tangent.  ``nan`` where the path is shorter.
    """
    d = np.r_[0.0, np.cumsum(np.linalg.norm(np.diff(p, axis=0), axis=1))]
    out = np.full(len(sigma), np.nan)
    for k, s in enumerate(sigma):
        if s > d[-1]:
            continue
        r = X.arc_metrics(p[d <= s])["radius"]
        if r > 1e-9:
            out[k] = np.degrees(s / r)
    return out


def measured_profiles(trials):
    """Mean measured turn profile of each set's first advance."""
    out = {}
    for es in X.EXPERIMENTS:
        if es.steps[0].kind != "T":
            continue
        rows = [turn_profile(trials[label].segment(0))
                for label, _ in es.trials if trials[label].segments]
        out[es.key] = np.nanmean(np.array(rows), axis=0)
    return out


def _translation_residual(robot, trials, n_per_step=21, advance=X.SEQUENTIAL,
                          profiles=None):
    """RMS turn-angle mismatch over the first advance of each set [deg].

    Scored on the accumulated turn rather than on the registration residual.
    A rigid registration absorbs most of a radius error -- an R = 50 mm and an
    R = 63 mm arc over 35 mm differ by only ~0.6 mm of sagitta, so that
    objective is nearly flat and identifies the geometry poorly.  The turn
    profile is the same information without the rigid transform free to hide
    it, is linear in curvature, and is the quantity that actually separates
    candidate tube shapes.

    Step 1 only: it is the one step whose waypoints need no prior state, so no
    rotation sweep has to be simulated, and the six sets already span the
    informative cases -- the inner tube alone at 35 mm (`set2`) pins the inner
    tube, and both tubes at 17.5 / 20 / 35 mm pin the outer tube.
    """
    profiles = profiles if profiles is not None else measured_profiles(trials)
    err, count = 0.0, 0
    for es in X.EXPERIMENTS:
        if es.steps[0].kind != "T":
            continue
        p_model, slices, _ = X.model_trajectory(robot, es, n_per_step, advance, steps=[0])
        model = turn_profile(p_model[slices[0]])
        meas = profiles[es.key]
        ok = np.isfinite(model) & np.isfinite(meas)
        if not np.any(ok):
            continue
        err += float(np.sum((model[ok] - meas[ok]) ** 2))
        count += int(ok.sum())
    return float(np.sqrt(err / max(count, 1)))


def fit_geometry(trials, x0=(46.0, 20.0, 42.0), verbose=True):
    """Identify (R_outer, S_outer, R_inner) from the translation segments.

    Nelder-Mead on three parameters.  Only the advances are scored, so the
    rotation behaviour stays an out-of-sample prediction.
    """
    from scipy.optimize import minimize

    history = []
    profiles = measured_profiles(trials)

    def obj(x):
        R_o, S_o, R_i = x
        if not (20.0 < R_o < 400.0 and 0.0 <= S_o < 40.0 and 20.0 < R_i < 400.0):
            return 1e3
        try:
            r = _translation_residual(build_robot(R_o, S_o, R_i), trials,
                                      profiles=profiles)
        except Exception:
            return 1e3
        history.append((x.copy(), r))
        if verbose:
            print(f"    R_outer={R_o:6.2f}  S_outer={S_o:5.2f}  R_inner={R_i:6.2f}"
                  f"  ->  turn residual {r:.4f} deg", flush=True)
        return r

    res = minimize(obj, np.array(x0, dtype=float), method="Nelder-Mead",
                   options=dict(xatol=0.05, fatol=1e-3, maxiter=120))
    return dict(R_outer=float(res.x[0]), S_outer=float(res.x[1]),
                R_inner=float(res.x[2])), float(res.fun), history


# --- per-trial comparison ---------------------------------------------------

def compare_trial(tr, es, p_model, slices, n_per_seg=120):
    """Register one trial against the model and collect metrics."""
    meas, mod, n_seg = X.correspond(tr, p_model, slices, n_per_seg)
    R, t, rms = X.kabsch(mod, meas)
    mod_reg = mod @ R.T + t

    dev = np.linalg.norm(meas - mod_reg, axis=1)
    per_step = []
    for k in range(n_seg):
        sl = slice(k * n_per_seg, (k + 1) * n_per_seg)
        mm = X.arc_metrics(tr.segment(k))
        md = X.arc_metrics(p_model[slices[k]])
        per_step.append(dict(
            step=k + 1, kind=es.steps[k].kind, moves=es.steps[k].moves,
            meas=mm, model=md,
            dev_mean=float(dev[sl].mean()), dev_max=float(dev[sl].max()),
        ))
    return dict(meas=meas, model=mod_reg, dev=dev, rms=rms,
                n_seg=n_seg, per_step=per_step, n_per_seg=n_per_seg)


# --- figures ----------------------------------------------------------------

def fig_trial(tr, es, cmp_spec, cmp_fit, path):
    vs.apply()
    fig = plt.figure(figsize=(13.2, 7.4))
    gs = fig.add_gridspec(2, 3, height_ratios=[1.35, 1.0], hspace=0.42, wspace=0.30)

    best = cmp_fit if cmp_fit is not None else cmp_spec
    meas, mod = best["meas"], best["model"]
    allp = np.vstack([meas, mod])

    ax3d = fig.add_subplot(gs[0, 0], projection="3d")
    ax3d.plot(meas[:, 0], meas[:, 1], meas[:, 2], color=vs.SERIES[0], lw=2.0,
              label="measured (NDI)")
    ax3d.plot(mod[:, 0], mod[:, 1], mod[:, 2], color=vs.SERIES[1], lw=2.0,
              label="model (fitted)")
    ax3d.set_title("registered 3D tip path")
    ax3d.set_xlabel("x (mm)")
    ax3d.set_ylabel("y (mm)")
    ax3d.set_zlabel("z (mm)")
    vs.tidy3d(ax3d)
    vs.equal_aspect_3d(ax3d, allp)
    ax3d.view_init(elev=22, azim=-60)

    # Project onto the measured path's own best-fit plane: for these
    # trajectories that is the plane the motion actually happens in, so it
    # shows the disagreement rather than an arbitrary axis-aligned slice.
    c, u_ax, v_ax, normal, _ = X.ndi.fit_plane(meas)
    ax_p = fig.add_subplot(gs[0, 1])
    for p, colour, lab in ((meas, vs.SERIES[0], "measured"), (mod, vs.SERIES[1], "model")):
        u, v = X.ndi.project_plane(p, c, u_ax, v_ax)
        ax_p.plot(u, v, color=colour, lw=2.0, label=lab)
        ax_p.plot(u[0], v[0], "o", ms=5, color=colour, mec=vs.SURFACE, mew=1.2)
    ax_p.set_title("in-plane view")
    ax_p.set_xlabel("in-plane u (mm)")
    ax_p.set_ylabel("in-plane v (mm)")
    ax_p.set_aspect("equal", adjustable="datalim")
    vs.tidy(ax_p)

    ax_o = fig.add_subplot(gs[0, 2])
    for p, colour, lab in ((meas, vs.SERIES[0], "measured"), (mod, vs.SERIES[1], "model")):
        u, _ = X.ndi.project_plane(p, c, u_ax, v_ax)
        w = (p - c) @ normal
        ax_o.plot(u, w, color=colour, lw=2.0, label=lab)
    ax_o.set_title("out-of-plane view (note the finer w scale)")
    ax_o.set_xlabel("in-plane u (mm)")
    ax_o.set_ylabel("out-of-plane w (mm)")
    vs.tidy(ax_o)
    # Both projections are onto the measured path's own best-fit plane, so the
    # in-plane panel keeps a true 1:1 aspect and the out-of-plane panel does
    # not -- w is a fraction of a millimetre against tens in u.

    ax_d = fig.add_subplot(gs[1, :2])
    prog = np.arange(len(best["dev"]))
    for cmp_, colour, lab in ((cmp_spec, vs.NEUTRAL, "measured tubes (R = 57 mm)"),
                              (cmp_fit, vs.SERIES[2], "identified (57 out / 39.5 in)")):
        if cmp_ is None:
            continue
        ax_d.plot(prog, cmp_["dev"], color=colour, lw=2.0, label=lab)
    ax_d.set_xlim(0, len(prog) - 1)
    # Headroom so the step banner and the legend each get their own band and
    # neither lands on the curves.
    top = max(c["dev"].max() for c in (cmp_spec, cmp_fit) if c is not None)
    ax_d.set_ylim(0, top * 1.45)
    for k in range(1, best["n_seg"]):
        ax_d.axvline(k * best["n_per_seg"], color=vs.GRID, lw=1.2, zorder=0)
    for k in range(best["n_seg"]):
        mid = (k + 0.5) * best["n_per_seg"] / len(prog)
        st = es.steps[k]
        ax_d.annotate(f"{st.kind}{k + 1}  {st.moves}", xy=(mid, 0.98),
                      xycoords="axes fraction", ha="center", va="top",
                      color=vs.INK_2, fontsize=8.5)
    ax_d.set_title("model-to-measurement deviation after rigid registration")
    ax_d.set_xlabel("progress through the commanded sequence (samples, equal per step)")
    ax_d.set_ylabel("deviation (mm)")
    ax_d.legend(loc="upper left", bbox_to_anchor=(0.0, 0.88))
    vs.tidy(ax_d)

    ax_t = fig.add_subplot(gs[1, 2])
    ax_t.axis("off")
    rows = [f"{'step':<11s}{'quantity':<9s}{'meas':>8s}{'model':>8s}"]
    for ps in best["per_step"]:
        tag = f"{ps['kind']}{ps['step']} {ps['moves']}"[:10]
        if ps["kind"] == "T":
            rows.append(f"{tag:<11s}{'length':<9s}{ps['meas']['length']:8.1f}"
                        f"{ps['model']['length']:8.1f}")
            rows.append(f"{'':<11s}{'arc R':<9s}{ps['meas']['radius']:8.1f}"
                        f"{ps['model']['radius']:8.1f}")
        else:
            rows.append(f"{tag:<11s}{'swept':<9s}{ps['meas']['swept']:8.1f}"
                        f"{ps['model']['swept']:8.1f}")
            rows.append(f"{'':<11s}{'circ R':<9s}{ps['meas']['radius']:8.1f}"
                        f"{ps['model']['radius']:8.1f}")
        rows.append(f"{'':<11s}{'dev mean':<9s}{ps['dev_mean']:8.2f}{'':>8s}")
    rows.append("")
    if cmp_spec is not None:
        rows.append(f"{'RMS as specified':<20s}{cmp_spec['rms']:8.2f} mm")
    if cmp_fit is not None:
        rows.append(f"{'RMS fitted':<20s}{cmp_fit['rms']:8.2f} mm")
    ax_t.text(0.0, 1.0, "\n".join(rows), va="top", ha="left", fontsize=8.4,
              family="monospace", color=vs.INK, transform=ax_t.transAxes)
    ax_t.set_title("metrics (mm, deg)")

    handles, labels = ax_p.get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper right", ncol=2, bbox_to_anchor=(0.99, 1.02))
    fig.suptitle(f"{tr.label} - {es.name}   [{es.config}]", x=0.005, ha="left",
                 fontsize=12.5, color=vs.INK, fontweight="semibold")
    vs.caption(fig, "Cosserat model (Rucker et al. 2010) vs NDI tip tracking, rigidly "
                    "registered (rotation + translation, no scaling).")
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def fig_curvature_profile(trials, robots, path):
    """What the advances say about the tubes' real curvature distribution."""
    vs.apply()
    sig = np.array([12, 16, 20, 25, 30, 35])

    def nested(p):
        d = np.r_[0.0, np.cumsum(np.linalg.norm(np.diff(p, axis=0), axis=1))]
        return np.array([X.arc_metrics(p[d <= s])["radius"] if s <= d[-1] else np.nan
                         for s in sig])

    groups = {"both tubes advancing": [], "inner tube alone": []}
    for es in X.EXPERIMENTS:
        for k, st in enumerate(es.steps):
            if st.kind != "T":
                continue
            key = "both tubes advancing" if "OTT" in st.moves else "inner tube alone"
            for label, _ in es.trials:
                tr = trials[label]
                if len(tr.segments) > k:
                    groups[key].append(nested(tr.segment(k)))
    if not all(groups.values()):
        return False

    fig, axes = plt.subplots(1, 2, figsize=(11.4, 4.3))
    for ax, (title, data) in zip(axes, groups.items()):
        arr = np.array(data)
        mean = np.nanmean(arr, axis=0)
        lo, hi = np.nanmin(arr, axis=0), np.nanmax(arr, axis=0)
        ax.fill_between(sig, lo, hi, color=vs.SERIES[0], alpha=0.16, lw=0)
        ax.plot(sig, mean, color=vs.SERIES[0], lw=2.0, marker="o", ms=5,
                mec=vs.SURFACE, mew=1.2, label=f"measured (n={len(arr)})")
        for (name, robot), colour, style in zip(
                robots.items(), (vs.NEUTRAL, vs.SERIES[1]), ((0, (4, 3)), "-")):
            es_ref = [e for e in X.EXPERIMENTS
                      if e.key == ("set3a" if "both" in title else "set2")][0]
            p, sl, _ = X.model_trajectory(robot, es_ref, 41, steps=[0])
            ax.plot(sig, nested(p[sl[0]]), color=colour, lw=2.0, ls=style, label=name)
        ax.set_title(title)
        ax.set_xlabel("advance completed (mm)")
        ax.set_ylabel("mean radius of the tip path (mm)")
        ax.set_ylim(0, 110)
        vs.tidy(ax)
    axes[0].legend(loc="upper left")
    fig.suptitle("Radius of curvature the tip actually traces as the tubes advance",
                 x=0.005, ha="left", fontsize=12, color=vs.INK, fontweight="semibold")
    vs.caption(fig, "Nested sub-arc circle fits within each advance, so the trend is "
                    "within-trial and free of any cross-session confound.  A tube of "
                    "uniform curvature would give a flat line.")
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    fig.savefig(path)
    plt.close(fig)
    return True


def fig_rotation_delivery(records, path):
    """Commanded vs delivered rotation: the study's dominant error source."""
    vs.apply()
    fig, axes = plt.subplots(1, 2, figsize=(11.6, 4.4))
    ax_b, ax_s = axes

    order = ["set2", "set3a", "set3b", "set4b", "set4a"]
    labels, meas, spec, fit, cmds = [], [], [], [], []
    for key in order:
        rs = [r for r in records if r["set"] == key and r["kind"] == "R"]
        if not rs:
            continue
        es = [e for e in X.EXPERIMENTS if e.key == key][0]
        st = [s for s in es.steps if s.kind == "R"][0]
        cmd = max(abs(st.otr), abs(st.itr))
        labels.append(f"{key}\n{st.moves} {cmd:.0f}°")
        cmds.append(cmd)
        meas.append(np.mean([r["meas_swept"] for r in rs]) / cmd * 100)
        spec.append(np.mean([r["spec_swept"] for r in rs]) / cmd * 100)
        fit.append(np.mean([r["fit_swept"] for r in rs]) / cmd * 100)

    if not labels:  # --sets restricted the run to sets without a rotation step
        plt.close(fig)
        return False

    xs = np.arange(len(labels))
    w = 0.26
    for off, vals, colour, lab in ((-w, meas, vs.SERIES[0], "measured"),
                                   (0.0, spec, vs.NEUTRAL, "model, measured R = 57 mm"),
                                   (w, fit, vs.SERIES[1], "model, identified geometry")):
        ax_b.bar(xs + off, vals, width=w - 0.03, color=colour, label=lab, zorder=3)
    ax_b.axhline(100.0, color=vs.INK_2, lw=1.0, ls=(0, (4, 3)), zorder=2)
    ax_b.set_xticks(xs)
    ax_b.set_xticklabels(labels, fontsize=8)
    ax_b.set_ylabel("tip rotation delivered (% of commanded)")
    ax_b.set_title("rotation reaching the tip")
    ax_b.set_ylim(0, 128)  # room for the legend above the 100 % line
    ax_b.legend(loc="upper center", ncol=3, fontsize=8, bbox_to_anchor=(0.5, 1.0))
    vs.tidy(ax_b)

    ax_s.scatter(spec, meas, s=54, color=vs.NEUTRAL, zorder=3,
                 edgecolors=vs.SURFACE, linewidths=1.2, label="measured R = 57 mm")
    ax_s.scatter(fit, meas, s=54, color=vs.SERIES[1], zorder=4,
                 edgecolors=vs.SURFACE, linewidths=1.2, label="identified (57 out / 39.5 in)")
    lim = [0, max(max(meas), max(spec), max(fit)) * 1.12 + 5]
    ax_s.plot(lim, lim, color=vs.INK_2, lw=1.0, ls=(0, (4, 3)), zorder=2)
    for x_, y_, lab in zip(fit, meas, labels):
        ax_s.annotate(lab.split("\n")[0], xy=(x_, y_), xytext=(7, -3),
                      textcoords="offset points", fontsize=8.5, color=vs.INK_2)
    ax_s.set_xlim(lim)
    ax_s.set_ylim(lim)
    ax_s.set_xlabel("model prediction (% of commanded)")
    ax_s.set_ylabel("measured (% of commanded)")
    ax_s.set_title("prediction vs measurement (dashed = perfect)")
    ax_s.legend(loc="lower right")
    vs.tidy(ax_s)

    fig.suptitle("Rotation delivered to the tip: the dominant error source in the study",
                 x=0.005, ha="left", fontsize=12, color=vs.INK, fontweight="semibold")
    vs.caption(fig, "Rotation is an out-of-sample prediction: only the translation "
                    "segments were used to identify the model's geometry.")
    fig.tight_layout(rect=(0, 0, 1, 0.91))
    fig.savefig(path)
    plt.close(fig)
    return True


def fig_summary_error(records, path):
    vs.apply()
    fig, ax = plt.subplots(figsize=(9.2, 4.2))
    keys = [e.key for e in X.EXPERIMENTS]
    xs = np.arange(len(keys))
    w = 0.34
    for off, field, colour, lab in ((-w / 2, "spec_rms", vs.NEUTRAL, "measured tubes (R = 57 mm, no lead-in)"),
                                    (w / 2, "fit_rms", vs.SERIES[1], "identified (57 out / 39.5 in)")):
        vals, errs = [], []
        for k in keys:
            v = [r[field] for r in records if r["set"] == k and r["step"] == 1]
            vals.append(np.mean(v) if v else np.nan)
            errs.append(np.std(v) if v else 0.0)
        ax.bar(xs + off, vals, width=w - 0.04, yerr=errs, color=colour, label=lab,
               zorder=3, error_kw=dict(ecolor=vs.INK_2, lw=1.0, capsize=3))
    ax.set_xticks(xs)
    ax.set_xticklabels(keys)
    ax.set_ylabel("whole-trial registration RMS (mm)")
    ax.set_title("Model-to-measurement agreement per experiment set")
    ax.legend(loc="upper left")
    vs.tidy(ax)
    vs.caption(fig, "Bars are the mean over the 3 trials in each set; whiskers are the "
                    "trial-to-trial standard deviation.")
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


# --- driver -----------------------------------------------------------------

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fit", action="store_true",
                    help="re-identify the tube geometry from the translation segments")
    ap.add_argument("--sets", nargs="+", default=None, help="restrict to these set keys")
    ap.add_argument("--n-per-step", type=int, default=41)
    ap.add_argument("--advance", choices=[X.SEQUENTIAL, X.SIMULTANEOUS],
                    default=X.SEQUENTIAL)
    ap.add_argument("--out-dir", default=OUT_DIR)
    args = ap.parse_args(argv)

    sets = [e for e in X.EXPERIMENTS if args.sets is None or e.key in args.sets]
    os.makedirs(args.out_dir, exist_ok=True)

    print("loading NDI trials ...")
    trials = {}
    for es in X.EXPERIMENTS:
        for label, fn in es.trials:
            trials[label] = X.load_trial(label, fn, es.key)
    print(f"  {len(trials)} trials, "
          f"{sum(t.p.shape[0] for t in trials.values())} samples")

    params = dict(FITTED)
    if args.fit:
        print("\nidentifying tube geometry from translation segments only ...")
        params, best, _ = fit_geometry(trials)
        print(f"  -> R_outer={params['R_outer']:.2f} mm  "
              f"S_outer={params['S_outer']:.2f} mm  "
              f"R_inner={params['R_inner']:.2f} mm   residual {best:.3f} mm")

    robot_spec = build_robot(**AS_SPECIFIED)
    robot_fit = build_robot(**params)
    print(f"\nmeasured tubes: R_outer={AS_SPECIFIED['R_outer']:.1f}  "
          f"S_outer={AS_SPECIFIED['S_outer']:.1f}  R_inner={AS_SPECIFIED['R_inner']:.1f}")
    print(f"fitted       : R_outer={params['R_outer']:.1f}  "
          f"S_outer={params['S_outer']:.1f}  R_inner={params['R_inner']:.1f}")

    records = []
    print("\nper-trial comparison")
    print(f"  {'trial':8s}{'set':7s}{'RMS spec':>10s}{'RMS fit':>9s}   steps")
    for es in sets:
        p_spec, sl_spec, _ = X.model_trajectory(robot_spec, es, args.n_per_step, args.advance)
        p_fit, sl_fit, _ = X.model_trajectory(robot_fit, es, args.n_per_step, args.advance)
        for label, _ in es.trials:
            tr = trials[label]
            c_spec = compare_trial(tr, es, p_spec, sl_spec)
            c_fit = compare_trial(tr, es, p_fit, sl_fit)
            fname = label.replace(" ", "").lower() + ".png"
            fig_trial(tr, es, c_spec, c_fit, os.path.join(args.out_dir, fname))
            bits = " ".join(f"{p['kind']}{p['step']}:{p['dev_mean']:.1f}"
                            for p in c_fit["per_step"])
            print(f"  {label:8s}{es.key:7s}{c_spec['rms']:10.2f}{c_fit['rms']:9.2f}   {bits}")
            for ps_s, ps_f in zip(c_spec["per_step"], c_fit["per_step"]):
                records.append(dict(
                    set=es.key, set_name=es.name, trial=label,
                    step=ps_f["step"], kind=ps_f["kind"], moves=ps_f["moves"],
                    spec_rms=c_spec["rms"], fit_rms=c_fit["rms"],
                    meas_length=ps_f["meas"]["length"], model_length=ps_f["model"]["length"],
                    meas_radius=ps_f["meas"]["radius"], model_radius=ps_f["model"]["radius"],
                    spec_radius=ps_s["model"]["radius"],
                    meas_swept=ps_f["meas"]["swept"], fit_swept=ps_f["model"]["swept"],
                    spec_swept=ps_s["model"]["swept"],
                    meas_circle_resid=ps_f["meas"]["circle_resid"],
                    dev_mean_spec=ps_s["dev_mean"], dev_mean_fit=ps_f["dev_mean"],
                    dev_max_fit=ps_f["dev_max"],
                ))

    csv_path = os.path.join(args.out_dir, "metrics.csv")
    with open(csv_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(records[0].keys()))
        w.writeheader()
        for r in records:
            w.writerow({k: (f"{v:.4f}" if isinstance(v, float) else v)
                        for k, v in r.items()})

    print("\nwriting summary figures ...")
    made = {
        "curvature_profile.png": fig_curvature_profile(
            trials, {"model, measured R = 57 mm": robot_spec,
                     "model, identified geometry": robot_fit},
            os.path.join(args.out_dir, "curvature_profile.png")),
        "rotation_delivery.png": fig_rotation_delivery(
            records, os.path.join(args.out_dir, "rotation_delivery.png")),
        "summary_error.png": True,
    }
    fig_summary_error(records, os.path.join(args.out_dir, "summary_error.png"))
    for name, ok in made.items():
        if not ok:
            print(f"  skipped {name} (needs the full set of experiments)")

    spec_rms = np.mean([r["spec_rms"] for r in records if r["step"] == 1])
    fit_rms = np.mean([r["fit_rms"] for r in records if r["step"] == 1])
    print(f"\nmean whole-trial RMS:  as specified {spec_rms:.2f} mm   "
          f"fitted {fit_rms:.2f} mm")
    print(f"wrote {csv_path} and figures into {args.out_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
