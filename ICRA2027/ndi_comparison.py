#!/usr/bin/env python3
"""
ndi_comparison.py -- cross-experiment comparison and workspace analysis for the
4-DoF concentric tube robot free-space tests.

Where ndi_tip_analysis.py reports each experiment SET on its own, this script
compares all 18 individual trials against each other and characterises the
sampled workspace.  Outputs go to figures/comparison/ :

    per_experiment_errors.jpg       translation + rotation error, every trial
    per_experiment_geometry.jpg     tip-path radius and plane RMS, every trial
    per_experiment_stability.jpg    sensor noise, post-motion settling, drift
    per_experiment_timing.jpg       duration and speed profile summary
    repeatability.jpg               endpoint and whole-trajectory repeatability
    closure_error.jpg               360 deg loop-closure error
    speed_profiles.jpg              speed vs time, all trials, threshold shown
    workspace_3d.jpg                sampled workspace envelope + convex hull
    workspace_projections.jpg       XY / XZ / YZ envelope projections
    workspace_metrics.jpg           reach, bounding box and hull volume per set
    metric_heatmap.jpg              z-scored metric matrix, outliers at a glance

plus comparison_metrics.csv with one fully-populated row per trial.

Usage
-----
    python3 ndi_comparison.py
    python3 ndi_comparison.py --out-dir /tmp/cmp

Measurement conventions that differ from ndi_tip_analysis.py
------------------------------------------------------------
Repeatability is evaluated at END OF MOTION, not end of recording.  Trials keep
recording for a few seconds after the last commanded move, and the tip is still
settling during that dwell.  Including it charges post-motion creep to
repeatability: set2 scores 0.739 mm at end-of-record but 0.079 mm at end of
motion, because Exp 5 alone drifts 1.76 mm during a 10.6 s trailing dwell.
That settling is reported separately as its own metric.
"""

import argparse
import csv
import os
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
from matplotlib.lines import Line2D

from ndi_tip_analysis import (EXPERIMENT_SETS, analyse_trial, apply_style,
                              fit_plane, project_plane, fit_circle_2d,
                              path_length, speed_profile, _save, _finish,
                              _equal_3d, DEG)

try:
    from scipy.spatial import ConvexHull
except ImportError:
    ConvexHull = None


# One colour per experiment set (Wong palette, colour-blind safe).
SET_COLORS = {'set1': '#0072B2', 'set2': '#D55E00', 'set3a': '#009E73',
              'set3b': '#CC79A7', 'set4a': '#E69F00', 'set4b': '#56B4E9'}


# ─── Metrics ─────────────────────────────────────────────────────────────────
def resample_arclength(P, n=200):
    """Resample a path to n points evenly spaced in normalised arc length."""
    d = np.concatenate([[0.0], np.cumsum(np.linalg.norm(np.diff(P, axis=0), axis=1))])
    if d[-1] <= 0:
        return np.repeat(P[:1], n, axis=0)
    s = d / d[-1]
    q = np.linspace(0.0, 1.0, n)
    return np.column_stack([np.interp(q, s, P[:, k]) for k in range(3)])


def motion_only(r):
    """Concatenate just the moving spans of a trial."""
    return np.vstack([r['Ps'][a:b + 1] for a, b in r['moving']])


def collect(args):
    """Analyse every trial and build the per-experiment metric table."""
    per_set, rows = {}, []
    for es in EXPERIMENT_SETS:
        runs = [analyse_trial(l, f, es['steps'], args) for l, f in es['trials']]
        per_set[es['key']] = (es, runs)

        # Whole-trajectory repeatability: resample each trial's motion to a
        # common arc-length parameterisation, then measure spread about the
        # set's mean path.  Endpoint spread alone hides mid-path divergence.
        R = np.stack([resample_arclength(motion_only(r)) for r in runs])
        mean_path = R.mean(axis=0)
        traj_rms = [float(np.sqrt(np.mean(np.linalg.norm(x - mean_path, axis=1) ** 2)))
                    for x in R]
        ends = np.array([r['Ps'][r['moving'][-1][1]] for r in runs])
        end_dev = np.linalg.norm(ends - ends.mean(axis=0), axis=1)

        for r, trms, edev in zip(runs, traj_rms, end_dev):
            t, Ps = r['t'], r['Ps']

            # Tip-path radius over the first (bending) segment.
            a, b = r['moving'][0]
            S = Ps[a:b + 1]
            c, ua, va, _, _ = fit_plane(S)
            u, v = project_plane(S, c, ua, va)
            _, _, radius, rad_resid = fit_circle_2d(u, v)
            ang = np.unwrap(np.arctan2(v - fit_circle_2d(u, v)[1],
                                       u - fit_circle_2d(u, v)[0]))
            subtend = float(np.degrees(abs(ang[-1] - ang[0])))

            # Sensor noise from the pre-motion dwell: the robot is at rest and
            # not yet settling, so this is instrument noise rather than creep.
            lead = [s for s in r['stationary'] if s[0] == 0]
            if lead and (t[lead[0][1]] - t[lead[0][0]]) > 1.0:
                la, lb = lead[0]
                seg = r['P'][la:lb + 1]
                noise = float(np.std(np.linalg.norm(seg - seg.mean(axis=0), axis=1)))
            else:
                noise = float('nan')

            # Post-motion settling during the trailing dwell.
            tail = [s for s in r['stationary'] if s[0] > r['moving'][-1][1]]
            if tail:
                ta, tb = tail[0]
                settle = float(np.linalg.norm(Ps[tb] - Ps[ta]))
                settle_s = float(t[tb] - t[ta])
            else:
                settle, settle_s = float('nan'), float('nan')

            # Loop closure: only meaningful for a full 360 deg revolution,
            # where the tip is commanded back to where the rotation began.
            closure = float('nan')
            for i, (kind, cmd) in enumerate(es['steps']):
                if kind == 'R' and abs(cmd - 360.0) < 1e-6 and i < len(r['moving']):
                    ra_, rb_ = r['moving'][i]
                    closure = float(np.linalg.norm(Ps[rb_] - Ps[ra_]))

            t_err = [abs(m['error']) for m in r['measures'] if m['kind'] == 'T']
            r_err = [abs(m['error']) for m in r['measures'] if m['kind'] == 'R']
            spd = r['speed'][np.concatenate([np.arange(a2, b2 + 1)
                                             for a2, b2 in r['moving']])]

            rows.append(dict(
                set=es['key'], set_name=es['name'], trial=r['label'],
                n_samples=len(t), duration_s=float(t[-1]),
                n_segments=len(r['moving']),
                motion_time_s=float(sum(t[b2] - t[a2] for a2, b2 in r['moving'])),
                trans_err_mm=float(np.mean(t_err)) if t_err else float('nan'),
                trans_err_max_mm=float(np.max(t_err)) if t_err else float('nan'),
                rot_err_deg=float(np.mean(r_err)) if r_err else float('nan'),
                tip_radius_mm=radius, radius_resid_mm=rad_resid,
                arc_subtend_deg=subtend,
                plane_rms_mm=r['plane_rms'],
                noise_sigma_mm=noise,
                settle_mm=settle, settle_dwell_s=settle_s,
                closure_mm=closure,
                end_repeat_mm=float(edev), traj_repeat_rms_mm=trms,
                peak_speed_mms=float(np.percentile(spd, 99)),
                mean_speed_mms=float(np.mean(spd)),
                max_reach_mm=float(np.linalg.norm(Ps, axis=1).max()),
            ))
    return per_set, rows


# ─── Small plot helpers ──────────────────────────────────────────────────────
def _bar_by_trial(ax, rows, key, ylabel, title, annot='%.2f'):
    """One bar per trial, coloured by set, with set-boundary separators."""
    vals = [r[key] for r in rows]
    ok = [i for i, v in enumerate(vals) if not np.isnan(v)]
    if not ok:
        ax.axis('off')
        return
    idx = np.arange(len(ok))
    ax.bar(idx, [vals[i] for i in ok],
           color=[SET_COLORS[rows[i]['set']] for i in ok],
           edgecolor='0.25', lw=0.5)
    for j, i in enumerate(ok):
        ax.annotate(annot % vals[i], xy=(j, vals[i]), xytext=(0, 3),
                    textcoords='offset points', ha='center', fontsize=6)
    ax.set_xticks(idx)
    ax.set_xticklabels([rows[i]['trial'].replace('Exp ', '') for i in ok],
                       fontsize=7)
    ax.set_xlabel('Experiment number', fontsize=8)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.margins(y=0.16)
    ax.grid(axis='x', visible=False)
    _finish(ax)


def _set_legend(fig, keys, ncol=6, y=0.5):
    handles = [Line2D([0], [0], marker='s', ls='none', ms=7,
                      color=SET_COLORS[k], label=k) for k in keys]
    fig.legend(handles=handles, loc='lower center', ncol=ncol, fontsize=8,
               frameon=False, bbox_to_anchor=(0.5, y))


# ─── Figures ─────────────────────────────────────────────────────────────────
def fig_errors(rows, out):
    fig, axes = plt.subplots(2, 1, figsize=(9.5, 7.4))
    _bar_by_trial(axes[0], rows, 'trans_err_mm', 'Mean |error| (mm)',
                  'Translation Error by Experiment')
    _bar_by_trial(axes[1], rows, 'rot_err_deg', 'Mean |error| (deg)',
                  'Rotation Error by Experiment (rotation trials only)', '%.1f')
    _set_legend(fig, list(SET_COLORS), y=0.0)
    fig.tight_layout(rect=(0, 0.045, 1, 1))
    return _save(fig, out, 'per_experiment_errors.jpg')


def fig_geometry(rows, out):
    fig, axes = plt.subplots(2, 1, figsize=(9.5, 7.4))
    _bar_by_trial(axes[0], rows, 'tip_radius_mm', 'Radius (mm)',
                  'Tip-Path Radius over the Extension Segment '
                  '(nominal tube RoC = 50 mm, dashed)')
    axes[0].axhline(50.0, color='0.25', ls='--', lw=1.0, zorder=5)
    _bar_by_trial(axes[1], rows, 'plane_rms_mm', 'Plane RMS (mm)',
                  'Planarity of the Bending Phase (first 10 s of motion)',
                  '%.3f')
    _set_legend(fig, list(SET_COLORS), y=0.0)
    fig.tight_layout(rect=(0, 0.045, 1, 1))
    return _save(fig, out, 'per_experiment_geometry.jpg')


def fig_stability(rows, out):
    fig, axes = plt.subplots(2, 1, figsize=(9.5, 7.4))
    _bar_by_trial(axes[0], rows, 'noise_sigma_mm', r'$\sigma$ (mm)',
                  'Tracker Noise during the Pre-Motion Dwell', '%.3f')
    _bar_by_trial(axes[1], rows, 'settle_mm', 'Settling (mm)',
                  'Post-Motion Settling during the Trailing Dwell', '%.3f')
    _set_legend(fig, list(SET_COLORS), y=0.0)
    fig.tight_layout(rect=(0, 0.045, 1, 1))
    return _save(fig, out, 'per_experiment_stability.jpg')


def fig_timing(rows, out):
    fig, axes = plt.subplots(2, 1, figsize=(9.5, 7.4))
    _bar_by_trial(axes[0], rows, 'motion_time_s', 'Motion time (s)',
                  'Commanded Motion Duration by Experiment', '%.1f')
    _bar_by_trial(axes[1], rows, 'peak_speed_mms', 'Speed (mm/s)',
                  'Peak Tip Speed (99th percentile while moving)', '%.1f')
    _set_legend(fig, list(SET_COLORS), y=0.0)
    fig.tight_layout(rect=(0, 0.045, 1, 1))
    return _save(fig, out, 'per_experiment_timing.jpg')


def fig_repeatability(rows, out):
    fig, axes = plt.subplots(2, 1, figsize=(9.5, 7.4))
    _bar_by_trial(axes[0], rows, 'end_repeat_mm', 'Deviation (mm)',
                  'Endpoint Repeatability (deviation from set mean, '
                  'measured at end of motion)', '%.3f')
    _bar_by_trial(axes[1], rows, 'traj_repeat_rms_mm', 'RMS deviation (mm)',
                  'Whole-Trajectory Repeatability (RMS deviation from the '
                  'set mean path)', '%.3f')
    _set_legend(fig, list(SET_COLORS), y=0.0)
    fig.tight_layout(rect=(0, 0.045, 1, 1))
    return _save(fig, out, 'repeatability.jpg')


def fig_closure(rows, out):
    fig, ax = plt.subplots(figsize=(7.6, 4.2))
    _bar_by_trial(ax, rows, 'closure_mm', 'Closure error (mm)',
                  r'Loop-Closure Error after a Full 360$^\circ$ Revolution'
                  '\n(tip displacement between start and end of the '
                  'rotation segment)', '%.2f')
    _set_legend(fig, ['set2', 'set3a', 'set3b'], ncol=3, y=0.0)
    fig.tight_layout(rect=(0, 0.07, 1, 1))
    return _save(fig, out, 'closure_error.jpg')


def fig_itr_coupling(per_set, out):
    """Fraction of commanded rotation that reaches the tip, vs outer extension.

    This is the study's key mechanism result: inner-tube rotation is
    progressively lost as the outer tube is advanced, while outer-tube
    rotation is delivered in full regardless.
    """
    # (set, outer-tube extension at the time of the rotation, which tube turns)
    itr = [('set2', 0.0), ('set4b', 17.5), ('set4a', 35.0)]
    otr = [('set3a', 35.0), ('set3b', 35.0)]

    def achieved(key):
        es, runs = per_set[key]
        out_ = []
        for r in runs:
            m = [m for m in r['measures'] if m['kind'] == 'R'][0]
            out_.append(100.0 * m['measured'] / m['commanded'])
        return np.array(out_)

    fig, axes = plt.subplots(1, 2, figsize=(11.4, 4.6))

    ax = axes[0]
    x = [o for _, o in itr]
    y = [achieved(k).mean() for k, _ in itr]
    e = [achieved(k).std() for k, _ in itr]
    ax.errorbar(x, y, yerr=e, marker='o', ms=8, lw=1.6, capsize=4,
                color='#D55E00', label='Inner tube rotates (ITR)', zorder=3)
    for (k, o), v in zip(itr, y):
        cmd = per_set[k][0]['steps'][[i for i, (kk, _) in
                                      enumerate(per_set[k][0]['steps'])
                                      if kk == 'R'][0]][1]
        ax.annotate('%s\n(cmd %.0f%s)' % (k, cmd, DEG), xy=(o, v),
                    xytext=(6, -16), textcoords='offset points', fontsize=7.5)
    xo = [o for _, o in otr]
    yo = [achieved(k).mean() for k, _ in otr]
    ax.plot(xo, yo, 's', ms=8, color='#009E73',
            label='Outer tube rotates (OTR)', zorder=3)
    ax.annotate('set3a / set3b\n(cmd 360%s)' % DEG, xy=(xo[0], yo[0]),
                xytext=(-16, -26), textcoords='offset points', fontsize=7.5)
    ax.axhline(100, color='0.4', ls='--', lw=0.9)
    ax.annotate('ideal', xy=(0.99, 100), xycoords=('axes fraction', 'data'),
                xytext=(0, 4), textcoords='offset points', ha='right',
                fontsize=7, color='0.4')
    ax.set_xlabel('Outer-tube extension OTT at time of rotation (mm)')
    ax.set_ylabel('Commanded rotation reaching the tip (%)')
    ax.set_title('Rotation Transmission vs. Outer-Tube Extension')
    ax.set_ylim(0, 115)
    ax.set_xlim(-4, 42)
    ax.legend(fontsize=8, loc='center left')
    _finish(ax)

    # Angular rate rules out a speed-related dynamic lag: the poorly-tracked
    # cases were commanded more slowly, not faster.
    ax = axes[1]
    keys, rate, ach = [], [], []
    for key in ['set2', 'set3a', 'set3b', 'set4a', 'set4b']:
        es, runs = per_set[key]
        rr = []
        for r in runs:
            i = [j for j, (k2, _) in enumerate(es['steps']) if k2 == 'R'][0]
            a, b = r['moving'][i]
            m = [m for m in r['measures'] if m['step'] == i + 1][0]
            rr.append(m['measured'] / (r['t'][b] - r['t'][a]))
        keys.append(key); rate.append(np.mean(rr)); ach.append(achieved(key).mean())
    ax.scatter(rate, ach, s=90, c=[SET_COLORS[k] for k in keys],
               edgecolor='0.25', zorder=3)
    # set3a/set3b sit almost on top of each other; fan the labels out.
    offsets = {'set2': (8, -3), 'set3a': (-4, 11), 'set3b': (-6, -15),
               'set4a': (8, -3), 'set4b': (8, -3)}
    for k, x2, y2 in zip(keys, rate, ach):
        ax.annotate(k, xy=(x2, y2), xytext=offsets.get(k, (7, -3)),
                    textcoords='offset points', fontsize=8)
    ax.set_xlabel('Mean angular rate of the rotation segment (deg/s)')
    ax.set_ylabel('Commanded rotation reaching the tip (%)')
    ax.set_title('Transmission vs. Rotation Speed\n'
                 '(slower commands track worse, so speed is not the cause)',
                 fontsize=9.5)
    ax.set_ylim(0, 115)
    _finish(ax)

    fig.suptitle('Where Commanded Rotation Is Lost', fontsize=11.5)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    return _save(fig, out, 'rotation_transmission.jpg')


def fig_speed_profiles(per_set, args, out):
    n = len(per_set)
    fig, axes = plt.subplots(2, 3, figsize=(13.5, 7.0))
    axes = axes.ravel()
    for ax, (key, (es, runs)) in zip(axes, per_set.items()):
        for r, ls in zip(runs, ['-', '--', ':']):
            ax.plot(r['t'], r['speed'], ls, color=SET_COLORS[key], lw=0.9,
                    alpha=0.85, label=r['label'])
        ax.axhline(args.speed_thresh, color='0.3', lw=0.8, ls='-.')
        ax.annotate('threshold %.1f mm/s' % args.speed_thresh,
                    xy=(0.98, args.speed_thresh), xycoords=('axes fraction', 'data'),
                    xytext=(0, 4), textcoords='offset points',
                    ha='right', fontsize=6.5, color='0.3')
        ax.set_title('%s -- %s' % (key, es['name']), fontsize=8.5)
        ax.set_xlabel('Time (s)', fontsize=8)
        ax.set_ylabel('Speed (mm/s)', fontsize=8)
        ax.legend(fontsize=6.5, loc='upper right')
        _finish(ax)
    for ax in axes[len(per_set):]:
        ax.axis('off')
    fig.suptitle('Tip Speed Profiles -- Motion Segmentation Basis', fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    return _save(fig, out, 'speed_profiles.jpg')


def fig_workspace_3d(per_set, out):
    fig = plt.figure(figsize=(12.6, 6.0))
    ax = fig.add_subplot(121, projection='3d')
    allP = []
    for key, (es, runs) in per_set.items():
        P = np.vstack([r['Ps'] for r in runs])
        allP.append(P)
        ax.plot(P[:, 0], P[:, 1], P[:, 2], '.', ms=0.6, alpha=0.5,
                color=SET_COLORS[key], label=key)
    A = np.vstack(allP)
    ax.set_xlabel(r'$\Delta X$ (mm)'); ax.set_ylabel(r'$\Delta Y$ (mm)')
    ax.set_zlabel(r'$\Delta Z$ (mm)')
    ax.set_title('Sampled Tip Positions, All 18 Trials', fontsize=9.5)
    _equal_3d(ax, A); ax.view_init(elev=22, azim=-58); ax.tick_params(labelsize=6.5)
    ax.legend(fontsize=7, loc='upper left', markerscale=14)

    ax2 = fig.add_subplot(122, projection='3d')
    for key, (es, runs) in per_set.items():
        P = np.vstack([r['Ps'] for r in runs])
        ax2.plot(P[:, 0], P[:, 1], P[:, 2], '-', lw=0.5, alpha=0.35,
                 color=SET_COLORS[key])
    if ConvexHull is not None:
        from mpl_toolkits.mplot3d.art3d import Poly3DCollection
        hull = ConvexHull(A)
        # A wireframe of ~2800 simplices reads as noise; a translucent surface
        # shows the envelope shape instead.
        ax2.add_collection3d(Poly3DCollection(
            [A[s] for s in hull.simplices], facecolor='#7BA7C7', alpha=0.16,
            edgecolor='0.5', linewidths=0.15))
        ax2.set_title('Convex Hull of Sampled Envelope\n'
                      r'volume = %.0f mm$^3$, %d faces'
                      % (hull.volume, len(hull.simplices)), fontsize=9.5)
    else:
        ax2.set_title('Convex hull unavailable (scipy missing)', fontsize=9.5)
    ax2.set_xlabel(r'$\Delta X$ (mm)'); ax2.set_ylabel(r'$\Delta Y$ (mm)')
    ax2.set_zlabel(r'$\Delta Z$ (mm)')
    _equal_3d(ax2, A); ax2.view_init(elev=22, azim=-58); ax2.tick_params(labelsize=6.5)

    fig.suptitle('Sampled Workspace Envelope  (envelope of 18 recorded '
                 'trajectories, not the full reachable workspace)', fontsize=10.5)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    return _save(fig, out, 'workspace_3d.jpg')


def fig_workspace_projections(per_set, out):
    fig, axes = plt.subplots(1, 3, figsize=(13.2, 4.6))
    for ax, (i, j, li, lj) in zip(axes, [(0, 1, 'X', 'Y'), (0, 2, 'X', 'Z'),
                                         (1, 2, 'Y', 'Z')]):
        for key, (es, runs) in per_set.items():
            for r in runs:
                P = r['Ps']
                ax.plot(P[:, i], P[:, j], '-', lw=0.8, alpha=0.8,
                        color=SET_COLORS[key])
        ax.set_xlabel(r'$\Delta %s$ (mm)' % li)
        ax.set_ylabel(r'$\Delta %s$ (mm)' % lj)
        ax.set_title('%s%s Projection' % (li, lj))
        ax.set_aspect('equal', adjustable='datalim')
        _finish(ax)
    _set_legend(fig, list(SET_COLORS), y=0.0)
    fig.suptitle('Workspace Envelope Projections, All Trials', fontsize=11)
    fig.tight_layout(rect=(0, 0.07, 1, 0.94))
    return _save(fig, out, 'workspace_projections.jpg')


def fig_workspace_metrics(per_set, out):
    keys = list(per_set)
    reach, vol, bbox = [], [], []
    for key in keys:
        P = np.vstack([r['Ps'] for r in per_set[key][1]])
        reach.append(float(np.linalg.norm(P, axis=1).max()))
        bbox.append(np.ptp(P, axis=0))
        vol.append(ConvexHull(P).volume if ConvexHull is not None else np.nan)
    idx = np.arange(len(keys))
    fig, axes = plt.subplots(1, 3, figsize=(13.2, 4.3))

    axes[0].bar(idx, reach, 0.6, color=[SET_COLORS[k] for k in keys],
                edgecolor='0.25', lw=0.5)
    for i, v in enumerate(reach):
        axes[0].annotate('%.1f' % v, xy=(i, v), xytext=(0, 3),
                         textcoords='offset points', ha='center', fontsize=7)
    axes[0].set_ylabel('Max reach from start (mm)')
    axes[0].set_title('Maximum Tip Excursion')

    bbox = np.array(bbox)
    w = 0.26
    for k, (lab, col) in enumerate(zip(['X', 'Y', 'Z'],
                                       ['#0072B2', '#D55E00', '#009E73'])):
        axes[1].bar(idx + (k - 1) * w, bbox[:, k], w, label=r'$\Delta %s$' % lab,
                    color=col, edgecolor='0.25', lw=0.4)
    axes[1].set_ylabel('Extent (mm)')
    axes[1].set_title('Bounding-Box Extent')
    axes[1].legend(fontsize=7.5)

    axes[2].bar(idx, vol, 0.6, color=[SET_COLORS[k] for k in keys],
                edgecolor='0.25', lw=0.5)
    for i, v in enumerate(vol):
        axes[2].annotate('%.0f' % v, xy=(i, v), xytext=(0, 3),
                         textcoords='offset points', ha='center', fontsize=7)
    axes[2].set_ylabel(r'Convex hull volume (mm$^3$)')
    axes[2].set_title('Swept Volume\n(near-zero = planar trajectory)',
                      fontsize=9.5)

    for ax in axes:
        ax.set_xticks(idx)
        ax.set_xticklabels(keys, fontsize=8)
        ax.grid(axis='x', visible=False)
        _finish(ax)
    fig.suptitle('Workspace Metrics by Experiment Set', fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    return _save(fig, out, 'workspace_metrics.jpg')


def fig_heatmap(rows, out):
    """z-scored metric matrix: makes per-trial outliers pop out at a glance."""
    metrics = [('trans_err_mm', 'Trans. err'), ('rot_err_deg', 'Rot. err'),
               ('tip_radius_mm', 'Tip radius'), ('plane_rms_mm', 'Plane RMS'),
               ('noise_sigma_mm', 'Noise'), ('settle_mm', 'Settling'),
               ('end_repeat_mm', 'End repeat'),
               ('traj_repeat_rms_mm', 'Traj repeat'),
               ('closure_mm', 'Closure'), ('peak_speed_mms', 'Peak speed')]
    M = np.full((len(rows), len(metrics)), np.nan)
    for j, (k, _) in enumerate(metrics):
        col = np.array([r[k] for r in rows], dtype=float)
        good = ~np.isnan(col)
        if good.sum() > 1 and np.std(col[good]) > 0:
            M[good, j] = (col[good] - np.mean(col[good])) / np.std(col[good])
        elif good.sum():
            M[good, j] = 0.0

    # The global style turns grids on, which imshow and colorbar both warn
    # about; disable it for this figure rather than per-artist.
    with plt.rc_context({'axes.grid': False}):
        return _heatmap_figure(rows, metrics, M, out)


def _heatmap_figure(rows, metrics, M, out):
    fig, ax = plt.subplots(figsize=(10.0, 8.2))
    vmax = float(np.nanmax(np.abs(M)))
    im = ax.imshow(M, cmap='RdBu_r', vmin=-vmax, vmax=vmax, aspect='auto')
    ax.set_xticks(range(len(metrics)))
    ax.set_xticklabels([n for _, n in metrics], rotation=40, ha='right',
                       fontsize=8)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels(['%s  %s' % (r['set'], r['trial']) for r in rows],
                       fontsize=7.5)
    for i, r in enumerate(rows):
        for j, (k, _) in enumerate(metrics):
            v = r[k]
            if np.isnan(v):
                ax.text(j, i, '--', ha='center', va='center', fontsize=6,
                        color='0.6')
            else:
                fmt = '%.2f' if abs(v) < 100 else '%.0f'
                ax.text(j, i, fmt % v, ha='center', va='center', fontsize=6,
                        color='0.15' if abs(M[i, j]) < 1.6 else 'white')
    ax.set_xticks(np.arange(-0.5, len(metrics), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(rows), 1), minor=True)
    ax.grid(which='minor', color='white', lw=0.8)
    ax.grid(which='major', visible=False)
    ax.tick_params(which='minor', length=0)
    fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02,
                 label='z-score within metric (cell text = raw value)')
    ax.set_title('Per-Experiment Metric Matrix\n'
                 'colour = deviation from the column mean; '
                 'saturated cells are outliers', fontsize=10.5)
    fig.tight_layout()
    return _save(fig, out, 'metric_heatmap.jpg')


# ─── Driver ──────────────────────────────────────────────────────────────────
FIELDS = ['set', 'set_name', 'trial', 'n_samples', 'duration_s',
          'motion_time_s', 'n_segments', 'trans_err_mm', 'trans_err_max_mm',
          'rot_err_deg', 'tip_radius_mm', 'radius_resid_mm', 'arc_subtend_deg',
          'plane_rms_mm', 'noise_sigma_mm', 'settle_mm', 'settle_dwell_s',
          'closure_mm', 'end_repeat_mm', 'traj_repeat_rms_mm',
          'peak_speed_mms', 'mean_speed_mms', 'max_reach_mm']


def write_csv(rows, out):
    path = os.path.join(out, 'comparison_metrics.csv')
    with open(path, 'w', newline='') as fh:
        wr = csv.DictWriter(fh, fieldnames=FIELDS)
        wr.writeheader()
        for r in rows:
            wr.writerow({k: ('' if isinstance(r[k], float) and np.isnan(r[k])
                             else (('%.4f' % r[k]) if isinstance(r[k], float)
                                   else r[k]))
                         for k in FIELDS})
    return path


def summarise(rows):
    def col(k, sub=None):
        v = [r[k] for r in (sub or rows) if not np.isnan(r[k])]
        return np.array(v)
    print('\n=== Per-experiment comparison (18 trials) ===')
    hdr = ('%-6s %-7s %7s %7s %8s %7s %7s %8s %8s'
           % ('set', 'trial', 'transE', 'rotE', 'radius', 'noise', 'settle',
              'endRep', 'trajRep'))
    print(hdr); print('-' * len(hdr))
    for r in rows:
        f = lambda v, p='%7.3f': ('%7s' % '--') if np.isnan(v) else (p % v)
        print('%-6s %-7s %s %s %s %s %s %s %s'
              % (r['set'], r['trial'], f(r['trans_err_mm']),
                 f(r['rot_err_deg'], '%7.1f'), f(r['tip_radius_mm'], '%8.2f'),
                 f(r['noise_sigma_mm']), f(r['settle_mm']),
                 f(r['end_repeat_mm'], '%8.3f'),
                 f(r['traj_repeat_rms_mm'], '%8.3f')))
    print('-' * len(hdr))
    print('noise sigma      : %.4f +/- %.4f mm (max %.4f)'
          % (col('noise_sigma_mm').mean(), col('noise_sigma_mm').std(),
             col('noise_sigma_mm').max()))
    s = col('settle_mm')
    print('post-motion settle: median %.3f mm, max %.3f mm' % (np.median(s), s.max()))
    print('endpoint repeat  : mean %.3f mm, max %.3f mm'
          % (col('end_repeat_mm').mean(), col('end_repeat_mm').max()))
    print('traj repeat RMS  : mean %.3f mm, max %.3f mm'
          % (col('traj_repeat_rms_mm').mean(), col('traj_repeat_rms_mm').max()))
    print('closure (360 deg): mean %.2f mm, max %.2f mm'
          % (col('closure_mm').mean(), col('closure_mm').max()))


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    here = os.path.dirname(os.path.abspath(__file__))
    ap.add_argument('--data-dir', default=here)
    ap.add_argument('--out-dir', default=None,
                    help='default: <data-dir>/figures/comparison')
    ap.add_argument('--plane-window', type=float, default=10.0)
    ap.add_argument('--speed-thresh', type=float, default=0.3)
    ap.add_argument('--min-move-s', type=float, default=1.0)
    ap.add_argument('--min-pause-s', type=float, default=0.8)
    ap.add_argument('--smooth-s', type=float, default=0.5)
    ap.add_argument('--arc-step-s', type=float, default=0.2)
    ap.add_argument('--jump-mm', type=float, default=2.0)
    args = ap.parse_args()

    os.chdir(args.data_dir)
    if args.out_dir is None:
        args.out_dir = os.path.join(args.data_dir, 'figures', 'comparison')
    os.makedirs(args.out_dir, exist_ok=True)
    apply_style()

    if ConvexHull is None:
        print('warning: scipy unavailable -- hull volumes will be blank')

    per_set, rows = collect(args)
    out = args.out_dir

    fig_errors(rows, out)
    fig_geometry(rows, out)
    fig_stability(rows, out)
    fig_timing(rows, out)
    fig_repeatability(rows, out)
    fig_closure(rows, out)
    fig_itr_coupling(per_set, out)
    fig_speed_profiles(per_set, args, out)
    fig_workspace_3d(per_set, out)
    fig_workspace_projections(per_set, out)
    fig_workspace_metrics(per_set, out)
    fig_heatmap(rows, out)

    summarise(rows)
    print('\nwrote %s' % write_csv(rows, out))
    print("wrote 12 figures to %s/" % out)
    return 0


if __name__ == '__main__':
    sys.exit(main())
