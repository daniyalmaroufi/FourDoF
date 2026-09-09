#!/usr/bin/env python3
"""
ndi_tip_analysis.py -- NDI tip-tracking analysis for the 4-DoF concentric
tube robot (ICRA 2027 free-space tests).

Runs every experiment set defined in exp_info.md and writes one JPG per
result into figures/<set_name>/ :

    trajectory_3d.jpg               full 3D tip trajectory, all trials
    projection_xy.jpg               top view
    projection_xz.jpg               side view
    plane_fit.jpg                   plane fitted to the first N s of motion
    position_vs_time.jpg            full record, motion/pause segments marked
    position_vs_time_truncated.jpg  stationary spans removed
    trajectory_3d_segments.jpg      one 3D subplot per motion segment
    kinematic_error.jpg             commanded vs measured translation/rotation

plus, at the top level, summary_kinematic_error.jpg and kinematic_results.csv.

Usage
-----
    python3 ndi_tip_analysis.py                  # everything, into figures/
    python3 ndi_tip_analysis.py --sets set1 set2
    python3 ndi_tip_analysis.py --plane-window 5

Method notes
------------
* Positions are converted to millimetres and referenced to each trial's first
  sample, so every trace starts at the origin.

* Spike rejection is discontinuity-based, not magnitude-based: a sample is
  dropped only if it is more than --jump-mm from BOTH neighbours.  A
  median/MAD filter on absolute position would flag the far end of a real
  trajectory as an outlier and silently truncate the motion.

* Motion segmentation thresholds a smoothed speed profile.  The measured
  noise floor is ~0.02-0.09 mm/s against 2.8-5.9 mm/s while moving, so the
  0.3 mm/s default sits comfortably between the two.

* Translation of a segment is the tip's PATH LENGTH: as a pre-curved tube is
  advanced, the tip travels along the tube's own curve, so arc length -- not
  straight-line displacement -- corresponds to the commanded advance.  It is
  accumulated on the smoothed path in ~0.2 s steps so per-sample jitter does
  not integrate into a spurious length.

* Rotation of a segment is measured by fitting a plane, then a circle, to the
  segment and unwrapping the angle about the circle centre.  The circle-fit
  residual is reported alongside every angle so a poor fit is visible rather
  than hidden behind a confident-looking number.
"""

import argparse
import csv
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
from matplotlib.lines import Line2D


# ─── Style ───────────────────────────────────────────────────────────────────
# Computer Modern is metrically the same design as Latin Modern; matplotlib
# ships it, so figures match an IEEE LaTeX template without needing a TeX
# install.  Set text.usetex below if you do have one and want true LM.
def apply_style():
    plt.rcParams.update({
        # 'text.usetex': True,
        # 'text.latex.preamble': r'\usepackage{lmodern}',
        'font.family':                 'serif',
        'font.serif':                  ['cmr10', 'STIXGeneral', 'DejaVu Serif'],
        'mathtext.fontset':            'cm',
        'axes.formatter.use_mathtext': True,
        'font.size':                   9,
        'axes.labelsize':              9,
        'axes.titlesize':              10,
        'xtick.labelsize':             8,
        'ytick.labelsize':             8,
        'legend.fontsize':             8,
        'legend.framealpha':           0.88,
        'legend.edgecolor':            '0.75',
        'figure.dpi':                  100,
        'savefig.dpi':                 300,
        'savefig.bbox':                'tight',
        'savefig.pad_inches':          0.06,
        'axes.linewidth':              0.65,
        'grid.linewidth':              0.35,
        'grid.alpha':                  0.45,
        'axes.grid':                   True,
    })


# The two report-hero figures get their text set at a multiple of the base
# sizes so they stay legible in a paper column or on a slide.  What matters for
# apparent size is font_scale / canvas_scale: enlarging both equally just
# renders the same figure bigger.  Overview is at 4.5/1.0 = 4.5x apparent,
# summary at 3.0/2.0 = 1.5x.
#
# At 4.5x apparent, text only fits if the panels are wide relative to the
# label lengths: a title fits when chars <= ~78 * panel_inches / (4.6 * scale
# ratio), which is why the overview uses a 2-column grid with short titles and
# one shared legend rather than 3 columns with six large per-panel legends.
OVERVIEW_FONT_SCALE = 4.5
OVERVIEW_CANVAS_SCALE = 1.0
SUMMARY_FONT_SCALE = 3.0
SUMMARY_CANVAS_SCALE = 2.0


def _scaled_font_rc(scale):
    """rcParams overrides scaling every text size by `scale`."""
    return {
        'font.size':        9 * scale,
        'axes.labelsize':   9 * scale,
        'axes.titlesize':  10 * scale,
        'xtick.labelsize':  8 * scale,
        'ytick.labelsize':  8 * scale,
        'legend.fontsize':  8 * scale,
    }


# Wong (2011) colour-blind-safe palette; also separable in greyscale print.
COLORS = ['#0072B2', '#D55E00', '#009E73']
# Repeat trials often coincide to ~0.1 mm, so taper the widths: without this
# only the last-drawn trial is visible.
LWS = [2.1, 1.35, 0.75]
DEG = r'$^\circ$'


# ─── Experiment definitions (from exp_info.md) ───────────────────────────────
# steps: ordered commanded motions, ('T', mm) translation or ('R', deg) rotation.
# Exp 4 is excluded ("Random Exp Don't count").  Where a run was recorded
# twice the aborted/empty file is skipped -- see `note` on each set.
EXPERIMENT_SETS = [
    dict(
        key='set1', name='1st Free Space Test',
        config='OTT = 20 mm, ITT = 20 mm  |  OTR = 0' + DEG + ', ITR = 0' + DEG,
        steps=[('T', 20.0)],
        trials=[('Exp 1', 'exp1_20260811_215059.csv'),
                ('Exp 2', 'exp2_20260811_215200.csv'),
                ('Exp 3', 'exp3_20260811_215253.csv')],
        note='',
    ),
    dict(
        key='set2', name='2nd Free Space Test',
        config='OTT = 0 mm, ITT = 35 mm  |  OTR = 0' + DEG + ', ITR = 360' + DEG,
        steps=[('T', 35.0), ('R', 360.0)],
        trials=[('Exp 5', 'exp5_20260811_222950.csv'),
                ('Exp 6', 'exp6_20260811_223714.csv'),
                ('Exp 7', 'exp7_20260811_224018.csv')],
        note='No outer tube movement.',
    ),
    dict(
        key='set3a', name='3rd Free Space Test (Top Right)',
        config='OTT = 35 mm, ITT = 35 mm  |  OTR = 360' + DEG + ', ITR = 0' + DEG,
        steps=[('T', 35.0), ('R', 360.0)],
        trials=[('Exp 8',  'exp8_20260811_230128.csv'),
                ('Exp 9',  'exp9_20260811_230810.csv'),
                ('Exp 10', 'exp10_20260811_231038.csv')],
        note='exp8_...225935.csv is a 5.6 s aborted run (0.05 mm net); '
             'exp8_...230128.csv used instead.',
    ),
    dict(
        key='set3b', name='3rd Free Space Test (Bottom Left)',
        config='OTT = 35 mm, ITT = 35 mm  |  OTR = 360' + DEG + ', ITR = 360' + DEG,
        steps=[('T', 35.0), ('R', 360.0)],
        trials=[('Exp 11', 'exp11_20260811_232209.csv'),
                ('Exp 12', 'exp12_20260811_232456.csv'),
                ('Exp 13', 'exp13_20260811_232803.csv')],
        note='',
    ),
    dict(
        key='set4a', name='4th Free Space Test (Bottom Mid)',
        config='OTT = 35 mm, ITT = 35 mm  |  ITR = 180' + DEG + '  |  ITT += 35 mm',
        steps=[('T', 35.0), ('R', 180.0), ('T', 35.0)],
        trials=[('Exp 14', 'exp14_20260811_233834.csv'),
                ('Exp 15', 'exp15_20260811_234353.csv'),
                ('Exp 16', 'exp16_20260811_235022.csv')],
        note='',
    ),
    dict(
        key='set4b', name='4th Free Space Test (Bottom Right)',
        config='OTT = 17.5 mm, ITT = 17.5 mm  |  ITR = 180' + DEG
               + '  |  ITT += 17.5 mm',
        steps=[('T', 17.5), ('R', 180.0), ('T', 17.5)],
        trials=[('Exp 17', 'exp17_20260812_001514.csv'),
                ('Exp 18', 'exp18_20260812_010331.csv'),
                ('Exp 19', 'exp19_20260812_010534.csv')],
        note='exp18_...010107.csv contains only a header; '
             'exp18_...010331.csv used instead.  '
             '"exp17 ... copy.csv" is byte-identical to exp17.',
    ),
]


# ─── Loading and conditioning ────────────────────────────────────────────────
def load_csv(path):
    """Read a /tip_xyz recording.  Returns (t [s], P [n x 3] mm)."""
    t, xs, ys, zs = [], [], [], []
    with open(path) as fh:
        reader = csv.DictReader(fh)
        missing = {'elapsed', 'x', 'y', 'z'} - set(reader.fieldnames or [])
        if missing:
            raise SystemExit('%s missing column(s): %s'
                             % (path, ', '.join(sorted(missing))))
        for row in reader:
            t.append(float(row['elapsed']))
            xs.append(float(row['x']))
            ys.append(float(row['y']))
            zs.append(float(row['z']))
    if len(t) < 2:
        raise SystemExit('%s contains no usable samples' % path)
    return np.asarray(t), np.column_stack([xs, ys, zs]) * 1000.0


def drop_spikes(t, P, jump_mm=2.0):
    """Remove tracker dropouts: samples far from BOTH neighbours."""
    keep = np.ones(len(t), dtype=bool)
    if len(t) >= 3:
        d_prev = np.linalg.norm(P[1:-1] - P[:-2], axis=1)
        d_next = np.linalg.norm(P[2:] - P[1:-1], axis=1)
        keep[1:-1] = ~((d_prev > jump_mm) & (d_next > jump_mm))
    return t[keep], P[keep], int((~keep).sum())


def smooth(t, P, win_s=0.5):
    """Centred moving average with a window given in seconds."""
    dt = float(np.median(np.diff(t)))
    w = max(3, int(round(win_s / dt)) | 1)          # force odd
    out = np.empty_like(P)
    pad = w // 2
    for k in range(P.shape[1]):
        padded = np.pad(P[:, k], (pad, pad), mode='edge')
        out[:, k] = np.convolve(padded, np.ones(w) / w, mode='valid')
    return out


def speed_profile(t, P, win_s=0.3):
    """Centred-difference speed in mm/s, smoothed over +/- win_s."""
    dt = float(np.median(np.diff(t)))
    w = max(1, int(round(win_s / dt)))
    n = len(t)
    i0 = np.clip(np.arange(n) - w, 0, n - 1)
    i1 = np.clip(np.arange(n) + w, 0, n - 1)
    return np.linalg.norm(P[i1] - P[i0], axis=1) / np.maximum(t[i1] - t[i0], 1e-9)


# ─── Motion segmentation ─────────────────────────────────────────────────────
def _runs(mask):
    """Contiguous runs of equal value as (start, end_inclusive, value)."""
    out, start = [], 0
    for i in range(1, len(mask)):
        if mask[i] != mask[i - 1]:
            out.append((start, i - 1, bool(mask[start])))
            start = i
    out.append((start, len(mask) - 1, bool(mask[start])))
    return out


def _flip_short(mask, t, value, min_s):
    """Flip runs of `value` shorter than min_s (de-fragments the mask)."""
    out = mask.copy()
    for a, b, v in _runs(mask):
        if v == value and (t[b] - t[a]) < min_s:
            out[a:b + 1] = not value
    return out


def segment_motion(t, P, thresh=0.3, min_move_s=1.0, min_pause_s=0.8):
    """Split a trial into moving / stationary spans.

    Returns (moving_segments, stationary_segments, speed), each segment an
    inclusive (start_idx, end_idx) pair.
    """
    v = speed_profile(t, P)
    moving = v > thresh
    # Drop noise blips, bridge momentary dips mid-motion, drop blips again.
    moving = _flip_short(moving, t, True, min_move_s)
    moving = _flip_short(moving, t, False, min_pause_s)
    moving = _flip_short(moving, t, True, min_move_s)
    mov = [(a, b) for a, b, val in _runs(moving) if val]
    sta = [(a, b) for a, b, val in _runs(moving) if not val]
    return mov, sta, v


# ─── Geometry ────────────────────────────────────────────────────────────────
def fit_plane(P):
    """SVD best-fit plane.  Returns (centroid, u_axis, v_axis, normal, rms)."""
    c = P.mean(axis=0)
    _, _, vt = np.linalg.svd(P - c, full_matrices=False)
    u_ax, v_ax, normal = vt[0], vt[1], vt[2]
    rms = float(np.sqrt(np.mean(((P - c) @ normal) ** 2)))
    return c, u_ax, v_ax, normal, rms


def project_plane(P, centroid, u_ax, v_ax):
    d = P - centroid
    return d @ u_ax, d @ v_ax


def fit_circle_2d(u, v):
    """Algebraic (Kasa) circle fit.  Returns (cu, cv, radius, rms_residual)."""
    A = np.column_stack([2 * u, 2 * v, np.ones(len(u))])
    sol, *_ = np.linalg.lstsq(A, u ** 2 + v ** 2, rcond=None)
    cu, cv = float(sol[0]), float(sol[1])
    r = float(np.sqrt(max(sol[2] + cu ** 2 + cv ** 2, 0.0)))
    resid = float(np.sqrt(np.mean((np.hypot(u - cu, v - cv) - r) ** 2)))
    return cu, cv, r, resid


def path_length(t, P, step_s=0.2):
    """Arc length accumulated in ~step_s hops.

    Summing every sample would integrate per-sample jitter into tens of mm of
    phantom length over a long record; hopping keeps real curvature while
    averaging the noise away.
    """
    dt = float(np.median(np.diff(t)))
    k = max(1, int(round(step_s / dt)))
    Q = P[::k]
    if len(Q) < 2:
        return 0.0
    return float(np.linalg.norm(np.diff(Q, axis=0), axis=1).sum())


def measure_rotation(P_seg):
    """Swept angle (deg) of a segment about its fitted circle centre."""
    c, u_ax, v_ax, _, plane_rms = fit_plane(P_seg)
    u, v = project_plane(P_seg, c, u_ax, v_ax)
    cu, cv, r, resid = fit_circle_2d(u, v)
    ang = np.unwrap(np.arctan2(v - cv, u - cu))
    return float(np.degrees(abs(ang[-1] - ang[0]))), r, resid, plane_rms


# ─── Per-trial analysis ──────────────────────────────────────────────────────
def analyse_trial(label, path, steps, args):
    t, P, n_spikes = drop_spikes(*load_csv(path), jump_mm=args.jump_mm)
    P = P - P[0]
    t = t - t[0]
    Ps = smooth(t, P, args.smooth_s)

    mov, sta, v = segment_motion(t, Ps, args.speed_thresh,
                                 args.min_move_s, args.min_pause_s)

    # Plane fit over the first `plane_window` seconds of motion.  Measuring
    # from motion onset rather than t=0 matters: several trials idle for 5-10 s
    # before moving, and a window of pure dwell yields a degenerate plane.
    t_start = t[mov[0][0]] if mov else t[0]
    win = (t >= t_start) & (t <= t_start + args.plane_window)
    if win.sum() >= 3:
        c, u_ax, v_ax, normal, plane_rms = fit_plane(Ps[win])
    else:
        c, u_ax, v_ax, normal, plane_rms = fit_plane(Ps)
    pu, pv = project_plane(Ps[win], c, u_ax, v_ax)

    # Pair detected segments with commanded steps, in order.
    measures = []
    n_pair = min(len(mov), len(steps))
    for i in range(n_pair):
        a, b = mov[i]
        kind, commanded = steps[i]
        seg = Ps[a:b + 1]
        if kind == 'T':
            measured = path_length(t[a:b + 1], seg, args.arc_step_s)
            radius = resid = float('nan')
        else:
            measured, radius, resid, _ = measure_rotation(seg)
        measures.append(dict(step=i + 1, kind=kind, commanded=commanded,
                             measured=measured, error=measured - commanded,
                             radius=radius, resid=resid,
                             t0=t[a], t1=t[b]))

    return dict(label=label, t=t, P=P, Ps=Ps, speed=v,
                moving=mov, stationary=sta, n_spikes=n_spikes,
                plane=(c, u_ax, v_ax, normal, plane_rms),
                plane_uv=(pu, pv), plane_rms=plane_rms,
                plane_t0=t_start, measures=measures,
                seg_mismatch=(len(mov) != len(steps)))


# ─── Plot helpers ────────────────────────────────────────────────────────────
def _finish(ax):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)


def _save(fig, out_dir, name):
    path = os.path.join(out_dir, name)
    fig.savefig(path, dpi=300, bbox_inches='tight',
                pil_kwargs={'quality': 95, 'optimize': True})
    plt.close(fig)
    return path


def _equal_3d(ax, P):
    half = max(float(np.ptp(P, axis=0).max()), 1e-6) / 2.0
    mid = (P.max(axis=0) + P.min(axis=0)) / 2.0
    ax.set_xlim(mid[0] - half, mid[0] + half)
    ax.set_ylim(mid[1] - half, mid[1] + half)
    ax.set_zlim(mid[2] - half, mid[2] + half)
    try:
        ax.set_box_aspect((1, 1, 1))
    except Exception:
        pass
    for pane in (ax.xaxis.pane, ax.yaxis.pane, ax.zaxis.pane):
        pane.fill = False
        pane.set_edgecolor('0.85')


def _marker_handles():
    return [Line2D([0], [0], ls='none', marker='o', ms=6, mfc='w', mec='k',
                   mew=1.2, label='Start'),
            Line2D([0], [0], ls='none', marker='s', ms=5, color='k',
                   label='End')]


# ─── Figures ─────────────────────────────────────────────────────────────────
def fig_trajectory_3d(runs, es, out_dir):
    fig = plt.figure(figsize=(5.6, 5.0))
    ax = fig.add_subplot(111, projection='3d')
    for r, c, w in zip(runs, COLORS, LWS):
        P = r['Ps']
        ax.plot(P[:, 0], P[:, 1], P[:, 2], color=c, lw=w, alpha=0.9,
                label=r['label'])
        ax.scatter(*P[0],  marker='o', s=55, facecolors='white',
                   edgecolors=c, linewidths=1.6, zorder=6)
        ax.scatter(*P[-1], marker='s', s=50, color=c, zorder=6)
    ax.set_xlabel(r'$\Delta X$ (mm)')
    ax.set_ylabel(r'$\Delta Y$ (mm)')
    ax.set_zlabel(r'$\Delta Z$ (mm)')
    ax.set_title('%s\n3D Tip Trajectory' % es['name'])
    _equal_3d(ax, np.vstack([r['Ps'] for r in runs]))
    ax.view_init(elev=22, azim=-58)
    ax.tick_params(labelsize=7)
    handles = [Line2D([0], [0], color=c, lw=w, label=r['label'])
               for r, c, w in zip(runs, COLORS, LWS)] + _marker_handles()
    ax.legend(handles=handles, fontsize=7, loc='upper left')
    return _save(fig, out_dir, 'trajectory_3d.jpg')


def fig_projection(runs, es, out_dir, i, j, li, lj, name):
    fig, ax = plt.subplots(figsize=(5.2, 4.4))
    for r, c, w in zip(runs, COLORS, LWS):
        P = r['Ps']
        ax.plot(P[:, i], P[:, j], color=c, lw=w, alpha=0.9, label=r['label'])
        ax.plot(P[0, i],  P[0, j],  'o', ms=7, mfc='white', mec=c, mew=1.5, zorder=5)
        ax.plot(P[-1, i], P[-1, j], 's', ms=6, color=c, zorder=5)
    ax.set_xlabel(r'$\Delta %s$ (mm)' % li)
    ax.set_ylabel(r'$\Delta %s$ (mm)' % lj)
    ax.set_title('%s\n%s%s Projection' % (es['name'], li, lj))
    ax.set_aspect('equal', adjustable='datalim')
    ax.legend(handles=[Line2D([0], [0], color=c, lw=w, label=r['label'])
                       for r, c, w in zip(runs, COLORS, LWS)] + _marker_handles(),
              fontsize=7, loc='best')
    _finish(ax)
    return _save(fig, out_dir, name)


def fig_plane_fit(runs, es, out_dir, window):
    fig, ax = plt.subplots(figsize=(5.2, 4.4))
    for r, c, w in zip(runs, COLORS, LWS):
        u, v = r['plane_uv']
        ax.plot(u, v, color=c, lw=w, alpha=0.9,
                label='%s (RMS = %.3f mm)' % (r['label'], r['plane_rms']))
        ax.plot(u[0],  v[0],  'o', ms=7, mfc='white', mec=c, mew=1.5, zorder=5)
        ax.plot(u[-1], v[-1], 's', ms=6, color=c, zorder=5)
    ax.set_xlabel(r'$u$ (mm)')
    ax.set_ylabel(r'$v$ (mm)')
    ax.set_title('%s\nBest-Fit Plane, First %.0f s of Motion'
                 % (es['name'], window))
    ax.set_aspect('equal', adjustable='datalim')
    ax.legend(fontsize=7, loc='best')
    _finish(ax)
    return _save(fig, out_dir, 'plane_fit.jpg')


def _draw_axes_vs_time(ax, t, P, c, w):
    ax.plot(t, P[:, 0], '-',  color=c, lw=w * 0.8, alpha=0.9)
    ax.plot(t, P[:, 1], '--', color=c, lw=w * 0.8, alpha=0.9)
    ax.plot(t, P[:, 2], ':',  color=c, lw=w * 0.8, alpha=0.9)


_AXIS_HANDLES = [
    Line2D([0], [0], color='0.35', ls='-',  lw=1, label=r'$\Delta X$'),
    Line2D([0], [0], color='0.35', ls='--', lw=1, label=r'$\Delta Y$'),
    Line2D([0], [0], color='0.35', ls=':',  lw=1, label=r'$\Delta Z$'),
]


def fig_position_vs_time(runs, es, out_dir):
    """Full record; motion spans shaded, boundaries marked, pauses annotated."""
    fig, axes = plt.subplots(len(runs), 1, figsize=(7.4, 2.7 * len(runs)))
    axes = np.atleast_1d(axes)
    for ax, r, c, w in zip(axes, runs, COLORS, LWS):
        t, P = r['t'], r['Ps']
        _draw_axes_vs_time(ax, t, P, c, 1.6)
        for k, (a, b) in enumerate(r['moving']):
            ax.axvspan(t[a], t[b], color='#0072B2', alpha=0.10, zorder=0)
            ax.axvline(t[a], color='0.35', lw=0.7, ls='--', zorder=1)
            ax.axvline(t[b], color='0.35', lw=0.7, ls='--', zorder=1)
            ax.text((t[a] + t[b]) / 2, 0.965, 'S%d' % (k + 1),
                    transform=ax.get_xaxis_transform(), ha='center', va='top',
                    fontsize=7.5,
                    bbox=dict(boxstyle='round,pad=0.22', fc='white',
                              ec='0.7', lw=0.5, alpha=0.9))
        for a, b in r['stationary']:
            if t[b] - t[a] > 0.5:
                ax.text((t[a] + t[b]) / 2, 0.03, 'pause\n%.1f s' % (t[b] - t[a]),
                        transform=ax.get_xaxis_transform(), ha='center',
                        va='bottom', fontsize=6.5, color='0.4')
        ax.axhline(0, color='0.55', lw=0.45, ls='--', zorder=0)
        ax.set_ylabel('%s\nDispl. (mm)' % r['label'])
        _finish(ax)
    axes[-1].set_xlabel('Elapsed Time (s)')
    fig.suptitle('%s\nPosition vs. Time (shaded = motion)' % es['name'])
    # One shared legend: a per-axes legend sits on top of the S1 badge.
    fig.legend(handles=_AXIS_HANDLES, fontsize=8, ncol=3,
               loc='upper center', bbox_to_anchor=(0.5, 0.945), frameon=False)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    return _save(fig, out_dir, 'position_vs_time.jpg')


def fig_position_vs_time_truncated(runs, es, out_dir):
    """Stationary spans removed; motion segments concatenated end to end."""
    fig, axes = plt.subplots(len(runs), 1, figsize=(7.4, 2.7 * len(runs)))
    axes = np.atleast_1d(axes)
    for ax, r, c, w in zip(axes, runs, COLORS, LWS):
        t, P = r['t'], r['Ps']
        offset, bounds, tt, PP = 0.0, [], [], []
        for a, b in r['moving']:
            seg_t = t[a:b + 1] - t[a] + offset
            tt.append(seg_t)
            PP.append(P[a:b + 1])
            offset = seg_t[-1]
            bounds.append(offset)
        if not tt:
            continue
        tt = np.concatenate(tt)
        PP = np.vstack(PP)
        _draw_axes_vs_time(ax, tt, PP, c, 1.6)
        start = 0.0
        for k, edge in enumerate(bounds):
            if k:
                ax.axvline(start, color='0.35', lw=0.8, ls='--')
            ax.text((start + edge) / 2, 0.965, 'S%d' % (k + 1),
                    transform=ax.get_xaxis_transform(), ha='center', va='top',
                    fontsize=7.5,
                    bbox=dict(boxstyle='round,pad=0.22', fc='white',
                              ec='0.7', lw=0.5, alpha=0.9))
            start = edge
        ax.axhline(0, color='0.55', lw=0.45, ls='--', zorder=0)
        ax.set_ylabel('%s\nDispl. (mm)' % r['label'])
        _finish(ax)
    axes[-1].set_xlabel('Motion Time (s), stationary spans removed')
    fig.suptitle('%s\nTruncated Position vs. Time' % es['name'])
    fig.legend(handles=_AXIS_HANDLES, fontsize=8, ncol=3,
               loc='upper center', bbox_to_anchor=(0.5, 0.945), frameon=False)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    return _save(fig, out_dir, 'position_vs_time_truncated.jpg')


def fig_segments_3d(runs, es, out_dir):
    """One 3D subplot per motion segment, all trials overlaid."""
    n = min(len(r['moving']) for r in runs)
    if n == 0:
        return None
    fig = plt.figure(figsize=(4.7 * n, 4.9))
    for k in range(n):
        ax = fig.add_subplot(1, n, k + 1, projection='3d')
        chunks = []
        for r, c, w in zip(runs, COLORS, LWS):
            a, b = r['moving'][k]
            S = r['Ps'][a:b + 1]
            chunks.append(S)
            ax.plot(S[:, 0], S[:, 1], S[:, 2], color=c, lw=w, alpha=0.9,
                    label=r['label'])
            ax.scatter(*S[0],  marker='o', s=45, facecolors='white',
                       edgecolors=c, linewidths=1.4, zorder=6)
            ax.scatter(*S[-1], marker='s', s=40, color=c, zorder=6)
        kind, commanded = (es['steps'][k] if k < len(es['steps'])
                           else ('?', float('nan')))
        desc = ('Translation %.4g mm' % commanded if kind == 'T'
                else 'Rotation %.4g%s' % (commanded, DEG) if kind == 'R'
                else 'Unmatched')
        ax.set_title('Segment %d\n%s' % (k + 1, desc), fontsize=9)
        ax.set_xlabel(r'$\Delta X$ (mm)', fontsize=8)
        ax.set_ylabel(r'$\Delta Y$ (mm)', fontsize=8)
        ax.set_zlabel(r'$\Delta Z$ (mm)', fontsize=8)
        _equal_3d(ax, np.vstack(chunks))
        ax.view_init(elev=22, azim=-58)
        ax.tick_params(labelsize=6.5)
        if k == 0:
            ax.legend(fontsize=7, loc='upper left')
    fig.suptitle('%s  --  Trajectory by Motion Segment' % es['name'],
                 fontsize=11, y=0.99)
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    return _save(fig, out_dir, 'trajectory_3d_segments.jpg')


def _grouped_error_axis(ax, entries, unit, title):
    """Grouped commanded-vs-measured bars with spread across trials."""
    idx = np.arange(len(entries))
    cmd = [e['commanded'] for e in entries]
    mean = [e['mean'] for e in entries]
    std = [e['std'] for e in entries]
    ax.bar(idx - 0.2, cmd, 0.38, color='#BBBBBB', edgecolor='0.35',
           lw=0.6, label='Commanded')
    ax.bar(idx + 0.2, mean, 0.38, yerr=std, capsize=3.5, color='#0072B2',
           edgecolor='0.25', lw=0.6, label='Measured', zorder=3,
           error_kw=dict(lw=0.9, ecolor='0.2'))
    for i, e in enumerate(entries):
        ax.annotate('%+.2f' % (e['mean'] - e['commanded']),
                    xy=(i + 0.2, e['mean'] + e['std']),
                    xytext=(0, 4), textcoords='offset points',
                    ha='center', fontsize=7.5)
    ax.set_xticks(idx)
    ax.set_xticklabels([e['name'] for e in entries])
    ax.set_ylabel(unit)
    ax.set_title(title)
    # Headroom so the legend clears the tallest bar and its annotation.
    ax.set_ylim(0, max(max(cmd), max(m + s for m, s in zip(mean, std))) * 1.22)
    ax.legend(fontsize=7.5, loc='upper left')
    ax.grid(axis='x', visible=False)
    _finish(ax)


def fig_kinematic_error(runs, es, out_dir):
    """Commanded vs measured translation / rotation, aggregated over trials."""
    groups = {'T': [], 'R': []}
    for i, (kind, commanded) in enumerate(es['steps']):
        vals = [m['measured'] for r in runs for m in r['measures']
                if m['step'] == i + 1]
        if not vals:
            continue
        groups[kind].append(dict(name='S%d' % (i + 1), commanded=commanded,
                                 mean=float(np.mean(vals)),
                                 std=float(np.std(vals))))
    panels = [(k, v) for k, v in groups.items() if v]
    if not panels:
        return None
    fig, axes = plt.subplots(1, len(panels),
                             figsize=(4.6 * len(panels), 4.2))
    axes = np.atleast_1d(axes)
    for ax, (kind, entries) in zip(axes, panels):
        if kind == 'T':
            _grouped_error_axis(ax, entries, 'Translation (mm)',
                                'Tip Translation (arc length)')
        else:
            _grouped_error_axis(ax, entries, 'Rotation (deg)',
                                'Tip Rotation (swept angle)')
    fig.suptitle('%s  --  Commanded vs. Measured' % es['name'], fontsize=11)
    fig.tight_layout()
    return _save(fig, out_dir, 'kinematic_error.jpg')


def fig_overview_3d(all_results, out_dir, ncols=2, fs=OVERVIEW_FONT_SCALE,
                    cs=OVERVIEW_CANVAS_SCALE):
    """One 3D trajectory panel per experiment set, for the report figure."""
    with plt.rc_context(_scaled_font_rc(fs)):
        return _overview_3d_figure(all_results, out_dir, ncols, fs, cs)


def _overview_3d_figure(all_results, out_dir, ncols, fs, cs):
    n = len(all_results)
    ncols = min(ncols, n)
    nrows = int(np.ceil(n / ncols))
    # Panels are wide (7.4in) relative to their height so the enlarged titles
    # and 3D axis labels have somewhere to go.
    fig = plt.figure(figsize=(7.4 * ncols * cs, 7.6 * nrows * cs))
    for k, (es, runs) in enumerate(all_results):
        ax = fig.add_subplot(nrows, ncols, k + 1, projection='3d')
        for r, c, w in zip(runs, COLORS, LWS):
            P = r['Ps']
            ax.plot(P[:, 0], P[:, 1], P[:, 2], color=c, lw=w * 1.7 * cs,
                    alpha=0.9, label=r['label'])
            ax.scatter(*P[0],  marker='o', s=150 * cs, facecolors='white',
                       edgecolors=c, linewidths=2.6 * cs, zorder=6)
            ax.scatter(*P[-1], marker='s', s=135 * cs, color=c, zorder=6)
        # Compact form ("T35mm" not "T 35 mm"): the three-step sequences are
        # 26 chars spelled out, which collides across the column gutter.
        steps = ' + '.join(
            ('T%.4gmm' % v) if kind == 'T' else ('R%.4g%s' % (v, DEG))
            for kind, v in es['steps'])
        # Short titles: the full set names do not fit at this text size, so the
        # set key plus the trial range carries the identification and the
        # descriptive names stay in the caption.
        nums = [t[0].replace('Exp ', '') for t in es['trials']]
        ax.set_title('(%s) %s   Exp %s-%s\n%s'
                     % ('abcdef'[k], es['key'], nums[0], nums[-1], steps),
                     fontsize=8.5 * fs, pad=12 * fs)
        ax.set_xlabel(r'$\Delta X$ (mm)', fontsize=7.5 * fs, labelpad=8 * fs)
        ax.set_ylabel(r'$\Delta Y$ (mm)', fontsize=7.5 * fs, labelpad=8 * fs)
        ax.set_zlabel(r'$\Delta Z$ (mm)', fontsize=7.5 * fs, labelpad=7 * fs)
        _equal_3d(ax, np.vstack([r['Ps'] for r in runs]))
        ax.view_init(elev=22, azim=-58)
        # Three ticks per axis: more than that collides at this text size.
        for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
            axis.set_major_locator(plt.MaxNLocator(3))
        ax.tick_params(labelsize=6 * fs, pad=3 * fs)

    # One shared legend rather than six: the colour order is the same in every
    # panel, so it encodes trial index, and six enlarged legend boxes would
    # cover the trajectories they describe.
    handles = [Line2D([0], [0], color=c, lw=2.4 * fs,
                      label='Trial %d' % (i + 1))
               for i, c in enumerate(COLORS)]
    handles += [Line2D([0], [0], ls='none', marker='o', ms=2.6 * fs, mfc='w',
                       mec='k', mew=0.5 * fs, label='Start'),
                Line2D([0], [0], ls='none', marker='s', ms=2.2 * fs,
                       color='k', label='End')]
    fig.legend(handles=handles, loc='lower center', ncol=5,
               fontsize=7.5 * fs, frameon=False, bbox_to_anchor=(0.5, 0.006))
    fig.suptitle('Tip Trajectories -- All Experiment Sets',
                 fontsize=10.5 * fs, y=0.995)
    # Bottom margin has to clear the last row's X/Y labels *and* the legend.
    fig.tight_layout(rect=(0, 0.075, 1, 0.955))
    # tight_layout cannot measure 3D axis labels -- they are drawn from the
    # projection and extend outside the axes bbox -- so at this text size the
    # inter-row gap has to be set by hand or each title lands on the row
    # above's X/Y labels.
    fig.subplots_adjust(hspace=0.95, wspace=0.20)
    return _save(fig, out_dir, 'overview_trajectories_3d.jpg')


def fig_summary(all_results, out_dir, fs=SUMMARY_FONT_SCALE,
                cs=SUMMARY_CANVAS_SCALE):
    """Absolute error per set, translation and rotation side by side."""
    with plt.rc_context(_scaled_font_rc(fs)):
        return _summary_figure(all_results, out_dir, fs, cs)


def _summary_figure(all_results, out_dir, fs, cs):
    fig, axes = plt.subplots(1, 2, figsize=(10.4 * cs, 4.3 * cs))
    for ax, kind, unit, title in (
            (axes[0], 'T', 'Absolute Error (mm)', 'Translation Error'),
            (axes[1], 'R', 'Absolute Error (deg)', 'Rotation Error')):
        names, means, stds = [], [], []
        for es, runs in all_results:
            errs = [abs(m['error']) for r in runs for m in r['measures']
                    if m['kind'] == kind]
            if errs:
                names.append(es['key'])
                means.append(float(np.mean(errs)))
                stds.append(float(np.std(errs)))
        if not names:
            ax.axis('off')
            continue
        idx = np.arange(len(names))
        ax.bar(idx, means, 0.6, yerr=stds, capsize=4 * cs, color='#0072B2',
               edgecolor='0.25', lw=0.6 * cs,
               error_kw=dict(lw=0.9 * cs, ecolor='0.2'))
        for i, (m, s) in enumerate(zip(means, stds)):
            ax.annotate('%.2f' % m, xy=(i, m + s), xytext=(0, 4 * fs),
                        textcoords='offset points', ha='center',
                        fontsize=7.5 * fs)
        ax.set_xticks(idx)
        ax.set_xticklabels(names)
        ax.set_ylabel(unit)
        ax.set_title(title)
        # Headroom for the enlarged value labels above the error bars.
        ax.set_ylim(0, max(m + s for m, s in zip(means, stds)) * 1.18)
        ax.grid(axis='x', visible=False)
        _finish(ax)
    fig.suptitle('Kinematic Accuracy Summary (mean $\\pm$ s.d. over trials)',
                 fontsize=11 * fs)
    fig.tight_layout()
    return _save(fig, out_dir, 'summary_kinematic_error.jpg')


# ─── Driver ──────────────────────────────────────────────────────────────────
def run_set(es, args):
    out_dir = os.path.join(args.out_dir, es['key'])
    os.makedirs(out_dir, exist_ok=True)

    runs = []
    for label, fname in es['trials']:
        path = os.path.join(args.data_dir, fname)
        if not os.path.exists(path):
            print('  ! missing %s -- skipped' % fname)
            continue
        runs.append(analyse_trial(label, path, es['steps'], args))
    if not runs:
        print('  ! no usable trials')
        return None

    print('\n%s  (%s)' % (es['name'], es['key']))
    if es['note']:
        print('  note: %s' % es['note'])
    for r in runs:
        flag = '  << segment count != commanded steps' if r['seg_mismatch'] else ''
        print('  %-7s %d samples, %d segment(s), plane RMS %.3f mm%s'
              % (r['label'], len(r['t']), len(r['moving']), r['plane_rms'], flag))
        for m in r['measures']:
            if m['kind'] == 'T':
                print('     S%d translation: %6.2f / %6.2f mm  (err %+6.2f mm)'
                      % (m['step'], m['measured'], m['commanded'], m['error']))
            else:
                print('     S%d rotation   : %6.1f / %6.1f deg (err %+6.1f deg)'
                      '  r=%.2f mm, fit resid %.3f mm'
                      % (m['step'], m['measured'], m['commanded'], m['error'],
                         m['radius'], m['resid']))

    fig_trajectory_3d(runs, es, out_dir)
    fig_projection(runs, es, out_dir, 0, 1, 'X', 'Y', 'projection_xy.jpg')
    fig_projection(runs, es, out_dir, 0, 2, 'X', 'Z', 'projection_xz.jpg')
    fig_plane_fit(runs, es, out_dir, args.plane_window)
    fig_position_vs_time(runs, es, out_dir)
    fig_position_vs_time_truncated(runs, es, out_dir)
    fig_segments_3d(runs, es, out_dir)
    fig_kinematic_error(runs, es, out_dir)
    print('  -> %s/' % out_dir)
    return runs


def write_csv(all_results, out_dir):
    path = os.path.join(out_dir, 'kinematic_results.csv')
    with open(path, 'w', newline='') as fh:
        wr = csv.writer(fh)
        wr.writerow(['set', 'set_name', 'trial', 'segment', 'type',
                     'commanded', 'measured', 'error', 'abs_error',
                     'circle_radius_mm', 'circle_resid_mm',
                     't_start_s', 't_end_s'])
        for es, runs in all_results:
            for r in runs:
                for m in r['measures']:
                    wr.writerow([
                        es['key'], es['name'], r['label'], m['step'],
                        'translation' if m['kind'] == 'T' else 'rotation',
                        '%.4f' % m['commanded'], '%.4f' % m['measured'],
                        '%+.4f' % m['error'], '%.4f' % abs(m['error']),
                        '' if np.isnan(m['radius']) else '%.4f' % m['radius'],
                        '' if np.isnan(m['resid']) else '%.4f' % m['resid'],
                        '%.3f' % m['t0'], '%.3f' % m['t1'],
                    ])
    return path


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--data-dir', default=os.path.dirname(os.path.abspath(__file__)),
                    help='directory holding the exp*.csv recordings')
    ap.add_argument('--out-dir', default=None,
                    help='output directory (default: <data-dir>/figures)')
    ap.add_argument('--sets', nargs='*', default=None,
                    help='only these set keys (e.g. set1 set4a)')
    ap.add_argument('--plane-window', type=float, default=10.0,
                    help='seconds of motion used for the plane fit (default: 10)')
    ap.add_argument('--speed-thresh', type=float, default=0.3,
                    help='moving/stationary speed threshold, mm/s (default: 0.3)')
    ap.add_argument('--min-move-s', type=float, default=1.0,
                    help='shortest accepted motion segment, s (default: 1.0)')
    ap.add_argument('--min-pause-s', type=float, default=0.8,
                    help='shortest accepted pause, s (default: 0.8).  Exp 19 '
                         'pauses for only ~0.9 s between steps; 0.6-0.8 '
                         'resolves all 18 trials, 1.0+ merges its first two.')
    ap.add_argument('--smooth-s', type=float, default=0.5,
                    help='position smoothing window, s (default: 0.5)')
    ap.add_argument('--arc-step-s', type=float, default=0.2,
                    help='arc-length integration step, s (default: 0.2)')
    ap.add_argument('--jump-mm', type=float, default=2.0,
                    help='spike rejection threshold, mm (default: 2.0)')
    args = ap.parse_args()

    if args.out_dir is None:
        args.out_dir = os.path.join(args.data_dir, 'figures')
    os.makedirs(args.out_dir, exist_ok=True)
    apply_style()

    wanted = EXPERIMENT_SETS
    if args.sets:
        wanted = [e for e in EXPERIMENT_SETS if e['key'] in args.sets]
        if not wanted:
            raise SystemExit('no matching sets; available: %s'
                             % ', '.join(e['key'] for e in EXPERIMENT_SETS))

    all_results = []
    for es in wanted:
        runs = run_set(es, args)
        if runs:
            all_results.append((es, runs))

    if all_results:
        fig_summary(all_results, args.out_dir)
        fig_overview_3d(all_results, args.out_dir)
        print('\nwrote %s' % write_csv(all_results, args.out_dir))
        print('wrote %s/summary_kinematic_error.jpg' % args.out_dir)
        print('wrote %s/overview_trajectories_3d.jpg' % args.out_dir)
    return 0


if __name__ == '__main__':
    sys.exit(main())
