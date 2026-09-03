#!/usr/bin/env python3
"""
ndi_plot_xyz.py -- plot a CSV recorded from the /tip_xyz topic.

    python3 ndi_plot_xyz.py run.csv
    python3 ndi_plot_xyz.py run.csv --save run.png
    python3 ndi_plot_xyz.py run.csv --no-show --save run.png     # headless
    python3 ndi_plot_xyz.py run.csv --no-3d                      # ZY plane instead
    python3 ndi_plot_xyz.py run.csv --no-outlier-removal         # keep every sample
    python3 ndi_plot_xyz.py run.csv --no-plane-fit               # skip plane projection

Outlier samples (robust per-axis modified z-score, default cutoff 3.5) are
dropped before plotting unless --no-outlier-removal is given.  A best-fit
plane (PCA/SVD through the centroid) is also fit to the (outlier-free)
points and the points are projected onto that plane's own U/V axes, shown
as an extra 2D subplot, unless --no-plane-fit is given.

Needs only numpy and matplotlib -- no ROS.  Expects the header written by
ndi_record_xyz.py: stamp,elapsed,x,y,z,frame_id
"""

import argparse
import csv
import sys


def load(path):
    """Read the CSV into parallel lists.  Returns (elapsed, x, y, z, frame_id)."""
    elapsed, xs, ys, zs = [], [], [], []
    frame_id = ""
    with open(path, "r") as handle:
        reader = csv.DictReader(handle)
        missing = {"elapsed", "x", "y", "z"} - set(reader.fieldnames or [])
        if missing:
            raise SystemExit("%s is missing column(s): %s"
                             % (path, ", ".join(sorted(missing))))
        for row in reader:
            elapsed.append(float(row["elapsed"]))
            xs.append(float(row["x"]))
            ys.append(float(row["y"]))
            zs.append(float(row["z"]))
            frame_id = row.get("frame_id", frame_id) or frame_id
    if not elapsed:
        raise SystemExit("%s contains no samples" % path)
    return elapsed, xs, ys, zs, frame_id


def remove_outliers(t, x, y, z, thresh=3.5):
    """Drop samples whose x/y/z is a robust outlier on any axis.

    Uses the modified z-score (median + MAD) per axis rather than mean/std,
    since a handful of tracker glitches can otherwise skew the mean enough
    to hide themselves.  Returns (t, x, y, z, n_removed).
    """
    import numpy as np

    mask = np.ones(len(t), dtype=bool)
    for series in (x, y, z):
        median = np.median(series)
        mad = np.median(np.abs(series - median))
        if mad == 0:
            continue
        modified_z = 0.6745 * (series - median) / mad
        mask &= np.abs(modified_z) <= thresh
    n_removed = int((~mask).sum())
    return t[mask], x[mask], y[mask], z[mask], n_removed


def fit_plane(x, y, z):
    """Least-squares best-fit plane through the points, via PCA/SVD.

    Returns (centroid, u_axis, v_axis, normal, rms_residual): u_axis/v_axis
    span the plane (the two directions of largest variance) and normal is
    the plane's normal (the direction of smallest variance).
    """
    import numpy as np

    points = np.column_stack((x, y, z))
    centroid = points.mean(axis=0)
    centered = points - centroid
    _, _, vt = np.linalg.svd(centered, full_matrices=False)
    u_axis, v_axis, normal = vt[0], vt[1], vt[2]
    rms_residual = float(np.sqrt(np.mean((centered @ normal) ** 2)))
    return centroid, u_axis, v_axis, normal, rms_residual


def project_to_plane(x, y, z, centroid, u_axis, v_axis):
    """Project 3D points onto a plane's (u, v) coordinate system."""
    import numpy as np

    centered = np.column_stack((x, y, z)) - centroid
    return centered @ u_axis, centered @ v_axis


def plot(path, units="m", save=None, show=True, connect=True, three_d=True,
         remove_outlier_points=True, outlier_thresh=3.5, plane_fit=True):
    import matplotlib
    if not show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    if three_d:
        # Before matplotlib 3.2 the "3d" projection is only registered as a
        # side effect of importing mplot3d, so add_subplot(projection="3d")
        # raises ValueError without this.  Ubuntu 20.04 ships 3.1.2.
        try:
            from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
        except ImportError:
            print("mpl_toolkits.mplot3d unavailable; showing the ZY plane "
                  "instead of the 3D trajectory")
            three_d = False

    elapsed, xs, ys, zs, frame_id = load(path)
    t = np.asarray(elapsed)
    x, y, z = np.asarray(xs), np.asarray(ys), np.asarray(zs)

    n_removed = 0
    if remove_outlier_points:
        t, x, y, z, n_removed = remove_outliers(t, x, y, z, thresh=outlier_thresh)
        if n_removed:
            print("removed %d outlier sample(s) (modified z-score > %.1f)"
                  % (n_removed, outlier_thresh))
        if len(t) < 2:
            raise SystemExit("only %d sample(s) left after outlier removal; "
                              "try --no-outlier-removal or a higher --outlier-thresh"
                              % len(t))

    ncols = 3 if plane_fit else 2
    fig = plt.figure(figsize=(6 * ncols, 9.5))
    fig.suptitle("%s  --  %d samples over %.1f s  (frame: %s)%s"
                 % (path, len(t), t[-1] - t[0], frame_id or "unknown",
                    "  [%d outliers removed]" % n_removed if n_removed else ""))

    ax_xz = fig.add_subplot(2, ncols, 1)
    ax_xy = fig.add_subplot(2, ncols, 2)
    ax_c = None
    if three_d:
        try:
            ax_c = fig.add_subplot(2, ncols, ncols + 1, projection="3d")
        except ValueError as exc:
            print("3D projection unavailable (%s); showing the ZY plane "
                  "instead" % exc)
            three_d = False
    if ax_c is None:
        ax_c = fig.add_subplot(2, ncols, ncols + 1)
    ax_t = fig.add_subplot(2, ncols, ncols + 2)
    ax_plane = fig.add_subplot(2, ncols, 3) if plane_fit else None

    # flat projections, coloured by time so direction of travel is visible
    planes = [(ax_xz, x, z, "X", "Z"), (ax_xy, x, y, "X", "Y")]
    if not three_d:
        planes.append((ax_c, z, y, "Z", "Y"))
    for ax, horiz, vert, hlabel, vlabel in planes:
        if connect:
            ax.plot(horiz, vert, "-", color="0.8", linewidth=0.8, zorder=1)
        dots = ax.scatter(horiz, vert, c=t, cmap="viridis", s=12, zorder=2)
        ax.set_xlabel("%s (%s)" % (hlabel, units))
        ax.set_ylabel("%s (%s)" % (vlabel, units))
        ax.set_title("%s%s plane" % (hlabel, vlabel))
        ax.set_aspect("equal", adjustable="datalim")
        ax.grid(True, linewidth=0.3, alpha=0.6)
        ax.plot(horiz[0], vert[0], "o", mfc="none", mec="green",
                ms=11, mew=1.6, zorder=3, label="start")
        ax.plot(horiz[-1], vert[-1], "s", mfc="none", mec="red",
                ms=11, mew=1.6, zorder=3, label="end")
    ax_xz.legend(loc="best", fontsize=8)
    fig.colorbar(dots, ax=ax_xy, label="time (s)", fraction=0.046)

    if three_d:
        if connect:
            ax_c.plot(x, y, z, "-", color="0.8", linewidth=0.8)
        ax_c.scatter(x, y, z, c=t, cmap="viridis", s=10, depthshade=False)
        ax_c.plot([x[0]], [y[0]], [z[0]], "o", mfc="none", mec="green",
                  ms=11, mew=1.6)
        ax_c.plot([x[-1]], [y[-1]], [z[-1]], "s", mfc="none", mec="red",
                  ms=11, mew=1.6)
        ax_c.set_xlabel("X (%s)" % units)
        ax_c.set_ylabel("Y (%s)" % units)
        ax_c.set_zlabel("Z (%s)" % units)
        ax_c.set_title("3D trajectory")
        # equal scale on all three axes: expand each about its midpoint to a
        # common span, otherwise a shallow path looks dramatic
        half = max(max(float(np.ptp(v)) for v in (x, y, z)), 1e-6) / 2.0
        for setlim, series in ((ax_c.set_xlim, x), (ax_c.set_ylim, y),
                               (ax_c.set_zlim, z)):
            mid = (float(series.max()) + float(series.min())) / 2.0
            setlim(mid - half, mid + half)
        try:
            ax_c.set_box_aspect((1.0, 1.0, 1.0), zoom=1.25)  # matplotlib >= 3.6
        except TypeError:
            try:
                ax_c.set_box_aspect((1.0, 1.0, 1.0))         # matplotlib >= 3.3
            except AttributeError:
                pass
        except AttributeError:
            pass
        ax_c.view_init(elev=22, azim=-60)

    rms_residual = None
    if plane_fit:
        centroid, u_axis, v_axis, normal, rms_residual = fit_plane(x, y, z)
        u, v = project_to_plane(x, y, z, centroid, u_axis, v_axis)
        if connect:
            ax_plane.plot(u, v, "-", color="0.8", linewidth=0.8, zorder=1)
        plane_dots = ax_plane.scatter(u, v, c=t, cmap="viridis", s=12, zorder=2)
        ax_plane.plot(u[0], v[0], "o", mfc="none", mec="green",
                      ms=11, mew=1.6, zorder=3, label="start")
        ax_plane.plot(u[-1], v[-1], "s", mfc="none", mec="red",
                      ms=11, mew=1.6, zorder=3, label="end")
        ax_plane.set_xlabel("U (%s)" % units)
        ax_plane.set_ylabel("V (%s)" % units)
        ax_plane.set_title("best-fit plane projection (rms %.4f %s)"
                            % (rms_residual, units))
        ax_plane.set_aspect("equal", adjustable="datalim")
        ax_plane.grid(True, linewidth=0.3, alpha=0.6)
        ax_plane.legend(loc="best", fontsize=8)
        fig.colorbar(plane_dots, ax=ax_plane, label="time (s)", fraction=0.046)

    # each axis against time, to spot dropouts and drift.  These are plotted
    # relative to their own mean: the absolute offsets differ by a metre or
    # more, which would flatten the small variations we actually care about.
    ax = ax_t
    for series, label in ((x, "x"), (y, "y"), (z, "z")):
        ax.plot(t, series - series.mean(), linewidth=1.0,
                label="%s (mean %+.4f)" % (label, series.mean()))
    ax.axhline(0.0, color="0.6", linewidth=0.6, zorder=0)
    ax.set_xlabel("time (s)")
    ax.set_ylabel("deviation from mean (%s)" % units)
    ax.set_title("position vs time, mean removed")
    ax.grid(True, linewidth=0.3, alpha=0.6)
    ax.legend(loc="best", fontsize=8)

    try:
        fig.tight_layout()
    except Exception:
        pass

    span = [float(np.ptp(v)) for v in (x, y, z)]
    print("samples : %d over %.2f s (%.1f Hz average)"
          % (len(t), t[-1] - t[0], (len(t) - 1) / (t[-1] - t[0]) if t[-1] > t[0] else 0.0))
    print("mean    : x=%.4f  y=%.4f  z=%.4f %s" % (x.mean(), y.mean(), z.mean(), units))
    print("std dev : x=%.4f  y=%.4f  z=%.4f %s" % (x.std(), y.std(), z.std(), units))
    print("range   : x=%.4f  y=%.4f  z=%.4f %s" % (span[0], span[1], span[2], units))
    if plane_fit:
        print("plane   : normal=[%.4f, %.4f, %.4f]  rms residual=%.4f %s"
              % (normal[0], normal[1], normal[2], rms_residual, units))

    if save:
        fig.savefig(save, dpi=150)
        print("saved plot to %s" % save)
    if show:
        plt.show()
    else:
        plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("csv", help="CSV file written by ndi_record_xyz.py")
    parser.add_argument("--units", default="m", help="axis label only (default: m)")
    parser.add_argument("--save", default=None, help="also write the figure to this path")
    parser.add_argument("--no-show", dest="show", action="store_false",
                        help="do not open a window (use with --save on a headless box)")
    parser.add_argument("--no-connect", dest="connect", action="store_false",
                        help="scatter only, do not join consecutive samples")
    parser.add_argument("--no-3d", dest="three_d", action="store_false",
                        help="show the ZY plane instead of the 3D trajectory")
    parser.add_argument("--no-outlier-removal", dest="remove_outlier_points",
                        action="store_false",
                        help="keep all samples, do not drop outliers before plotting")
    parser.add_argument("--outlier-thresh", type=float, default=3.5,
                        help="modified z-score cutoff per axis for outlier "
                             "removal (default: 3.5; lower = stricter)")
    parser.add_argument("--no-plane-fit", dest="plane_fit", action="store_false",
                        help="do not fit/plot the best-fit plane projection")
    args = parser.parse_args()
    plot(args.csv, units=args.units, save=args.save,
         show=args.show, connect=args.connect, three_d=args.three_d,
         remove_outlier_points=args.remove_outlier_points,
         outlier_thresh=args.outlier_thresh, plane_fit=args.plane_fit)
    return 0


if __name__ == "__main__":
    sys.exit(main())