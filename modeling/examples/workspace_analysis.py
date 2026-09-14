#!/usr/bin/env python3
"""
workspace_analysis.py -- 3D robot workspace analysis for the CT-SDR.

Analyzes and visualizes the 3D reachable workspace of the 4-DoF Concentric
Tube Steerable Drilling Robot (CT-SDR) across three operational scenarios:
    1) Inner tube alone (stroke 0 to 83 mm, 360 deg rotation)
    2) Outer tube alone (stroke 0 to 46 mm, 360 deg rotation)
    3) Both tubes together (coupled Cosserat rod mechanics, 4 DoFs)

Generates four publication-quality figures in figures/ (300 DPI, cmr10 font):
    - workspace_inner_alone.png
    - workspace_outer_alone.png
    - workspace_both_tubes.png (with extreme poses in meridian cross-section & 3D cutaway)
    - workspace_comparison.png

Usage:
    ./venv/bin/python3 modeling/examples/workspace_analysis.py
    ./venv/bin/python3 modeling/examples/workspace_analysis.py --quick
"""

from __future__ import annotations

import argparse
import os
import sys
import time
from typing import Dict, List, Tuple

import matplotlib
if not getattr(matplotlib, "_interactive_mode", False):
    try:
        matplotlib.use("Agg")
    except Exception:
        pass
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np
from scipy.spatial import ConvexHull, Delaunay
import collections

# Ensure modeling root is in path
MODELING_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if MODELING_DIR not in sys.path:
    sys.path.insert(0, MODELING_DIR)

import examples.vizstyle as vs
from ctr import CTSDR, Joints
from ctr.tube import Tube
import yaml

FIGURES_DIR = os.path.join(MODELING_DIR, "figures")
CONFIG_PATH = os.path.join(MODELING_DIR, "config", "ct_sdr.yaml")
DATA_CACHE_PATH = os.path.join(FIGURES_DIR, "workspace_data.npz")


def compute_alpha_shape_boundary(
    points: np.ndarray,
    alpha: float = 9.0,
) -> Tuple[np.ndarray, np.ndarray, float, float]:
    """Compute alpha shape (concave hull) for 2D points to strictly eliminate unreachable space.
    
    Returns:
        boundary_rz: (K, 2) ordered perimeter loop
        valid_tri: (M, 3) indices of valid Delaunay triangles
        area_2d_mm2: sum of valid triangle areas (exact non-convex area)
        volume_mm3: revolved volume via Pappus's theorem on valid triangles
    """
    tri = Delaunay(points)
    pa = points[tri.simplices[:, 0]]
    pb = points[tri.simplices[:, 1]]
    pc = points[tri.simplices[:, 2]]
    
    a_len = np.hypot(pb[:, 0] - pa[:, 0], pb[:, 1] - pa[:, 1])
    b_len = np.hypot(pc[:, 0] - pb[:, 0], pc[:, 1] - pb[:, 1])
    c_len = np.hypot(pa[:, 0] - pc[:, 0], pa[:, 1] - pc[:, 1])
    
    s_len = (a_len + b_len + c_len) / 2.0
    tri_areas = np.sqrt(np.maximum(s_len * (s_len - a_len) * (s_len - b_len) * (s_len - c_len), 1e-12))
    circum_r = (a_len * b_len * c_len) / (4.0 * tri_areas)
    
    max_edges = np.maximum(a_len, np.maximum(b_len, c_len))
    valid = (circum_r < alpha) & (max_edges < 2.2 * alpha)
    valid_tri = tri.simplices[valid]
    valid_areas = tri_areas[valid]
    
    area_2d_mm2 = float(np.sum(valid_areas))
    centroids_r = (pa[valid, 0] + pb[valid, 0] + pc[valid, 0]) / 3.0
    volume_mm3 = float(2.0 * np.pi * np.sum(centroids_r * valid_areas))
    
    edge_counts = collections.defaultdict(int)
    for simp in valid_tri:
        for i in range(3):
            edge = tuple(sorted([simp[i], simp[(i+1)%3]]))
            edge_counts[edge] += 1
    boundary_edges = [edge for edge, count in edge_counts.items() if count == 1]
    
    adj = collections.defaultdict(list)
    for u, v in boundary_edges:
        adj[u].append(v)
        adj[v].append(u)
        
    visited = set()
    cycles = []
    for start_node in adj:
        if start_node in visited:
            continue
        cycle = [start_node]
        visited.add(start_node)
        curr = start_node
        while True:
            nxt_cands = [n for n in adj[curr] if n not in visited]
            if not nxt_cands:
                break
            nxt = nxt_cands[0]
            cycle.append(nxt)
            visited.add(nxt)
            curr = nxt
        if len(cycle) >= 3:
            cycle.append(cycle[0])
            cycles.append(cycle)
            
    if not cycles:
        hull = ConvexHull(points)
        boundary_rz = points[np.append(hull.vertices, hull.vertices[0])]
    else:
        best_cycle = max(cycles, key=lambda c: len(c))
        boundary_rz = points[best_cycle]
        
    return boundary_rz, valid_tri, area_2d_mm2, volume_mm3


# ---------------------------------------------------------------------------
# Kinematic Models
# ---------------------------------------------------------------------------

def single_tube_backbone(
    s_arr: np.ndarray,
    radius_of_curvature: float,
    curved_length: float,
    rotation_rad: float = 0.0,
) -> np.ndarray:
    """Compute 3D backbone curve p(s) for a single pre-curved tube in free space.
    
    Coordinate convention: s = 0 is guide exit, initial tangent along +z.
    The tube is straight over s in [0, s_straight] and circularly curved
    over s in [s_straight, s_straight + curved_length].
    Rotation rotates the bending plane about the z-axis.
    """
    R = radius_of_curvature
    pts = []
    
    s_max = s_arr[-1] if len(s_arr) > 0 else 0.0
    s_straight = max(0.0, s_max - curved_length)
    
    for s in s_arr:
        if s <= s_straight:
            x = 0.0
            y = 0.0
            z = s
        else:
            s_curv = s - s_straight
            y_unrot = -R * (1.0 - np.cos(s_curv / R))
            z_unrot = s_straight + R * np.sin(s_curv / R)
            
            x = -y_unrot * np.sin(rotation_rad)
            y = y_unrot * np.cos(rotation_rad)
            z = z_unrot
            
        pts.append([x, y, z])
        
    return np.array(pts)


def compute_single_tube_workspace(
    max_stroke_mm: float,
    radius_of_curvature_mm: float,
    curved_length_mm: float,
    n_stroke: int = 100,
    n_rot: int = 72,
) -> Dict[str, np.ndarray]:
    """Compute the 2D meridian profile and 3D surface mesh for a single tube."""
    strokes = np.linspace(0.0, max_stroke_mm, n_stroke)
    thetas = np.linspace(0.0, 2.0 * np.pi, n_rot)
    
    R = radius_of_curvature_mm
    
    r_profile = []
    z_profile = []
    psi_deg = []
    
    for L in strokes:
        s_straight = max(0.0, L - curved_length_mm)
        s_curv = L - s_straight
        r = R * (1.0 - np.cos(s_curv / R))
        z = s_straight + R * np.sin(s_curv / R)
        psi = np.degrees(s_curv / R)
        r_profile.append(r)
        z_profile.append(z)
        psi_deg.append(psi)
        
    r_profile = np.array(r_profile)
    z_profile = np.array(z_profile)
    psi_deg = np.array(psi_deg)
    
    L_grid, TH_grid = np.meshgrid(strokes, thetas)
    
    S_straight_grid = np.maximum(0.0, L_grid - curved_length_mm)
    S_curv_grid = L_grid - S_straight_grid
    
    R_grid = R * (1.0 - np.cos(S_curv_grid / R))
    Z_grid = S_straight_grid + R * np.sin(S_curv_grid / R)
    X_grid = R_grid * np.sin(TH_grid)
    Y_grid = -R_grid * np.cos(TH_grid)
    
    dr_dL = np.gradient(r_profile, strokes)
    dz_dL = np.gradient(z_profile, strokes)
    ds_arc = np.sqrt(dr_dL**2 + dz_dL**2)
    try:
        from scipy.integrate import trapezoid
        surface_area_mm2 = float(trapezoid(2.0 * np.pi * r_profile * ds_arc, strokes))
    except ImportError:
        surface_area_mm2 = float(np.trapezoid(2.0 * np.pi * r_profile * ds_arc, strokes))
    
    return {
        "strokes_mm": strokes,
        "thetas_rad": thetas,
        "r_profile_mm": r_profile,
        "z_profile_mm": z_profile,
        "psi_deg": psi_deg,
        "X_grid_mm": X_grid,
        "Y_grid_mm": Y_grid,
        "Z_grid_mm": Z_grid,
        "max_r_mm": float(r_profile[-1]),
        "max_z_mm": float(z_profile.max()),
        "max_psi_deg": float(psi_deg[-1]),
        "surface_area_mm2": surface_area_mm2,
    }


def compute_coupled_workspace(
    robot: CTSDR,
    max_ott_mm: float = 46.0,
    max_itt_mm: float = 83.0,
    n_ott: int = 10,
    n_itt: int = 14,
    n_angles: int = 19,
    n_rot_rev: int = 72,
    quick: bool = False,
    recompute: bool = False,
) -> Dict[str, any]:
    """Compute the 3D reachable workspace of both tubes interacting via Cosserat mechanics.
    
    Exploits rotational axisymmetry about +z: solving on (q_ott, q_itt, Delta_alpha)
    yields the complete meridian cross-section D_{rz} = {(r, z)}, which is then
    revolved 360 deg to form the full 3D solid volume.
    """
    rz_points = []
    configurations = []
    cached_ready = False
    
    if not recompute and os.path.exists(DATA_CACHE_PATH):
        try:
            cached = np.load(DATA_CACHE_PATH, allow_pickle=True)
            cached_rz = cached["rz_points"]
            if cached_rz[:, 1].max() > 75.0 and "valid_triangles" in cached:
                print(f"Loading coupled workspace points from cache: {DATA_CACHE_PATH} ...")
                rz_pts = cached_rz
                rz_points = rz_pts.tolist()
                valid_tri = cached["valid_triangles"]
                boundary_rz = cached["boundary_rz"]
                area_2d_mm2 = float(cached["area_2d_mm2"])
                volume_mm3 = float(cached["volume_mm3"])
                configs = cached["configs"] if "configs" in cached else np.array([])
                cached_ready = True
        except Exception:
            cached_ready = False
            
    if not cached_ready:
        if quick:
            n_ott = max(6, n_ott // 2)
            n_itt = max(8, n_itt // 2)
            n_angles = max(7, n_angles // 2)
            
        print(f"Sampling coupled workspace (bidirectional): {n_ott} OTT x {n_itt} ITT x {n_angles} dAlpha ...")
        t0 = time.time()
        
        ott_values = np.linspace(2.0, max_ott_mm, n_ott)
        delta_angles_deg = np.linspace(0.0, 180.0, n_angles)
        
        total_combinations = n_ott * n_itt
        combo_idx = 0
        
        for ott in ott_values:
            itt_values = np.linspace(ott, max_itt_mm, n_itt)
            for itt in itt_values:
                combo_idx += 1
                # Forward sweep (0 -> 180 deg) tracks wound torsional branch
                path_fwd = [Joints(ott=ott, itt=itt, otr=0.0, itr=a) for a in delta_angles_deg]
                try:
                    sols_fwd = robot.solve_path(path_fwd)
                    for a_deg, sol in zip(delta_angles_deg, sols_fwd):
                        tip = sol.tip_position * 1e3  # in mm
                        r = np.hypot(tip[0], tip[1])
                        z = tip[2]
                        rz_points.append([r, z])
                        configurations.append([ott, itt, a_deg, r, z])
                except Exception:
                    pass
                
                # Reverse sweep (180 -> 0 deg) tracks antagonistic unwound branch
                path_rev = [Joints(ott=ott, itt=itt, otr=0.0, itr=a) for a in delta_angles_deg[::-1]]
                try:
                    sols_rev = robot.solve_path(path_rev)
                    for a_deg, sol in zip(delta_angles_deg[::-1], sols_rev):
                        tip = sol.tip_position * 1e3
                        r = np.hypot(tip[0], tip[1])
                        z = tip[2]
                        rz_points.append([r, z])
                except Exception:
                    pass
                    
                if combo_idx % 20 == 0 or combo_idx == total_combinations:
                    pct = (combo_idx / total_combinations) * 100.0
                    sys.stdout.write(f"\r  Progress: {pct:5.1f}% ({combo_idx}/{total_combinations} sweeps)")
                    sys.stdout.flush()
                    
        sys.stdout.write("\n")
        dt = time.time() - t0
        print(f"Coupled Cosserat solves completed in {dt:.2f} s ({len(rz_points)} tip samples)")
        
        rz_pts = np.array(rz_points)
        configs = np.array(configurations)
    
    # -----------------------------------------------------------------------
    # Compute Key Landmark Extreme Poses with High Fidelity
    # -----------------------------------------------------------------------
    print("Computing landmark extreme poses for meridian cross-section...")
    extreme_poses = []
    
    # Pose 1: Aligned (0 deg) at Max Stroke (46 mm / 83 mm) -> Outermost lateral reach
    sol_0 = robot.solve(Joints(ott=46.0, itt=83.0, otr=0.0, itr=0.0))
    p1 = sol_0.p * 1e3
    tip1 = sol_0.tip_position * 1e3
    extreme_poses.append({
        "id": "1",
        "name": r"Aligned ($\Delta\alpha = 0^\circ$)",
        "desc": r"Max Lateral Reach: $q_{\mathrm{OTT}}=46, q_{\mathrm{ITT}}=83$ mm",
        "ott": 46.0,
        "itt": 83.0,
        "dalpha": 0.0,
        "p": p1,
        "tip": tip1,
        "r": float(np.hypot(tip1[0], tip1[1])),
        "z": float(tip1[2]),
        "color": "#2a78d6",  # Series blue
    })
    
    # Pose 3: S-Shape Antagonistic Inflection (ott=46, itt=83 mm, otr=0 deg, itr=180 deg)
    sol_s = robot.solve(Joints(ott=46.0, itt=83.0, otr=0.0, itr=180.0))
    p3 = sol_s.p * 1e3
    tip3 = sol_s.tip_position * 1e3
    extreme_poses.append({
        "id": "3",
        "name": r"S-Shape Inflection ($180^\circ$)",
        "desc": r"Antagonistic S-Curve: $q_{\mathrm{OTT}}=46, q_{\mathrm{ITT}}=83$ mm",
        "ott": 46.0,
        "itt": 83.0,
        "dalpha": 180.0,
        "p": p3,
        "tip": tip3,
        "r": float(np.hypot(tip3[0], tip3[1])),
        "z": float(tip3[2]),
        "color": "#8e24aa",  # Purple
    })
    
    # Pose 4: Intermediate Aligned Extension (ott=46, itt=65 mm, otr=0 deg, itr=0 deg)
    sol_4 = robot.solve(Joints(ott=46.0, itt=65.0, otr=0.0, itr=0.0))
    p4 = sol_4.p * 1e3
    tip4 = sol_4.tip_position * 1e3
    extreme_poses.append({
        "id": "4",
        "name": r"Aligned Intermediate ($0^\circ$)",
        "desc": r"Aligned Extension: $q_{\mathrm{OTT}}=46, q_{\mathrm{ITT}}=65$ mm",
        "ott": 46.0,
        "itt": 65.0,
        "dalpha": 0.0,
        "p": p4,
        "tip": tip4,
        "r": float(np.hypot(tip4[0], tip4[1])),
        "z": float(tip4[2]),
        "color": "#d32f2f",  # Crimson red
    })
    
    # Pose 5: Outer Flush Boundary (46/46 mm, 0 deg)
    sol_flush = robot.solve(Joints(ott=46.0, itt=46.0, otr=0.0, itr=0.0))
    p5 = sol_flush.p * 1e3
    tip5 = sol_flush.tip_position * 1e3
    extreme_poses.append({
        "id": "5",
        "name": r"Outer Flush Boundary ($46/46$)",
        "desc": r"Outer Max Reach: $r = 17.6, z = 41.2$ mm",
        "ott": 46.0,
        "itt": 46.0,
        "dalpha": 0.0,
        "p": p5,
        "tip": tip5,
        "r": float(np.hypot(tip5[0], tip5[1])),
        "z": float(tip5[2]),
        "color": "#00897b",  # Teal
    })
    
    # Include extreme pose tips in the boundary set
    for ep in extreme_poses:
        rz_points.append([ep["r"], ep["z"]])
    rz_pts = np.array(rz_points)
    
    # Compute true non-convex reachable boundary and volume via alpha shape
    if not cached_ready:
        boundary_rz, valid_tri, area_2d_mm2, volume_mm3 = compute_alpha_shape_boundary(rz_pts, alpha=9.0)
    
    phi_rev = np.linspace(0.0, 2.0 * np.pi, n_rot_rev)
    R_b = boundary_rz[:, 0]
    Z_b = boundary_rz[:, 1]
    
    R_b_grid, PHI_grid = np.meshgrid(R_b, phi_rev)
    Z_b_grid, _ = np.meshgrid(Z_b, phi_rev)
    X_b_grid = R_b_grid * np.cos(PHI_grid)
    Y_b_grid = R_b_grid * np.sin(PHI_grid)
    
    angles_3d = np.linspace(0.0, 2.0 * np.pi, 24, endpoint=False)
    pts_3d = []
    for r, z in rz_pts[::3]:
        for phi in angles_3d:
            pts_3d.append([r * np.cos(phi), r * np.sin(phi), z])
    pts_3d = np.array(pts_3d)
    
    data = {
        "rz_points": rz_pts,
        "configs": configs,
        "boundary_rz": boundary_rz,
        "valid_triangles": valid_tri,
        "area_2d_mm2": float(area_2d_mm2),
        "volume_mm3": float(volume_mm3),
        "min_r_mm": float(rz_pts[:, 0].min()),
        "max_r_mm": float(rz_pts[:, 0].max()),
        "min_z_mm": float(rz_pts[:, 1].min()),
        "max_z_mm": float(rz_pts[:, 1].max()),
        "X_b_grid": X_b_grid,
        "Y_b_grid": Y_b_grid,
        "Z_b_grid": Z_b_grid,
        "pts_3d": pts_3d,
        "extreme_poses": extreme_poses,
    }
    
    # Save cache for interactive viewer
    if not cached_ready:
        try:
            np.savez_compressed(
                DATA_CACHE_PATH,
                rz_points=rz_pts,
                boundary_rz=boundary_rz,
                valid_triangles=valid_tri,
                area_2d_mm2=area_2d_mm2,
                volume_mm3=volume_mm3,
                X_b_grid=X_b_grid,
                Y_b_grid=Y_b_grid,
                Z_b_grid=Z_b_grid,
                pts_3d=pts_3d,
                configs=configs,
            )
            print(f"Saved workspace cache to {DATA_CACHE_PATH}")
        except Exception as e:
            print(f"Notice: Cache save skipped ({e})")
        
    return data


# ---------------------------------------------------------------------------
# Visualization Helpers
# ---------------------------------------------------------------------------

def draw_guide_tube_3d(ax, length_mm: float = 20.0, radius_mm: float = 2.5) -> None:
    """Render a stylized stainless-steel rigid guide tube entering at z <= 0."""
    z_cyl = np.linspace(-length_mm, 0.0, 10)
    theta_cyl = np.linspace(0.0, 2.0 * np.pi, 24)
    Z_cyl, TH_cyl = np.meshgrid(z_cyl, theta_cyl)
    X_cyl = radius_mm * np.cos(TH_cyl)
    Y_cyl = radius_mm * np.sin(TH_cyl)
    
    ax.plot_surface(X_cyl, Y_cyl, Z_cyl, color="#7d7b75", alpha=0.65,
                    edgecolor="none", shade=True, zorder=2)
    ax.plot(radius_mm * np.cos(theta_cyl), radius_mm * np.sin(theta_cyl),
            np.zeros_like(theta_cyl), color="#3a3834", lw=1.8, zorder=3)


# ---------------------------------------------------------------------------
# Figure 1: Inner Tube Alone Workspace
# ---------------------------------------------------------------------------

def plot_inner_alone(data_in: Dict, out_path: str) -> None:
    """Create 3D and multi-view figures for the inner tube alone workspace."""
    vs.apply(use_cmr=True)
    fig = plt.figure(figsize=(10.2, 4.0))
    
    # Panel 1: 3D Perspective View
    ax3d = fig.add_subplot(1, 3, 1, projection="3d")
    draw_guide_tube_3d(ax3d, length_mm=25.0, radius_mm=2.2)
    
    surf = ax3d.plot_surface(
        data_in["X_grid_mm"], data_in["Y_grid_mm"], data_in["Z_grid_mm"],
        cmap="Blues", alpha=0.55, edgecolor="none", shade=True, antialiased=True
    )
    
    sample_strokes = [20.0, 40.0, 60.0, 83.0]
    sample_angles = [0.0, np.pi/2, np.pi, 3*np.pi/2]
    colors_bb = [vs.RAMP[0], vs.RAMP[1], vs.RAMP[2], vs.RAMP[4]]
    
    for L, c in zip(sample_strokes, colors_bb):
        s_vals = np.linspace(0.0, L, 50)
        bb0 = single_tube_backbone(s_vals, 57.0, 89.54, rotation_rad=0.0)
        ax3d.plot(bb0[:, 0], bb0[:, 1], bb0[:, 2], color=c, lw=2.0,
                  label=rf"$L = {L:.0f}$ mm")
        for th in sample_angles[1:]:
            bb_th = single_tube_backbone(s_vals, 57.0, 89.54, rotation_rad=th)
            ax3d.plot(bb_th[:, 0], bb_th[:, 1], bb_th[:, 2], color=c, lw=1.2,
                      ls="--", alpha=0.45)
            
    ring_idx = -1
    ax3d.plot(data_in["X_grid_mm"][:, ring_idx],
              data_in["Y_grid_mm"][:, ring_idx],
              data_in["Z_grid_mm"][:, ring_idx],
              color="#0d47a1", lw=2.4, ls="-", label=r"Tip Boundary ($83$ mm)")
              
    ax3d.set_title("3D Reachable Surface (Inner Alone)", pad=12)
    ax3d.set_xlabel(r"$x$ (mm)")
    ax3d.set_ylabel(r"$y$ (mm)")
    ax3d.set_zlabel(r"$z$ (mm)")
    vs.tidy3d(ax3d)
    ax3d.view_init(elev=26, azim=-55)
    
    all_pts = np.column_stack([data_in["X_grid_mm"].ravel(),
                               data_in["Y_grid_mm"].ravel(),
                               data_in["Z_grid_mm"].ravel()])
    vs.equal_aspect_3d(ax3d, all_pts)
    
    # Panel 2: 2D Meridian (r, z) Profile
    ax_rz = fig.add_subplot(1, 3, 2)
    ax_rz.plot(data_in["r_profile_mm"], data_in["z_profile_mm"],
               color=vs.SERIES[0], lw=2.6, label=r"Reachable Tip Arc $s \in [0, 83]$ mm")
    ax_rz.plot([0, 0], [0, -25], color="#7d7b75", lw=4.5, label=r"Guide Cannula ($s \leq 0$)")
    
    r_max = data_in["max_r_mm"]
    z_end = data_in["z_profile_mm"][-1]
    ax_rz.scatter([r_max], [z_end], color="#0d47a1", s=60, zorder=5)
    ax_rz.annotate(
        rf"Tip ($83$ mm):" + "\n" +
        rf"$r = {r_max:.1f}$ mm" + "\n" +
        rf"$z = {z_end:.1f}$ mm" + "\n" +
        rf"$\psi = {data_in['max_psi_deg']:.1f}^\circ$",
        xy=(r_max, z_end), xytext=(-65, -28), textcoords="offset points",
        arrowprops=dict(arrowstyle="->", color=vs.INK_2, lw=1.0),
        fontsize=8.5, color=vs.INK_2, bbox=dict(boxstyle="round,pad=0.3", fc=vs.SURFACE, ec=vs.GRID)
    )
    
    R = 57.0
    th_circ = np.linspace(0, data_in["max_psi_deg"] * np.pi / 180, 50)
    xc = R * (1.0 - np.cos(th_circ))
    zc = R * np.sin(th_circ)
    ax_rz.plot(xc, zc, color=vs.NEUTRAL, ls=":", lw=1.2, label=rf"Arc Center ($R = {R:.0f}$ mm)")
    
    ax_rz.set_title(r"Meridian Profile ($r$-$z$ Plane)", pad=12)
    ax_rz.set_xlabel(r"Radial Reach $r = \sqrt{x^2+y^2}$ (mm)")
    ax_rz.set_ylabel(r"Axial Position $z$ (mm)")
    ax_rz.set_xlim(-5, 60)
    ax_rz.set_ylim(-10, 70)
    ax_rz.set_aspect("equal", adjustable="datalim")
    vs.tidy(ax_rz)
    ax_rz.legend(loc="lower right", frameon=False, fontsize=8.5)
    
    # Panel 3: 2D Top View (x-y Projection)
    ax_xy = fig.add_subplot(1, 3, 3)
    th_full = np.linspace(0, 2*np.pi, 100)
    
    for L, c in zip(sample_strokes, colors_bb):
        idx = int(np.argmin(np.abs(data_in["strokes_mm"] - L)))
        r_L = data_in["r_profile_mm"][idx]
        ax_xy.plot(r_L * np.cos(th_full), r_L * np.sin(th_full),
                   color=c, lw=1.8, label=rf"$L = {L:.0f}$ mm ($r = {r_L:.1f}$ mm)")
        
    ax_xy.scatter([0], [0], color="#3a3834", s=50, marker="o", label=r"Guide Origin ($s=0$)")
    ax_xy.set_title(r"Top View ($x$-$y$ Footprint)", pad=12)
    ax_xy.set_xlabel(r"$x$ (mm)")
    ax_xy.set_ylabel(r"$y$ (mm)")
    ax_xy.set_aspect("equal", adjustable="datalim")
    vs.tidy(ax_xy)
    ax_xy.legend(loc="upper right", frameon=False, fontsize=8.5)
    
    fig.suptitle(
        r"CT-SDR Workspace: Inner Tube Alone $\quad (q_{\mathrm{ITT}} \in [0, 83]\ \mathrm{mm},\ R = 57.0\ \mathrm{mm},\ \theta \in [0, 360^\circ])$",
        x=0.01, ha="left", fontsize=12, color=vs.INK
    )
    vs.caption(
        fig,
        rf"2-DoF workspace forms a 2D surface of revolution. Max radial reach: {r_max:.1f} mm, "
        rf"axial depth: {z_end:.1f} mm. Total reachable surface area: {data_in['surface_area_mm2']:.0f} mm$^2$. "
        rf"Curvature: $\kappa = 17.54$ 1/m."
    )
    fig.tight_layout(rect=(0, 0.05, 1, 0.95), w_pad=2.5)
    fig.savefig(out_path, dpi=300)
    plt.close(fig)
    print(f"Wrote {out_path}")


# ---------------------------------------------------------------------------
# Figure 2: Outer Tube Alone Workspace
# ---------------------------------------------------------------------------

def plot_outer_alone(data_out: Dict, out_path: str) -> None:
    """Create 3D and multi-view figures for the outer tube alone workspace."""
    vs.apply(use_cmr=True)
    fig = plt.figure(figsize=(10.2, 4.0))
    
    # Panel 1: 3D Perspective View
    ax3d = fig.add_subplot(1, 3, 1, projection="3d")
    draw_guide_tube_3d(ax3d, length_mm=20.0, radius_mm=3.0)
    
    surf = ax3d.plot_surface(
        data_out["X_grid_mm"], data_out["Y_grid_mm"], data_out["Z_grid_mm"],
        cmap="Oranges", alpha=0.6, edgecolor="none", shade=True, antialiased=True
    )
    
    sample_strokes = [10.0, 20.0, 35.0, 46.0]
    sample_angles = [0.0, np.pi/2, np.pi, 3*np.pi/2]
    colors_bb = ["#f7c297", "#f39958", "#eb6834", "#b84518"]
    
    for L, c in zip(sample_strokes, colors_bb):
        s_vals = np.linspace(0.0, L, 40)
        bb0 = single_tube_backbone(s_vals, 57.0, 78.54, rotation_rad=0.0)
        ax3d.plot(bb0[:, 0], bb0[:, 1], bb0[:, 2], color=c, lw=2.2,
                  label=rf"$L = {L:.0f}$ mm")
        for th in sample_angles[1:]:
            bb_th = single_tube_backbone(s_vals, 57.0, 78.54, rotation_rad=th)
            ax3d.plot(bb_th[:, 0], bb_th[:, 1], bb_th[:, 2], color=c, lw=1.2,
                      ls="--", alpha=0.45)
            
    ax3d.plot(data_out["X_grid_mm"][:, -1],
              data_out["Y_grid_mm"][:, -1],
              data_out["Z_grid_mm"][:, -1],
              color="#b84518", lw=2.4, ls="-", label=r"Tip Boundary ($46$ mm)")
              
    ax3d.set_title("3D Reachable Surface (Outer Alone)", pad=12)
    ax3d.set_xlabel(r"$x$ (mm)")
    ax3d.set_ylabel(r"$y$ (mm)")
    ax3d.set_zlabel(r"$z$ (mm)")
    vs.tidy3d(ax3d)
    ax3d.view_init(elev=26, azim=-55)
    
    all_pts = np.column_stack([data_out["X_grid_mm"].ravel(),
                               data_out["Y_grid_mm"].ravel(),
                               data_out["Z_grid_mm"].ravel()])
    vs.equal_aspect_3d(ax3d, all_pts)
    
    # Panel 2: 2D Meridian (r, z) Profile
    ax_rz = fig.add_subplot(1, 3, 2)
    ax_rz.plot(data_out["r_profile_mm"], data_out["z_profile_mm"],
               color=vs.SERIES[1], lw=2.6, label=r"Reachable Tip Arc $s \in [0, 46]$ mm")
    ax_rz.plot([0, 0], [0, -20], color="#7d7b75", lw=5.0, label=r"Guide Cannula ($s \leq 0$)")
    
    r_max = data_out["max_r_mm"]
    z_end = data_out["z_profile_mm"][-1]
    ax_rz.scatter([r_max], [z_end], color="#bf360c", s=60, zorder=5)
    ax_rz.annotate(
        rf"Tip ($46$ mm):" + "\n" +
        rf"$r = {r_max:.1f}$ mm" + "\n" +
        rf"$z = {z_end:.1f}$ mm" + "\n" +
        rf"$\psi = {data_out['max_psi_deg']:.1f}^\circ$",
        xy=(r_max, z_end), xytext=(-65, -28), textcoords="offset points",
        arrowprops=dict(arrowstyle="->", color=vs.INK_2, lw=1.0),
        fontsize=8.5, color=vs.INK_2, bbox=dict(boxstyle="round,pad=0.3", fc=vs.SURFACE, ec=vs.GRID)
    )
    
    R = 57.0
    th_circ = np.linspace(0, data_out["max_psi_deg"] * np.pi / 180, 50)
    xc = R * (1.0 - np.cos(th_circ))
    zc = R * np.sin(th_circ)
    ax_rz.plot(xc, zc, color=vs.NEUTRAL, ls=":", lw=1.2, label=rf"Arc Center ($R = {R:.0f}$ mm)")
    
    ax_rz.set_title(r"Meridian Profile ($r$-$z$ Plane)", pad=12)
    ax_rz.set_xlabel(r"Radial Reach $r = \sqrt{x^2+y^2}$ (mm)")
    ax_rz.set_ylabel(r"Axial Position $z$ (mm)")
    ax_rz.set_xlim(-5, 60)
    ax_rz.set_ylim(-10, 70)
    ax_rz.set_aspect("equal", adjustable="datalim")
    vs.tidy(ax_rz)
    ax_rz.legend(loc="lower right", frameon=False, fontsize=8.5)
    
    # Panel 3: 2D Top View (x-y Projection)
    ax_xy = fig.add_subplot(1, 3, 3)
    th_full = np.linspace(0, 2*np.pi, 100)
    
    for L, c in zip(sample_strokes, colors_bb):
        idx = int(np.argmin(np.abs(data_out["strokes_mm"] - L)))
        r_L = data_out["r_profile_mm"][idx]
        ax_xy.plot(r_L * np.cos(th_full), r_L * np.sin(th_full),
                   color=c, lw=1.8, label=rf"$L = {L:.0f}$ mm ($r = {r_L:.1f}$ mm)")
        
    ax_xy.scatter([0], [0], color="#3a3834", s=50, marker="o", label=r"Guide Origin ($s=0$)")
    ax_xy.set_title(r"Top View ($x$-$y$ Footprint)", pad=12)
    ax_xy.set_xlabel(r"$x$ (mm)")
    ax_xy.set_ylabel(r"$y$ (mm)")
    ax_xy.set_aspect("equal", adjustable="datalim")
    vs.tidy(ax_xy)
    ax_xy.legend(loc="upper right", frameon=False, fontsize=8.5)
    
    fig.suptitle(
        r"CT-SDR Workspace: Outer Tube Alone $\quad (q_{\mathrm{OTT}} \in [0, 46]\ \mathrm{mm},\ R = 57.0\ \mathrm{mm},\ \theta \in [0, 360^\circ])$",
        x=0.01, ha="left", fontsize=12, color=vs.INK
    )
    vs.caption(
        fig,
        rf"2-DoF workspace forms a 2D surface of revolution. Max radial reach: {r_max:.1f} mm, "
        rf"axial depth: {z_end:.1f} mm. Reachable surface area: {data_out['surface_area_mm2']:.0f} mm$^2$. "
        rf"Outer tube OD = 3.6 mm, wall = 0.25 mm."
    )
    fig.tight_layout(rect=(0, 0.05, 1, 0.95), w_pad=2.5)
    fig.savefig(out_path, dpi=300)
    plt.close(fig)
    print(f"Wrote {out_path}")


# ---------------------------------------------------------------------------
# Figure 3: Both Tubes Together (Coupled Cosserat Mechanics)
# ---------------------------------------------------------------------------

def plot_both_together(data_both: Dict, out_path: str) -> None:
    """Create 4-panel publication visualization for both tubes combined with landmark poses."""
    vs.apply(use_cmr=True)
    fig = plt.figure(figsize=(10.5, 7.2))
    
    extreme_poses = data_both["extreme_poses"]
    
    # -----------------------------------------------------------------------
    # Subplot 1: 3D Reachable Volume Envelope & Dense Point Cloud
    # -----------------------------------------------------------------------
    ax3d_1 = fig.add_subplot(2, 2, 1, projection="3d")
    draw_guide_tube_3d(ax3d_1, length_mm=20.0, radius_mm=2.5)
    
    ax3d_1.plot_surface(
        data_both["X_b_grid"], data_both["Y_b_grid"], data_both["Z_b_grid"],
        color="#1baf7a", alpha=0.30, edgecolor="none", shade=True
    )
    pts = data_both["pts_3d"]
    sc = ax3d_1.scatter(
        pts[:, 0], pts[:, 1], pts[:, 2],
        c=np.hypot(pts[:, 0], pts[:, 1]), cmap="viridis", s=2.5, alpha=0.35,
        depthshade=True
    )
    ax3d_1.set_title("(a) 3D Reachable Volume Envelope & Point Cloud", pad=10)
    ax3d_1.set_xlabel(r"$x$ (mm)")
    ax3d_1.set_ylabel(r"$y$ (mm)")
    ax3d_1.set_zlabel(r"$z$ (mm)")
    vs.tidy3d(ax3d_1)
    vs.equal_aspect_3d(ax3d_1, pts)
    ax3d_1.view_init(elev=24, azim=-50)
    
    # -----------------------------------------------------------------------
    # Subplot 2: 3D Cutaway with Overlaid Landmark Poses
    # -----------------------------------------------------------------------
    ax3d_2 = fig.add_subplot(2, 2, 2, projection="3d")
    draw_guide_tube_3d(ax3d_2, length_mm=20.0, radius_mm=2.5)
    
    X_b = data_both["X_b_grid"]
    Y_b = data_both["Y_b_grid"]
    Z_b = data_both["Z_b_grid"]
    n_rev = X_b.shape[0]
    n_half = n_rev // 2 + 1
    ax3d_2.plot_surface(
        X_b[:n_half, :], Y_b[:n_half, :], Z_b[:n_half, :],
        color="#1baf7a", alpha=0.18, edgecolor="none", shade=True
    )
    
    # Draw landmark robot backbones in 3D
    for ep in extreme_poses:
        p = ep["p"]
        c = ep["color"]
        ott = ep["ott"]
        s_approx = np.linspace(0, ep["itt"], len(p))
        split_idx = int(np.argmin(np.abs(s_approx - ott)))
        lbl = rf"{ep['name']}"
        
        # Outer tube (thicker)
        ax3d_2.plot(p[:split_idx+1, 0], p[:split_idx+1, 1], p[:split_idx+1, 2],
                    color=c, lw=3.2, solid_capstyle="round")
        # Inner tube extension (thinner)
        if split_idx < len(p) - 1:
            ax3d_2.plot(p[split_idx:, 0], p[split_idx:, 1], p[split_idx:, 2],
                        color=c, lw=1.8, ls="-")
            
        # Tip marker
        ax3d_2.scatter([ep["tip"][0]], [ep["tip"][1]], [ep["tip"][2]],
                       color=c, s=40, edgecolor="#ffffff", lw=1.0, zorder=6,
                       label=rf"Pose {ep['id']}: {lbl}")
                       
    ax3d_2.set_title("(b) 3D Cutaway with Landmark Boundary Poses", pad=10)
    ax3d_2.set_xlabel(r"$x$ (mm)")
    ax3d_2.set_ylabel(r"$y$ (mm)")
    ax3d_2.set_zlabel(r"$z$ (mm)")
    vs.tidy3d(ax3d_2)
    vs.equal_aspect_3d(ax3d_2, pts)
    ax3d_2.view_init(elev=20, azim=-40)
    ax3d_2.legend(loc="upper left", frameon=False, fontsize=7.5)
    
    # -----------------------------------------------------------------------
    # Subplot 3: 2D Meridian (r, z) Cross-Section with Overlaid Robot Poses
    # -----------------------------------------------------------------------
    ax_rz = fig.add_subplot(2, 2, 3)
    rz = data_both["rz_points"]
    b_rz = data_both["boundary_rz"]
    valid_tri = data_both.get("valid_triangles", None)
    
    # Fill reachable area strictly without overfill using valid alpha-shape triangles
    if valid_tri is not None and len(valid_tri) > 0:
        from matplotlib.collections import PolyCollection
        tri_polys = rz[valid_tri]
        poly_coll = PolyCollection(
            tri_polys, facecolors="#1baf7a", edgecolors="none", alpha=0.22,
            label=rf"Reachable Area $\mathcal{{D}}_{{rz}}$ (${data_both['area_2d_mm2']:.0f}\ \mathrm{{mm}}^2$)"
        )
        ax_rz.add_collection(poly_coll)
    else:
        ax_rz.fill(b_rz[:, 0], b_rz[:, 1], color="#1baf7a", alpha=0.22,
                   label=rf"Reachable Area $\mathcal{{D}}_{{rz}}$ (${data_both['area_2d_mm2']:.0f}\ \mathrm{{mm}}^2$)")
                   
    ax_rz.plot(b_rz[:, 0], b_rz[:, 1], color="#1baf7a", lw=1.8, label=r"Volume Boundary Envelope")
    
    # Sample point cloud in (r, z)
    ax_rz.scatter(rz[:, 0], rz[:, 1], c=rz[:, 0], cmap="viridis", s=4.5, alpha=0.35,
                  edgecolor="none", zorder=2)
                  
    # Guide cannula
    ax_rz.plot([0, 0], [0, -25], color="#7d7b75", lw=5.0, label=r"Rigid Cannula ($s \leq 0$)")
    
    # Badge offsets for landmark planar poses 1, 3, 4, 5
    badge_offsets = {
        "1": (10, -4),
        "3": (8, 6),
        "4": (8, -10),
        "5": (-14, -8),
    }
    
    for ep in extreme_poses:
        p = ep["p"]
        c = ep["color"]
        ott = ep["ott"]
        
        # Compute r(s) = sqrt(x(s)^2 + y(s)^2) and z(s) along the backbone
        r_backbone = np.hypot(p[:, 0], p[:, 1])
        z_backbone = p[:, 2]
        
        s_approx = np.linspace(0, ep["itt"], len(p))
        split_idx = int(np.argmin(np.abs(s_approx - ott)))
        
        # Outer tube (thicker)
        ax_rz.plot(r_backbone[:split_idx+1], z_backbone[:split_idx+1],
                   color=c, lw=3.2, solid_capstyle="round", zorder=4)
        # Inner tube extension (thinner)
        if split_idx < len(p) - 1:
            ax_rz.plot(r_backbone[split_idx:], z_backbone[split_idx:],
                       color=c, lw=1.8, ls="-", zorder=4)
            
        # Distal tip marker
        ax_rz.scatter([ep["r"]], [ep["z"]], color=c, s=55, edgecolor="#ffffff",
                      lw=1.2, zorder=6)
                      
        # Clean numbered badge (1, 3, 4, 5)
        dx, dy = badge_offsets.get(ep["id"], (8, 8))
        ax_rz.annotate(
            rf"$\mathbf{{{ep['id']}}}$",
            xy=(ep["r"], ep["z"]), xytext=(dx, dy), textcoords="offset points",
            fontsize=8.5, color="#ffffff",
            bbox=dict(boxstyle="circle,pad=0.20", fc=c, ec="#ffffff", lw=1.1),
            zorder=8
        )
        
    ax_rz.set_title(r"(c) Meridian Cross-Section $\mathcal{D}_{rz}$", pad=10)
    ax_rz.set_xlabel(r"Radial Reach $r = \sqrt{x^2+y^2}$ (mm)")
    ax_rz.set_ylabel(r"Axial Reach $z$ (mm)")
    ax_rz.set_xlim(-4, 58)
    ax_rz.set_ylim(-15, 88)
    ax_rz.set_aspect("equal", adjustable="box")
    vs.tidy(ax_rz)
    ax_rz.legend(loc="lower right", frameon=False, fontsize=7.5)
    
    # -----------------------------------------------------------------------
    # Subplot 4: Radial & Axial Reach Modulation vs Relative Roll
    # -----------------------------------------------------------------------
    ax_w = fig.add_subplot(2, 2, 4)
    cfg = data_both["configs"]
    if len(cfg) > 0:
        unique_otts = np.unique(cfg[:, 0])
        sel_otts = [unique_otts[0], unique_otts[len(unique_otts)//2], unique_otts[-1]]
        
        colors_sweep = vs.ramp(len(sel_otts))
        for ott_val, col in zip(sel_otts, colors_sweep):
            mask_ott = np.abs(cfg[:, 0] - ott_val) < 1e-3
            itt_max = cfg[mask_ott, 1].max()
            mask_pair = mask_ott & (np.abs(cfg[:, 1] - itt_max) < 1e-3)
            
            d_alphas = cfg[mask_pair, 2]
            r_vals = cfg[mask_pair, 3]
            z_vals = cfg[mask_pair, 4]
            
            lbl = rf"$q_{{\mathrm{{OTT}}}} = {ott_val:.0f}, q_{{\mathrm{{ITT}}}} = {itt_max:.0f}$ mm"
            ax_w.plot(d_alphas, r_vals, color=col, lw=1.8, marker="o", ms=3.5, label=lbl + r" ($r$)")
            ax_w.plot(d_alphas, z_vals, color=col, lw=1.4, ls="--", marker="s", ms=3.0, label=lbl + r" ($z$)")
            
    ax_w.set_title(r"(d) Tip Coordinates vs Roll $\Delta\alpha$", pad=10)
    ax_w.set_xlabel(r"Commanded Relative Roll $\Delta\alpha = q_{\mathrm{ITR}} - q_{\mathrm{OTR}}$ (deg)")
    ax_w.set_ylabel(r"Tip Coordinates $r, z$ (mm)")
    ax_w.set_xlim(-5, 365)
    vs.tidy(ax_w)
    ax_w.legend(loc="best", frameon=False, fontsize=7.2)
    
    fig.suptitle(
        "CT-SDR Workspace: Both Tubes Combined (Cosserat Mechanics, 4 DoFs)",
        x=0.01, ha="left", fontsize=11.5, color=vs.INK
    )
    vs.caption(
        fig,
        rf"Coupled 4-DoF mechanics forms an annular 3D workspace volume of {data_both['volume_mm3']:,.0f} mm$^3$ "
        rf"({data_both['volume_mm3']*1e-3:.1f} cm$^3$), swept by the non-convex 2D meridian cross-section $\mathcal{{D}}_{{rz}}$ "
        rf"({data_both['area_2d_mm2']:.0f} mm$^2$). Numbered badges denote planar boundary configurations in the meridian section: "
        rf"1 (Aligned maximum lateral reach), 3 (S-shape antagonistic inflection), 4 (Aligned intermediate extension), "
        rf"and 5 (Outer flush boundary). All landmark poses lie strictly in the meridian cutting plane ($x \equiv 0$)."
    )
    fig.tight_layout(rect=(0, 0.05, 0.98, 0.96), w_pad=2.8, h_pad=2.0)
    fig.savefig(out_path, dpi=300)
    plt.close(fig)
    print(f"Wrote {out_path}")


# ---------------------------------------------------------------------------
# Figure 4: Comparative Workspace Analysis
# ---------------------------------------------------------------------------

def plot_comparison(
    data_in: Dict,
    data_out: Dict,
    data_both: Dict,
    out_path: str,
) -> None:
    """Create comprehensive comparative figure overlaying all three cases."""
    vs.apply(use_cmr=True)
    fig = plt.figure(figsize=(10.5, 4.2))
    
    # Subplot 1: 3D Overlay Comparison
    ax3d = fig.add_subplot(1, 3, 1, projection="3d")
    draw_guide_tube_3d(ax3d, length_mm=25.0, radius_mm=3.2)
    
    ax3d.plot_surface(
        data_out["X_grid_mm"], data_out["Y_grid_mm"], data_out["Z_grid_mm"],
        color=vs.SERIES[1], alpha=0.35, edgecolor="none", shade=True
    )
    ax3d.plot_surface(
        data_in["X_grid_mm"], data_in["Y_grid_mm"], data_in["Z_grid_mm"],
        color=vs.SERIES[0], alpha=0.22, edgecolor="none", shade=True
    )
    ax3d.plot_wireframe(
        data_both["X_b_grid"], data_both["Y_b_grid"], data_both["Z_b_grid"],
        color="#1baf7a", alpha=0.45, lw=0.6, rstride=3, cstride=3
    )
    
    p_out = plt.Line2D([0], [0], color=vs.SERIES[1], lw=3, label=r"Outer Alone ($46$ mm)")
    p_in = plt.Line2D([0], [0], color=vs.SERIES[0], lw=3, label=r"Inner Alone ($83$ mm)")
    p_both = plt.Line2D([0], [0], color="#1baf7a", lw=3, label=r"Both Together (Solid Volume)")
    
    ax3d.set_title("(a) 3D Workspace Overlay", pad=12)
    ax3d.set_xlabel(r"$x$ (mm)")
    ax3d.set_ylabel(r"$y$ (mm)")
    ax3d.set_zlabel(r"$z$ (mm)")
    vs.tidy3d(ax3d)
    ax3d.view_init(elev=24, azim=-55)
    vs.equal_aspect_3d(ax3d, data_both["pts_3d"])
    ax3d.legend(handles=[p_out, p_in, p_both], loc="upper left", frameon=False, fontsize=8.5)
    
    # Subplot 2: Meridian (r, z) Profile Comparison
    ax_rz = fig.add_subplot(1, 3, 2)
    
    b_rz = data_both["boundary_rz"]
    ax_rz.fill(b_rz[:, 0], b_rz[:, 1], color="#1baf7a", alpha=0.22,
               label=rf"Both Together: 2D Slice (${data_both['area_2d_mm2']:.0f}\ \mathrm{{mm}}^2$)")
    ax_rz.plot(b_rz[:, 0], b_rz[:, 1], color="#1baf7a", lw=2.0)
    
    ax_rz.plot(data_out["r_profile_mm"], data_out["z_profile_mm"],
               color=vs.SERIES[1], lw=2.8, ls="--",
               label=rf"Outer Alone: 1D Arc ($r_{{max}} = {data_out['max_r_mm']:.1f}$ mm)")
               
    ax_rz.plot(data_in["r_profile_mm"], data_in["z_profile_mm"],
               color=vs.SERIES[0], lw=2.6, ls="-",
               label=rf"Inner Alone: 1D Arc ($r_{{max}} = {data_in['max_r_mm']:.1f}$ mm)")
               
    ax_rz.plot([0, 0], [0, -25], color="#7d7b75", lw=5.0, label=r"Guide ($s \leq 0$)")
    
    ax_rz.set_title(r"(b) Meridian Profile Comparison", pad=10)
    ax_rz.set_xlabel(r"Radial Distance $r$ (mm)")
    ax_rz.set_ylabel(r"Axial Distance $z$ (mm)")
    ax_rz.set_xlim(-5, 62)
    ax_rz.set_ylim(-15, 88)
    ax_rz.set_aspect("equal", adjustable="box")
    vs.tidy(ax_rz)
    ax_rz.legend(loc="lower right", frameon=False, fontsize=8.0)
    
    # Subplot 3: Metrics Table & Bar Charts
    ax_bar = fig.add_subplot(1, 3, 3)
    
    labels = ["Outer Alone\n(2-DoF)", "Inner Alone\n(2-DoF)", "Both Combined\n(4-DoF)"]
    r_maxes = [data_out["max_r_mm"], data_in["max_r_mm"], data_both["max_r_mm"]]
    z_maxes = [data_out["max_z_mm"], data_in["max_z_mm"], data_both["max_z_mm"]]
    
    x = np.arange(len(labels))
    width = 0.35
    
    rects1 = ax_bar.bar(x - width/2, r_maxes, width, label=r"Max Radial Reach $r_{\max}$ (mm)",
                        color="#5598e7", edgecolor=vs.SPINE)
    rects2 = ax_bar.bar(x + width/2, z_maxes, width, label=r"Max Axial Reach $z_{\max}$ (mm)",
                        color="#1c5cab", edgecolor=vs.SPINE)
                        
    for rect in rects1:
        h = rect.get_height()
        ax_bar.annotate(rf"${h:.1f}$", xy=(rect.get_x() + rect.get_width()/2, h),
                        xytext=(0, 3), textcoords="offset points", ha="center",
                        fontsize=8.5, color=vs.INK_2)
    for rect in rects2:
        h = rect.get_height()
        ax_bar.annotate(rf"${h:.1f}$", xy=(rect.get_x() + rect.get_width()/2, h),
                        xytext=(0, 3), textcoords="offset points", ha="center",
                        fontsize=8.5, color=vs.INK_2)
                        
    ax_bar.set_title("(c) Reach & Volume Metrics", pad=10)
    ax_bar.set_ylabel(r"Reach Distance (mm)")
    ax_bar.set_xticks(x)
    ax_bar.set_xticklabels(labels)
    ax_bar.set_ylim(0, 100)
    vs.tidy(ax_bar)
    ax_bar.legend(loc="upper left", frameon=False, fontsize=8.5)
    
    vol_text = (
        "Reachable 3D Volumes:\n" +
        rf"- Outer alone: 0 mm$^3$ (2D shell, {data_out['surface_area_mm2']:.0f} mm$^2$)" + "\n" +
        rf"- Inner alone: 0 mm$^3$ (2D shell, {data_in['surface_area_mm2']:.0f} mm$^2$)" + "\n" +
        rf"- Both combined: {data_both['volume_mm3']:,.0f} mm$^3$ ({data_both['volume_mm3']*1e-3:.1f} cm$^3$)"
    )
    ax_bar.text(
        0.5, 0.52, vol_text, transform=ax_bar.transAxes,
        fontsize=8.0, va="center", ha="center", color=vs.INK_2,
        bbox=dict(boxstyle="round,pad=0.5", fc=vs.SURFACE, ec=vs.GRID, lw=1.0)
    )
    
    fig.suptitle(
        "CT-SDR Workspace Comparison: Single Tubes vs 4-DoF Coupled Robot",
        x=0.01, ha="left", fontsize=12, color=vs.INK
    )
    vs.caption(
        fig,
        rf"Comparison confirms dimension expansion: single tubes trace hollow 2D surfaces of revolution, "
        rf"while concentric coupling fills a solid 3D annular volume of {data_both['volume_mm3']*1e-3:.1f} cm$^3$, expanding reach and interior dexterity."
    )
    fig.tight_layout(rect=(0, 0.05, 1, 0.95), w_pad=2.8)
    fig.savefig(out_path, dpi=300)
    plt.close(fig)
    print(f"Wrote {out_path}")


# ---------------------------------------------------------------------------
# Main Routine
# ---------------------------------------------------------------------------

def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--config", default=CONFIG_PATH, help="Path to ct_sdr.yaml")
    parser.add_argument("--inner-max", type=float, default=83.0, help="Max inner stroke (mm)")
    parser.add_argument("--outer-max", type=float, default=46.0, help="Max outer stroke (mm)")
    parser.add_argument("--out-dir", default=FIGURES_DIR, help="Directory to save figures")
    parser.add_argument("--quick", action="store_true", help="Use coarser grid for fast run")
    parser.add_argument("--recompute", action="store_true", help="Force full recomputation of coupled workspace")
    args = parser.parse_args(argv)
    
    os.makedirs(args.out_dir, exist_ok=True)
    
    with open(args.config, "r") as f:
        cfg = yaml.safe_load(f)
        
    E = float(cfg["material"]["youngs_modulus_gpa"]) * 1e9
    nu = float(cfg["material"]["poisson_ratio"])
    rho = float(cfg["material"].get("density_kg_m3", 6450.0))
    
    t_out = cfg["tubes"]["outer"]
    od_out = float(t_out["outer_diameter_mm"]) * 1e-3
    wall_out = float(t_out["wall_thickness_mm"]) * 1e-3
    len_out = float(t_out["length_mm"]) * 1e-3
    roc_out = float(t_out["radius_of_curvature_mm"]) * 1e-3
    tube_outer = Tube(
        name="outer", outer_diameter=od_out, inner_diameter=od_out - 2*wall_out,
        length=len_out, curved_length=float(t_out["curved_length_mm"]) * 1e-3,
        curvature=1.0 / roc_out, youngs_modulus=E, poisson_ratio=nu, density=rho
    )
    
    t_in = cfg["tubes"]["inner"]
    od_in = float(t_in["outer_diameter_mm"]) * 1e-3
    wall_in = float(t_in["wall_thickness_mm"]) * 1e-3
    len_in = float(t_in["length_mm"]) * 1e-3
    roc_in = float(t_in["radius_of_curvature_mm"]) * 1e-3
    curved_in = max(float(t_in["curved_length_mm"]), args.inner_max + 1.0) * 1e-3
    
    tube_inner = Tube(
        name="inner", outer_diameter=od_in, inner_diameter=od_in - 2*wall_in,
        length=len_in, curved_length=curved_in,
        curvature=1.0 / roc_in, youngs_modulus=E, poisson_ratio=nu, density=rho
    )
    
    robot = CTSDR(tube_outer, tube_inner)
    
    print("=" * 70)
    print("CT-SDR 3D ROBOT WORKSPACE ANALYSIS (Publication Paper Style: cmr10)")
    print(f"Inner Tube Stroke: 0.0 to {args.inner_max:.1f} mm  (R = {roc_in*1e3:.1f} mm)")
    print(f"Outer Tube Stroke: 0.0 to {args.outer_max:.1f} mm  (R = {roc_out*1e3:.1f} mm)")
    print(f"Rotations: Continuous 360 deg")
    print("=" * 70)
    
    # 1. Inner Tube Alone
    print("\n[1/3] Computing Inner Tube Alone Workspace...")
    data_in = compute_single_tube_workspace(
        max_stroke_mm=args.inner_max,
        radius_of_curvature_mm=roc_in * 1e3,
        curved_length_mm=curved_in * 1e3,
    )
    fig1_path = os.path.join(args.out_dir, "workspace_inner_alone.png")
    plot_inner_alone(data_in, fig1_path)
    
    # 2. Outer Tube Alone
    print("\n[2/3] Computing Outer Tube Alone Workspace...")
    data_out = compute_single_tube_workspace(
        max_stroke_mm=args.outer_max,
        radius_of_curvature_mm=roc_out * 1e3,
        curved_length_mm=float(t_out["curved_length_mm"]),
    )
    fig2_path = os.path.join(args.out_dir, "workspace_outer_alone.png")
    plot_outer_alone(data_out, fig2_path)
    
    # 3. Both Tubes Together (Coupled Cosserat Mechanics)
    print("\n[3/3] Computing Coupled Robot Workspace & Extreme Poses...")
    data_both = compute_coupled_workspace(
        robot=robot,
        max_ott_mm=args.outer_max,
        max_itt_mm=args.inner_max,
        quick=args.quick,
        recompute=args.recompute,
    )
    fig3_path = os.path.join(args.out_dir, "workspace_both_tubes.png")
    plot_both_together(data_both, fig3_path)
    
    # 4. Comparative Analysis
    print("\n[4/4] Generating Comparative Workspace Figure...")
    fig4_path = os.path.join(args.out_dir, "workspace_comparison.png")
    plot_comparison(data_in, data_out, data_both, fig4_path)
    
    print("\n" + "=" * 70)
    print("WORKSPACE ANALYSIS SUMMARY:")
    print(f"  • Inner Alone (83 mm stroke):  r_max = {data_in['max_r_mm']:.2f} mm,  z_max = {data_in['max_z_mm']:.2f} mm,  Area = {data_in['surface_area_mm2']:.1f} mm^2")
    print(f"  • Outer Alone (46 mm stroke):  r_max = {data_out['max_r_mm']:.2f} mm,  z_max = {data_out['max_z_mm']:.2f} mm,  Area = {data_out['surface_area_mm2']:.1f} mm^2")
    print(f"  • Both Combined (4-DoF):       r_max = {data_both['max_r_mm']:.2f} mm,  z_max = {data_both['max_z_mm']:.2f} mm,  Volume = {data_both['volume_mm3']:,.1f} mm^3 ({data_both['volume_mm3']*1e-3:.2f} cm^3)")
    print("=" * 70)


if __name__ == "__main__":
    main()
