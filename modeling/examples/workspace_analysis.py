#!/usr/bin/env python3
"""
workspace_analysis.py -- 3D robot workspace analysis for the CT-SDR.

Analyzes and visualizes the 3D reachable workspace of the 4-DoF Concentric
Tube Steerable Drilling Robot (CT-SDR) across three operational scenarios:
    1) Inner tube alone (stroke 0 to 83 mm, 360 deg rotation)
    2) Outer tube alone (stroke 0 to 46 mm, 360 deg rotation)
    3) Both tubes together (coupled Cosserat rod mechanics, 4 DoFs)

Generates four publication-quality figures in figures/:
    - workspace_inner_alone.png
    - workspace_outer_alone.png
    - workspace_both_tubes.png
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
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np
from scipy.spatial import ConvexHull, Delaunay

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
    
    # If deployed length s_max > curved_length, proximal straight section emerges:
    s_max = s_arr[-1] if len(s_arr) > 0 else 0.0
    s_straight = max(0.0, s_max - curved_length)
    
    for s in s_arr:
        if s <= s_straight:
            # In straight proximal run
            x = 0.0
            y = 0.0
            z = s
        else:
            # In curved distal run
            s_curv = s - s_straight
            # Unrotated bend in y-z plane bending toward -y
            y_unrot = -R * (1.0 - np.cos(s_curv / R))
            z_unrot = s_straight + R * np.sin(s_curv / R)
            x_unrot = 0.0
            
            # Rotate by rotation_rad about +z
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
    
    # Meridian profile: r(L) and z(L) for rotation = 0
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
    
    # 3D surface mesh: Grid of (stroke, theta)
    L_grid, TH_grid = np.meshgrid(strokes, thetas)
    
    S_straight_grid = np.maximum(0.0, L_grid - curved_length_mm)
    S_curv_grid = L_grid - S_straight_grid
    
    R_grid = R * (1.0 - np.cos(S_curv_grid / R))
    Z_grid = S_straight_grid + R * np.sin(S_curv_grid / R)
    X_grid = R_grid * np.sin(TH_grid)
    Y_grid = -R_grid * np.cos(TH_grid)
    
    # Surface area via numerical integration: int 2*pi*r*dL
    dL = np.gradient(strokes)
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
) -> Dict[str, np.ndarray]:
    """Compute the 3D reachable workspace of both tubes interacting via Cosserat mechanics.
    
    Exploits rotational axisymmetry about +z: solving on (q_ott, q_itt, Delta_alpha)
    yields the complete meridian cross-section D_{rz} = {(r, z)}, which is then
    revolved 360 deg to form the full 3D solid volume.
    """
    if quick:
        n_ott = max(6, n_ott // 2)
        n_itt = max(8, n_itt // 2)
        n_angles = max(7, n_angles // 2)
        
    print(f"Sampling coupled workspace: {n_ott} OTT x {n_itt} ITT x {n_angles} dAlpha ...")
    t0 = time.time()
    
    ott_values = np.linspace(2.0, max_ott_mm, n_ott)
    delta_angles_deg = np.linspace(0.0, 180.0, n_angles)
    
    rz_points = []
    configurations = []
    backbone_samples = []
    
    total_combinations = n_ott * n_itt
    combo_idx = 0
    
    for ott in ott_values:
        itt_values = np.linspace(ott, max_itt_mm, n_itt)
        for itt in itt_values:
            combo_idx += 1
            path = [Joints(ott=ott, itt=itt, otr=0.0, itr=a) for a in delta_angles_deg]
            try:
                sols = robot.solve_path(path)
                for a_deg, sol in zip(delta_angles_deg, sols):
                    tip = sol.tip_position * 1e3  # in mm
                    r = np.hypot(tip[0], tip[1])
                    z = tip[2]
                    rz_points.append([r, z])
                    configurations.append([ott, itt, a_deg, r, z])
                    
                # Save a few sample backbones for plotting
                if (abs(ott - ott_values[len(ott_values)//2]) < 3.0 and
                    abs(itt - itt_values[-1]) < 3.0):
                    for a_deg, sol in zip(delta_angles_deg, sols):
                        if a_deg in (0.0, 90.0, 180.0):
                            backbone_samples.append({
                                "ott": ott,
                                "itt": itt,
                                "dalpha": a_deg,
                                "p": sol.p * 1e3,
                            })
            except Exception as e:
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
    
    # Compute convex hull / boundary in (r, z) plane
    hull = ConvexHull(rz_pts)
    hull_rz = rz_pts[hull.vertices]
    
    # Sort hull vertices counterclockwise
    center = np.mean(hull_rz, axis=0)
    angles_hull = np.arctan2(hull_rz[:, 1] - center[1], hull_rz[:, 0] - center[0])
    sort_idx = np.argsort(angles_hull)
    boundary_rz = hull_rz[sort_idx]
    # Close polygon
    boundary_rz = np.vstack([boundary_rz, boundary_rz[0]])
    
    # 2D cross-section area via Green's theorem / polygon shoelace
    x_h = boundary_rz[:, 0]
    y_h = boundary_rz[:, 1]
    area_2d_mm2 = 0.5 * np.abs(np.dot(x_h[:-1], y_h[1:]) - np.dot(x_h[1:], y_h[:-1]))
    
    # 3D volume by Pappus's centroid theorem: V = 2 * pi * r_centroid * Area
    # Polygon centroid r_c:
    factor = (x_h[:-1] * y_h[1:] - x_h[1:] * y_h[:-1])
    r_centroid = (1.0 / (6.0 * area_2d_mm2)) * np.abs(np.sum((x_h[:-1] + x_h[1:]) * factor))
    volume_mm3 = 2.0 * np.pi * r_centroid * area_2d_mm2
    
    # Revolve boundary points to create 3D bounding mesh
    phi_rev = np.linspace(0.0, 2.0 * np.pi, n_rot_rev)
    R_b = boundary_rz[:, 0]
    Z_b = boundary_rz[:, 1]
    
    R_b_grid, PHI_grid = np.meshgrid(R_b, phi_rev)
    Z_b_grid, _ = np.meshgrid(Z_b, phi_rev)
    X_b_grid = R_b_grid * np.cos(PHI_grid)
    Y_b_grid = R_b_grid * np.sin(PHI_grid)
    
    # Dense point cloud in 3D by revolving all sample rz points
    # (sample subset of angles for clean visualization)
    angles_3d = np.linspace(0.0, 2.0 * np.pi, 24, endpoint=False)
    pts_3d = []
    for r, z in rz_pts[::3]:  # sub-sample for plot performance
        for phi in angles_3d:
            pts_3d.append([r * np.cos(phi), r * np.sin(phi), z])
    pts_3d = np.array(pts_3d)
    
    return {
        "rz_points": rz_pts,
        "configs": configs,
        "boundary_rz": boundary_rz,
        "area_2d_mm2": float(area_2d_mm2),
        "r_centroid_mm": float(r_centroid),
        "volume_mm3": float(volume_mm3),
        "min_r_mm": float(rz_pts[:, 0].min()),
        "max_r_mm": float(rz_pts[:, 0].max()),
        "min_z_mm": float(rz_pts[:, 1].min()),
        "max_z_mm": float(rz_pts[:, 1].max()),
        "X_b_grid": X_b_grid,
        "Y_b_grid": Y_b_grid,
        "Z_b_grid": Z_b_grid,
        "pts_3d": pts_3d,
        "backbone_samples": backbone_samples,
    }


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
    
    ax.plot_surface(X_cyl, Y_cyl, Z_cyl, color="#7d7b75", alpha=0.6,
                    edgecolor="none", shade=True, zorder=2)
    # Distal collar ring at z = 0
    ax.plot(radius_mm * np.cos(theta_cyl), radius_mm * np.sin(theta_cyl),
            np.zeros_like(theta_cyl), color="#3a3834", lw=2.0, zorder=3)


# ---------------------------------------------------------------------------
# Figure 1: Inner Tube Alone Workspace
# ---------------------------------------------------------------------------

def plot_inner_alone(data_in: Dict, out_path: str) -> None:
    """Create 3D and multi-view figures for the inner tube alone workspace."""
    vs.apply()
    fig = plt.figure(figsize=(14.0, 5.2))
    
    # Panel 1: 3D Perspective View
    ax3d = fig.add_subplot(1, 3, 1, projection="3d")
    draw_guide_tube_3d(ax3d, length_mm=25.0, radius_mm=2.2)
    
    # 3D Surface of Revolution
    surf = ax3d.plot_surface(
        data_in["X_grid_mm"], data_in["Y_grid_mm"], data_in["Z_grid_mm"],
        cmap="Blues", alpha=0.55, edgecolor="none", shade=True, antialiased=True
    )
    
    # Sample backbone curves
    sample_strokes = [20.0, 40.0, 60.0, 83.0]
    sample_angles = [0.0, np.pi/2, np.pi, 3*np.pi/2]
    colors_bb = [vs.RAMP[0], vs.RAMP[1], vs.RAMP[2], vs.RAMP[4]]
    
    for L, c in zip(sample_strokes, colors_bb):
        s_vals = np.linspace(0.0, L, 50)
        # Backbone along 0 deg
        bb0 = single_tube_backbone(s_vals, 57.0, 89.54, rotation_rad=0.0)
        ax3d.plot(bb0[:, 0], bb0[:, 1], bb0[:, 2], color=c, lw=2.0,
                  label=f"L = {L:.0f} mm")
        # Trace other quadrant backbones lightly
        for th in sample_angles[1:]:
            bb_th = single_tube_backbone(s_vals, 57.0, 89.54, rotation_rad=th)
            ax3d.plot(bb_th[:, 0], bb_th[:, 1], bb_th[:, 2], color=c, lw=1.2,
                      ls="--", alpha=0.5)
            
    # Maximum reach perimeter ring at L = 83 mm
    ring_idx = -1
    ax3d.plot(data_in["X_grid_mm"][:, ring_idx],
              data_in["Y_grid_mm"][:, ring_idx],
              data_in["Z_grid_mm"][:, ring_idx],
              color="#0d47a1", lw=2.5, ls="-", label="Tip Boundary (83 mm)")
              
    ax3d.set_title("3D Reachable Surface (Inner Alone)")
    ax3d.set_xlabel("x (mm)")
    ax3d.set_ylabel("y (mm)")
    ax3d.set_zlabel("z (mm)")
    vs.tidy3d(ax3d)
    ax3d.view_init(elev=26, azim=-55)
    
    # Equal aspect
    all_pts = np.column_stack([data_in["X_grid_mm"].ravel(),
                               data_in["Y_grid_mm"].ravel(),
                               data_in["Z_grid_mm"].ravel()])
    vs.equal_aspect_3d(ax3d, all_pts)
    
    # Panel 2: 2D Meridian (r, z) Profile
    ax_rz = fig.add_subplot(1, 3, 2)
    ax_rz.plot(data_in["r_profile_mm"], data_in["z_profile_mm"],
               color=vs.SERIES[0], lw=2.6, label="Reachable Tip Arc")
    ax_rz.plot([0, 0], [0, -25], color="#7d7b75", lw=4.0, label="Guide cannula")
    
    # Mark notable coordinates
    r_max = data_in["max_r_mm"]
    z_end = data_in["z_profile_mm"][-1]
    ax_rz.scatter([r_max], [z_end], color="#0d47a1", s=60, zorder=5)
    ax_rz.annotate(
        f"Tip (83 mm):\n$r = {r_max:.1f}$ mm\n$z = {z_end:.1f}$ mm\n$\\psi = {data_in['max_psi_deg']:.1f}^\\circ$",
        xy=(r_max, z_end), xytext=(-55, -25), textcoords="offset points",
        arrowprops=dict(arrowstyle="->", color=vs.INK_2, lw=1.0),
        fontsize=8.5, color=vs.INK_2, bbox=dict(boxstyle="round,pad=0.3", fc=vs.SURFACE, ec=vs.GRID)
    )
    
    # Radius of curvature circle indicator
    R = 57.0
    th_circ = np.linspace(0, data_in["max_psi_deg"] * np.pi / 180, 50)
    xc = R * (1.0 - np.cos(th_circ))
    zc = R * np.sin(th_circ)
    ax_rz.plot(xc, zc, color=vs.NEUTRAL, ls=":", lw=1.2)
    
    ax_rz.set_title("Meridian Profile (r-z Plane)")
    ax_rz.set_xlabel("Radial Distance $r = \\sqrt{x^2+y^2}$ (mm)")
    ax_rz.set_ylabel("Axial Position $z$ (mm)")
    ax_rz.set_xlim(-5, 60)
    ax_rz.set_ylim(-10, 70)
    ax_rz.set_aspect("equal", adjustable="datalim")
    vs.tidy(ax_rz)
    ax_rz.legend(loc="lower right", frameon=False, fontsize=8)
    
    # Panel 3: 2D Top View (x-y Projection)
    ax_xy = fig.add_subplot(1, 3, 3)
    th_full = np.linspace(0, 2*np.pi, 100)
    
    for L, c in zip(sample_strokes, colors_bb):
        idx = int(np.argmin(np.abs(data_in["strokes_mm"] - L)))
        r_L = data_in["r_profile_mm"][idx]
        ax_xy.plot(r_L * np.cos(th_full), r_L * np.sin(th_full),
                   color=c, lw=1.8, label=f"L = {L:.0f} mm (r = {r_L:.1f} mm)")
        
    ax_xy.scatter([0], [0], color="#3a3834", s=50, marker="o", label="Guide exit (s=0)")
    ax_xy.set_title("Top View (x-y Footprint)")
    ax_xy.set_xlabel("x (mm)")
    ax_xy.set_ylabel("y (mm)")
    ax_xy.set_aspect("equal", adjustable="datalim")
    vs.tidy(ax_xy)
    ax_xy.legend(loc="upper right", frameon=False, fontsize=8)
    
    fig.suptitle(
        f"CT-SDR Workspace: Inner Tube Alone   (Stroke = 0–83 mm, R = 57.0 mm, $\\theta \\in [0, 360^\\circ]$)",
        x=0.01, ha="left", fontsize=12, fontweight="semibold", color=vs.INK
    )
    vs.caption(
        fig,
        f"2-DoF workspace forms a 2D surface of revolution in $\\mathbb{{R}}^3$. Max radial reach: {r_max:.1f} mm, "
        f"axial depth: {z_end:.1f} mm. Total reachable surface area: {data_in['surface_area_mm2']:.0f} mm$^2$. "
        f"Curvature: $\\kappa = 17.54$ m$^{{-1}}$."
    )
    fig.tight_layout(rect=(0, 0.02, 1, 0.94))
    fig.savefig(out_path)
    plt.close(fig)
    print(f"Wrote {out_path}")


# ---------------------------------------------------------------------------
# Figure 2: Outer Tube Alone Workspace
# ---------------------------------------------------------------------------

def plot_outer_alone(data_out: Dict, out_path: str) -> None:
    """Create 3D and multi-view figures for the outer tube alone workspace."""
    vs.apply()
    fig = plt.figure(figsize=(14.0, 5.2))
    
    # Panel 1: 3D Perspective View
    ax3d = fig.add_subplot(1, 3, 1, projection="3d")
    draw_guide_tube_3d(ax3d, length_mm=20.0, radius_mm=3.0)
    
    # 3D Surface of Revolution (Amber/Orange)
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
                  label=f"L = {L:.0f} mm")
        for th in sample_angles[1:]:
            bb_th = single_tube_backbone(s_vals, 57.0, 78.54, rotation_rad=th)
            ax3d.plot(bb_th[:, 0], bb_th[:, 1], bb_th[:, 2], color=c, lw=1.2,
                      ls="--", alpha=0.5)
            
    # Tip boundary ring at L = 46 mm
    ax3d.plot(data_out["X_grid_mm"][:, -1],
              data_out["Y_grid_mm"][:, -1],
              data_out["Z_grid_mm"][:, -1],
              color="#b84518", lw=2.5, ls="-", label="Tip Boundary (46 mm)")
              
    ax3d.set_title("3D Reachable Surface (Outer Alone)")
    ax3d.set_xlabel("x (mm)")
    ax3d.set_ylabel("y (mm)")
    ax3d.set_zlabel("z (mm)")
    vs.tidy3d(ax3d)
    ax3d.view_init(elev=26, azim=-55)
    
    all_pts = np.column_stack([data_out["X_grid_mm"].ravel(),
                               data_out["Y_grid_mm"].ravel(),
                               data_out["Z_grid_mm"].ravel()])
    vs.equal_aspect_3d(ax3d, all_pts)
    
    # Panel 2: 2D Meridian (r, z) Profile
    ax_rz = fig.add_subplot(1, 3, 2)
    ax_rz.plot(data_out["r_profile_mm"], data_out["z_profile_mm"],
               color=vs.SERIES[1], lw=2.6, label="Reachable Tip Arc")
    ax_rz.plot([0, 0], [0, -20], color="#7d7b75", lw=5.0, label="Guide cannula")
    
    r_max = data_out["max_r_mm"]
    z_end = data_out["z_profile_mm"][-1]
    ax_rz.scatter([r_max], [z_end], color="#b84518", s=60, zorder=5)
    ax_rz.annotate(
        f"Tip (46 mm):\n$r = {r_max:.1f}$ mm\n$z = {z_end:.1f}$ mm\n$\\psi = {data_out['max_psi_deg']:.1f}^\\circ$",
        xy=(r_max, z_end), xytext=(-50, -25), textcoords="offset points",
        arrowprops=dict(arrowstyle="->", color=vs.INK_2, lw=1.0),
        fontsize=8.5, color=vs.INK_2, bbox=dict(boxstyle="round,pad=0.3", fc=vs.SURFACE, ec=vs.GRID)
    )
    
    ax_rz.set_title("Meridian Profile (r-z Plane)")
    ax_rz.set_xlabel("Radial Distance $r = \\sqrt{x^2+y^2}$ (mm)")
    ax_rz.set_ylabel("Axial Position $z$ (mm)")
    ax_rz.set_xlim(-3, 30)
    ax_rz.set_ylim(-10, 55)
    ax_rz.set_aspect("equal", adjustable="datalim")
    vs.tidy(ax_rz)
    ax_rz.legend(loc="lower right", frameon=False, fontsize=8)
    
    # Panel 3: 2D Top View (x-y Projection)
    ax_xy = fig.add_subplot(1, 3, 3)
    th_full = np.linspace(0, 2*np.pi, 100)
    
    for L, c in zip(sample_strokes, colors_bb):
        idx = int(np.argmin(np.abs(data_out["strokes_mm"] - L)))
        r_L = data_out["r_profile_mm"][idx]
        ax_xy.plot(r_L * np.cos(th_full), r_L * np.sin(th_full),
                   color=c, lw=1.8, label=f"L = {L:.0f} mm (r = {r_L:.1f} mm)")
        
    ax_xy.scatter([0], [0], color="#3a3834", s=50, marker="o", label="Guide exit (s=0)")
    ax_xy.set_title("Top View (x-y Footprint)")
    ax_xy.set_xlabel("x (mm)")
    ax_xy.set_ylabel("y (mm)")
    ax_xy.set_aspect("equal", adjustable="datalim")
    vs.tidy(ax_xy)
    ax_xy.legend(loc="upper right", frameon=False, fontsize=8)
    
    fig.suptitle(
        f"CT-SDR Workspace: Outer Tube Alone   (Stroke = 0–46 mm, R = 57.0 mm, $\\theta \\in [0, 360^\\circ]$)",
        x=0.01, ha="left", fontsize=12, fontweight="semibold", color=vs.INK
    )
    vs.caption(
        fig,
        f"2-DoF workspace forms a 2D surface of revolution. Max radial reach: {r_max:.1f} mm, "
        f"axial depth: {z_end:.1f} mm. Reachable surface area: {data_out['surface_area_mm2']:.0f} mm$^2$. "
        f"Outer tube OD = 3.6 mm, wall = 0.25 mm."
    )
    fig.tight_layout(rect=(0, 0.02, 1, 0.94))
    fig.savefig(out_path)
    plt.close(fig)
    print(f"Wrote {out_path}")


# ---------------------------------------------------------------------------
# Figure 3: Both Tubes Together (Coupled Cosserat Mechanics)
# ---------------------------------------------------------------------------

def plot_both_together(data_both: Dict, out_path: str) -> None:
    """Create 4-panel comprehensive visualization for both tubes combined."""
    vs.apply()
    fig = plt.figure(figsize=(15.0, 10.0))
    
    # Subplot 1: 3D Reachable Volume (Outer Envelope + Point Cloud)
    ax3d_1 = fig.add_subplot(2, 2, 1, projection="3d")
    draw_guide_tube_3d(ax3d_1, length_mm=25.0, radius_mm=3.2)
    
    # Bounding shell
    ax3d_1.plot_surface(
        data_both["X_b_grid"], data_both["Y_b_grid"], data_both["Z_b_grid"],
        color="#1baf7a", alpha=0.35, edgecolor="none", shade=True
    )
    # Scatter points of reachable tip positions
    pts = data_both["pts_3d"]
    sc = ax3d_1.scatter(
        pts[:, 0], pts[:, 1], pts[:, 2],
        c=np.hypot(pts[:, 0], pts[:, 1]), cmap="viridis", s=3, alpha=0.45,
        depthshade=True
    )
    ax3d_1.set_title("3D Reachable Volume Envelope")
    ax3d_1.set_xlabel("x (mm)")
    ax3d_1.set_ylabel("y (mm)")
    ax3d_1.set_zlabel("z (mm)")
    vs.tidy3d(ax3d_1)
    vs.equal_aspect_3d(ax3d_1, pts)
    ax3d_1.view_init(elev=24, azim=-50)
    
    # Subplot 2: 3D Cutaway with Representative Backbones
    ax3d_2 = fig.add_subplot(2, 2, 2, projection="3d")
    draw_guide_tube_3d(ax3d_2, length_mm=25.0, radius_mm=3.2)
    
    # Half envelope cutaway (y <= 0)
    mask_cut = data_both["Y_b_grid"] <= 1e-6
    X_cut = np.where(mask_cut, data_both["X_b_grid"], np.nan)
    Y_cut = np.where(mask_cut, data_both["Y_b_grid"], np.nan)
    Z_cut = np.where(mask_cut, data_both["Z_b_grid"], np.nan)
    
    ax3d_2.plot_surface(X_cut, Y_cut, Z_cut, color="#1baf7a", alpha=0.3,
                        edgecolor="none", shade=True)
                        
    # Plot representative backbones inside volume
    colors_bb = {"0.0": "#2a78d6", "90.0": "#eb6834", "180.0": "#8e24aa"}
    labels_bb = {
        "0.0": "Aligned ($\\Delta\\alpha=0^\\circ$): Max Bend",
        "90.0": "Orthogonal ($\\Delta\\alpha=90^\\circ$): 3D Curve",
        "180.0": "Antagonistic ($\\Delta\\alpha=180^\\circ$): Min Bend",
    }
    
    for bb in data_both["backbone_samples"]:
        k = f"{bb['dalpha']:.1f}"
        if k in colors_bb:
            p = bb["p"]
            c = colors_bb[k]
            lbl = labels_bb[k]
            ax3d_2.plot(p[:, 0], p[:, 1], p[:, 2], color=c, lw=2.6, label=lbl)
            ax3d_2.scatter([p[-1, 0]], [p[-1, 1]], [p[-1, 2]], color=c, s=50, zorder=6)
            
    ax3d_2.set_title("Interior Cutaway & Sample Backbones")
    ax3d_2.set_xlabel("x (mm)")
    ax3d_2.set_ylabel("y (mm)")
    ax3d_2.set_zlabel("z (mm)")
    vs.tidy3d(ax3d_2)
    vs.equal_aspect_3d(ax3d_2, pts)
    ax3d_2.view_init(elev=20, azim=-40)
    ax3d_2.legend(loc="upper left", frameon=False, fontsize=8)
    
    # Subplot 3: 2D Meridian (r, z) Reachable Cross-Section
    ax_rz = fig.add_subplot(2, 2, 3)
    rz = data_both["rz_points"]
    b_rz = data_both["boundary_rz"]
    
    # Shaded polygon for 2D cross-section
    ax_rz.fill(b_rz[:, 0], b_rz[:, 1], color="#1baf7a", alpha=0.25,
               label=f"Reachable Area $\\mathcal{{D}}_{{rz}}$ ({data_both['area_2d_mm2']:.0f} mm$^2$)")
    ax_rz.plot(b_rz[:, 0], b_rz[:, 1], color="#1baf7a", lw=2.2, label="Envelope Boundary")
    
    # Scatter of internal configurations
    ax_rz.scatter(rz[:, 0], rz[:, 1], c=rz[:, 0], cmap="viridis", s=8, alpha=0.6,
                  edgecolor="none", zorder=3)
                  
    # Guide cannula representation
    ax_rz.plot([0, 0], [0, -25], color="#7d7b75", lw=6.0, label="Rigid Guide (s=0)")
    
    # Mark centroid of meridian cross-section
    rc = data_both["r_centroid_mm"]
    zc = (data_both["min_z_mm"] + data_both["max_z_mm"]) / 2.0
    ax_rz.scatter([rc], [zc], color="#104281", marker="x", s=70, lw=2.0,
                  label=f"Centroid ($r_c = {rc:.1f}$ mm)")
                  
    ax_rz.set_title("Meridian Cross-Section $\\mathcal{D}_{rz}$ (r-z Plane)")
    ax_rz.set_xlabel("Radial Reach $r = \\sqrt{x^2+y^2}$ (mm)")
    ax_rz.set_ylabel("Axial Reach $z$ (mm)")
    ax_rz.set_xlim(-3, 60)
    ax_rz.set_ylim(-10, 75)
    ax_rz.set_aspect("equal", adjustable="datalim")
    vs.tidy(ax_rz)
    ax_rz.legend(loc="lower right", frameon=False, fontsize=8)
    
    # Subplot 4: Relative Angle Impact on Radial Extension & Windup
    ax_w = fig.add_subplot(2, 2, 4)
    cfg = data_both["configs"]
    # Group by (ott, itt) pairs and plot r vs delta_alpha
    unique_otts = np.unique(cfg[:, 0])
    sel_otts = [unique_otts[0], unique_otts[len(unique_otts)//2], unique_otts[-1]]
    
    colors_sweep = vs.ramp(len(sel_otts))
    for ott_val, col in zip(sel_otts, colors_sweep):
        mask_ott = np.abs(cfg[:, 0] - ott_val) < 1e-3
        itt_max = cfg[mask_ott, 1].max()
        mask_pair = mask_ott & (np.abs(cfg[:, 1] - itt_max) < 1e-3)
        
        d_alphas = cfg[mask_pair, 2]
        r_vals = cfg[mask_pair, 3]
        
        lbl = f"OTT = {ott_val:.0f} mm, ITT = {itt_max:.0f} mm"
        ax_w.plot(d_alphas, r_vals, color=col, lw=2.2, marker="o", ms=4, label=lbl)
        
    ax_w.set_title("Radial Modulation via Relative Roll $\\Delta\\alpha$")
    ax_w.set_xlabel("Commanded Relative Roll $\\Delta\\alpha = q_{\\mathrm{ITR}} - q_{\\mathrm{OTR}}$ (deg)")
    ax_w.set_ylabel("Distal Tip Radius $r$ (mm)")
    ax_w.set_xticks(np.arange(0, 181, 30))
    vs.tidy(ax_w)
    ax_w.legend(loc="upper right", frameon=False, fontsize=8.5)
    
    fig.suptitle(
        "CT-SDR Workspace: Both Tubes Combined (Cosserat Rod Model, 4 DoFs)",
        x=0.01, ha="left", fontsize=13, fontweight="semibold", color=vs.INK
    )
    vs.caption(
        fig,
        f"4-DoF coupled mechanics: OTT $\\in [0, 46]$ mm, ITT $\\in [0, 83]$ mm, OTR/ITR $\\in [0, 360^\\circ]$. "
        f"Total 3D workspace volume: $V = {data_both['volume_mm3']:,.0f}$ mm$^3$ ({data_both['volume_mm3']*1e-3:.1f} cm$^3$). "
        f"Radial span: [{data_both['min_r_mm']:.1f}, {data_both['max_r_mm']:.1f}] mm, "
        f"Axial span: [{data_both['min_z_mm']:.1f}, {data_both['max_z_mm']:.1f}] mm."
    )
    fig.tight_layout(rect=(0, 0.02, 1, 0.95))
    fig.savefig(out_path)
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
    vs.apply()
    fig = plt.figure(figsize=(15.0, 5.4))
    
    # Subplot 1: 3D Overlay Comparison
    ax3d = fig.add_subplot(1, 3, 1, projection="3d")
    draw_guide_tube_3d(ax3d, length_mm=25.0, radius_mm=3.2)
    
    # 1. Outer alone (Amber)
    ax3d.plot_surface(
        data_out["X_grid_mm"], data_out["Y_grid_mm"], data_out["Z_grid_mm"],
        color=vs.SERIES[1], alpha=0.35, edgecolor="none", shade=True
    )
    # 2. Inner alone (Blue)
    ax3d.plot_surface(
        data_in["X_grid_mm"], data_in["Y_grid_mm"], data_in["Z_grid_mm"],
        color=vs.SERIES[0], alpha=0.25, edgecolor="none", shade=True
    )
    # 3. Both together bounding boundary (Green wireframe/shell)
    ax3d.plot_wireframe(
        data_both["X_b_grid"], data_both["Y_b_grid"], data_both["Z_b_grid"],
        color="#1baf7a", alpha=0.45, lw=0.6, rstride=3, cstride=3
    )
    
    # Fake proxies for 3D legend
    p_out = plt.Line2D([0], [0], color=vs.SERIES[1], lw=3, label="Outer Alone (46 mm)")
    p_in = plt.Line2D([0], [0], color=vs.SERIES[0], lw=3, label="Inner Alone (83 mm)")
    p_both = plt.Line2D([0], [0], color="#1baf7a", lw=3, label="Both Together (Volume)")
    
    ax3d.set_title("3D Workspace Overlay")
    ax3d.set_xlabel("x (mm)")
    ax3d.set_ylabel("y (mm)")
    ax3d.set_zlabel("z (mm)")
    vs.tidy3d(ax3d)
    ax3d.view_init(elev=24, azim=-55)
    vs.equal_aspect_3d(ax3d, data_both["pts_3d"])
    ax3d.legend(handles=[p_out, p_in, p_both], loc="upper left", frameon=False, fontsize=8)
    
    # Subplot 2: Meridian (r, z) Profile Comparison
    ax_rz = fig.add_subplot(1, 3, 2)
    
    # Shaded both together
    b_rz = data_both["boundary_rz"]
    ax_rz.fill(b_rz[:, 0], b_rz[:, 1], color="#1baf7a", alpha=0.22,
               label=f"Both Together: 2D Cross-Section ({data_both['area_2d_mm2']:.0f} mm$^2$)")
    ax_rz.plot(b_rz[:, 0], b_rz[:, 1], color="#1baf7a", lw=2.0)
    
    # Outer alone arc
    ax_rz.plot(data_out["r_profile_mm"], data_out["z_profile_mm"],
               color=vs.SERIES[1], lw=2.6, ls="--",
               label=f"Outer Alone: 1D Arc ($r_{{max}} = {data_out['max_r_mm']:.1f}$ mm)")
               
    # Inner alone arc
    ax_rz.plot(data_in["r_profile_mm"], data_in["z_profile_mm"],
               color=vs.SERIES[0], lw=2.6, ls="-",
               label=f"Inner Alone: 1D Arc ($r_{{max}} = {data_in['max_r_mm']:.1f}$ mm)")
               
    ax_rz.plot([0, 0], [0, -25], color="#7d7b75", lw=5.0, label="Guide (s=0)")
    
    ax_rz.set_title("Meridian Cross-Section Comparison (r-z Plane)")
    ax_rz.set_xlabel("Radial Distance $r$ (mm)")
    ax_rz.set_ylabel("Axial Distance $z$ (mm)")
    ax_rz.set_xlim(-3, 60)
    ax_rz.set_ylim(-10, 75)
    ax_rz.set_aspect("equal", adjustable="datalim")
    vs.tidy(ax_rz)
    ax_rz.legend(loc="lower right", frameon=False, fontsize=8)
    
    # Subplot 3: Metrics Table & Bar Charts
    ax_bar = fig.add_subplot(1, 3, 3)
    
    labels = ["Outer Alone\n(2-DoF)", "Inner Alone\n(2-DoF)", "Both Combined\n(4-DoF)"]
    r_maxes = [data_out["max_r_mm"], data_in["max_r_mm"], data_both["max_r_mm"]]
    z_maxes = [data_out["max_z_mm"], data_in["max_z_mm"], data_both["max_z_mm"]]
    
    x = np.arange(len(labels))
    width = 0.35
    
    rects1 = ax_bar.bar(x - width/2, r_maxes, width, label="Max Radial Reach $r_{max}$ (mm)",
                        color="#5598e7", edgecolor=vs.SPINE)
    rects2 = ax_bar.bar(x + width/2, z_maxes, width, label="Max Axial Reach $z_{max}$ (mm)",
                        color="#1c5cab", edgecolor=vs.SPINE)
                        
    # Direct data labels
    for rect in rects1:
        h = rect.get_height()
        ax_bar.annotate(f"{h:.1f}", xy=(rect.get_x() + rect.get_width()/2, h),
                        xytext=(0, 3), textcoords="offset points", ha="center",
                        fontsize=8.5, color=vs.INK_2)
    for rect in rects2:
        h = rect.get_height()
        ax_bar.annotate(f"{h:.1f}", xy=(rect.get_x() + rect.get_width()/2, h),
                        xytext=(0, 3), textcoords="offset points", ha="center",
                        fontsize=8.5, color=vs.INK_2)
                        
    ax_bar.set_title("Reach Metrics & Volume Comparison")
    ax_bar.set_ylabel("Reach Distance (mm)")
    ax_bar.set_xticks(x)
    ax_bar.set_xticklabels(labels)
    ax_bar.set_ylim(0, 85)
    vs.tidy(ax_bar)
    ax_bar.legend(loc="upper left", frameon=False, fontsize=8.5)
    
    # Volume callout box inside bar chart
    vol_text = (
        "Reachable 3D Volumes:\n"
        f"• Outer alone: 0 mm$^3$ (2D surface, {data_out['surface_area_mm2']:.0f} mm$^2$)\n"
        f"• Inner alone: 0 mm$^3$ (2D surface, {data_in['surface_area_mm2']:.0f} mm$^2$)\n"
        f"• Both combined: {data_both['volume_mm3']:,.0f} mm$^3$ ({data_both['volume_mm3']*1e-3:.1f} cm$^3$)"
    )
    ax_bar.text(
        0.5, 0.45, vol_text, transform=ax_bar.transAxes,
        fontsize=8.0, va="center", ha="center", color=vs.INK_2,
        bbox=dict(boxstyle="round,pad=0.5", fc=vs.SURFACE, ec=vs.GRID, lw=1.0)
    )
    
    fig.suptitle(
        "CT-SDR Workspace Comparison: Single Tubes vs 4-DoF Coupled Robot",
        x=0.01, ha="left", fontsize=13, fontweight="semibold", color=vs.INK
    )
    vs.caption(
        fig,
        f"Comparison confirms dimension expansion: single tubes trace hollow 2D surfaces of revolution, "
        f"while concentric coupling fills a solid 3D annular volume of {data_both['volume_mm3']*1e-3:.1f} cm$^3$, expanding reach and interior dexterity."
    )
    fig.tight_layout(rect=(0, 0.02, 1, 0.95))
    fig.savefig(out_path)
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
    args = parser.parse_args(argv)
    
    os.makedirs(args.out_dir, exist_ok=True)
    
    # Load config and adjust inner tube curved length if needed so 83 mm is valid
    with open(args.config, "r") as f:
        cfg = yaml.safe_load(f)
        
    E = float(cfg["material"]["youngs_modulus_gpa"]) * 1e9
    nu = float(cfg["material"]["poisson_ratio"])
    rho = float(cfg["material"].get("density_kg_m3", 6450.0))
    
    # Outer tube
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
    
    # Inner tube: ensure curved length covers the 83 mm stroke (e.g. 90 deg arc at R=57 mm is 89.54 mm)
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
    print("CT-SDR 3D ROBOT WORKSPACE ANALYSIS")
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
    print("\n[3/3] Computing Coupled Robot Workspace...")
    data_both = compute_coupled_workspace(
        robot=robot,
        max_ott_mm=args.outer_max,
        max_itt_mm=args.inner_max,
        quick=args.quick,
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
