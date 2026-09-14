#!/usr/bin/env python3
"""
interactive_workspace_3d.py -- Interactive 3D CT-SDR Robot Workspace Viewer.

Opens an interactive Matplotlib 3D window allowing the user to rotate, pan,
and zoom through the 3D workspace of the 4-DoF Concentric Tube Steerable
Drilling Robot (CT-SDR).

Displays:
    - 3D Bounding Volume Envelope (translucent emerald green shell)
    - 3D Reachable Tip Point Cloud
    - Rigid Stainless-Steel Guide Cannula
    - 3D Robot Backbones for 6 Extreme Poses (Aligned, Orthogonal, Antagonistic,
      Pre-Snap Max Axial, Outer Flush, Base Emergence)

Interactive Controls:
    - Left Click + Drag: Rotate 3D orientation (azimuth & elevation)
    - Right Click + Drag / Scroll: Zoom in / out
    - Middle Click + Drag: Pan
    - Key '1': Top View (x-y)
    - Key '2': Front View (y-z)
    - Key '3': Side View (x-z)
    - Key '4': Isometric 3D View
    - Key 'r': Reset view to default
    - Key 's': Save current view as 'interactive_view.png' (300 DPI)

Usage:
    ./venv/bin/python3 modeling/examples/interactive_workspace_3d.py
"""

from __future__ import annotations

import os
import sys
import numpy as np

# Ensure interactive GUI backend is used (e.g. MacOSX or TkAgg)
import matplotlib
matplotlib._interactive_mode = True
for gui_backend in ["MacOSX", "TkAgg", "QtAgg", "Qt5Agg"]:
    try:
        matplotlib.use(gui_backend, force=True)
        break
    except Exception:
        pass

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

MODELING_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if MODELING_DIR not in sys.path:
    sys.path.insert(0, MODELING_DIR)

import examples.vizstyle as vs
from ctr import CTSDR, Joints
from ctr.tube import Tube
import yaml

DATA_CACHE_PATH = os.path.join(MODELING_DIR, "figures", "workspace_data.npz")
CONFIG_PATH = os.path.join(MODELING_DIR, "config", "ct_sdr.yaml")


def load_or_compute_workspace_data():
    """Load cached workspace data or compute on the fly."""
    if os.path.exists(DATA_CACHE_PATH):
        print(f"Loading workspace data from cache: {DATA_CACHE_PATH} ...")
        cached = np.load(DATA_CACHE_PATH)
        return {
            "rz_points": cached["rz_points"],
            "boundary_rz": cached["boundary_rz"],
            "X_b_grid": cached["X_b_grid"],
            "Y_b_grid": cached["Y_b_grid"],
            "Z_b_grid": cached["Z_b_grid"],
            "pts_3d": cached["pts_3d"],
            "volume_mm3": float(cached["volume_mm3"]),
        }
    
    print("No cache found. Running workspace computation...")
    from examples.workspace_analysis import compute_coupled_workspace
    with open(CONFIG_PATH, "r") as f:
        cfg = yaml.safe_load(f)
        
    E = float(cfg["material"]["youngs_modulus_gpa"]) * 1e9
    nu = float(cfg["material"]["poisson_ratio"])
    rho = float(cfg["material"].get("density_kg_m3", 6450.0))
    
    t_out = cfg["tubes"]["outer"]
    od_out = float(t_out["outer_diameter_mm"]) * 1e-3
    wall_out = float(t_out["wall_thickness_mm"]) * 1e-3
    tube_outer = Tube(
        name="outer", outer_diameter=od_out, inner_diameter=od_out - 2*wall_out,
        length=float(t_out["length_mm"]) * 1e-3, curved_length=float(t_out["curved_length_mm"]) * 1e-3,
        curvature=1.0 / (float(t_out["radius_of_curvature_mm"]) * 1e-3),
        youngs_modulus=E, poisson_ratio=nu, density=rho
    )
    
    t_in = cfg["tubes"]["inner"]
    od_in = float(t_in["outer_diameter_mm"]) * 1e-3
    wall_in = float(t_in["wall_thickness_mm"]) * 1e-3
    tube_inner = Tube(
        name="inner", outer_diameter=od_in, inner_diameter=od_in - 2*wall_in,
        length=float(t_in["length_mm"]) * 1e-3, curved_length=89.54 * 1e-3,
        curvature=1.0 / (float(t_in["radius_of_curvature_mm"]) * 1e-3),
        youngs_modulus=E, poisson_ratio=nu, density=rho
    )
    robot = CTSDR(tube_outer, tube_inner)
    return compute_coupled_workspace(robot, quick=True)


def compute_extreme_poses_for_viewer(robot: CTSDR):
    """Compute the 4 landmark planar robot poses (all in x=0 plane)."""
    poses = []
    
    # 1. Aligned Max Lateral (0 deg, 46/83 mm)
    s1 = robot.solve(Joints(ott=46.0, itt=83.0, otr=0.0, itr=0.0))
    t1 = s1.tip_position * 1e3
    poses.append({
        "id": "1",
        "name": r"Pose 1: Aligned Max Lateral (46/83)",
        "ott": 46.0,
        "itt": 83.0,
        "p": s1.p * 1e3,
        "tip": t1,
        "r": float(np.hypot(t1[0], t1[1])),
        "z": float(t1[2]),
        "color": "#2a78d6",  # Blue
    })
    
    # 3. S-Shape Antagonistic Inflection (180 deg, 46/83 mm)
    s3 = robot.solve(Joints(ott=46.0, itt=83.0, otr=0.0, itr=180.0))
    t3 = s3.tip_position * 1e3
    poses.append({
        "id": "3",
        "name": r"Pose 3: S-Shape Inflection (46/83, 180$^\circ$)",
        "ott": 46.0,
        "itt": 83.0,
        "p": s3.p * 1e3,
        "tip": t3,
        "r": float(np.hypot(t3[0], t3[1])),
        "z": float(t3[2]),
        "color": "#8e24aa",  # Purple
    })
    
    # 4. Aligned Intermediate Extension (0 deg, 46/65 mm)
    s4 = robot.solve(Joints(ott=46.0, itt=65.0, otr=0.0, itr=0.0))
    t4 = s4.tip_position * 1e3
    poses.append({
        "id": "4",
        "name": r"Pose 4: Aligned Intermediate (46/65)",
        "ott": 46.0,
        "itt": 65.0,
        "p": s4.p * 1e3,
        "tip": t4,
        "r": float(np.hypot(t4[0], t4[1])),
        "z": float(t4[2]),
        "color": "#d32f2f",  # Crimson
    })
    
    # 5. Outer Flush Boundary (0 deg, 46/46 mm)
    s5 = robot.solve(Joints(ott=46.0, itt=46.0, otr=0.0, itr=0.0))
    t5 = s5.tip_position * 1e3
    poses.append({
        "id": "5",
        "name": r"Pose 5: Outer Flush (46/46)",
        "ott": 46.0,
        "itt": 46.0,
        "p": s5.p * 1e3,
        "tip": t5,
        "r": float(np.hypot(t5[0], t5[1])),
        "z": float(t5[2]),
        "color": "#00897b",  # Teal
    })
    
    return poses


def draw_guide_tube_3d(ax, length_mm: float = 25.0, radius_mm: float = 3.2) -> None:
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


def main():
    print("=" * 70)
    print("CT-SDR INTERACTIVE 3D WORKSPACE VIEWER")
    print("=" * 70)
    print("Controls:")
    print("  • Click & Drag Left Mouse Button: Rotate 3D view")
    print("  • Right Click / Scroll: Zoom in and out")
    print("  • Middle Click: Pan")
    print("  • Press '1': Top View (x-y)")
    print("  • Press '2': Front View (y-z)")
    print("  • Press '3': Side View (x-z)")
    print("  • Press '4': Isometric 3D View")
    print("  • Press 'r': Reset view")
    print("  • Press 's': Save 3D screenshot to 'interactive_view.png'")
    print("=" * 70)
    
    # Apply paper styling
    vs.apply(use_cmr=True)
    
    # Load workspace data
    data = load_or_compute_workspace_data()
    
    # Build robot and compute extreme poses
    with open(CONFIG_PATH, "r") as f:
        cfg = yaml.safe_load(f)
    E = float(cfg["material"]["youngs_modulus_gpa"]) * 1e9
    nu = float(cfg["material"]["poisson_ratio"])
    rho = float(cfg["material"].get("density_kg_m3", 6450.0))
    t_out = cfg["tubes"]["outer"]
    od_out = float(t_out["outer_diameter_mm"]) * 1e-3
    wall_out = float(t_out["wall_thickness_mm"]) * 1e-3
    tube_outer = Tube("outer", od_out, od_out - 2*wall_out, float(t_out["length_mm"])*1e-3,
                      float(t_out["curved_length_mm"])*1e-3, 1.0/(float(t_out["radius_of_curvature_mm"])*1e-3),
                      E, nu, rho)
    t_in = cfg["tubes"]["inner"]
    od_in = float(t_in["outer_diameter_mm"]) * 1e-3
    wall_in = float(t_in["wall_thickness_mm"]) * 1e-3
    tube_inner = Tube("inner", od_in, od_in - 2*wall_in, float(t_in["length_mm"])*1e-3,
                      89.54*1e-3, 1.0/(float(t_in["radius_of_curvature_mm"])*1e-3),
                      E, nu, rho)
    robot = CTSDR(tube_outer, tube_inner)
    extreme_poses = compute_extreme_poses_for_viewer(robot)
    
    # Create interactive figure
    fig = plt.figure(figsize=(12.0, 9.0))
    ax = fig.add_subplot(1, 1, 1, projection="3d")
    
    # 1. Guide cannula
    draw_guide_tube_3d(ax, length_mm=25.0, radius_mm=3.2)
    
    # 2. 3D Bounding Volume Shell (Translucent Emerald Green)
    shell = ax.plot_surface(
        data["X_b_grid"], data["Y_b_grid"], data["Z_b_grid"],
        color="#1baf7a", alpha=0.28, edgecolor="none", shade=True
    )
    
    # 3. 3D Scatter Point Cloud
    pts = data["pts_3d"]
    # sc = ax.scatter(
    #     pts[:, 0], pts[:, 1], pts[:, 2],
    #     c=np.hypot(pts[:, 0], pts[:, 1]), cmap="viridis", s=4, alpha=0.35,
    #     depthshade=True
    # )
    
    # 4. Extreme 3D Robot Poses
    for ep in extreme_poses:
        p = ep["p"]
        c = ep["color"]
        ott = ep["ott"]
        s_approx = np.linspace(0, ep["itt"], len(p))
        split_idx = int(np.argmin(np.abs(s_approx - ott)))
        
        # Outer tube (thicker)
        ax.plot(p[:split_idx+1, 0], p[:split_idx+1, 1], p[:split_idx+1, 2],
                color=c, lw=4.5, solid_capstyle="round")
        # Inner tube extension (thinner)
        if split_idx < len(p) - 1:
            ax.plot(p[split_idx:, 0], p[split_idx:, 1], p[split_idx:, 2],
                    color=c, lw=2.4, ls="-")
            
        # Distal tip marker with white outline
        ax.scatter([ep["tip"][0]], [ep["tip"][1]], [ep["tip"][2]],
                   color=c, s=70, edgecolor="#ffffff", lw=1.5, zorder=6,
                   label=rf"Pose {ep['id']}: {ep['name']}")
                   
    ax.set_title(
        "CT-SDR Interactive 3D Workspace Viewer (Cosserat Mechanics)\n"
        "Left-click & drag to rotate | Right-click/scroll to zoom | Keys: 1 (Top), 2 (Front), 3 (Side), 4 (Iso), R (Reset)",
        fontsize=11, pad=15
    )
    ax.set_xlabel(r"$x$ (mm)", fontsize=10, labelpad=8)
    ax.set_ylabel(r"$y$ (mm)", fontsize=10, labelpad=8)
    ax.set_zlabel(r"$z$ (mm)", fontsize=10, labelpad=8)
    
    vs.tidy3d(ax)
    vs.equal_aspect_3d(ax, pts)
    ax.view_init(elev=24, azim=-50)
    # ax.legend(loc="upper left", frameon=False, fontsize=8.5)
    
    # Keyboard interaction handler
    def on_key(event):
        if event.key == "1":
            print("Switching to Top View (x-y)...")
            ax.view_init(elev=90, azim=-90)
            fig.canvas.draw_idle()
        elif event.key == "2":
            print("Switching to Front View (y-z)...")
            ax.view_init(elev=0, azim=-90)
            fig.canvas.draw_idle()
        elif event.key == "3":
            print("Switching to Side View (x-z)...")
            ax.view_init(elev=0, azim=0)
            fig.canvas.draw_idle()
        elif event.key == "4":
            print("Switching to Isometric View...")
            ax.view_init(elev=35.264, azim=-45)
            fig.canvas.draw_idle()
        elif event.key == "r":
            print("Resetting view...")
            ax.view_init(elev=24, azim=-50)
            fig.canvas.draw_idle()
        elif event.key == "s":
            out_snap = os.path.join(MODELING_DIR, "figures", "interactive_view.png")
            fig.savefig(out_snap, dpi=300)
            print(f"Saved 3D screenshot to {out_snap}")
            
    fig.canvas.mpl_connect("key_press_event", on_key)
    
    plt.tight_layout()
    print("Opening interactive 3D window (using backend:", plt.get_backend(), ")...")
    plt.show(block=False)

    while plt.get_fignums():
        plt.pause(0.1)


if __name__ == "__main__":
    main()
