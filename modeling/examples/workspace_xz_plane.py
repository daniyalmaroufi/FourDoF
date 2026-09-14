#!/usr/bin/env python3
"""
workspace_xz_plane.py -- Bilateral Planar Cross-Section Analysis for CT-SDR.

Visualizes the complete planar cross-section of the 4-DoF Concentric Tube
Steerable Drilling Robot (CT-SDR) workspace in the x-z Cartesian plane (y = 0),
corresponding to the diameter slice across theta = 0 deg and theta = 180 deg.

Features:
    - Seamless bilateral reachable cross-section with no interior center seam.
    - Demonstrates both boundary limits and interior reach/dexterity.
    - 7 curated landmark configurations (4 boundary anchors + 3 interior poses).
    - Publication typography (cmr10) scaled up 50% for maximum paper readability.
    - Panel (a): 3D Workspace Volume Envelope with x-z Sectioning Plane (y = 0)
    - Panel (b): Complete Bilateral 2D x-z Planar Cross-Section with Overlaid Poses

Usage:
    ./venv/bin/python3 modeling/examples/workspace_xz_plane.py
"""

from __future__ import annotations

import argparse
import os
import sys
from typing import Dict, List, Tuple

import matplotlib
if not getattr(matplotlib, "_interactive_mode", False):
    try:
        matplotlib.use("Agg")
    except Exception:
        pass
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
import numpy as np
import yaml

# Ensure modeling root is in path
MODELING_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if MODELING_DIR not in sys.path:
    sys.path.insert(0, MODELING_DIR)

import examples.vizstyle as vs
from ctr import CTSDR, Joints
from ctr.tube import Tube
from examples.workspace_analysis import draw_guide_tube_3d

FIGURES_DIR = os.path.join(MODELING_DIR, "figures")
CONFIG_PATH = os.path.join(MODELING_DIR, "config", "ct_sdr.yaml")
DATA_CACHE_PATH = os.path.join(FIGURES_DIR, "workspace_data.npz")


def build_robot_from_yaml(config_path: str, max_inner_stroke: float = 83.0) -> CTSDR:
    """Build CTSDR instance from YAML configuration."""
    with open(config_path, "r") as f:
        cfg = yaml.safe_load(f)

    E = float(cfg["material"]["youngs_modulus_gpa"]) * 1e9
    nu = float(cfg["material"]["poisson_ratio"])
    rho = float(cfg["material"].get("density_kg_m3", 6450.0))

    t_out = cfg["tubes"]["outer"]
    od_out = float(t_out["outer_diameter_mm"]) * 1e-3
    wall_out = float(t_out["wall_thickness_mm"]) * 1e-3
    tube_outer = Tube(
        name="outer", outer_diameter=od_out, inner_diameter=od_out - 2 * wall_out,
        length=float(t_out["length_mm"]) * 1e-3, curved_length=float(t_out["curved_length_mm"]) * 1e-3,
        curvature=1.0 / (float(t_out["radius_of_curvature_mm"]) * 1e-3),
        youngs_modulus=E, poisson_ratio=nu, density=rho
    )

    t_in = cfg["tubes"]["inner"]
    od_in = float(t_in["outer_diameter_mm"]) * 1e-3
    wall_in = float(t_in["wall_thickness_mm"]) * 1e-3
    curved_in = max(float(t_in["curved_length_mm"]), max_inner_stroke + 1.0) * 1e-3
    tube_inner = Tube(
        name="inner", outer_diameter=od_in, inner_diameter=od_in - 2 * wall_in,
        length=float(t_in["length_mm"]) * 1e-3, curved_length=curved_in,
        curvature=1.0 / (float(t_in["radius_of_curvature_mm"]) * 1e-3),
        youngs_modulus=E, poisson_ratio=nu, density=rho
    )

    return CTSDR(tube_outer, tube_inner)


def compute_xz_poses(robot: CTSDR) -> List[Dict]:
    """Compute 7 curated planar configurations (4 boundary anchors + 3 interior dexterity poses).
    
    All backbones lie strictly in the x-z plane (y = 0).
    """
    poses = []

    # ---------------- Boundary Anchors ----------------
    # 1. Aligned Maximum Lateral Reach (+x): otr=90, itr=90
    sol_1 = robot.solve(Joints(ott=46.0, itt=83.0, otr=90.0, itr=90.0))
    poses.append({
        "id": "1",
        "name": r"Aligned Max Reach ($+x$)",
        "desc": r"Outer Lateral Boundary: $q_{\mathrm{OTT}}=46, q_{\mathrm{ITT}}=83\ \mathrm{mm}$",
        "ott": 46.0,
        "itt": 83.0,
        "p": sol_1.p * 1e3,
        "tip": sol_1.tip_position * 1e3,
        "color": "#2a78d6",  # Blue
        "badge_offset": (12, -5),
    })

    # 2. Aligned Symmetric Lateral Reach (-x): otr=270, itr=270
    sol_2 = robot.solve(Joints(ott=46.0, itt=83.0, otr=270.0, itr=270.0))
    poses.append({
        "id": "2",
        "name": r"Symmetric Max Reach ($-x$)",
        "desc": r"Opposite Lateral Boundary: $q_{\mathrm{OTT}}=46, q_{\mathrm{ITT}}=83\ \mathrm{mm}$",
        "ott": 46.0,
        "itt": 83.0,
        "p": sol_2.p * 1e3,
        "tip": sol_2.tip_position * 1e3,
        "color": "#2a78d6",  # Blue
        "badge_offset": (-18, -5),
    })

    # 3. S-Shape Antagonistic Inflection Crossing Centerline (+x -> -x): otr=90, itr=270
    sol_3 = robot.solve(Joints(ott=29.0, itt=83.0, otr=90.0, itr=270.0))
    poses.append({
        "id": "3",
        "name": r"S-Shape Centerline Crossing",
        "desc": r"Inflected S-curve crossing $x=0$: $q_{\mathrm{OTT}}=29, q_{\mathrm{ITT}}=83\ \mathrm{mm}$",
        "ott": 29.0,
        "itt": 83.0,
        "p": sol_3.p * 1e3,
        "tip": sol_3.tip_position * 1e3,
        "color": "#8e24aa",  # Purple
        "badge_offset": (-18, 9),
    })

    # 4. S-Shape Full-Stroke Top Reach (+x Lobe): otr=90, itr=270
    sol_4 = robot.solve(Joints(ott=46.0, itt=83.0, otr=90.0, itr=270.0))
    poses.append({
        "id": "4",
        "name": r"S-Shape Top Axial Reach",
        "desc": r"Full-Stroke S-curve: $q_{\mathrm{OTT}}=46, q_{\mathrm{ITT}}=83\ \mathrm{mm}$",
        "ott": 46.0,
        "itt": 83.0,
        "p": sol_4.p * 1e3,
        "tip": sol_4.tip_position * 1e3,
        "color": "#e65100",  # Dark Orange
        "badge_offset": (12, 8),
    })

    # ---------------- Interior Dexterity Poses ----------------
    # 5. Deep Right Interior Reach: otr=270, itr=90 (antagonistic)
    sol_5 = robot.solve(Joints(ott=15.0, itt=75.0, otr=270.0, itr=90.0))
    poses.append({
        "id": "5",
        "name": r"Right Interior Reach",
        "desc": r"Deep Interior: $q_{\mathrm{OTT}}=15, q_{\mathrm{ITT}}=75\ \mathrm{mm}, x=+20.3\ \mathrm{mm}$",
        "ott": 15.0,
        "itt": 75.0,
        "p": sol_5.p * 1e3,
        "tip": sol_5.tip_position * 1e3,
        "color": "#00897b",  # Teal
        "badge_offset": (12, 2),
    })

    # 6. Deep Left Interior Reach: otr=90, itr=270 (antagonistic)
    sol_6 = robot.solve(Joints(ott=15.0, itt=75.0, otr=90.0, itr=270.0))
    poses.append({
        "id": "6",
        "name": r"Left Interior Reach",
        "desc": r"Deep Interior: $q_{\mathrm{OTT}}=15, q_{\mathrm{ITT}}=75\ \mathrm{mm}, x=-20.3\ \mathrm{mm}$",
        "ott": 15.0,
        "itt": 75.0,
        "p": sol_6.p * 1e3,
        "tip": sol_6.tip_position * 1e3,
        "color": "#d97706",  # Amber/Ochre
        "badge_offset": (-18, 2),
    })

    # 7. Central-Axis Deep Interior Reach: otr=90, itr=270 (antagonistic)
    sol_7 = robot.solve(Joints(ott=25.0, itt=60.0, otr=90.0, itr=270.0))
    poses.append({
        "id": "7",
        "name": r"Central Axis Interior",
        "desc": r"Centerline Reach ($x \approx 0$): $q_{\mathrm{OTT}}=25, q_{\mathrm{ITT}}=60\ \mathrm{mm}, z=59.1\ \mathrm{mm}$",
        "ott": 25.0,
        "itt": 60.0,
        "p": sol_7.p * 1e3,
        "tip": sol_7.tip_position * 1e3,
        "color": "#dc2626",  # Crimson
        "badge_offset": (12, -6),
    })

    return poses


def plot_xz_workspace(
    data_both: Dict,
    poses_xz: List[Dict],
    out_path: str,
) -> None:
    """Generate 2-panel publication figure for x-z planar cross-section with 50% larger typography."""
    vs.apply(use_cmr=True)

    # Configure 50% larger typography globally
    plt.rc("font", size=12.5)
    plt.rc("axes", titlesize=14.5, labelsize=14.0)
    plt.rc("xtick", labelsize=12.0)
    plt.rc("ytick", labelsize=12.0)
    plt.rc("legend", fontsize=11.0)
    plt.rc("figure", titlesize=16.0)

    fig = plt.figure(figsize=(12.0, 5.8))

    # -----------------------------------------------------------------------
    # Subplot 1: 3D Workspace with x-z Sectioning Plane (y = 0)
    # -----------------------------------------------------------------------
    ax1 = fig.add_subplot(1, 2, 1, projection="3d")
    draw_guide_tube_3d(ax1, length_mm=20.0, radius_mm=2.5)

    # 3D Revolved Volume Surface
    ax1.plot_surface(
        data_both["X_b_grid"], data_both["Y_b_grid"], data_both["Z_b_grid"],
        color="#1baf7a", alpha=0.18, edgecolor="none", shade=True
    )

    # Translucent Sectioning Plane at y = 0
    x_plane = np.linspace(-55, 55, 10)
    z_plane = np.linspace(0, 85, 10)
    X_p, Z_p = np.meshgrid(x_plane, z_plane)
    Y_p = np.zeros_like(X_p)
    ax1.plot_surface(X_p, Y_p, Z_p, color="#2a78d6", alpha=0.10, edgecolor="none")
    ax1.plot([-55, 55, 55, -55, -55], [0, 0, 0, 0, 0], [0, 0, 85, 85, 0],
             color="#2a78d6", lw=1.2, ls="--", alpha=0.65)

    # Plot planar robot backbones in 3D
    for ep in poses_xz:
        p = ep["p"]
        c = ep["color"]
        ott = ep["ott"]
        s_approx = np.linspace(0, ep["itt"], len(p))
        split_idx = int(np.argmin(np.abs(s_approx - ott)))

        # Outer tube (thicker)
        ax1.plot(p[:split_idx+1, 0], p[:split_idx+1, 1], p[:split_idx+1, 2],
                 color=c, lw=2.6, solid_capstyle="round")
        # Inner tube extension (thinner)
        if split_idx < len(p) - 1:
            ax1.plot(p[split_idx:, 0], p[split_idx:, 1], p[split_idx:, 2],
                     color=c, lw=1.5, ls="-")
        # Distal tip marker
        ax1.scatter([p[-1, 0]], [p[-1, 1]], [p[-1, 2]],
                    color=c, s=34, edgecolor="#ffffff", lw=0.9, zorder=6)

    ax1.set_title(r"(a) 3D Workspace with $x$-$z$ Cutting Plane ($y = 0$)", pad=12, fontsize=14.5)
    ax1.set_xlabel(r"$x$ (mm)", labelpad=8, fontsize=13.5)
    ax1.set_ylabel(r"$y$ (mm)", labelpad=8, fontsize=13.5)
    ax1.set_zlabel(r"$z$ (mm)", labelpad=8, fontsize=13.5)
    ax1.tick_params(axis="both", which="major", labelsize=11.5)
    vs.tidy3d(ax1)
    vs.equal_aspect_3d(ax1, data_both["pts_3d"])
    ax1.view_init(elev=22, azim=-55)

    # -----------------------------------------------------------------------
    # Subplot 2: Complete Bilateral 2D x-z Planar Cross-Section
    # -----------------------------------------------------------------------
    ax2 = fig.add_subplot(1, 2, 2)

    rz = data_both["rz_points"]
    b = data_both["boundary_rz"]

    # Construct continuous seamless bilateral boundary without interior seam
    min_z_idx = int(np.argmin(b[:, 1]))
    max_z_idx = int(np.argmax(b[:, 1]))

    outer_indices = list(range(min_z_idx, -1, -1)) + list(range(len(b) - 1, max_z_idx - 1, -1))
    outer_r = b[outer_indices, 0]
    outer_z = b[outer_indices, 1]

    top_x = 0.0
    top_z = float(outer_z[-1])
    bottom_x = 0.0
    bottom_z = float(outer_z[0])

    full_x = np.concatenate([[bottom_x], outer_r, [top_x], -outer_r[::-1], [bottom_x]])
    full_z = np.concatenate([[bottom_z], outer_z, [top_z], outer_z[::-1], [bottom_z]])

    # Fill bilateral reachable area seamlessly
    ax2.fill(full_x, full_z, color="#1baf7a", alpha=0.22, zorder=2,
             label=r"Reachable Slice $\mathcal{D}_{xz}$ (" + f"{2 * data_both['area_2d_mm2']:.0f} " + r"$\mathrm{mm}^2$)")

    # Scatter reachable points in both lobes
    ax2.scatter(rz[:, 0], rz[:, 1], s=2.2, color="#1baf7a", alpha=0.22, zorder=3)
    ax2.scatter(-rz[:, 0], rz[:, 1], s=2.2, color="#1baf7a", alpha=0.22, zorder=3)

    # Seamless bilateral boundary envelope
    ax2.plot(full_x, full_z, color="#0f8a5f", lw=2.2, zorder=4, label=r"Volume Boundary Envelope")

    # Rigid cannula at base
    ax2.plot([0, 0], [0, -25], color="#7d7b75", lw=6.0, solid_capstyle="butt",
             zorder=5, label=r"Rigid Cannula ($s \leq 0$)")
    ax2.axvline(0, color="#999999", ls=":", lw=0.9, alpha=0.6, zorder=1)

    # Draw planar robot backbones in 2D x-z plane
    for ep in poses_xz:
        p = ep["p"]
        c = ep["color"]
        ott = ep["ott"]
        s_approx = np.linspace(0, ep["itt"], len(p))
        split_idx = int(np.argmin(np.abs(s_approx - ott)))

        # Outer tube (thicker)
        ax2.plot(p[:split_idx+1, 0], p[:split_idx+1, 2], color=c, lw=2.8,
                 solid_capstyle="round", zorder=6)
        # Inner tube extension (thinner)
        if split_idx < len(p) - 1:
            ax2.plot(p[split_idx:, 0], p[split_idx:, 2], color=c, lw=1.6,
                     ls="-", zorder=6)

        # Distal tip marker
        ax2.scatter([p[-1, 0]], [p[-1, 2]], color=c, s=52, edgecolor="#ffffff",
                    lw=1.2, zorder=7)

        # Clean numbered circular badge (50% larger)
        dx, dy = ep["badge_offset"]
        ax2.annotate(
            rf"$\mathbf{{{ep['id']}}}$",
            xy=(p[-1, 0], p[-1, 2]), xytext=(dx, dy), textcoords="offset points",
            fontsize=11.0, color="#ffffff",
            bbox=dict(boxstyle="circle,pad=0.22", fc=c, ec="#ffffff", lw=1.2),
            zorder=8
        )

    ax2.set_title(r"(b) Planar Reachable Cross-Section ($x$-$z$ Plane, $y = 0$)", pad=12, fontsize=14.5)
    ax2.set_xlabel(r"$x$ Lateral Reach (mm)", labelpad=6, fontsize=14.0)
    ax2.set_ylabel(r"$z$ Axial Reach (mm)", labelpad=6, fontsize=14.0)
    ax2.set_xlim(-58, 58)
    ax2.set_ylim(-15, 88)
    ax2.set_aspect("equal", adjustable="box")
    ax2.tick_params(axis="both", which="major", labelsize=12.0)
    vs.tidy(ax2)
    ax2.legend(loc="lower left", frameon=False, fontsize=10.5)

    fig.suptitle(
        r"CT-SDR Workspace: Planar Cross-Section in $x$-$z$ Plane ($\theta = 0^\circ / 180^\circ$)",
        x=0.01, ha="left", fontsize=15.5, color=vs.INK
    )
    vs.caption(
        fig,
        rf"Bilateral planar cross-section of the 4-DoF CT-SDR workspace in the $x$-$z$ plane ($y = 0$), "
        rf"spanning both lateral lobes (total cross-sectional area: {2 * data_both['area_2d_mm2']:.0f} mm$^2$). "
        rf"Numbered badges highlight boundary and interior configurations: "
        rf"1 (+x maximum lateral reach, 50.5 mm), 2 (-x symmetric lateral reach), "
        rf"3 (inflected S-curve crossing centerline x = 0 at z = 70 mm to reach x = -6.2 mm), "
        rf"4 (full-stroke S-curve reaching top boundary at z = 80.6 mm), "
        rf"5 (deep right interior reach at x = +20.3 mm, z = 68.1 mm), "
        rf"6 (deep left interior reach at x = -20.3 mm, z = 68.1 mm), "
        rf"and 7 (central-axis deep interior reach at x = +0.6 mm, z = 59.1 mm)."
    )
    fig.tight_layout(rect=(0, 0.06, 0.98, 0.95), w_pad=3.2)
    fig.savefig(out_path, dpi=300)
    plt.close(fig)
    print(f"Wrote {out_path}")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default=CONFIG_PATH, help="Path to ct_sdr.yaml")
    parser.add_argument("--cache", default=DATA_CACHE_PATH, help="Path to workspace_data.npz")
    parser.add_argument("--out-dir", default=FIGURES_DIR, help="Directory to save figures")
    args = parser.parse_args(argv)

    os.makedirs(args.out_dir, exist_ok=True)
    out_png = os.path.join(args.out_dir, "workspace_xz_plane.png")

    print(f"Loading workspace cache from: {args.cache} ...")
    if not os.path.exists(args.cache):
        print(f"Error: Cache {args.cache} not found! Run workspace_analysis.py first.")
        sys.exit(1)

    data_both = np.load(args.cache, allow_pickle=True)

    print("Building robot model and solving planar landmark poses in x-z plane...")
    robot = build_robot_from_yaml(args.config)
    poses_xz = compute_xz_poses(robot)

    print("Generating publication figure: workspace_xz_plane.png ...")
    plot_xz_workspace(data_both, poses_xz, out_png)
    print("Done!")


if __name__ == "__main__":
    main()
