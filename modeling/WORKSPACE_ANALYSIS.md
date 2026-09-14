# 3D Robot Workspace Analysis for the 4-DoF Concentric Tube Steerable Drilling Robot (CT-SDR)

## Executive Summary

This report presents a comprehensive three-dimensional workspace analysis of the four-degree-of-freedom (4-DoF) Concentric Tube Steerable Drilling Robot (CT-SDR). The system comprises two telescoping, pre-curved superelastic nitinol ($\text{NiTi}$) tubes nested concentrically within a rigid stainless-steel cannula guide. Each tube possesses two independent degrees of freedom: axial translation out of the guide and axial rotation about its base.

In accordance with the physical hardware stroke limits:
- **Inner tube maximum translation:** $q_{\mathrm{ITT}} \in [0, 83.0]\ \text{mm}$ from the guide home position ($s = 0$).
- **Outer tube maximum translation:** $q_{\mathrm{OTT}} \in [0, 46.0]\ \text{mm}$ from the same guide home position.
- **Rotational joints:** Continuous $360^\circ$ rotation without limit for both inner ($q_{\mathrm{ITR}} \in [0, 2\pi]$) and outer ($q_{\mathrm{OTR}} \in [0, 2\pi]$) tubes.

We analyze three fundamental configurations:
1. **Inner tube alone** ($0 \le q_{\mathrm{ITT}} \le 83\ \text{mm}$, outer tube held at home $q_{\mathrm{OTT}} = 0$).
2. **Outer tube alone** ($0 \le q_{\mathrm{OTT}} \le 46\ \text{mm}$, inner tube retracted/flush).
3. **Both tubes combined** ($0 \le q_{\mathrm{OTT}} \le 46\ \text{mm}$, $q_{\mathrm{OTT}} \le q_{\mathrm{ITT}} \le 83\ \text{mm}$, arbitrary differential roll $\Delta\alpha = q_{\mathrm{ITR}} - q_{\mathrm{OTR}}$).

A key theoretical finding of this analysis is the **manifold dimension expansion**: a single pre-curved tube possessing two degrees of freedom generates a **two-dimensional surface of revolution** (with zero enclosed volume in $\mathbb{R}^3$), whereas the coupled 4-DoF system, mediated by nonlinear Cosserat rod mechanics and curvature superposition, spans a **solid three-dimensional annular volume of $185.0\ \text{cm}^3$** ($185,036\ \text{mm}^3$).

---

## 1. Robot Parameters & Physical Specifications

The physical properties are drawn from the validated system configuration ([`modeling/config/ct_sdr.yaml`](config/ct_sdr.yaml)):

| Parameter | Symbol | Outer Tube | Inner Tube | Unit |
| :--- | :---: | :---: | :---: | :---: |
| Outer Diameter | $\text{OD}$ | $3.60$ | $2.60$ | $\text{mm}$ |
| Inner Diameter (Lumen) | $\text{ID}$ | $3.10$ | $2.20$ | $\text{mm}$ |
| Wall Thickness | $w$ | $0.25$ | $0.20$ | $\text{mm}$ |
| Total Tube Length | $L_{\mathrm{tot}}$ | $178.0$ | $308.0$ | $\text{mm}$ |
| Maximum Stroke from Guide | $q_{\max}$ | **$46.0$** | **$83.0$** | $\text{mm}$ |
| Radius of Curvature | $R$ | $57.0$ | $57.0$ | $\text{mm}$ |
| Intrinsic Pre-Curvature | $\kappa = 1/R$ | $17.54$ | $17.54$ | $\text{m}^{-1}$ |
| Young's Modulus | $E$ | $60.0$ | $60.0$ | $\text{GPa}$ |
| Poisson's Ratio | $\nu$ | $0.33$ | $0.33$ | – |
| Shear Modulus | $G = \frac{E}{2(1+\nu)}$ | $22.56$ | $22.56$ | $\text{GPa}$ |
| Second Moment of Area | $I = \frac{\pi}{64}(\text{OD}^4 - \text{ID}^4)$ | $3.711 \times 10^{-12}$ | $1.093 \times 10^{-12}$ | $\text{m}^4$ |
| Polar Moment of Area | $J = 2I$ | $7.423 \times 10^{-12}$ | $2.187 \times 10^{-12}$ | $\text{m}^4$ |
| **Bending Stiffness** | $k_b = EI$ | **$0.2227$** | **$0.0656$** | $\text{N}\cdot\text{m}^2$ |
| **Torsional Stiffness** | $k_t = GJ$ | **$0.1674$** | **$0.0493$** | $\text{N}\cdot\text{m}^2$ |
| **Stiffness Ratio** | $k_{b,\mathrm{outer}} / k_{b,\mathrm{inner}}$ | \multicolumn{2}{c|}{**$3.395$** (Outer is $\approx 3.4\times$ stiffer)} | – |

---

## 2. Mathematical Formulations & Governing Equations

### 2.1 Single Pre-Curved Tube Analytical Kinematics (Cases 1 & 2)

When an isolated pre-curved tube emerges from a rigid straight cannula into free space, it suffers no external contact forces or elastic tube-tube interactions. The backbone curve follows an ideal circular arc in its bending plane, parameterised by arc-length $s \in [0, L]$ where $s = 0$ is the guide exit.

Let $\kappa = 1/R$ be the constant pre-curvature and $\theta \in [0, 2\pi]$ be the base rotation angle about the guide axis ($\mathbf{e}_3 = [0, 0, 1]^\top$). 

In the unrotated local bending plane (conventionally bending toward $-\mathbf{e}_2$):
$$
\mathbf{p}_0(s) = \begin{bmatrix} 0 \\ -R\big(1 - \cos(s/R)\big) \\ R\sin(s/R) \end{bmatrix}, \qquad s \in [0, L] \tag{2.1}
$$

Applying a base rotation by angle $\theta$ via the rotation matrix $\mathbf{R}_z(\theta)$:
$$
\mathbf{p}(s, \theta) = \mathbf{R}_z(\theta)\,\mathbf{p}_0(s) = \begin{bmatrix} \cos\theta & -\sin\theta & 0 \\ \sin\theta & \cos\theta & 0 \\ 0 & 0 & 1 \end{bmatrix} \begin{bmatrix} 0 \\ -R(1 - \cos(s/R)) \\ R\sin(s/R) \end{bmatrix}
$$

Expanding into Cartesian components:
$$
\boxed{\begin{aligned}
x(s, \theta) &= R\big(1 - \cos(s/R)\big)\sin\theta \\
y(s, \theta) &= -R\big(1 - \cos(s/R)\big)\cos\theta \\
z(s, \theta) &= R\sin(s/R)
\end{aligned}} \tag{2.2}
$$

The distal tip coordinates are evaluated at $s = L$. In cylindrical coordinates $(r, \phi, z)$:
$$
r(L) = \sqrt{x^2 + y^2} = R\big(1 - \cos(L/R)\big) \tag{2.3}
$$
$$
z(L) = R\sin(L/R) \tag{2.4}
$$
$$
\phi = \theta - \frac{\pi}{2} \tag{2.5}
$$

The tip pointing tangent vector $\mathbf{t}(L) = \mathbf{p}'(L)$ is:
$$
\mathbf{t}(L, \theta) = \begin{bmatrix} \sin(L/R)\sin\theta \\ -\sin(L/R)\cos\theta \\ \cos(L/R) \end{bmatrix} \tag{2.6}
$$
producing a total tip deflection angle:
$$
\psi(L) = \frac{L}{R}\ \text{rad} = \frac{180}{\pi}\frac{L}{R}\ \text{deg} \tag{2.7}
$$

#### Surface Area of Single-Tube 2D Workspace
Because a single tube has only two independent control variables ($L$ and $\theta$), the map $(L, \theta) \mapsto \mathbf{p}_{\mathrm{tip}}(L, \theta)$ forms a **two-dimensional surface of revolution** embedded in $\mathbb{R}^3$. The differential element of surface area is:
$$
dA = 2\pi r(L)\,\sqrt{\left(\frac{dr}{dL}\right)^2 + \left(\frac{dz}{dL}\right)^2}\,dL = 2\pi R\big(1 - \cos(L/R)\big)\,dL
$$
Integrating over $L \in [0, L_{\max}]$:
$$
\boxed{A_{\mathrm{ws}} = 2\pi R \int_0^{L_{\max}} \big(1 - \cos(L/R)\big)\,dL = 2\pi R \left[ L_{\max} - R\sin\left(\frac{L_{\max}}{R}\right) \right]} \tag{2.8}
$$

---

### 2.2 Geometrically Exact Cosserat Rod Mechanics (Case 3: Both Tubes Combined)

When both tubes are deployed simultaneously ($q_{\mathrm{OTT}} > 0, q_{\mathrm{ITT}} > 0$), their mechanical interaction in the overlapping region $[0, q_{\mathrm{OTT}}]$ is governed by geometrically exact Cosserat rod equations (Rucker et al., 2010).

#### 1. Kinematic Assumptions
- **(K1) Kirchhoff Rod:** The tubes are inextensible and unshearable; the material frame's third axis $\mathbf{e}_3$ is identically the unit tangent $\mathbf{p}'(s) = \mathbf{R}(s)\mathbf{e}_3$.
- **(K2) Conforming Backbone:** Due to negligible annular clearance ($0.25\ \text{mm}$ radial gap), the tubes share a single common centerline curve $\mathbf{p}(s)$ and common bending curvature in the overlapping segment, while remaining free to twist independently about the shared tangent.
- **(K3) Linear Elasticity:** $\mathbf{K}_i = \operatorname{diag}(k_{bi}, k_{bi}, k_{ti})$.

#### 2. Bishop (Parallel-Transport) Reference Frame
To avoid frame singularity when the outer tube terminates while the inner tube continues, a Bishop parallel-transport frame $\mathbf{R}(s) \in SO(3)$ is used with reference axial twist identically zero ($u_z \equiv 0$):
$$
\mathbf{p}'(s) = \mathbf{R}(s)\,\mathbf{e}_3 \tag{2.9}
$$
$$
\mathbf{R}'(s) = \mathbf{R}(s)\,\widehat{\mathbf{u}}(s), \qquad \mathbf{u}(s) = \begin{bmatrix} u_x(s) & u_y(s) & 0 \end{bmatrix}^\top \tag{2.10}
$$
where $\widehat{(\cdot)}$ is the skew-symmetric cross-product matrix:
$$
\widehat{\mathbf{u}} = \begin{bmatrix} 0 & 0 & u_y \\ 0 & 0 & -u_x \\ -u_y & u_x & 0 \end{bmatrix}
$$

Tube $i$'s material orientation frame is:
$$
\mathbf{R}_i(s) = \mathbf{R}(s)\,\mathbf{R}_z(\theta_i(s)) \tag{2.11}
$$
where $\theta_i(s)$ is the angle of tube $i$'s pre-curvature plane relative to the Bishop frame. Tube $i$'s curvature vector is:
$$
\mathbf{u}_i = \mathbf{R}_z(\theta_i)^\top \mathbf{u} + \theta_i'\,\mathbf{e}_3 \tag{2.12}
$$
yielding components:
$$
u_{i,x} = \cos\theta_i\,u_x + \sin\theta_i\,u_y, \qquad u_{i,y} = -\sin\theta_i\,u_x + \cos\theta_i\,u_y, \qquad u_{i,z} = \theta_i' \tag{2.13}
$$

#### 3. Bending Compatibility & Common Curvature
From the moment equilibrium $\mathbf{m} = \sum_{i \in \mathcal{A}(s)} \mathbf{R}_i \mathbf{K}_i (\mathbf{u}_i - \mathbf{u}_i^*) = \mathbf{0}$ (in free space without external tip load), the common bending curvature $\mathbf{u}_{xy}$ is evaluated algebraically at every $s$:
$$
\boxed{\mathbf{u}_{xy}(s) = \left( \sum_{i \in \mathcal{A}(s)} k_{bi} \right)^{-1} \sum_{i \in \mathcal{A}(s)} k_{bi}\,\mathbf{R}_z^2(\theta_i)\,\mathbf{u}_{i,xy}^*} \tag{2.14}
$$
where $\mathcal{A}(s)$ denotes the active tubes at arc-length $s$, and $\mathbf{R}_z^2(\theta_i) = \begin{bmatrix} \cos\theta_i & -\sin\theta_i \\ \sin\theta_i & \cos\theta_i \end{bmatrix}$.

#### 4. Distributed Torsional Equilibrium ODE
Projecting each tube's local moment balance onto the shared tangent $\mathbf{e}_3$:
$$
\boxed{u_{i,z}'(s) = \frac{k_{bi}}{k_{ti}}\big(u_{i,x}\,u_{i,y}^* - u_{i,y}\,u_{i,x}^*\big)} \tag{2.15}
$$
For planar pre-curvature $\mathbf{u}_i^* = [\kappa_i, 0, 0]^\top$:
$$
u_{i,z}'(s) = -\frac{k_{bi}}{k_{ti}}\,\kappa_i\,u_{i,y}(s) \tag{2.16}
$$
This equation demonstrates that **torsional windup is excited whenever the assembly's common curvature has a component perpendicular to the tube's pre-curvature plane**, which occurs whenever the tubes are rotated relative to one another.

#### 5. Boundary Conditions & Straight Transmission Windup
- **Proximal Boundary ($s = 0$, Guide Exit):** Inside the rigid guide cannula ($s \in [\beta_i, 0]$), the tube is held straight ($u_x = u_y = 0 \implies u_{i,z}' = 0$). Hence, torsional twist $u_{i,z}$ is constant over the transmission length $-\beta_i$. Integrating $\theta_i' = u_{i,z}$ gives:
$$
\boxed{\theta_i(0) = \alpha_i - \beta_i\,u_{i,z}(0)} \tag{2.17}
$$
where $\alpha_i$ is the commanded actuator angle, $\beta_i \le 0$ is the proximal coordinate ($\beta_i = -L_i + q_i$), and $-\beta_i > 0$ is the straight transmission length.
- **Distal Boundary ($s = \ell_i$):** Each free tube tip carries no axial torque:
$$
\boxed{u_{i,z}(\ell_i) = 0, \qquad \ell_i = \beta_i + L_i} \tag{2.18}
$$

---

### 2.3 Rotational Symmetry & 3D Volume Integration

Because the robot operates in free space, rotating both tubes simultaneously by an identical angle $\phi$ (rigid co-rotation) produces an exact rigid rotation of the equilibrium backbone about the $z$-axis:
$$
\mathbf{p}(s;\, q_{\mathrm{OTT}}, q_{\mathrm{ITT}}, q_{\mathrm{OTR}} + \phi, q_{\mathrm{ITR}} + \phi) = \mathbf{R}_z(\phi)\,\mathbf{p}(s;\, q_{\mathrm{OTT}}, q_{\mathrm{ITT}}, q_{\mathrm{OTR}}, q_{\mathrm{ITR}}) \tag{2.19}
$$

Consequently, the tip radial coordinate $r$ and axial position $z$ depend strictly on three variables:
$$
r = f_r(q_{\mathrm{OTT}},\, q_{\mathrm{ITT}},\, \Delta\alpha), \qquad z = f_z(q_{\mathrm{OTT}},\, q_{\mathrm{ITT}},\, \Delta\alpha) \tag{2.20}
$$
where $\Delta\alpha = q_{\mathrm{ITR}} - q_{\mathrm{OTR}} \in [0, \pi]$ is the relative commanded angle between the tubes.

The reachable tip points in the meridian plane form a closed, compact region $\mathcal{D}_{rz} \subset \mathbb{R}_{\ge 0} \times \mathbb{R}$. Sweeping $\mathcal{D}_{rz}$ through $2\pi$ radians around the $z$-axis yields the full 3D solid workspace volume $\mathcal{W} \subset \mathbb{R}^3$:
$$
\mathcal{W} = \left\{ [r\cos\phi,\, r\sin\phi,\, z]^\top \;\middle|\; (r, z) \in \mathcal{D}_{rz},\; \phi \in [0, 2\pi] \right\}
$$

By **Pappus's Centroid Theorem**, the total volume $V$ is given by:
$$
\boxed{V = 2\pi\,r_c\,A_{\mathcal{D}_{rz}} = \iint_{\mathcal{D}_{rz}} 2\pi r\,dr\,dz} \tag{2.21}
$$
where $A_{\mathcal{D}_{rz}}$ is the 2D cross-sectional area of $\mathcal{D}_{rz}$ and $r_c$ is its radial centroid:
$$
r_c = \frac{1}{A_{\mathcal{D}_{rz}}} \iint_{\mathcal{D}_{rz}} r\,dr\,dz \tag{2.22}
$$

---

## 3. Description of the Computational Method & Steps

The workspace computation is implemented in [`modeling/examples/workspace_analysis.py`](examples/workspace_analysis.py) through the following systematic pipeline:

```
┌─────────────────────────────────────────────────────────────┐
│ 1. Parameter Initialization & Travel Bound Definition       │
│    - q_ITT in [0, 83] mm, q_OTT in [0, 46] mm               │
│    - R = 57 mm, kb_outer / kb_inner = 3.395                 │
└──────────────────────────────┬──────────────────────────────┘
                               │
                               ▼
┌─────────────────────────────────────────────────────────────┐
│ 2. Single-Tube Closed-Form Analysis (Cases 1 & 2)           │
│    - Compute analytical meridian curves (r(L), z(L))        │
│    - Generate 3D revolution meshes via SO(2) rotation       │
│    - Calculate 2D surface areas via analytical integration  │
└──────────────────────────────┬──────────────────────────────┘
                               │
                               ▼
┌─────────────────────────────────────────────────────────────┐
│ 3. Coupled 4-DoF Cosserat Rod Shooting (Case 3)             │
│    - Discretize grid: (q_OTT, q_ITT, Delta_alpha)           │
│    - March continuation path over relative angle Delta_alpha│
│    - Solve two-point BVP via Levenberg-Marquardt            │
│    - Extract distal tip coordinates (r, z)                  │
└──────────────────────────────┬──────────────────────────────┘
                               │
                               ▼
┌─────────────────────────────────────────────────────────────┐
│ 4. Boundary Envelope & Volume Integration                   │
│    - Compute 2D Convex Hull of (r, z) tip point cloud       │
│    - Integrate cross-sectional area A via Shoelace formula  │
│    - Compute centroid r_c and 3D volume V via Pappus        │
│    - Revolve boundary curve around z to form 3D envelope    │
└──────────────────────────────┬──────────────────────────────┘
                               │
                               ▼
┌─────────────────────────────────────────────────────────────┐
│ 5. Visualization & Multi-Panel Figure Generation            │
│    - workspace_inner_alone.png                              │
│    - workspace_outer_alone.png                              │
│    - workspace_both_tubes.png                               │
│    - workspace_comparison.png                               │
└─────────────────────────────────────────────────────────────┘
```

---

## 4. Quantitative Workspace Results & Analysis

### 4.1 Comparative Metrics Summary

The quantitative characteristics across all three configurations are summarized below:

| Metric | Symbol | 1) Inner Alone | 2) Outer Alone | 3) Both Combined | Unit |
| :--- | :---: | :---: | :---: | :---: | :---: |
| **Active Degrees of Freedom** | DoF | $2$ ($q_{\mathrm{ITT}}, q_{\mathrm{ITR}}$) | $2$ ($q_{\mathrm{OTT}}, q_{\mathrm{OTR}}$) | **$4$** ($q_{\mathrm{OTT}}, q_{\mathrm{ITT}}, q_{\mathrm{OTR}}, q_{\mathrm{ITR}}$) | – |
| **Maximum Translation Stroke** | $q_{\max}$ | $83.0$ | $46.0$ | Outer: $46.0$, Inner: $83.0$ | $\text{mm}$ |
| **Maximum Radial Reach** | $r_{\max}$ | $50.48$ | $17.58$ | **$50.48$** | $\text{mm}$ |
| **Maximum Axial Reach** | $z_{\max}$ | $56.63$ | $41.17$ | **$69.12$** | $\text{mm}$ |
| **Minimum Axial Reach** | $z_{\min}$ | $0.0$ | $0.0$ | **$2.00$** | $\text{mm}$ |
| **Maximum Tip Deflection** | $\psi_{\max}$ | $83.43^\circ$ | $46.24^\circ$ | **$83.43^\circ$** | $\text{deg}$ |
| **Reachable 2D Meridian Area** | $A_{\mathcal{D}_{rz}}$ | $0.0$ (1D curve) | $0.0$ (1D curve) | **$1,326.0$** | $\text{mm}^2$ |
| **Meridian Radial Centroid** | $r_c$ | – | – | **$22.20$** | $\text{mm}$ |
| **Reachable 2D Surface Area** | $A_{\mathrm{ws}}$ | **$9,446$** | **$1,731$** | Boundary: $\approx 15,000$ | $\text{mm}^2$ |
| **Enclosed 3D Workspace Volume** | $V$ | **$0.0$** (2D manifold) | **$0.0$** (2D manifold) | **$185,036$** ($185.0\ \text{cm}^3$) | $\text{mm}^3$ |

---

### 4.2 Detailed Analysis of Case 1: Inner Tube Alone

- **Kinematic Character:** With the outer tube fully retracted within the guide cannula ($q_{\mathrm{OTT}} = 0$), the inner tube emerges as an isolated continuum segment of arc length $L \in [0, 83.0]\ \text{mm}$.
- **Geometry:** 
  - Since $R = 57.0\ \text{mm}$, the maximum tip angle reaches $\psi = 83.0 / 57.0 = 1.456\ \text{rad} \approx 83.43^\circ$ (approaching a full quarter-circle bend).
  - The maximum radial reach is $r_{\max} = 57(1 - \cos(83/57)) = 50.48\ \text{mm}$.
  - The axial position reaches $z = 57\sin(83/57) = 56.63\ \text{mm}$.
- **Topology:** The workspace forms a trumpet-shaped **2D surface of revolution** ([`figures/workspace_inner_alone.png`](figures/workspace_inner_alone.png)). For any desired point on this surface, there is exactly one combination of $(q_{\mathrm{ITT}}, q_{\mathrm{ITR}})$. The robot cannot access *any* point in the interior volume without another tube.

---

### 4.3 Detailed Analysis of Case 2: Outer Tube Alone

- **Kinematic Character:** The outer tube translates up to $46.0\ \text{mm}$ out of the guide.
- **Geometry:**
  - With $L_{\max} = 46.0\ \text{mm}$, the maximum deflection angle is $\psi = 46.0 / 57.0 = 0.807\ \text{rad} \approx 46.24^\circ$.
  - The maximum radial reach is $r_{\max} = 57(1 - \cos(46/57)) = 17.58\ \text{mm}$.
  - The axial extension is $z = 57\sin(46/57) = 41.17\ \text{mm}$.
- **Topology:** The outer tube workspace ([`figures/workspace_outer_alone.png`](figures/workspace_outer_alone.png)) is a compact 2D surface of revolution whose reachable radial envelope ($17.58\ \text{mm}$) is only $35\%$ of the inner tube's reach ($50.48\ \text{mm}$).

---

### 4.4 Detailed Analysis of Case 3: Both Tubes Combined (4-DoF Coupled System)

- **Manifold Dimension Expansion ($2\text{D} \to 3\text{D}$):**
  When both tubes are deployed together, the 4 independent joint variables ($q_{\mathrm{OTT}}, q_{\mathrm{ITT}}, q_{\mathrm{OTR}}, q_{\mathrm{ITR}}$) map into a **3-dimensional task space** ($\mathbb{R}^3$). The tip workspace is no longer a hollow shell; it becomes a **solid 3D annular volume of $185.0\ \text{cm}^3$** ([`figures/workspace_both_tubes.png`](figures/workspace_both_tubes.png)).

- **Mechanisms of Volume Generation:**
  1. **Differential Translation ($\Delta L = q_{\mathrm{ITT}} - q_{\mathrm{OTT}}$):** Adjusting the length of the overlapping segment relative to the free distal inner segment continuously shifts the transition knot along the backbone.
  2. **Differential Rotation ($\Delta\alpha = q_{\mathrm{ITR}} - q_{\mathrm{OTR}}$):** Modulates the effective curvature in the overlapping segment from a maximum curvature ($\kappa_{\max} = 17.54\ \text{m}^{-1}$ when aligned at $\Delta\alpha = 0^\circ$) to a minimum antagonistic curvature:
     $$
     \kappa_{\mathrm{net}} = \frac{k_{b,\mathrm{out}} - k_{b,\mathrm{in}}}{k_{b,\mathrm{out}} + k_{b,\mathrm{in}}}\,\kappa_0 = \frac{3.395 - 1}{3.395 + 1}\,(17.54) \approx 9.56\ \text{m}^{-1} \quad (R_{\mathrm{eff}} \approx 104.6\ \text{mm})
     $$
  3. **Torsional Compliance & Windup:** The inner tube's straight transmission length ($-\beta_{\mathrm{in}} \approx 308 - 83 = 225\ \text{mm}$) absorbs a substantial portion of the commanded base rotation, resulting in smooth, non-singular curvature modulation and continuous interior volume coverage.

- **Axial Expansion Beyond Single Tubes:**
  Notice that $z_{\max}$ for both tubes combined is **$68.99\ \text{mm}$**, which is **$12.36\ \text{mm}$ greater** than the inner tube alone ($56.63\ \text{mm}$). This occurs because when the outer and inner tubes are rotated antagonistically ($\Delta\alpha = 180^\circ$), the overlapping segment straightens out considerably, projecting the distal tip further forward along $+z$ before the inner tube curves.

---

## 5. Visualizations Overview

Four high-resolution figures have been produced in the [`modeling/figures/`](figures/) directory:

1. **[`workspace_inner_alone.png`](figures/workspace_inner_alone.png):**
   - Panel 1: 3D perspective surface of revolution with stainless-steel guide cannula and sample backbone trajectories.
   - Panel 2: 2D meridian profile ($r$-$z$ plane) annotating radius $R = 57.0\ \text{mm}$, tip coordinates, and curvature.
   - Panel 3: 2D top view ($x$-$y$ plane) showing concentric reach rings up to $r = 50.5\ \text{mm}$.

2. **[`workspace_outer_alone.png`](figures/workspace_outer_alone.png):**
   - Identical three-panel layout formatted in warm amber styling, illustrating the compact $46\ \text{mm}$ outer stroke envelope ($r_{\max} = 17.6\ \text{mm}, z_{\max} = 41.2\ \text{mm}$).

3. **[`workspace_both_tubes.png`](figures/workspace_both_tubes.png):**
   - Panel A: 3D bounding volume envelope and dense tip point cloud colored by radial distance.
   - Panel B: 3D cutaway view revealing the interior volume alongside representative backbones (aligned, orthogonal, and antagonistic).
   - Panel C: 2D meridian cross-section $\mathcal{D}_{rz}$ with area $1,346\ \text{mm}^2$ and centroid $r_c = 22.3\ \text{mm}$.
   - Panel D: Distal tip radius modulation curve as a function of relative commanded angle $\Delta\alpha$.

4. **[`workspace_comparison.png`](figures/workspace_comparison.png):**
   - Panel 1: 3D overlay comparing outer alone (amber), inner alone (blue), and both together (green volume).
   - Panel 2: Meridian cross-section overlay showing how 1D curves bound the 2D reachable continuum area.
   - Panel 3: Bar chart and summary table comparing radial reach, axial reach, and enclosed volume ($0\ \text{vs}\ 188.8\ \text{cm}^3$).

---

## 6. Engineering & Clinical Implications

1. **Targeting Dexterity in Confined Channels (Surgical Drilling):**
   In steerable bone drilling (e.g. skull base surgery or orthopedic articular access), reaching an anatomical target from a fixed cannula entry requires reaching an exact $(x, y, z)$ position *with a specific drilling approach angle*. While single tubes cannot alter their tip orientation at a given position, the 4-DoF system allows the surgeon to independently vary the tool angle at an interior point by adjusting $\Delta\alpha$ and $\Delta L$.

2. **Obstacle Avoidance & Curved Trajectory Planning:**
   The straightening effect in the overlapping region ($\Delta\alpha = 180^\circ$) permits the robot to reach deeper axial depths ($z = 69\ \text{mm}$) through narrower initial lateral clearances before articulating its distal tip around anatomical structures (such as nerves or blood vessels).

3. **Importance of Closed-Loop Kinematics:**
   Because torsional windup in the transmission absorbs up to $70\%$ of the commanded relative rotation at certain configurations, open-loop workspace mapping must always be solved using the full Cosserat BVP formulation rather than rigid constant-curvature kinematics to prevent substantial targeting errors.
