# Cosserat-rod model of the CT-SDR

A geometrically exact mechanics model of the two-nitinol-tube concentric-tube
steerable drilling robot this repo controls, implementing

> D. C. Rucker, B. A. Jones and R. J. Webster III, **"A Geometrically Exact
> Model for Externally Loaded Concentric-Tube Continuum Robots,"** *IEEE
> Transactions on Robotics*, vol. 26, no. 5, pp. 769–780, 2010.
> [doi:10.1109/TRO.2010.2062570](https://doi.org/10.1109/TRO.2010.2062570) ·
> [PubMed 21566688](https://pubmed.ncbi.nlm.nih.gov/21566688/)

Given the four joint values the ROS1 stack already speaks — `OTT`, `ITT` (mm)
and `OTR`, `ITR` (deg) — it returns the full backbone shape `p(s)`, each tube's
roll and twist along it, and the internal force and moment, with or without an
external load at the tip.

This is deliberately **not** a constant-curvature kinematic model. Two effects
it captures that constant curvature cannot are both large on this robot:

* **Torsional windup.** The inner tube has ~270 mm of torsionally free length
  behind the guide when 35 mm is deployed. Commanding `ITR = 90°` turns the tip
  by about **14°** — the rest is wound into the tubes.
* **Snap-through.** Past a critical rotation the wound equilibrium branch folds
  away and the robot jumps discontinuously to another one. The model finds it;
  `figures/workspace_itr_sweep.png` shows where.

## Layout

```
modeling/
├── ctr/                       # the model
│   ├── tube.py                # geometry, material, pre-curvature, EI/GJ, strain & stress
│   ├── loads.py               # external tip wrench and distributed loads
│   ├── model.py               # the Cosserat ODEs + the shooting solver  <- the paper
│   └── robot.py               # CT-SDR: joints -> tube placement -> solve
├── config/ct_sdr.yaml         # the nitinol tube parameters (edit this, not the code)
├── examples/
│   ├── vizstyle.py            # shared figure style
│   ├── plot_shape.py          # one configuration: shape, windup, material utilisation
│   ├── plot_workspace.py      # ITR sweep: tip loci, windup and snap-through
│   └── plot_drilling_load.py  # tip deviation under a drilling reaction wrench
├── figures/                   # generated PNGs
└── tests/test_ctr.py          # 67 verification tests
```

## Quick start

Only `numpy`, `scipy`, `pyyaml` and `matplotlib` are needed — no ROS, no
hardware.

```python
from ctr import CTSDR, Joints

robot = CTSDR.from_yaml()
sol = robot.solve(Joints(ott=35, itt=35, otr=0, itr=180))

sol.tip_position     # (3,)   metres, in the guide frame
sol.tip_tangent      # (3,)   drilling direction
sol.p, sol.s         # (M,3) backbone and its arc lengths
sol.theta, sol.u_z   # (M,N) per-tube roll and twist; nan where a tube is absent
sol.n, sol.m         # (M,3) internal force and moment
```

```bash
cd modeling
python3 examples/plot_shape.py --ott 35 --itt 35 --itr 90
python3 examples/plot_workspace.py
python3 examples/plot_drilling_load.py
python3 tests/test_ctr.py
```

## The robot as modelled

Two pre-curved nitinol tubes inside a rigid stainless-steel guide, from the
design description:

| | Outer tube | Inner tube |
|---|---|---|
| Outer diameter | 3.6 mm | 2.6 mm |
| Wall thickness | 0.25 mm | 0.2 mm |
| Inner diameter (derived) | 3.1 mm | 2.2 mm |
| Total length | ≈ 178 mm | ≈ 308 mm |
| Radius of curvature | 50 mm | 50 mm |
| `EI` at E = 60 GPa | 0.2227 N·m² | 0.0656 N·m² |
| `GJ` at ν = 0.33 | 0.1674 N·m² | 0.0493 N·m² |

The outer tube is **3.4× stiffer in bending**, so where the two overlap the
combined shape sits much closer to the outer tube's pre-curvature than to a
simple average.

`s = 0` is the distal face of the guide, +z along the guide axis. With all
joints at zero the tube tips are flush with it, so `OTT` and `ITT` read out
deployed length directly.

### Two parameters you should check before trusting absolute numbers

Everything above except these two comes straight from the design description.

1. **`curved_length_mm` (78.54 mm, a placeholder).** How much of each tube's
   distal end is pre-curved was not specified. 78.54 mm is a 90° arc at
   R = 50 mm, which covers the largest advance in the ICRA2027 tests (35 mm
   outer, 70 mm inner) with margin. **Measure it and set it in
   `config/ct_sdr.yaml`.** It does not affect any configuration whose
   deployment stays inside the arc — `CTSDR.check()` raises if you leave it.
2. **`youngs_modulus_gpa` (60 GPa).** See the next section — this is an
   *effective* modulus to be fitted, not a handbook value.

### Nitinol is not linear, and the model is

The 50 mm pre-curvature alone puts the outer tube's surface at **3.6 % strain**
(`kappa * r = 20 / m × 1.8 mm`). Superelastic NiTi leaves its linear austenitic
branch around 1 % and spends the rest on the stress-induced martensite plateau,
where stress is nearly flat in strain. So the linear-elastic constitutive law
`K = diag(EI, EI, GJ)` — the paper's, and every standard concentric-tube
model's — is being applied outside the regime it describes.

This is the normal state of affairs for concentric-tube robots, and it is
workable, but it changes what `E` means: **treat it as an effective modulus
fitted to measured shapes at your operating curvature**, not as a material
constant. The 60 GPa default is the value most common in the literature and is
a starting point for that fit, nothing more. Values from 40 to 83 GPa all
appear for binary NiTi.

`CTSDR.material_limits(sol)` reports peak bending strain, torque and torsional
shear stress per tube, and flags configurations whose predicted shear passes a
conservative 250 MPa — past which the numbers are extrapolation:

```
material utilisation (linear-elastic model):
  outer  peak bending strain  3.58 %   torque +0.1936 N.m   peak shear   47.0 MPa
  inner  peak bending strain  2.59 %   torque -0.1936 N.m   peak shear  115.1 MPa
```

(The two torques being exactly equal and opposite is not a coincidence — it is
the axial moment balance, and the test suite checks it.)

## The model

### Assumptions

* Each tube is a **Kirchhoff rod**: inextensible, no transverse shear, so the
  material frame's z-axis is the backbone tangent and only bending and torsion
  store energy.
* The tubes are concentric with **negligible clearance**, so at every arc
  length the tubes present share one backbone curve and therefore one bending
  curvature. They remain free to twist relative to one another. (The real
  clearance is 0.25 mm per side; that is the assumption most likely to cost
  you accuracy at low deployment.)
* Tube-to-tube contact is **frictionless and transmits no axial torque**, so
  each tube's moment balance closes about its own z-axis independently.
* The guide is **rigid and straight**. Inside it the tubes cannot bend but are
  free to twist.

### State and equations

The reference frame `R(s)` is chosen torsion-free (a Bishop frame) rather than
tied to one particular tube — bookkeeping only, but it is what lets tubes start
and end anywhere along `s` without the reference tube disappearing when the
inner tube is advanced past the outer one. Tube `i`'s material frame is
`R_i = R Rz(theta_i)`.

Bending compatibility makes the common curvature algebraic, from the
constitutive law `m = sum_i R_i K_i (u_i - u_i*)`:

```
u_xy = ( sum_i k_bi )^-1 [ (R^T m)_xy + sum_i k_bi Rz(theta_i) u*_i,xy ]
```

and the integrated states are

```
p'       =  R e3
R'       =  R hat(u)
theta_i' =  u_i,z
u_i,z'   =  (k_bi / k_ti) (u_i,x u*_i,y  -  u_i,y u*_i,x)   -   (R e3)·l / k_ti
n'       = -f
m'       = -p' x n  -  l
```

with `u_i,xy = Rz(theta_i)^T u_xy`, `k_bi = E_i I_i`, `k_ti = G_i J_i`.
Unloaded, `n` and `m` are identically zero and the solver drops them.

### Boundary conditions

Shooting unknowns are the base twists `u_i,z(0)`, plus `n(0)` and `m(0)` when
loaded. The conditions are:

* **`u_i,z(l_i) = 0`** at each tube's distal end — a free end carries no
  torsional moment.
* **`n(L) = F_tip`, `m(L) = M_tip`** at the distal-most tip.
* **`theta_i(0) = alpha_i - beta_i u_i,z(0)`** — windup in the straight
  transmission between actuator and guide exit. Inside the guide the tubes
  cannot bend, so the torsion ODE's right-hand side vanishes and `u_i,z` is
  constant there. **On this robot this term dominates**: at `ITT = 35 mm` the
  inner tube has 273 mm of transmission against 35 mm of deployed length, and
  roughly 70° of a 90° `ITR` command is absorbed before the robot even starts.

One exception to the first condition. Torque **about the tangent** is the only
load the tubes do not share — they are free to twist relative to each other, so
it stays in whichever tube it was applied to. The tube carrying an applied
drill torque therefore ends at `k_ti u_i,z(l_i) = tau` instead of zero, while
every other tube's end stays free. `ExternalLoad.axial_torque_tube` selects it;
the default is the tube reaching the tip, which here is the inner tube — the
one the drill bit is on.

### Multiple equilibria

A wound concentric-tube robot has several torsional equilibria, and branches
are created and destroyed at fold points as the tubes rotate. Two consequences:

* **Use continuation for sweeps.** `CTSDR.solve_path()` marches each
  configuration from its neighbour's converged twists, keeping the solver on
  the branch the hardware actually follows instead of letting it jump. Passing
  `guess=sol.u_z0` does the same for a single step.
* **A fold is not a solver failure.** When the current branch folds away,
  `solve()` retries from a spread of physically scaled starts and reports how
  many it needed in `Solution.restarts`. Across 5 deployments × 73 rotations,
  3 of 365 configurations needed a restart and none failed.

## Verification

`tests/test_ctr.py` checks the model against cases with closed-form answers
rather than against itself:

| What | Check |
|---|---|
| Single pre-curved tube | traces an exact circular arc; every sample lies on the circle to 1 nm |
| Identical tubes, aligned | reduce to the single-tube arc |
| Identical tubes, opposed | cancel to a straight line, zero curvature |
| Dissimilar tubes, opposed | give the stiffness-weighted curvature `(k1 κ1 − k2 κ2)/(k1 + k2)` |
| Straight tube, small tip load | matches Euler–Bernoulli `δ = FL³/3EI` to 2 % |
| Straight tube, pure end moment | exact constant curvature `M/EI` |
| Distributed load | integrates to its resultant at the base |
| Axial moment balance | `Σ GJ_i u_i,z = (R e3)·m` pointwise — an exact identity |
| `GJ → ∞` | recovers the torsionally rigid constant-curvature model |
| Inextensibility | `|p'| = 1`; chord deficit only the expected O(h²) |
| Snap-through | full 360° sweeps solve at every deployment |

```
$ python3 tests/test_ctr.py
Ran 67 tests in 64s
OK
```

## Known limitations

* **Linear elasticity on a superelastic material** — see above. The single
  biggest source of absolute error; fit `E` before quoting numbers.
* **Zero clearance.** The 0.25 mm per-side gap between the tubes is ignored.
  Real tubes can lag each other slightly, which matters most at short overlap.
* **Frictionless.** No tube-to-tube or tube-to-guide friction, and no
  hysteresis. Real windup shows both, so measured rotation will lag the model's
  and will not retrace on the way back.
* **Statics only.** No dynamics, no drilling process model — the drilling
  reaction is an input wrench you supply, not something the model predicts.
* **The guide is perfectly rigid and straight**, and the transition at its exit
  is treated as an ideal clamp.

## Comparing against the NDI measurements

`plot_shape.py --csv out.csv` writes the modelled backbone in the same
millimetre convention `ICRA2027/ndi_tip_analysis.py` uses, so a modelled shape
can be laid over a measured one. The model's origin is the guide exit with +z
along the guide axis; the NDI traces are referenced to each trial's first
sample, so a rigid registration between the two frames is still needed before
the comparison means anything.

The most informative comparison to run first is the windup panel of
`plot_workspace.py` against the measured rotation in `exp_info.md`'s Exp 5–13
(360° commanded rotations): if the measured tip rotation is much closer to the
commanded value than the model says, the transmission is stiffer than modelled
— most likely because the tubes are supported over more of their retracted
length than the free-twist assumption allows.
