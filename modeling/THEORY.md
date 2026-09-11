# Mathematics of the CT-SDR Cosserat model

The equations solved by [`modeling/ctr/`](ctr/), their derivations, and the
numerical method. Implements

> D. C. Rucker, B. A. Jones and R. J. Webster III, "A Geometrically Exact Model
> for Externally Loaded Concentric-Tube Continuum Robots," *IEEE Transactions
> on Robotics* **26**(5):769–780, 2010.
> [doi:10.1109/TRO.2010.2062570](https://doi.org/10.1109/TRO.2010.2062570)

with two additions of our own: a straight lead-in in the pre-curvature profile
(§3.2) and Coulomb friction (§8).

Every equation below is cross-referenced to the code that implements it.
Numerical results and the comparison against hardware are in
[VALIDATION_REPORT.md](VALIDATION_REPORT.md); this file is the maths only.

---

## 1. Notation

| Symbol | Meaning | Units |
| :-- | :-- | :-- |
| $N$ | number of tubes; $i=0$ outermost, $i=N-1$ innermost | – |
| $s$ | arc length along the shared backbone, $s=0$ at the guide exit | m |
| $L$ | arc length of the distal-most tube tip | m |
| $\beta_i$ | arc-length coordinate of tube $i$'s proximal end, $\beta_i\le 0$ | m |
| $\ell_i=\beta_i+L_i$ | arc length of tube $i$'s distal tip | m |
| $\mathbf p(s)$ | backbone position | m |
| $\mathbf R(s)\in SO(3)$ | reference frame | – |
| $\theta_i(s)$ | roll of tube $i$ about the backbone, in the reference frame | rad |
| $\mathbf u(s)$ | curvature of the reference frame, in that frame | m⁻¹ |
| $\mathbf u_i(s)$ | curvature of tube $i$, in tube $i$'s own frame | m⁻¹ |
| $\mathbf u_i^*(\sigma)$ | pre-curvature of tube $i$ at its material coordinate $\sigma$ | m⁻¹ |
| $\mathbf n(s),\ \mathbf m(s)$ | internal force / moment of the assembly, global frame | N, N·m |
| $k_{bi}=E_iI_i$ | bending stiffness | N·m² |
| $k_{ti}=G_iJ_i$ | torsional stiffness | N·m² |
| $\alpha_i$ | commanded actuator roll of tube $i$ | rad |
| $\mathbf f,\ \mathbf l$ | distributed external force / moment per unit length | N/m, N·m/m |

$\widehat{(\cdot)}:\mathbb R^3\to\mathfrak{so}(3)$ is the hat map,
$\widehat{\mathbf a}\mathbf b=\mathbf a\times\mathbf b$;
$\mathbf R_z(\theta)$ is a rotation about $\mathbf e_3$; and
$\mathbf R_z^{2}(\theta)=\begin{bmatrix}\cos\theta&-\sin\theta\\ \sin\theta&\cos\theta\end{bmatrix}$
is its $2\times2$ block.

---

## 2. Kinematic assumptions

**(K1) Kirchhoff rod.** Each tube is inextensible and unshearable, so the
material frame's third axis is the unit tangent:

$$\mathbf p'(s)=\mathbf R(s)\,\mathbf e_3,\qquad \|\mathbf p'\|=1 .$$

Only bending and torsion store energy.

**(K2) Conforming backbone.** Tubes are concentric with negligible clearance,
so at every $s$ all tubes present share **one** curve $\mathbf p(s)$ and hence
one bending curvature. They remain free to twist relative to one another about
the common tangent.

**(K3) Frictionless axial contact.** Tube-to-tube contact transmits no torque
about the tangent, so each tube's moment balance closes about its own axis
independently. §8 relaxes this.

**(K4) Linear elasticity.** $\mathbf K_i=\mathrm{diag}(k_{bi},k_{bi},k_{ti})$,
with $J_i = 2I_i$ for a circular annulus.

### 2.1 Choice of reference frame

Rucker et al. tie $\mathbf R$ to tube 1's material frame. We instead take a
**Bishop (parallel-transport) frame**, defined by

$$u_z\equiv 0 .$$

This is a change of bookkeeping only — $\mathbf R$ still satisfies
$\mathbf p'=\mathbf R\mathbf e_3$ — but it removes the requirement that a
particular tube exist everywhere. On this robot the inner tube is advanced past
the outer one, so a frame tied to the outer tube would vanish part way along
$s$. Every tube is then treated symmetrically, at the cost of $\mathbf R$ no
longer being any tube's material frame.

Tube $i$'s material frame is

$$\mathbf R_i(s)=\mathbf R(s)\,\mathbf R_z(\theta_i(s)). \tag{2.1}$$

### 2.2 Curvature of an individual tube

Differentiating (2.1) with $\mathbf R'=\mathbf R\widehat{\mathbf u}$:

$$\widehat{\mathbf u_i}=\mathbf R_i^\mathsf T\mathbf R_i'
=\mathbf R_z^\mathsf T\widehat{\mathbf u}\mathbf R_z+\mathbf R_z^\mathsf T\mathbf R_z'
=\widehat{\mathbf R_z^\mathsf T\mathbf u}+\theta_i'\widehat{\mathbf e_3},$$

using $\mathbf R^\mathsf T\widehat{\mathbf a}\mathbf R=\widehat{\mathbf R^\mathsf T\mathbf a}$
and $\mathbf R_z^\mathsf T\mathbf R_z'=\theta_i'\widehat{\mathbf e_3}$. Hence

$$\boxed{\ \mathbf u_i=\mathbf R_z(\theta_i)^\mathsf T\mathbf u+\theta_i'\mathbf e_3\ } \tag{2.2}$$

and, componentwise, with $u_z=0$:

$$u_{i,x}=\ \ \cos\theta_i\,u_x+\sin\theta_i\,u_y,\qquad
u_{i,y}=-\sin\theta_i\,u_x+\cos\theta_i\,u_y,\qquad
u_{i,z}=\theta_i'. \tag{2.3}$$

The last of these is what makes the Bishop choice convenient: **$\theta_i'=u_{i,z}$
exactly**, with no reference-tube twist to subtract.

---

## 3. Constitutive law and pre-curvature

The moment carried by tube $i$, in the global frame, is

$$\mathbf m_i=\mathbf R_i\,\mathbf K_i\,(\mathbf u_i-\mathbf u_i^*),
\qquad \mathbf m=\sum_i \mathbf m_i . \tag{3.1}$$

### 3.1 Pre-curvature convention

$\mathbf u_i^*$ is expressed in tube $i$'s own frame and taken about its local
$x$-axis:

$$\mathbf u_i^*(\sigma)=\begin{cases}(\kappa_i,\,0,\,0)^\mathsf T & \sigma\in[\sigma_i^{\text{start}},\ \sigma_i^{\text{end}}]\\[2pt] \mathbf 0 & \text{otherwise,}\end{cases}
\qquad \sigma=s-\beta_i . \tag{3.2}$$

With $\mathbf u=(\kappa,0,0)$ and $\mathbf R(0)=\mathbf I$,
$\mathbf p''(0)=\widehat{\mathbf u}\mathbf e_3=\mathbf u\times\mathbf e_3=(0,-\kappa,0)$,
so a lone tube curves toward $-\mathbf e_2$ with centre of curvature at
$(0,-1/\kappa,0)$. Rolling the tube rotates this bending plane; that rotation
enters through $\theta_i$, never through $\mathbf u_i^*$.

### 3.2 Where the curve sits along the tube

Rucker et al. put the curve at the distal tip. We allow a straight **lead-in**
of length $S_i$ between the curve and the tip, so with total length $L_i$ and
curved length $C_i$:

$$\sigma_i^{\text{end}}=L_i-S_i,\qquad \sigma_i^{\text{start}}=L_i-S_i-C_i . \tag{3.3}$$

$S_i=0$ recovers the standard distally-curved tube. This matters because the
*emerged* part of a tube is its distal part: with $S_i>0$ the curvature sits
**proximal**, near the guide exit, with the straight run beyond it. That is the
opposite of the usual arrangement, and it is what the NDI advances appear to
show (VALIDATION_REPORT.md §6). The shipped config uses $S_i=0$.

*Code: [`ctr/tube.py`](ctr/tube.py) — `precurvature`, `curve_end`,
`straight_length`, `deployable_length`.*

---

## 4. Bending compatibility gives $\mathbf u$ algebraically

Rotate (3.1) into the reference frame and substitute (2.2):

$$\mathbf R^\mathsf T\mathbf m
=\sum_i \mathbf R_z(\theta_i)\mathbf K_i(\mathbf u_i-\mathbf u_i^*)
=\sum_i \mathbf R_z\mathbf K_i\mathbf R_z^\mathsf T\mathbf u
+\mathbf e_3\sum_i k_{ti}\theta_i'
-\sum_i \mathbf R_z\mathbf K_i\mathbf u_i^* .$$

Because the $xy$ block of $\mathbf K_i$ is isotropic ($k_{bi}\mathbf I_2$), it
commutes with $\mathbf R_z$:

$$\mathbf R_z(\theta)\,\mathbf K_i\,\mathbf R_z(\theta)^\mathsf T=\mathbf K_i .$$

Taking the $xy$ components (the $\mathbf e_3$ term drops out) and solving:

$$\boxed{\ \mathbf u_{xy}=\Big(\textstyle\sum_{i\in\mathcal A(s)} k_{bi}\Big)^{-1}
\Big[(\mathbf R^\mathsf T\mathbf m)_{xy}
+\sum_{i\in\mathcal A(s)} k_{bi}\,\mathbf R_z^{2}(\theta_i)\,\mathbf u^*_{i,xy}\Big],
\qquad u_z=0\ } \tag{4.1}$$

where $\mathcal A(s)$ is the set of tubes present at $s$. **$\mathbf u$ is not a
state** — it is evaluated pointwise from $\mathbf R$, $\mathbf m$ and the
$\theta_i$. Unloaded ($\mathbf m\equiv\mathbf 0$) this is the familiar
stiffness-weighted superposition of pre-curvatures.

*Code: [`ctr/model.py`](ctr/model.py) `CosseratModel._deriv`, the `acc` accumulator.*

---

## 5. Torsion: one ODE per tube

This is the step that makes the model more than constant-curvature kinematics.

Each tube separately obeys

$$\mathbf n_i'+\mathbf f_i=\mathbf 0,\qquad
\mathbf m_i'+\mathbf p'\times\mathbf n_i+\mathbf l_i=\mathbf 0, \tag{5.1}$$

where $\mathbf f_i,\mathbf l_i$ are the loads on tube $i$ *including* contact
from its neighbours. Write $\mathbf m_i=\mathbf R_i\mathbf M_i$ with
$\mathbf M_i=\mathbf K_i(\mathbf u_i-\mathbf u_i^*)$, so that

$$\mathbf m_i'=\mathbf R_i\big(\mathbf u_i\times\mathbf M_i+\mathbf M_i'\big).$$

Rotate (5.1) into tube $i$'s frame:

$$\mathbf u_i\times\mathbf M_i+\mathbf M_i'
+\mathbf R_i^\mathsf T(\mathbf p'\times\mathbf n_i)+\mathbf R_i^\mathsf T\mathbf l_i=\mathbf 0 .$$

**Key simplification.** All material frames share the tangent:
$\mathbf p'=\mathbf R\mathbf e_3=\mathbf R_i\mathbf R_z^\mathsf T\mathbf e_3=\mathbf R_i\mathbf e_3$,
so $\mathbf R_i^\mathsf T\mathbf p'=\mathbf e_3$ and

$$\mathbf R_i^\mathsf T(\mathbf p'\times\mathbf n_i)=\mathbf e_3\times(\mathbf R_i^\mathsf T\mathbf n_i),$$

which has **no $\mathbf e_3$ component**. Taking the $\mathbf e_3$ component of
the whole equation therefore eliminates $\mathbf n_i$ — the shear force, which
we never have to compute:

$$\underbrace{\mathbf e_3\cdot(\mathbf u_i\times\mathbf M_i)}_{k_{bi}(u_{i,y}u^*_{i,x}-u_{i,x}u^*_{i,y})}
+\underbrace{\mathbf e_3\cdot\mathbf M_i'}_{k_{ti}u_{i,z}'}
+\underbrace{\mathbf e_3\cdot(\mathbf R_i^\mathsf T\mathbf l_i)}_{\mathbf p'\cdot\mathbf l_i}=0 .$$

Hence

$$\boxed{\ u_{i,z}'=\frac{k_{bi}}{k_{ti}}\big(u_{i,x}u^*_{i,y}-u_{i,y}u^*_{i,x}\big)
-\frac{\mathbf p'\cdot\mathbf l_i}{k_{ti}}\ } \tag{5.2}$$

Under (K3) the contact part of $\mathbf l_i$ is normal to $\mathbf p'$, so it
drops out and only *external* axial torque survives. §8 puts friction back.

For a planar tube ($\mathbf u_i^*=(\kappa_i,0,0)$) the driving term is
$-\,(k_{bi}/k_{ti})\,\kappa_i\,u_{i,y}$: a tube twists only when the common
curvature has a component **perpendicular** to its own pre-curvature, i.e. only
when it is rotated relative to its neighbours.

---

## 6. Assembly equilibrium

Summing (5.1) over tubes, internal contact forces cancel pairwise and

$$\mathbf n'=-\mathbf f,\qquad \mathbf m'=-\mathbf p'\times\mathbf n-\mathbf l . \tag{6.1}$$

Unloaded with a free tip, $\mathbf n\equiv\mathbf m\equiv\mathbf 0$ and (4.1)
reduces to pure superposition; the solver then drops these six states.

---

## 7. The complete system

State vector, $\dim = 18+2N$:

$$\mathbf y=\big(\underbrace{\mathbf p}_{3},\ \underbrace{\mathrm{vec}\,\mathbf R}_{9},\ \underbrace{\theta_{0..N-1}}_{N},\ \underbrace{u_{0,z..N-1,z}}_{N},\ \underbrace{\mathbf n}_{3},\ \underbrace{\mathbf m}_{3}\big).$$

$$
\begin{aligned}
\mathbf p' &= \mathbf R\,\mathbf e_3\\
\mathbf R' &= \mathbf R\,\widehat{\mathbf u}, &&\mathbf u \text{ from (4.1)}\\
\theta_i' &= u_{i,z}\\
u_{i,z}' &= \tfrac{k_{bi}}{k_{ti}}\big(u_{i,x}u^*_{i,y}-u_{i,y}u^*_{i,x}\big)-\tfrac{\mathbf p'\cdot\mathbf l_i}{k_{ti}}\\
\mathbf n' &= -\mathbf f\\
\mathbf m' &= -\mathbf p'\times\mathbf n-\mathbf l
\end{aligned}
\tag{7.1}
$$

with $u_{i,x},u_{i,y}$ from (2.3). Tubes not present at $s$ have
$\theta_i'=u_{i,z}'=0$ and contribute nothing to (4.1).

---

## 8. Friction

Relaxing (K3). At a contact between adjacent tubes $a$ (outer) and $b$ (inner)
with normal force per unit length $w$, contact radius $r_b$ and coefficient
$\mu$, the friction torque per unit length is $\tau=\mu\,w\,r_b$, applied
equal and opposite:

$$u_{b,z}'\ \mathrel{-}=\ \frac{\varsigma\,\tau}{k_{tb}},\qquad
u_{a,z}'\ \mathrel{+}=\ \frac{\varsigma\,\tau}{k_{ta}}, \tag{8.1}$$

where $\varsigma=\mathrm{sgn}\,\dfrac{\mathrm d(\theta_b-\theta_a)}{\mathrm dt}$
is the sense of *relative angular velocity* as the manoeuvre proceeds — **not**
$\mathrm d u_z/\mathrm ds$, which is a spatial derivative and says nothing about
sliding. Because (8.1) is antisymmetric, the identity (11.1) below is preserved
exactly: friction moves torque between tubes, it never creates any.

**Normal force.** A rod held at curvature $\mathbf u$ against a natural
curvature $\mathbf u^*$ carries $\Delta M=k_b\|\mathbf u-\mathbf u^*\|$; holding
it there over a characteristic length $1/\kappa$ requires a transverse load of
order

$$w \simeq \kappa^2\,k_{b}\,\|\mathbf u_{xy}-\mathbf u^*_{xy}\|,\qquad \kappa=\|\mathbf u_{xy}\| . \tag{8.2}$$

This is an estimate good to about a factor of two, exposed as `force_scale` so
it can be overridden or calibrated. It correctly vanishes where a tube sits at
its natural curvature or where the assembly is straight.

**Guide contact.** A uniformly pre-curved tube held *straight* carries a
constant moment, so $\mathbf m_i'=\mathbf 0$ and it needs **no distributed
transverse load** — the reaction is concentrated at the guide's lip. (8.2)
therefore gives zero inside the guide, and the transmission's normal force is
supplied separately as $w_g$, estimated from the guide reacting the
straightening moment $k_{bi}\kappa_i$ as a couple over a grip length $\lambda$:

$$w_g\simeq \frac{k_{bi}\,\kappa_i}{\lambda^{2}} . \tag{8.3}$$

**Coefficients.** Dry nitinol–nitinol $\mu=0.35$ (adhesive, galls readily);
dry nitinol–stainless $\mu=0.25$; wetted, 0.15 / 0.12.

*Code: [`ctr/friction.py`](ctr/friction.py), and `CosseratModel._apply_friction`.*

---

## 9. Boundary conditions

### 9.1 Distal

A free tube end carries no torsional moment. If an external torque $\tau_i$
about the tangent is applied to tube $i$ at its tip (a drill bit),

$$k_{ti}\,u_{i,z}(\ell_i)=\tau_i . \tag{9.1}$$

Torque about the tangent is the **one** load the tubes do not share — they are
free to twist relative to each other, so it stays in the tube it was applied to
and every *other* tube still ends at $u_{i,z}=0$. At the distal-most tip,

$$\mathbf n(L)=\mathbf F_{\text{tip}},\qquad \mathbf m(L)=\mathbf M_{\text{tip}}. \tag{9.2}$$

### 9.2 Proximal: the transmission

Inside the rigid guide the tubes cannot bend, so $u_{i,x}=u_{i,y}=0$ and (5.2)
gives $u_{i,z}'=0$: **twist is constant along the transmission**. Integrating
$\theta_i'=u_{i,z}$ from $\beta_i$ to $0$,

$$\boxed{\ \theta_i(0)=\alpha_i-\beta_i\,u_{i,z}(0)\ } \tag{9.3}$$

($\beta_i\le 0$, so this *adds* windup.) On this robot the term dominates: at
ITT = 35 mm the inner tube has ~273 mm of torsionally free transmission against
35 mm of deployed length.

The actuator angle carries a **datum offset**,

$$\alpha_i=\varsigma_i\,(\text{commanded angle})+\phi_i , \tag{9.3a}$$

because the rotation joints have no absolute reference: homing zeroes whatever
position the tube currently sits at, so $\phi_i$ is whatever the pre-curvature
plane happened to be at homing. Only the difference
$\psi_0=\phi_{\text{in}}-\phi_{\text{out}}$ is physical — a common offset
rigidly rolls the whole robot about $\mathbf e_3$ — and it can differ between
sessions. Note that (9.3) then *absorbs* most of $\psi_0$: a clamp offset is
wound into the transmission just as a commanded rotation is, so the advances are
nearly blind to it and it has to be identified from the rotation segments.

With guide friction the torque is largest at the actuator and bleeds off
distally, $|T_i(s)|=|T_i(0)|+\tau_i^{f}|s|$, so integrating $T_i/k_{ti}$ adds a
quadratic term:

$$\theta_i(0)=\alpha_i-\beta_i u_{i,z}(0)
+\mathrm{sgn}\big(u_{i,z}(0)\big)\,\frac{\tau_i^{f}\,\beta_i^{2}}{2\,k_{ti}} . \tag{9.4}$$

*Code: `CosseratModel._pack_y0`, `_transmission_friction_twist`.*

---

## 10. Solution method

### 10.1 Shooting

(7.1) is an initial value problem, but the conditions are split between $s=0$
and $s=\ell_i,L$. We shoot. Let $\mathcal E=\{i:\ell_i>0\}$ be the emerged
tubes, $k=|\mathcal E|$. Unknowns:

$$\mathbf x=\big(\underbrace{\mathbf n(0),\ \mathbf m(0)}_{6,\ \text{loaded only}},\ \{u_{i,z}(0)\}_{i\in\mathcal E}\big)\in\mathbb R^{6+k}\ \text{ or }\ \mathbb R^{k}.$$

Initial state: $\mathbf p(0)=\mathbf 0$, $\mathbf R(0)=\mathbf I$,
$\theta_i(0)$ from (9.3)/(9.4). Integrate to $L$, then form

$$\mathbf F(\mathbf x)=\begin{pmatrix}
\big\{u_{i,z}(\ell_i)-\tau_i/k_{ti}\big\}_{i\in\mathcal E}\\
\mathbf n(L)-\mathbf F_{\text{tip}}\\
\mathbf m(L)-\mathbf M_{\text{tip}}
\end{pmatrix}=\mathbf 0 . \tag{10.1}$$

Square system, solved by Levenberg–Marquardt with a forward-difference
Jacobian. The residual mixes units — m⁻¹ in the torsion block, N and N·m in the
wrench blocks — and convergence is judged on the plain 2-norm of the whole
vector, so the threshold has to be read against the largest-magnitude block.

### 10.2 Segmentation

The right-hand side is only *piecewise* smooth: $\mathcal A(s)$ changes where a
tube starts or ends, and $\mathbf u_i^*$ jumps at each straight↔curved
transition. Breakpoints are placed at

$$\{0,\ L\}\ \cup\ \bigcup_i\{\beta_i,\ \beta_i+\sigma_i^{\text{start}},\ \beta_i+\sigma_i^{\text{end}},\ \ell_i\}$$

clipped to $(0,L)$ and de-duplicated at $10^{-9}$ m. Each segment is integrated
separately with a constant active set, and $\mathbf R$ is re-projected onto
$SO(3)$ at every boundary by

$$\mathbf R\leftarrow \mathbf U\mathbf V^\mathsf T,\qquad \mathbf R=\mathbf U\boldsymbol\Sigma\mathbf V^\mathsf T,$$

since RK45 does not preserve orthonormality.

### 10.3 Numerical parameters

| Parameter | Value | Why |
| :-- | :-- | :-- |
| integrator | RK45 | non-stiff; the system is a smooth ODE within a segment |
| `rtol`, `atol` | $10^{-9}$, $10^{-11}$ | must sit below the FD step (next row) |
| FD step `diff_step` | $10^{-6}$ (relative) | **must exceed the adaptive integrator's step-to-step jitter** ($\sim$`rtol`), or the Jacobian is dominated by quadrature noise. Three decades of margin. |
| LM tolerances | $10^{-14}$ | run to stagnation; convergence judged on $\|\mathbf F\|$ |
| residual threshold | $10^{-6}$ | $\sim$µrad of twist over these lengths |

---

## 11. Invariants

Two identities hold exactly and are asserted in the test suite.

**Axial moment balance.** Since $\mathbf R_z$ preserves the $\mathbf e_3$
component,

$$\mathbf p'\cdot\mathbf m=\mathbf e_3\cdot(\mathbf R^\mathsf T\mathbf m)
=\sum_i \mathbf e_3\cdot\mathbf K_i(\mathbf u_i-\mathbf u_i^*)
=\sum_i k_{ti}\,u_{i,z}. \tag{11.1}$$

Unloaded, $\sum_i k_{ti}u_{i,z}=0$ everywhere: with two tubes they carry equal
and opposite torque. Friction (8.1) preserves this by construction.

**Inextensibility.** $\|\mathbf p'\|=\|\mathbf R\mathbf e_3\|=1$ exactly.

---

## 12. Multiple equilibria, folds and continuation

The BVP does **not** have a unique solution, and this drives the whole solver
design.

For two unloaded tubes with $\mathbf u_i^*=(\kappa_i,0,0)$, write
$\psi=\theta_2-\theta_1$ for the relative roll. From (4.1),

$$u_{1,y}=\frac{k_{b2}\kappa_2\sin\psi}{k_{b1}+k_{b2}},\qquad
u_{2,y}=-\frac{k_{b1}\kappa_1\sin\psi}{k_{b1}+k_{b2}},$$

and substituting into (5.2) with $u^*_{i,y}=0$ gives $u_{i,z}'=-(k_{bi}/k_{ti})\kappa_i u_{i,y}$, so

$$\boxed{\ \psi''=u_{2,z}'-u_{1,z}'=C\sin\psi,\qquad
C=\frac{k_{b1}k_{b2}\,\kappa_1\kappa_2}{k_{b1}+k_{b2}}\Big(\frac1{k_{t1}}+\frac1{k_{t2}}\Big)>0\ } \tag{12.1}$$

An **inverted** pendulum: the aligned state $\psi=0$ minimises elastic energy
yet is a *hyperbolic* fixed point of the ODE, with solutions growing like
$e^{\sqrt C s}$. That is not a paradox — it is the generic Euler–Lagrange
structure of an energy minimiser — but it has three consequences.

For this robot $C=409\ \mathrm{m^{-2}}$, so $\sqrt C=20.2\ \mathrm{m^{-1}}$ and
the e-folding length is **49 mm** — comparable to the deployed lengths of
17.5–70 mm. The instability is therefore mild over one deployment and severe
over several, which is exactly the regime where shooting needs help:

1. **Shooting is ill-conditioned** for large $\sqrt C\,L$: small changes in
   $u_{i,z}(0)$ grow exponentially toward the tip.
2. **Multiple roots.** (12.1) with the BVP conditions admits several solutions
   once the tubes are wound up.
3. **Folds (snap-through).** As $\alpha$ varies, branches are created and
   destroyed at limit points. At 10 mm deployment the wound branch of this
   robot dies between ITR 190° and 200°, leaving exactly one root — the
   physical snap.

**Continuation.** Sweeps march the solution, seeding each solve with the
previous $u_{i,z}(0)$ (`CTSDR.solve_path`). This keeps the solver on the branch
the hardware follows rather than letting it jump.

**Restarts.** When a branch folds away, LM stalls in its basin. `solve` then
retries from a deterministic spread of starts scaled by

$$u_{i,z}^{(0)}\sim \frac{\pi}{|\beta_i|}\times\{0.15,0.35,0.6,1.0\}\times\{\pm1\}^{k},$$

i.e. the twist that would wind half a turn into tube $i$'s transmission — a
problem-scaled guess rather than an arbitrary constant. Sign patterns come
first because the tubes' torques balance and therefore usually have opposite
signs. Over 5 deployments × 73 rotations, 3 of 365 configurations needed a
restart and none failed.

> ⚠️ The same non-uniqueness bites friction: cold-starting a friction solve can
> land on a different branch (220° of delivered roll for a 90° command). Always
> continue from the frictionless solution.

---

## 13. Closed-form cases used for verification

These are the analytic answers the test suite checks against, not model output.

**Single pre-curved tube, deployed $d$.**

$$\mathbf p(d)=\Big(0,\ \tfrac{\cos\kappa d-1}{\kappa},\ \tfrac{\sin\kappa d}{\kappa}\Big),\qquad
\mathbf p'(d)=(0,-\sin\kappa d,\cos\kappa d). \tag{13.1}$$

The tip traces the arc of radius $1/\kappa$, so **the traced path radius equals
the tube's radius** — the identity used to read curvature off the NDI advances.

**With a lead-in $S$** (§3.2): an arc of length $d-S$ followed by a straight
tangent of length $S$.

**$n$ tubes, torsionally rigid, aligned deployment.** Constant curvature

$$\kappa_{\text{eff}}=\frac{\big\|\sum_i k_{bi}\kappa_i e^{\mathrm i\theta_i}\big\|}{\sum_i k_{bi}} . \tag{13.2}$$

Two identical tubes at $\psi=\pi$ cancel exactly; at $\psi=0$ they reproduce the
single-tube arc.

**Straight tube, small tip load.** Euler–Bernoulli cantilever
$\delta=FL^3/3EI$; and a pure end moment gives exact constant curvature $M/EI$.

**Material utilisation.** Peak surface bending strain and torsional shear:

$$\varepsilon=\kappa\,r_o,\qquad \tau=G\,u_z\,r_o . \tag{13.3}$$

At $\kappa=1/57$ mm⁻¹ and $r_o=1.8$ mm this is 3.2 % — several times the ~1 %
where binary NiTi leaves its linear austenitic branch, which is why $E$ here is
an *effective* modulus to be fitted, not a handbook value.

---

## 14. Post-processing used for validation

**Rigid registration (Kabsch).** Model and tracker frames are unrelated, so
shapes are compared after minimising $\|\mathbf A\mathbf R^\mathsf T+\mathbf t-\mathbf B\|_F$
over $\mathbf R\in SO(3),\mathbf t$. With centred $\mathbf A_c,\mathbf B_c$ and
$\mathbf H=\mathbf A_c^\mathsf T\mathbf B_c=\mathbf U\boldsymbol\Sigma\mathbf V^\mathsf T$,

$$\mathbf R=\mathbf V\,\mathrm{diag}\big(1,1,\det(\mathbf V\mathbf U^\mathsf T)\big)\,\mathbf U^\mathsf T,
\qquad \mathbf t=\bar{\mathbf B}-\mathbf R\bar{\mathbf A}. \tag{14.1}$$

The $\det$ correction forces a proper rotation. **This is why a handedness
error cannot be absorbed**: a mirrored path has no proper rotation onto its
original, which is how the inverted rotation sense was detected.

**Circle fit (Kåsa).** Project onto the SVD best-fit plane, then solve the
*linear* least squares

$$\begin{bmatrix}2u_j & 2v_j & 1\end{bmatrix}\begin{pmatrix}c_u\\c_v\\d\end{pmatrix}=u_j^2+v_j^2,
\qquad r=\sqrt{d+c_u^2+c_v^2}, \tag{14.2}$$

and unwrap $\operatorname{atan2}(v-c_v,\,u-c_u)$ for the swept angle.

**Turn profile.** $\Theta(\sigma)=\sigma/r(\sigma)$ from a circle fit to the
sub-arc $[0,\sigma]$. Linear in curvature and well-conditioned, unlike the
radius itself — which is why geometry is identified on $\Theta$ rather than on
the registration residual, whose sensitivity a rigid transform largely absorbs.

---

## 15. Code map

| Equation | Implementation |
| :-- | :-- |
| (2.2), (2.3) | `CosseratModel._deriv` |
| (3.2), (3.3) | `Tube.precurvature`, `Tube.curve_end` |
| (4.1) | `_deriv`, `acc` accumulator |
| (5.2), (7.1) | `_deriv` |
| (8.1)–(8.3) | `ctr/friction.py`, `_apply_friction` |
| (9.1), (9.2) | `solve`, `residual` |
| (9.3), (9.4) | `_pack_y0`, `_transmission_friction_twist` |
| (10.1) | `solve`, `least_squares(method="lm")` |
| §10.2 | `CosseratModel.segments`, `_integrate`, `orthonormalize` |
| §12 restarts | `_restart_points` |
| (11.1) | asserted in `TestAxialMomentBalance` |
| (13.1)–(13.3) | `tests/test_ctr.py`, `Tube.bending_strain`/`shear_stress` |
| (14.1), (14.2) | `validation/ndi_experiments.py` — `kabsch`, `arc_metrics` |

---

## 16. Assumptions, collected

Ordered by how much they are believed to cost on this robot.

1. **Linear elasticity** on a material at 3.2 % strain, well onto its
   superelastic plateau. $E$ is an effective modulus. Note that torsional
   windup is *invariant* to a uniform scaling of $E$ — both the driving torque
   and the transmission compliance scale linearly — so this does not explain
   windup discrepancies; it does set absolute moments and stresses.
2. **Zero clearance** (K2). The real gap is 0.25 mm per side, which lets the
   tubes lag each other; the shared-backbone assumption cannot represent it.
3. **Frictionless by default** (K3), optionally relaxed by §8 with a single
   sliding sense per commanded step — right on a monotonic leg, wrong where the
   local sliding reverses. No stick, no friction cone, no axial friction.
4. **Rigid, perfectly straight guide**, with an ideal clamp at its exit.
5. **Statics only.** No dynamics; a drilling reaction is an input wrench, not
   something the model predicts.
