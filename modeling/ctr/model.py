"""
model.py -- geometrically exact Cosserat-rod model of a concentric-tube robot.

Implements the formulation of

    D. C. Rucker, B. A. Jones and R. J. Webster III, "A Geometrically Exact
    Model for Externally Loaded Concentric-Tube Continuum Robots," IEEE
    Transactions on Robotics, vol. 26, no. 5, pp. 769-780, 2010.

Model assumptions (the paper's, restated for this robot)
--------------------------------------------------------
* Each tube is a Kirchhoff rod: inextensible, no transverse shear, so the
  material frame's z-axis is the backbone tangent.  Only bending and torsion
  store energy.
* The tubes are concentric with negligible clearance, so at every arc length
  the tubes present share **one** backbone curve ``p(s)`` and therefore one
  bending curvature.  They may, however, twist arbitrarily relative to one
  another about the common tangent.
* Tube-to-tube interaction is frictionless and transmits no axial torque, so
  each tube's own moment balance closes about its z-axis independently.
* Linear elasticity with ``K_i = diag(E_i I_i, E_i I_i, G_i J_i)``.

State and equations
-------------------
The reference frame ``R(s)`` is chosen torsion-free (a Bishop frame,
``u_ref,z = 0``) rather than tied to one particular tube.  That is a change of
bookkeeping only, and it is what lets tubes start and end anywhere along ``s``
without the reference tube disappearing.  Each tube's material frame is
``R_i = R R_z(theta_i)``, so ``theta_i`` is tube ``i``'s roll about the
backbone measured in the reference frame, and ``theta_i' = u_i,z``.

With ``k_bi = E_i I_i``, ``k_ti = G_i J_i`` and ``A(s)`` the set of tubes
present at ``s``, the common bending curvature follows algebraically from the
constitutive law ``m = sum_i R_i K_i (u_i - u_i*)``::

    u_ref,xy = (sum_{i in A} k_bi)^-1 [ (R^T m)_xy
                                        + sum_{i in A} k_bi Rz(theta_i) u*_i,xy ]

and the integrated states obey::

    p'       =  R e3
    R'       =  R hat(u_ref)
    theta_i' =  u_i,z
    u_i,z'   =  (k_bi / k_ti) (u_i,x u*_i,y - u_i,y u*_i,x)  -  (R e3) . l / k_ti
    n'       = -f
    m'       = -p' x n  -  l

where ``u_i,xy = Rz(theta_i)^T u_ref,xy``.  (The distributed-moment term drops
out of ``theta_i`` because ``R_i e3 = R e3`` for every tube -- all material
frames share the tangent.)

Boundary conditions and shooting
--------------------------------
Unknowns are the base torsions ``u_i,z(0)`` and, when the robot is loaded,
``n(0)`` and ``m(0)``.  Residuals are:

* ``u_i,z(l_i) = 0`` at each tube's own distal end -- a free tube end carries
  no torsional moment.  The one exception is the tube an external torque about
  the tangent is applied to (the drill bit's tube), whose end instead satisfies
  ``k_ti u_i,z(l_i) = tau``.  Torque about the tangent is the one load the
  tubes do *not* share: they are free to twist relative to each other, so it
  stays in the tube it was applied to.
* ``n(L) = F_tip`` and ``m(L) = M_tip`` at the distal-most tip.

The proximal roll boundary condition accounts for windup in the straight
transmission between the actuator and the robot base.  Inside the rigid
stainless-steel guide the tubes cannot bend, so the torsion ODE's right-hand
side vanishes and ``u_i,z`` is constant there, giving

    theta_i(0) = alpha_i - beta_i u_i,z(0)

for a tube whose proximal end sits at ``beta_i <= 0`` and whose actuator is
commanded to ``alpha_i``.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Optional, Sequence, Tuple

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares

from .loads import ExternalLoad
from .tube import Tube

E3 = np.array([0.0, 0.0, 1.0])

#: Arc lengths closer than this are treated as the same breakpoint [m].
KNOT_TOL = 1e-9


def hat(u: Sequence[float]) -> np.ndarray:
    """Skew-symmetric matrix with ``hat(u) v == cross(u, v)``."""
    return np.array(
        [
            [0.0, -u[2], u[1]],
            [u[2], 0.0, -u[0]],
            [-u[1], u[0], 0.0],
        ]
    )


def rot_z(theta: float) -> np.ndarray:
    """Rotation by ``theta`` about z."""
    c, s = np.cos(theta), np.sin(theta)
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])


def orthonormalize(R: np.ndarray) -> np.ndarray:
    """Nearest rotation matrix to ``R`` in the Frobenius sense."""
    u, _, vt = np.linalg.svd(R)
    Rn = u @ vt
    if np.linalg.det(Rn) < 0.0:  # pragma: no cover - integration never gets here
        u[:, -1] *= -1.0
        Rn = u @ vt
    return Rn


class ShootingError(RuntimeError):
    """Raised when the two-point boundary value problem does not converge."""


@dataclass
class Solution:
    """Solved shape of the robot.

    Attributes
    ----------
    s
        Arc-length samples [m], from the robot base ``s = 0`` (the distal face
        of the rigid guide) to the distal-most tube tip.
    p
        ``(M, 3)`` backbone positions [m] in the global frame.
    R
        ``(M, 3, 3)`` reference-frame orientations (Bishop frame, see module
        docstring).  ``R[k][:, 2]`` is the backbone tangent.
    theta
        ``(M, N)`` roll of each tube about the backbone [rad].  Entries for a
        tube that is not present at that ``s`` are ``nan``.
    u_z
        ``(M, N)`` torsional curvature of each tube [1/m], ``nan`` where absent.
    u
        ``(M, 3)`` common curvature vector in the reference frame [1/m].
    n, m
        ``(M, 3)`` internal force [N] and moment [N.m] in the global frame.
    active
        Per-sample tuple of tube indices present at that arc length.
    u_z0
        Converged base torsions, indexed like ``tubes``; ``nan`` for tubes that
        have not emerged.  Pass back as ``guess=`` to warm-start a nearby
        configuration.
    residual
        Final shooting residual vector; ``residual_norm`` is its 2-norm.
    """

    s: np.ndarray
    p: np.ndarray
    R: np.ndarray
    theta: np.ndarray
    u_z: np.ndarray
    u: np.ndarray
    n: np.ndarray
    m: np.ndarray
    active: List[Tuple[int, ...]]
    u_z0: np.ndarray
    residual: np.ndarray
    n_iter: int = 0
    #: Extra shooting starts needed; > 0 means the first branch folded away.
    restarts: int = 0
    betas: np.ndarray = field(default_factory=lambda: np.zeros(0))
    alphas: np.ndarray = field(default_factory=lambda: np.zeros(0))

    @property
    def tip_position(self) -> np.ndarray:
        return self.p[-1]

    @property
    def tip_rotation(self) -> np.ndarray:
        return self.R[-1]

    @property
    def tip_tangent(self) -> np.ndarray:
        return self.R[-1][:, 2]

    @property
    def arc_length(self) -> float:
        """Total deployed length of the backbone [m]."""
        return float(self.s[-1] - self.s[0])

    @property
    def residual_norm(self) -> float:
        return float(np.linalg.norm(self.residual))

    def tube_frame(self, k: int, i: int) -> np.ndarray:
        """Material frame ``R_i`` of tube ``i`` at sample ``k``."""
        return self.R[k] @ rot_z(self.theta[k, i])

    def position_at(self, s_query: float) -> np.ndarray:
        """Backbone position at arc length ``s_query`` [m], linearly interpolated.

        Needed to track a tube tip that is not the distal-most point: while the
        outer tube is advanced ahead of the inner one, the inner tube's tip --
        which is where the drill bit and the tracker marker sit -- lies part way
        along the backbone, not at its end.  Queries outside the deployed range
        clamp to the ends.
        """
        s = float(np.clip(s_query, self.s[0], self.s[-1]))
        return np.array([np.interp(s, self.s, self.p[:, k]) for k in range(3)])


class CosseratModel:
    """Geometrically exact model for a set of concentric tubes.

    Parameters
    ----------
    tubes
        Ordered outermost-to-innermost (the order only affects indexing).

    Notes
    -----
    Tube placement is given per solve as ``betas``: the arc-length coordinate
    of each tube's proximal end, which is zero or negative because that end
    sits inside the rigid guide.  Advancing a tube by ``d`` raises its
    ``beta`` by ``d``.  :class:`ctr.robot.CTSDR` builds these from the robot's
    OTT/ITT/OTR/ITR joint values.
    """

    def __init__(self, tubes: Sequence[Tube]):
        self.tubes: List[Tube] = list(tubes)
        if not self.tubes:
            raise ValueError("need at least one tube")
        self.n_tubes = len(self.tubes)
        self.k_b = np.array([t.bending_stiffness for t in self.tubes])
        self.k_t = np.array([t.torsional_stiffness for t in self.tubes])

        N = self.n_tubes
        self._i_p = slice(0, 3)
        self._i_R = slice(3, 12)
        self._i_th = slice(12, 12 + N)
        self._i_uz = slice(12 + N, 12 + 2 * N)
        self._i_n = slice(12 + 2 * N, 15 + 2 * N)
        self._i_m = slice(15 + 2 * N, 18 + 2 * N)
        self.state_dim = 18 + 2 * N

    # -- segmentation -------------------------------------------------------

    def segments(
        self, betas: Sequence[float]
    ) -> Tuple[List[Tuple[float, float, Tuple[int, ...]]], float]:
        """Split ``[0, L]`` where the tube set or the pre-curvature changes.

        Returns ``(segments, L)`` with each segment ``(s_start, s_end,
        active_tube_indices)``.  Breakpoints are placed at every tube start,
        every tube end and both straight-to-curved transitions of every tube,
        because the right-hand side is only piecewise smooth across those.
        """
        betas = np.asarray(betas, dtype=float)
        if betas.shape != (self.n_tubes,):
            raise ValueError(f"betas must have length {self.n_tubes}")

        ends = betas + np.array([t.length for t in self.tubes])
        L = float(np.max(ends))
        if L <= KNOT_TOL:
            raise ValueError(
                "no tube has emerged from the guide (all tube tips are at or "
                "behind s = 0); advance at least one translational joint"
            )
        emerged = ends > KNOT_TOL
        if np.any(betas[emerged] > KNOT_TOL):
            raise ValueError(
                "a tube's proximal end lies beyond the robot base (beta > 0), "
                "which the rigid-guide boundary condition does not model"
            )

        knots = [0.0, L]
        for i, t in enumerate(self.tubes):
            for x in (betas[i], betas[i] + t.straight_length,
                      betas[i] + t.curve_end, ends[i]):
                if KNOT_TOL < x < L - KNOT_TOL:
                    knots.append(float(x))
        knots.sort()

        merged = [knots[0]]
        for x in knots[1:]:
            if x - merged[-1] > KNOT_TOL:
                merged.append(x)

        segs = []
        for a, b in zip(merged[:-1], merged[1:]):
            mid = 0.5 * (a + b)
            act = tuple(
                i
                for i in range(self.n_tubes)
                if betas[i] - KNOT_TOL < mid < ends[i] + KNOT_TOL
            )
            if not act:  # pragma: no cover - impossible for beta <= 0
                raise ValueError(f"no tube covers arc length {mid:.6g} m")
            segs.append((a, b, act))
        return segs, L

    # -- right-hand side ----------------------------------------------------

    def _deriv(
        self,
        s: float,
        y: np.ndarray,
        active: Tuple[int, ...],
        betas: np.ndarray,
        load: ExternalLoad,
        want_l: bool,
        torque_tube: int,
    ) -> np.ndarray:
        R = y[self._i_R].reshape(3, 3)
        theta = y[self._i_th]
        u_z = y[self._i_uz]
        n_vec = y[self._i_n]
        m_vec = y[self._i_m]

        # Common bending curvature from the constitutive law (see docstring).
        acc = (R.T @ m_vec)[:2].copy()
        k_b_sum = 0.0
        star = {}
        for i in active:
            ux, uy = self.tubes[i].precurvature(s - betas[i])
            star[i] = (ux, uy)
            c, sn = np.cos(theta[i]), np.sin(theta[i])
            kb = self.k_b[i]
            acc[0] += kb * (c * ux - sn * uy)
            acc[1] += kb * (sn * ux + c * uy)
            k_b_sum += kb
        u = np.array([acc[0] / k_b_sum, acc[1] / k_b_sum, 0.0])

        dy = np.zeros_like(y)
        tangent = R @ E3
        dy[self._i_p] = tangent
        dy[self._i_R] = (R @ hat(u)).ravel()

        l_vec = load.l(s) if want_l else None
        l_axial = float(tangent @ l_vec) if want_l else 0.0

        dth = dy[self._i_th]
        duz = dy[self._i_uz]
        for i in active:
            c, sn = np.cos(theta[i]), np.sin(theta[i])
            u_ix = c * u[0] + sn * u[1]
            u_iy = -sn * u[0] + c * u[1]
            ux, uy = star[i]
            dth[i] = u_z[i]
            duz[i] = (self.k_b[i] / self.k_t[i]) * (u_ix * uy - u_iy * ux)
            # Only the tube the distributed moment acts on takes its axial
            # part; the others are free to twist past it.
            if want_l and i == torque_tube:
                duz[i] -= l_axial / self.k_t[i]

        dy[self._i_n] = -load.f(s)
        dy[self._i_m] = -np.cross(tangent, n_vec) - (l_vec if want_l else 0.0)
        return dy

    # -- forward integration ------------------------------------------------

    def _integrate(
        self,
        y0: np.ndarray,
        segs,
        betas: np.ndarray,
        load: ExternalLoad,
        want_l: bool,
        torque_tube: int,
        rtol: float,
        atol: float,
        n_points: Optional[int],
    ):
        """Integrate base-to-tip, recording each tube's torsion at its own tip.

        With ``n_points is None`` only the tip state and the tube-end torsions
        are returned (the shooting inner loop); otherwise the trajectory is
        densely sampled for a :class:`Solution`.
        """
        ends = betas + np.array([t.length for t in self.tubes])
        total = segs[-1][1] - segs[0][0]

        y = y0.copy()
        tip_u_z = np.full(self.n_tubes, np.nan)
        traj_s: List[np.ndarray] = []
        traj_y: List[np.ndarray] = []
        traj_act: List[Tuple[int, ...]] = []

        for a, b, act in segs:
            if n_points is None:
                t_eval = None
            else:
                k = max(2, int(round(n_points * (b - a) / total)) + 1)
                t_eval = np.linspace(a, b, k)

            sol = solve_ivp(
                self._deriv,
                (a, b),
                y,
                args=(act, betas, load, want_l, torque_tube),
                method="RK45",
                rtol=rtol,
                atol=atol,
                t_eval=t_eval,
                dense_output=False,
            )
            if not sol.success:  # pragma: no cover - stiff/ill-posed inputs
                raise ShootingError(f"integration failed on [{a:.4g}, {b:.4g}]: {sol.message}")

            y = sol.y[:, -1].copy()
            # Re-project the frame; RK45 does not preserve SO(3) exactly.
            y[self._i_R] = orthonormalize(y[self._i_R].reshape(3, 3)).ravel()

            if t_eval is not None:
                keep = slice(1, None) if traj_s else slice(None)
                traj_s.append(sol.t[keep])
                traj_y.append(sol.y[:, keep])
                traj_act.extend([act] * sol.t[keep].shape[0])

            for i in act:
                if abs(ends[i] - b) <= 1e-7:
                    tip_u_z[i] = y[self._i_uz][i]

        if n_points is None:
            return y, tip_u_z, None
        return y, tip_u_z, (np.concatenate(traj_s), np.concatenate(traj_y, axis=1), traj_act)

    # -- boundary value problem ---------------------------------------------

    def _pack_y0(
        self,
        x: np.ndarray,
        emerged: np.ndarray,
        betas: np.ndarray,
        alphas: np.ndarray,
        loaded: bool,
    ) -> np.ndarray:
        y0 = np.zeros(self.state_dim)
        y0[self._i_R] = np.eye(3).ravel()

        offset = 6 if loaded else 0
        u_z0 = np.zeros(self.n_tubes)
        u_z0[emerged] = x[offset:]
        y0[self._i_uz] = u_z0
        # Windup accumulated in the straight transmission, s in [beta_i, 0].
        y0[self._i_th] = alphas - betas * u_z0
        if loaded:
            y0[self._i_n] = x[0:3]
            y0[self._i_m] = x[3:6]
        return y0

    def _restart_points(
        self, x0: np.ndarray, emerged: np.ndarray, betas: np.ndarray, loaded: bool
    ):
        """Deterministic spread of shooting starts, ordered cheapest-first.

        Scaled by ``pi / |beta_i|`` -- the base torsion that would wind half a
        turn into tube ``i``'s transmission -- so the spread tracks the actual
        compliance rather than a fixed number.  Sign patterns come first
        because the tubes' torsions balance and therefore usually have opposite
        signs.
        """
        k = int(emerged.sum())
        off = 6 if loaded else 0
        scale = np.pi / np.maximum(np.abs(betas[emerged]), 1e-3)

        def emit(vec):
            x = x0.copy()
            x[off:] = vec
            return x

        if np.any(x0[off:]):
            yield emit(np.zeros(k))  # cold start, if the caller gave a guess
            yield emit(-x0[off:])  # the mirrored branch

        signs = [
            np.array([1.0 if (j >> i) & 1 else -1.0 for i in range(k)])
            for j in range(min(1 << k, 16))
        ]
        for frac in (0.15, 0.35, 0.6, 1.0):
            for sgn in signs:
                yield emit(frac * scale * sgn)

    def solve(
        self,
        betas: Sequence[float],
        alphas: Sequence[float],
        load: Optional[ExternalLoad] = None,
        guess: Optional[Sequence[float]] = None,
        n_points: int = 200,
        rtol: float = 1e-9,
        atol: float = 1e-11,
        tol: float = 1e-6,
        max_nfev: int = 400,
        strict: bool = True,
    ) -> Solution:
        """Solve for the robot's shape.

        Parameters
        ----------
        betas
            Arc-length position of each tube's proximal end [m], ``<= 0``.
        alphas
            Commanded actuator roll of each tube [rad], applied at the tube's
            proximal end.
        load
            External loading; ``None`` means unloaded.
        guess
            Base torsions ``u_i,z(0)`` to start from -- typically ``u_z0`` from
            a neighbouring solution.  ``nan`` entries fall back to zero.
        n_points
            Approximate number of arc-length samples in the returned solution.
        rtol, atol
            Tolerances passed to the RK45 integrator.
        tol
            Convergence threshold on the shooting residual norm.  The residual
            mixes units -- torsion components are 1/m, force components N and
            moment components N.m -- so the default 1e-6 corresponds to a
            twist error well below a microradian over these tube lengths.
        strict
            Raise :class:`ShootingError` if the residual exceeds ``tol``.
            With ``strict=False`` the unconverged solution is returned and the
            caller can inspect ``residual_norm``.
        """
        betas = np.asarray(betas, dtype=float)
        alphas = np.asarray(alphas, dtype=float)
        if alphas.shape != (self.n_tubes,):
            raise ValueError(f"alphas must have length {self.n_tubes}")
        load = load if load is not None else ExternalLoad()

        segs, L = self.segments(betas)
        ends = betas + np.array([t.length for t in self.tubes])
        emerged = ends > KNOT_TOL
        loaded = not load.is_zero
        want_l = load.has_distributed_moment
        n_unknown = int(emerged.sum()) + (6 if loaded else 0)

        # Which tube carries torque applied about the backbone tangent.
        # Default: whichever tube reaches the distal tip.
        if load.axial_torque_tube is None:
            torque_tube = int(np.argmax(ends))
        else:
            torque_tube = int(load.axial_torque_tube)
            if not (0 <= torque_tube < self.n_tubes) or not emerged[torque_tube]:
                raise ValueError(
                    f"axial_torque_tube={load.axial_torque_tube} is not an emerged tube"
                )

        x0 = np.zeros(n_unknown)
        if guess is not None:
            g = np.asarray(guess, dtype=float).reshape(-1)
            if g.size == self.n_tubes:
                g = g[emerged]
            if g.size == int(emerged.sum()):
                x0[(6 if loaded else 0):] = np.nan_to_num(g)
        if loaded:
            # A free tip carries exactly the applied wrench; with light
            # distributed loads that is already close to n(0), m(0).
            x0[0:3] = load.tip_force
            x0[3:6] = load.tip_moment

        def residual(x: np.ndarray) -> np.ndarray:
            y0 = self._pack_y0(x, emerged, betas, alphas, loaded)
            y_end, tip_u_z, _ = self._integrate(
                y0, segs, betas, load, want_l, torque_tube, rtol, atol, None
            )
            bc = tip_u_z.copy()
            if loaded and np.any(load.tip_moment):
                # The tube the tip wrench is attached to does not have a free
                # end: it carries the applied torque about the tangent.  Every
                # other tube still ends free, at u_z = 0.
                tangent = y_end[self._i_R].reshape(3, 3) @ E3
                tau = float(tangent @ load.tip_moment)
                bc[torque_tube] -= tau / self.k_t[torque_tube]
            r = list(bc[emerged])
            if loaded:
                r.extend(y_end[self._i_n] - load.tip_force)
                r.extend(y_end[self._i_m] - load.tip_moment)
            return np.asarray(r, dtype=float)

        if n_unknown == 0:  # pragma: no cover - segments() rejects this first
            raise ValueError("no emerged tubes to solve for")

        def run(start: np.ndarray):
            # The finite-difference step has to sit above the adaptive
            # integrator's own step-to-step jitter (~rtol) or the Jacobian is
            # dominated by quadrature noise; 1e-6 leaves three decades of margin.
            return least_squares(
                residual, start, method="lm", diff_step=1e-6,
                xtol=1e-14, ftol=1e-14, gtol=1e-14, max_nfev=max_nfev,
            )

        out = run(x0)
        restarts = 0
        if float(np.linalg.norm(out.fun)) > tol:
            # A wound concentric-tube robot has several torsional equilibria,
            # and branches are born and destroyed at fold points as the tubes
            # are rotated (the snap-through this design is prone to).  When the
            # start point's branch has just folded away, Levenberg-Marquardt
            # stalls in its basin instead of finding the surviving root, so
            # retry from a spread of physically scaled starts.
            for start in self._restart_points(x0, emerged, betas, loaded):
                restarts += 1
                trial = run(start)
                if float(np.linalg.norm(trial.fun)) < float(np.linalg.norm(out.fun)):
                    out = trial
                if float(np.linalg.norm(out.fun)) <= tol:
                    break

        res_norm = float(np.linalg.norm(out.fun))
        if strict and res_norm > tol:
            raise ShootingError(
                f"shooting did not converge: |residual| = {res_norm:.3e} > {tol:.3e} "
                f"after {restarts + 1} starts ({out.message}).  If this is a "
                "rotation sweep, step the actuator angles more finely and pass "
                "each solution's u_z0 as the next guess."
            )

        y0 = self._pack_y0(out.x, emerged, betas, alphas, loaded)
        _, _, traj = self._integrate(
            y0, segs, betas, load, want_l, torque_tube, rtol, atol, n_points
        )
        s_arr, y_arr, act = traj

        M = s_arr.size
        theta = y_arr[self._i_th].T.copy()
        u_z = y_arr[self._i_uz].T.copy()
        R_arr = y_arr[self._i_R].T.reshape(M, 3, 3)
        m_arr = y_arr[self._i_m].T.copy()

        u_arr = np.zeros((M, 3))
        for k in range(M):
            acc = (R_arr[k].T @ m_arr[k])[:2].copy()
            k_b_sum = 0.0
            for i in act[k]:
                ux, uy = self.tubes[i].precurvature(s_arr[k] - betas[i])
                c, sn = np.cos(theta[k, i]), np.sin(theta[k, i])
                acc[0] += self.k_b[i] * (c * ux - sn * uy)
                acc[1] += self.k_b[i] * (sn * ux + c * uy)
                k_b_sum += self.k_b[i]
            u_arr[k, :2] = acc / k_b_sum

        for k in range(M):
            absent = [i for i in range(self.n_tubes) if i not in act[k]]
            theta[k, absent] = np.nan
            u_z[k, absent] = np.nan

        u_z0 = np.full(self.n_tubes, np.nan)
        u_z0[emerged] = out.x[(6 if loaded else 0):]

        return Solution(
            s=s_arr,
            p=y_arr[self._i_p].T.copy(),
            R=R_arr,
            theta=theta,
            u_z=u_z,
            u=u_arr,
            n=y_arr[self._i_n].T.copy(),
            m=m_arr,
            active=act,
            u_z0=u_z0,
            residual=out.fun,
            n_iter=int(out.nfev),
            restarts=restarts,
            betas=betas,
            alphas=alphas,
        )
