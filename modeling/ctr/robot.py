"""
robot.py -- the CT-SDR itself: joint values in, backbone shape out.

Bridges the repo's joint space (OTT/ITT in mm, OTR/ITR in degrees -- the same
names and units the ROS1 nodes and ``dxl_control_4dof_cli.py`` use) to the SI
tube placement (``betas``, ``alphas``) that :class:`ctr.model.CosseratModel`
integrates.

Frame convention
----------------
``s = 0`` is the distal face of the rigid stainless-steel guide; +z points
along the guide axis, out of it.  With all four joints at zero the tube tips
sit at the origin, so ``OTT`` and ``ITT`` read out deployed length directly.
"""

from __future__ import annotations

import os
from dataclasses import dataclass, replace
from typing import Optional, Sequence, Tuple

import numpy as np
import yaml

from .friction import FrictionModel
from .loads import ExternalLoad
from .model import CosseratModel, Solution
from .tube import Tube

DEFAULT_CONFIG = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "config", "ct_sdr.yaml"
)

OUTER = 0
INNER = 1


@dataclass(frozen=True)
class Joints:
    """A CT-SDR joint-space configuration.

    Attributes
    ----------
    ott, itt
        Outer / inner tube translation [mm], positive out of the guide.
    otr, itr
        Outer / inner tube rotation [deg], right-handed about +z.
    """

    ott: float = 0.0
    itt: float = 0.0
    otr: float = 0.0
    itr: float = 0.0

    def as_tuple(self) -> Tuple[float, float, float, float]:
        return (self.ott, self.itt, self.otr, self.itr)

    def __str__(self) -> str:
        return (
            f"OTT={self.ott:.1f} mm  ITT={self.itt:.1f} mm  "
            f"OTR={self.otr:.1f} deg  ITR={self.itr:.1f} deg"
        )


class CTSDR:
    """Concentric-tube steerable drilling robot: two nitinol tubes in a guide.

    Parameters
    ----------
    outer, inner
        Pre-built :class:`~ctr.tube.Tube` objects.  Usually you want
        :meth:`from_yaml` instead.
    base_offsets
        Arc-length position of each tube's proximal end at zero translation
        [m].  Defaults to ``-length`` for each tube, i.e. tips flush with the
        guide exit.
    signs
        Multipliers applied to ``(ott, itt, otr, itr)``.
    roll_offsets
        Roll of each tube's pre-curvature plane when its rotation joint reads
        zero [rad].  **Not a nuisance parameter -- the rotation joints have no
        absolute reference.**  ``dxl_control_4dof_cli.py``'s ``home <joint>``
        zeroes *wherever the tube currently is*, so OTR = ITR = 0 means "homed
        here", not "pre-curvatures aligned".  The offset that matters is the
        difference ``roll_offsets[INNER] - roll_offsets[OUTER]``; it can differ
        between sessions and must be calibrated per experiment.
    """

    def __init__(
        self,
        outer: Tube,
        inner: Tube,
        base_offsets: Optional[Sequence[float]] = None,
        signs: Sequence[float] = (1.0, 1.0, 1.0, 1.0),
        roll_offsets: Optional[Sequence[float]] = None,
    ):
        self.tubes = [outer, inner]
        self.model = CosseratModel(self.tubes)
        if base_offsets is None:
            self.base_offsets = np.array([-outer.length, -inner.length])
        else:
            self.base_offsets = np.asarray(base_offsets, dtype=float)
        self.signs = np.asarray(signs, dtype=float)
        self.roll_offsets = (np.zeros(2) if roll_offsets is None
                             else np.asarray(roll_offsets, dtype=float))

    # -- construction -------------------------------------------------------

    @classmethod
    def from_yaml(cls, path: str = DEFAULT_CONFIG) -> "CTSDR":
        with open(path, "r") as fh:
            cfg = yaml.safe_load(fh)

        mat = cfg["material"]
        E = float(mat["youngs_modulus_gpa"]) * 1e9
        nu = float(mat["poisson_ratio"])
        rho = float(mat.get("density_kg_m3", 6450.0))

        tubes, offsets, rolls = [], [], []
        for name in ("outer", "inner"):
            t = cfg["tubes"][name]
            od = float(t["outer_diameter_mm"]) * 1e-3
            wall = float(t["wall_thickness_mm"]) * 1e-3
            length = float(t["length_mm"]) * 1e-3
            roc = float(t["radius_of_curvature_mm"]) * 1e-3
            tubes.append(
                Tube(
                    name=name,
                    outer_diameter=od,
                    inner_diameter=od - 2.0 * wall,
                    length=length,
                    curved_length=float(t["curved_length_mm"]) * 1e-3,
                    curvature=1.0 / roc,
                    youngs_modulus=E,
                    poisson_ratio=nu,
                    density=rho,
                    tip_straight_length=float(t.get("tip_straight_length_mm", 0.0)) * 1e-3,
                )
            )
            off = t.get("base_offset_mm")
            offsets.append(-length if off is None else float(off) * 1e-3)
            rolls.append(np.radians(float(t.get("roll_offset_deg", 0.0))))

        if tubes[INNER].outer_diameter > tubes[OUTER].inner_diameter:
            raise ValueError(
                f"inner tube OD {tubes[INNER].outer_diameter * 1e3:.2f} mm does not "
                f"fit in outer tube ID {tubes[OUTER].inner_diameter * 1e3:.2f} mm"
            )

        j = cfg.get("joints", {})
        signs = (
            float(j.get("ott_sign", 1.0)),
            float(j.get("itt_sign", 1.0)),
            float(j.get("otr_sign", 1.0)),
            float(j.get("itr_sign", 1.0)),
        )
        return cls(tubes[OUTER], tubes[INNER], base_offsets=offsets, signs=signs,
                   roll_offsets=rolls)

    # -- joint mapping ------------------------------------------------------

    def joints_to_model(self, joints: Joints) -> Tuple[np.ndarray, np.ndarray]:
        """Map ``Joints`` to ``(betas, alphas)`` in metres and radians."""
        ott, itt, otr, itr = joints.as_tuple()
        s = self.signs
        betas = self.base_offsets + np.array([s[0] * ott, s[1] * itt]) * 1e-3
        alphas = np.radians([s[2] * otr, s[3] * itr]) + self.roll_offsets
        return betas, alphas

    def deployed_lengths(self, joints: Joints) -> np.ndarray:
        """Length of each tube protruding from the guide [m] (clipped at 0)."""
        betas, _ = self.joints_to_model(joints)
        ends = betas + np.array([t.length for t in self.tubes])
        return np.maximum(ends, 0.0)

    def check(self, joints: Joints) -> None:
        """Raise if a configuration is outside what the model can represent.

        Catches the two failures that would otherwise produce a silently wrong
        answer: retracting past the guide exit, and advancing a tube so far
        that its straight proximal section has emerged (which the constant
        pre-curvature bookkeeping handles, but which almost always means the
        joint value is a mistake rather than a real configuration).
        """
        betas, _ = self.joints_to_model(joints)
        for i, t in enumerate(self.tubes):
            end = betas[i] + t.length
            if end > t.length + 1e-9:
                raise ValueError(f"{t.name} tube: commanded past its full length")
            if end > t.deployable_length + 1e-9:
                raise ValueError(
                    f"{t.name} tube: deployed {end * 1e3:.1f} mm exceeds its "
                    f"{t.deployable_length * 1e3:.1f} mm shaped section "
                    f"({t.curved_length * 1e3:.1f} mm curved + "
                    f"{t.tip_straight_length * 1e3:.1f} mm tip lead-in), so part of the "
                    "straight proximal section has emerged"
                )
        if np.max(betas + np.array([t.length for t in self.tubes])) <= 0.0:
            raise ValueError("nothing is deployed: advance OTT and/or ITT")

    # -- solving ------------------------------------------------------------

    def solve(
        self,
        joints: Joints,
        load: Optional[ExternalLoad] = None,
        friction: Optional[FrictionModel] = None,
        guess: Optional[Sequence[float]] = None,
        **kwargs,
    ) -> Solution:
        """Solve the boundary value problem for one joint configuration."""
        betas, alphas = self.joints_to_model(joints)
        return self.model.solve(betas, alphas, load=load, friction=friction,
                                guess=guess, **kwargs)

    def solve_path(
        self,
        path: Sequence[Joints],
        load: Optional[ExternalLoad] = None,
        friction: Optional[FrictionModel] = None,
        **kwargs,
    ):
        """Solve a sequence of configurations, warm-starting each from the last.

        Continuation matters here: a concentric-tube robot's torsional
        equations admit multiple equilibria once the tubes are wound up, and
        marching from a neighbouring solution keeps the solver on the branch
        the hardware actually follows instead of letting it jump.

        Friction, if given, has its sliding ``direction`` set from each step's
        own rotation increment rather than used as supplied.  That is what
        makes a forward sweep and the reverse of it different: run the same
        path in both orders and the gap between them is the modelled
        hysteresis.
        """
        guess = None
        out = []
        prev = None
        for j in path:
            # Nothing has moved yet at the first waypoint, so no contact has
            # a sliding sense: friction is off there by construction.
            step_friction = None if prev is None else friction
            if friction is not None and prev is not None:
                # Relative roll, through the configured signs: a co-rotation
                # leaves this unchanged and correctly yields no sliding.
                _, a_now = self.joints_to_model(j)
                _, a_prev = self.joints_to_model(prev)
                d = (a_now[INNER] - a_now[OUTER]) - (a_prev[INNER] - a_prev[OUTER])
                step_friction = replace(friction, direction=float(np.sign(d)))
            sol = self.solve(j, load=load, friction=step_friction,
                             guess=guess, **kwargs)
            guess = sol.u_z0
            prev = j
            out.append(sol)
        return out

    # -- loading helpers ----------------------------------------------------

    def gravity_load(
        self,
        joints: Joints,
        gravity: Sequence[float] = (0.0, 0.0, -9.81),
        tip_force: Sequence[float] = (0.0, 0.0, 0.0),
        tip_moment: Sequence[float] = (0.0, 0.0, 0.0),
    ) -> ExternalLoad:
        """Self-weight of the deployed tubes, plus an optional tip wrench.

        The distributed weight steps at each tube's tip, so it is built as a
        function of ``s`` from whichever tubes are present there.
        """
        betas, _ = self.joints_to_model(joints)
        ends = betas + np.array([t.length for t in self.tubes])
        g = np.asarray(gravity, dtype=float)
        mu = np.array([t.mass_per_length for t in self.tubes])

        def f(s: float) -> np.ndarray:
            present = (betas <= s) & (s <= ends)
            return g * float(mu[present].sum())

        return ExternalLoad(tip_force=tip_force, tip_moment=tip_moment, force_density=f)

    def drilling_load(
        self,
        thrust: float,
        lateral: Sequence[float] = (0.0, 0.0),
        torque: float = 0.0,
        solution: Optional[Solution] = None,
        joints: Optional[Joints] = None,
    ) -> ExternalLoad:
        """Build a tip wrench for a drilling reaction.

        The reaction is defined in the *tip* frame -- ``thrust`` opposes the
        feed direction along the tip tangent, ``lateral`` is the transverse
        cutting force, ``torque`` is the drill reaction about the tangent --
        and is rotated into the global frame that the model needs.  Because
        the tip frame depends on the deflected shape, pass either an unloaded
        ``solution`` or the ``joints`` to solve one from; for the modest
        deflections here one such pass is enough, and a second call using the
        loaded solution's tip frame gives the fixed point if needed.

        ``torque`` is left attributed to the tube that reaches the tip, which
        on this robot is the inner tube -- the one the drill bit is on.  Pass
        an :class:`~ctr.loads.ExternalLoad` with ``axial_torque_tube`` set
        explicitly if that is not true of your end effector.
        """
        if solution is None:
            if joints is None:
                raise ValueError("pass either solution= or joints=")
            solution = self.solve(joints)
        R_tip = solution.tip_rotation
        f_tip = np.array([lateral[0], lateral[1], -float(thrust)])
        m_tip = np.array([0.0, 0.0, -float(torque)])
        return ExternalLoad(tip_force=R_tip @ f_tip, tip_moment=R_tip @ m_tip)

    # -- reporting ----------------------------------------------------------

    def material_limits(self, solution: Solution, shear_limit: float = 250e6):
        """Peak strain and torsional shear stress per tube for a solution.

        Returns one dict per tube with ``bending_strain``, ``shear_stress``
        [Pa], ``torque`` [N.m] and ``shear_ok``.  Worth checking on any
        heavily wound configuration: the model is linear-elastic, so once the
        predicted shear passes nitinol's inelastic threshold its numbers are
        an extrapolation rather than a prediction.  ``shear_limit`` defaults to
        a conservative 250 MPa -- set it to your own material data.
        """
        out = []
        kappa = np.linalg.norm(solution.u[:, :2], axis=1)
        for i, t in enumerate(self.tubes):
            present = ~np.isnan(solution.u_z[:, i])
            if not np.any(present):
                out.append(None)
                continue
            col = solution.u_z[present, i]
            u_z_peak = float(col[np.argmax(np.abs(col))])  # signed peak
            tau = t.shear_stress(u_z_peak)
            out.append(
                {
                    "tube": t.name,
                    "bending_strain": t.bending_strain(float(kappa[present].max())),
                    "shear_stress": tau,
                    "torque": t.torsional_moment(u_z_peak),
                    "shear_ok": tau <= shear_limit,
                }
            )
        return out

    def summary(self) -> str:
        lines = ["CT-SDR: two pre-curved nitinol tubes in a rigid steel guide"]
        lines += ["  " + t.summary() for t in self.tubes]
        ratio = self.tubes[OUTER].bending_stiffness / self.tubes[INNER].bending_stiffness
        lines.append(f"  outer/inner bending stiffness ratio: {ratio:.2f}")
        return "\n".join(lines)
