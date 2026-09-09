"""
friction.py -- Coulomb friction between the tubes and against the guide.

The Cosserat model in :mod:`ctr.model` is frictionless: tube-to-tube contact is
assumed to transmit no axial torque, which is the standard assumption in the
concentric-tube literature and is what makes each tube's moment balance close
independently about its own axis.  This module relaxes that.

Where friction acts
-------------------
Two contacts matter on this robot, with different materials:

* **nitinol on nitinol** -- the inner tube inside the outer tube, wherever they
  overlap.  Dry NiTi-on-NiTi is adhesive and sits high, ``mu ~ 0.3-0.6``;
  :data:`MU_NITINOL_NITINOL` defaults to 0.35.
* **nitinol on stainless steel** -- the outer tube inside the rigid guide, and
  the inner tube against the guide wherever the outer tube is not between them.
  Lower, ``mu ~ 0.2-0.35``; :data:`MU_STEEL_NITINOL` defaults to 0.25.

Both are *dry* values, appropriate for a drilling robot running in air or bone.
Wetted with saline or blood they drop to roughly 0.1-0.2, which
:meth:`FrictionModel.lubricated` returns.  These are order-of-magnitude
literature values, not measurements of your tubes -- treat them as a band and
run :func:`ctr.friction.sensitivity` style sweeps rather than trusting one
number.

What it does to the torque balance
----------------------------------
For tube ``i`` rubbing on tube ``j`` with normal force per unit length ``w`` at
contact radius ``r``, the friction torque per unit length about the tube axis is

    tau = -d_ij * mu * w * r

with ``d_ij`` the sign of tube ``i``'s rotation *relative to* tube ``j``.  It
enters each tube's torsion ODE as an extra distributed axial moment, equal and
opposite on the two tubes, so total axial moment is still conserved.

This is a **kinetic** model: every contact is assumed to be sliding, with the
direction set by the sense of the commanded manoeuvre
(:attr:`FrictionModel.direction`).  That is right for a monotonic rotation and
wrong at a reversal, where parts of the tube stick before they slip.  Modelling
stick properly needs a complementarity solve; what this model does capture is
the *hysteresis* between a forward and a reverse sweep, which is the
measurable consequence.

Contact force
-------------
The normal force per unit length is estimated from curvature mismatch.  A rod
held at curvature ``u`` when its natural curvature is ``u*`` carries a bending
moment ``dM = k_b |u - u*|``; holding it there over a characteristic length
``1/kappa`` needs a transverse load of order

    w ~ kappa^2 * dM

which is the scaling used here (:func:`contact_force_per_length`).  It is an
estimate, good to a factor of about two, and it is deliberately exposed so it
can be overridden with a measured value -- see ``force_scale``.

Two things it correctly predicts as *small* on this robot:

* a uniformly pre-curved tube held straight inside the guide carries a constant
  moment, so it needs **no** distributed transverse load -- the reaction is
  concentrated at the guide's lip, not spread along it.  That is why the
  transmission's friction is modelled from a separate, explicit
  :attr:`guide_normal_force` rather than from curvature;
* in the deployed section ``kappa`` is only ~10-25 /m and the overlap is tens
  of millimetres, so the total normal force is a few newtons.

Whether that adds up to anything is an empirical question; see
``modeling/validation/friction_study.py``.

What is *not* modelled
----------------------
* **Axial friction.**  Only torque about the tube axis is resisted here.  A
  pure translation changes no relative roll, so ``direction`` is zero and this
  module contributes nothing at all -- correct for the torsional part, missing
  for the axial part.  Advancing a tube does slide it along the guide, and at
  the ``set2``-calibrated contact force that is ~40 N of resistance; against
  the inner tube's 90 kN axial stiffness it would cost only ~0.06 mm of tip
  lag, against a measured 0.66 mm translation undershoot.  So it is a real
  omission but a small one -- see VALIDATION_REPORT.md section 7.
* **The friction cone.**  Sliding axially and circumferentially at once splits
  one friction budget between the two directions, so torsional resistance
  should *fall* while a tube is translating.  Here the two are independent.
  The ICRA2027 trials translate and rotate in separate steps, so it does not
  bite; a helical or simultaneous motion would need it.
* **Stick.**  Every contact is assumed to be sliding -- see above.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Sequence

import numpy as np

#: Dry nitinol on nitinol.  Literature spread is wide (~0.3-0.6); NiTi-NiTi is
#: adhesive and galls readily, which is why concentric-tube designs often
#: passivate or coat the mating surfaces.
MU_NITINOL_NITINOL = 0.35

#: Dry nitinol on austenitic stainless steel, ~0.2-0.35.
MU_STEEL_NITINOL = 0.25

#: Wetted (saline / blood / irrigation) values for the same two pairs.
MU_NITINOL_NITINOL_WET = 0.15
MU_STEEL_NITINOL_WET = 0.12


def contact_force_per_length(kappa: float, k_b: float,
                             u_xy: Sequence[float],
                             ustar_xy: Sequence[float]) -> float:
    """Normal force per unit length between a tube and its neighbour [N/m].

    ``kappa`` is the magnitude of the common curvature, ``k_b`` the tube's
    bending stiffness, and ``u_xy`` / ``ustar_xy`` its actual and natural
    curvature in its own frame.  See the module docstring for the scaling; this
    returns ``kappa^2 * k_b * |u - u*|`` and is zero where a tube is at its
    natural curvature or where the assembly is straight.
    """
    d = np.asarray(u_xy, dtype=float) - np.asarray(ustar_xy, dtype=float)
    return float(kappa * kappa * k_b * np.linalg.norm(d))


@dataclass
class FrictionModel:
    """Coulomb friction between the tubes and against the guide.

    Parameters
    ----------
    mu_tube_tube, mu_guide
        Friction coefficients for the nitinol-nitinol and steel-nitinol
        contacts.
    direction
        Sense in which the **inner** tube of each contacting pair is rotating
        relative to the outer one as the manoeuvre proceeds: ``+1``, ``-1``, or
        ``0`` for no relative sliding (which disables friction).  Friction
        opposes it.  :meth:`ctr.robot.CTSDR.solve_path` sets it per step from
        the change in commanded *relative* roll, so a pure co-rotation or a
        pure translation correctly gets ``0`` and a reversal flips the sign --
        that flip is what produces hysteresis.
    guide_normal_force
        Normal force per unit length between a tube and the guide over the
        retracted length [N/m].  Not derived from curvature: a uniformly
        pre-curved tube held straight needs no distributed load (module
        docstring), so what is left is the guide's end reaction spread over
        its grip, plus tube weight and clearance-driven contact.  Default is
        :meth:`estimate_guide_force`.
    force_scale
        Multiplies the estimated contact force.  The estimate is good to about
        a factor of two, so sweeping this from ~0.5 to ~3 is the honest way to
        report a friction result.
    """

    mu_tube_tube: float = MU_NITINOL_NITINOL
    mu_guide: float = MU_STEEL_NITINOL
    direction: float = 1.0
    guide_normal_force: float = 0.0
    force_scale: float = 1.0

    @classmethod
    def dry(cls, direction: float = 1.0, **kw) -> "FrictionModel":
        return cls(MU_NITINOL_NITINOL, MU_STEEL_NITINOL, direction, **kw)

    @classmethod
    def lubricated(cls, direction: float = 1.0, **kw) -> "FrictionModel":
        """Saline / blood / irrigation wetted contacts."""
        return cls(MU_NITINOL_NITINOL_WET, MU_STEEL_NITINOL_WET, direction, **kw)

    @property
    def is_zero(self) -> bool:
        return (self.direction == 0.0
                or (self.mu_tube_tube == 0.0 and self.mu_guide == 0.0)
                or self.force_scale == 0.0)

    # -- deployed section ---------------------------------------------------

    def tube_tube_torque(self, kappa: float, k_b: float, u_xy, ustar_xy,
                         contact_radius: float) -> float:
        """Friction torque per unit length at a tube-tube contact [N.m/m].

        Returned as a positive magnitude; the caller applies the sign from the
        relative rotation sense.
        """
        w = contact_force_per_length(kappa, k_b, u_xy, ustar_xy) * self.force_scale
        return self.mu_tube_tube * w * contact_radius

    # -- retracted transmission ---------------------------------------------

    def transmission_torque(self, contact_radius: float) -> float:
        """Friction torque per unit length in the guide [N.m/m]."""
        return self.mu_guide * self.guide_normal_force * self.force_scale * contact_radius

    @staticmethod
    def estimate_guide_force(tube, grip_length: float) -> float:
        """Guide reaction per unit length while straightening a pre-curved tube.

        The tube's own moment while held straight is ``k_b * kappa*``.  The
        guide reacts it as a couple over ``grip_length``, giving a force
        ``k_b kappa* / grip_length`` at each end; spreading that over the grip
        gives the per-unit-length figure returned here.  It is a crude but
        bounded estimate -- the real distribution depends on the bore
        clearance, which this model does not resolve.
        """
        if grip_length <= 0.0:
            return 0.0
        return tube.bending_stiffness * tube.curvature / (grip_length ** 2)

    def __repr__(self) -> str:  # pragma: no cover - display only
        if self.is_zero:
            return "FrictionModel(off)"
        return (f"FrictionModel(mu_tt={self.mu_tube_tube:.2f}, "
                f"mu_guide={self.mu_guide:.2f}, dir={self.direction:+.0f}, "
                f"w_guide={self.guide_normal_force:.1f} N/m, "
                f"scale={self.force_scale:g})")


#: Shared frictionless instance.
NO_FRICTION = FrictionModel(direction=0.0)
