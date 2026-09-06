"""
loads.py -- external loading for the geometrically exact concentric-tube model.

The model of Rucker, Jones & Webster (2010) carries the internal force ``n(s)``
and internal moment ``m(s)`` as ODE states, driven by

    n' = -f          m' = -p' x n - l

with ``f`` a distributed force per unit length [N/m] and ``l`` a distributed
moment per unit length [N.m/m], both expressed in the **global** frame, and
closed by the tip boundary conditions ``n(L) = F_tip``, ``m(L) = M_tip``.

An :class:`ExternalLoad` holds exactly those four quantities.  The default is
the unloaded case, for which the solver drops ``n`` and ``m`` from the shooting
unknowns entirely (they are identically zero) and reduces to the classical
torsionally-compliant concentric-tube kinematics.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Optional, Union

import numpy as np

Density = Union[None, np.ndarray, Callable[[float], np.ndarray]]


def _as_vec3(v, name: str) -> np.ndarray:
    a = np.asarray(v, dtype=float).reshape(-1)
    if a.size != 3:
        raise ValueError(f"{name} must be a 3-vector, got shape {np.shape(v)}")
    return a


@dataclass
class ExternalLoad:
    """Point and distributed loads applied to the robot, in the global frame.

    Parameters
    ----------
    tip_force, tip_moment
        Wrench applied at the distal tip ``s = L`` [N], [N.m].
    force_density, moment_density
        Either a constant 3-vector or a callable ``f(s) -> (3,)`` giving the
        distributed load per unit arc length [N/m], [N.m/m].  ``None`` means
        zero.
    axial_torque_tube
        Index of the tube that physically carries any torque applied *about
        the backbone tangent*.  ``None`` means the tube that reaches the
        distal tip -- the inner tube on the CT-SDR, which is the one the drill
        bit is mounted to.

        This attribution is required, not cosmetic.  Forces and transverse
        moments are shared by every tube at a cross-section, because the tubes
        are forced onto one backbone.  Axial torque is not: the tubes are free
        to twist relative to each other, so torque about the tangent stays in
        whichever tube it was applied to, and each *other* tube's free end
        still satisfies ``u_z = 0``.

    Notes
    -----
    A drilling reaction is naturally expressed as a ``tip_force`` (thrust plus
    lateral cutting force) together with a ``tip_moment`` (the drill torque
    reaction about the tip tangent).
    """

    tip_force: np.ndarray = field(default_factory=lambda: np.zeros(3))
    tip_moment: np.ndarray = field(default_factory=lambda: np.zeros(3))
    force_density: Density = None
    moment_density: Density = None
    axial_torque_tube: Optional[int] = None

    def __post_init__(self) -> None:
        self.tip_force = _as_vec3(self.tip_force, "tip_force")
        self.tip_moment = _as_vec3(self.tip_moment, "tip_moment")
        if self.force_density is not None and not callable(self.force_density):
            self.force_density = _as_vec3(self.force_density, "force_density")
        if self.moment_density is not None and not callable(self.moment_density):
            self.moment_density = _as_vec3(self.moment_density, "moment_density")

    # -- evaluation ---------------------------------------------------------

    def f(self, s: float) -> np.ndarray:
        """Distributed force per unit length at arc length ``s`` [N/m]."""
        return self._eval(self.force_density, s)

    def l(self, s: float) -> np.ndarray:
        """Distributed moment per unit length at arc length ``s`` [N.m/m]."""
        return self._eval(self.moment_density, s)

    @staticmethod
    def _eval(d: Density, s: float) -> np.ndarray:
        if d is None:
            return _ZERO3
        if callable(d):
            return np.asarray(d(s), dtype=float).reshape(3)
        return d

    # -- introspection ------------------------------------------------------

    @property
    def is_zero(self) -> bool:
        """True when the robot is completely unloaded.

        Callable densities are assumed non-zero; pass ``None`` (the default)
        for the genuinely unloaded case so the solver can take the cheaper
        path.
        """
        if callable(self.force_density) or callable(self.moment_density):
            return False
        return not (
            np.any(self.tip_force)
            or np.any(self.tip_moment)
            or (self.force_density is not None and np.any(self.force_density))
            or (self.moment_density is not None and np.any(self.moment_density))
        )

    @property
    def has_distributed_moment(self) -> bool:
        if callable(self.moment_density):
            return True
        return self.moment_density is not None and bool(np.any(self.moment_density))

    def __repr__(self) -> str:  # pragma: no cover - display only
        if self.is_zero:
            return "ExternalLoad(unloaded)"
        return (
            f"ExternalLoad(tip_force={np.round(self.tip_force, 6).tolist()} N, "
            f"tip_moment={np.round(self.tip_moment, 8).tolist()} N.m, "
            f"distributed={'yes' if self.force_density is not None else 'no'})"
        )


_ZERO3 = np.zeros(3)

#: Convenient shared instance for the unloaded case.
NO_LOAD = ExternalLoad()
