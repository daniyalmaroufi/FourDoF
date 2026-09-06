"""
tube.py -- geometry, material and pre-curvature of a single concentric tube.

Everything here is SI (metres, pascals, 1/metres, radians).  The millimetre
numbers from the CT-SDR datasheet live in ``modeling/config/ct_sdr.yaml`` and
are converted once, on load, by :mod:`ctr.robot`.

Pre-curvature convention
------------------------
Each tube is straight over its proximal ``length - curved_length`` and
circularly pre-curved over its distal ``curved_length``.  The pre-curvature is
taken about the tube's own local x-axis::

    u*(sigma) = [kappa, 0, 0]^T      sigma >= straight_length
                [0,     0, 0]^T      otherwise

so a lone tube released from a rigid guide bends in the local y-z plane with
its centre of curvature on the -y side.  Rotating the tube about its axis
(the OTR/ITR joints) rotates that bending plane; that rotation enters the
model as ``theta``, not as a change of ``u*``.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Tuple

import numpy as np


@dataclass(frozen=True)
class Tube:
    """One pre-curved elastic tube of a concentric-tube robot.

    Parameters
    ----------
    name
        Label used in reports and plots, e.g. ``"outer"``.
    outer_diameter, inner_diameter
        Cross-section diameters [m].  Assumed constant along the tube.
    length
        Total tube length [m], including the proximal part that stays inside
        the actuation unit / rigid guide.
    curved_length
        Length of the distal pre-curved portion [m].  The remaining
        ``length - curved_length`` at the proximal end is straight.
    curvature
        Pre-curvature ``kappa`` of the curved portion [1/m] (the reciprocal of
        the radius of curvature).
    youngs_modulus, poisson_ratio
        Isotropic linear-elastic properties.  Nitinol in its austenitic,
        superelastic state is modelled as linear here, which is the standard
        assumption in the concentric-tube-robot literature and is reasonable
        for the small strains a 50 mm radius of curvature produces in these
        wall thicknesses.
    density
        Mass density [kg/m^3], used only when gravity loading is requested.
    """

    name: str
    outer_diameter: float
    inner_diameter: float
    length: float
    curved_length: float
    curvature: float
    youngs_modulus: float
    poisson_ratio: float
    density: float = 6450.0  # nitinol

    def __post_init__(self) -> None:
        if self.inner_diameter < 0.0 or self.outer_diameter <= self.inner_diameter:
            raise ValueError(
                f"{self.name}: need 0 <= inner_diameter < outer_diameter, got "
                f"{self.inner_diameter} and {self.outer_diameter}"
            )
        if self.length <= 0.0:
            raise ValueError(f"{self.name}: length must be positive")
        if not 0.0 <= self.curved_length <= self.length:
            raise ValueError(
                f"{self.name}: curved_length {self.curved_length} must lie in "
                f"[0, length={self.length}]"
            )
        if self.curvature < 0.0:
            raise ValueError(f"{self.name}: curvature must be non-negative")
        if self.youngs_modulus <= 0.0:
            raise ValueError(f"{self.name}: youngs_modulus must be positive")
        if not -1.0 < self.poisson_ratio < 0.5:
            raise ValueError(f"{self.name}: poisson_ratio out of physical range")

    # -- derived geometry ---------------------------------------------------

    @property
    def straight_length(self) -> float:
        """Length of the proximal straight portion [m]."""
        return self.length - self.curved_length

    @property
    def wall_thickness(self) -> float:
        return 0.5 * (self.outer_diameter - self.inner_diameter)

    @property
    def area(self) -> float:
        """Cross-sectional area [m^2]."""
        return 0.25 * np.pi * (self.outer_diameter ** 2 - self.inner_diameter ** 2)

    @property
    def second_moment(self) -> float:
        """Second moment of area ``I`` about a diameter [m^4]."""
        return (np.pi / 64.0) * (self.outer_diameter ** 4 - self.inner_diameter ** 4)

    @property
    def polar_moment(self) -> float:
        """Polar second moment ``J = 2I`` [m^4] (thin/thick circular annulus)."""
        return 2.0 * self.second_moment

    @property
    def shear_modulus(self) -> float:
        """``G = E / (2 (1 + nu))`` [Pa]."""
        return self.youngs_modulus / (2.0 * (1.0 + self.poisson_ratio))

    # -- stiffnesses --------------------------------------------------------

    @property
    def bending_stiffness(self) -> float:
        """``k_b = E I`` [N.m^2].  Isotropic in x and y for a round tube."""
        return self.youngs_modulus * self.second_moment

    @property
    def torsional_stiffness(self) -> float:
        """``k_t = G J`` [N.m^2]."""
        return self.shear_modulus * self.polar_moment

    @property
    def stiffness_diag(self) -> np.ndarray:
        """Diagonal of ``K = diag(EI, EI, GJ)``."""
        kb = self.bending_stiffness
        return np.array([kb, kb, self.torsional_stiffness])

    @property
    def arc_angle(self) -> float:
        """Total angle subtended by the pre-curved portion [rad]."""
        return self.curvature * self.curved_length

    @property
    def mass_per_length(self) -> float:
        """Linear mass density [kg/m]."""
        return self.density * self.area

    # -- pre-curvature ------------------------------------------------------

    def precurvature(self, sigma: float) -> Tuple[float, float]:
        """Pre-curvature ``(u*_x, u*_y)`` at material coordinate ``sigma``.

        ``sigma`` is measured from the tube's own proximal end.  Values
        outside ``[0, length]`` return zero so the caller never has to guard
        the segment edges.
        """
        if self.straight_length <= sigma <= self.length:
            return self.curvature, 0.0
        return 0.0, 0.0

    # -- material utilisation ----------------------------------------------

    def bending_strain(self, kappa: float) -> float:
        """Peak surface strain from a total curvature ``kappa`` [-].

        ``eps = kappa * r_outer``.  For reference, the 50 mm pre-curvature
        alone puts the 3.6 mm tube at 3.6 % -- far past the ~1 % where binary
        nitinol leaves its linear austenitic branch and onto the superelastic
        plateau.  That is normal for concentric-tube robots and is why the
        model's ``youngs_modulus`` is an *effective* modulus to be fitted
        against measured shapes, not a handbook value.
        """
        return abs(kappa) * 0.5 * self.outer_diameter

    def shear_stress(self, u_z: float) -> float:
        """Peak torsional shear stress from a twist rate ``u_z`` [Pa].

        ``tau = G * u_z * r_outer``.  Superelastic nitinol shears
        inelastically somewhere around 200-300 MPa, so a configuration that
        exceeds this is outside what a linear-elastic model describes.
        """
        return self.shear_modulus * abs(u_z) * 0.5 * self.outer_diameter

    def torsional_moment(self, u_z: float) -> float:
        """Torque carried at twist rate ``u_z`` [N.m]."""
        return self.torsional_stiffness * u_z

    def summary(self) -> str:
        mm = 1e3
        return (
            f"{self.name:<6s} OD {self.outer_diameter * mm:5.2f} mm  "
            f"ID {self.inner_diameter * mm:5.2f} mm  "
            f"wall {self.wall_thickness * mm:4.2f} mm  "
            f"L {self.length * mm:6.1f} mm  "
            f"curved {self.curved_length * mm:5.1f} mm "
            f"(R {1e3 / self.curvature:4.1f} mm, {np.degrees(self.arc_angle):5.1f} deg)  "
            f"EI {self.bending_stiffness:8.4f}  GJ {self.torsional_stiffness:8.4f} N.m^2"
        )
