"""
ctr -- geometrically exact Cosserat-rod model of the CT-SDR.

Implements Rucker, Jones & Webster III, "A Geometrically Exact Model for
Externally Loaded Concentric-Tube Continuum Robots", IEEE Transactions on
Robotics 26(5):769-780, 2010, specialised to this repo's two-nitinol-tube
steerable drilling robot.

    from ctr import CTSDR, Joints

    robot = CTSDR.from_yaml()
    sol = robot.solve(Joints(ott=35, itt=35, otr=0, itr=180))
    print(sol.tip_position * 1e3, "mm")
"""

from .friction import (
    MU_NITINOL_NITINOL, MU_STEEL_NITINOL, NO_FRICTION, FrictionModel,
)
from .loads import NO_LOAD, ExternalLoad
from .model import CosseratModel, ShootingError, Solution, hat, rot_z
from .robot import DEFAULT_CONFIG, INNER, OUTER, CTSDR, Joints
from .tube import Tube

__all__ = [
    "CTSDR",
    "CosseratModel",
    "ExternalLoad",
    "FrictionModel",
    "NO_FRICTION",
    "MU_NITINOL_NITINOL",
    "MU_STEEL_NITINOL",
    "Joints",
    "NO_LOAD",
    "ShootingError",
    "Solution",
    "Tube",
    "DEFAULT_CONFIG",
    "OUTER",
    "INNER",
    "hat",
    "rot_z",
]
