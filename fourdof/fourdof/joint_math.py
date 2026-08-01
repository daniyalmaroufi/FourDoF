"""Unit conversion helpers for the 4-DOF robot's joints.

Rotational joints (ITR, OTR) ride a 1:40 worm gear between the XL330 motor
and the tube. Translational joints (ITT, OTT) ride a lead screw between the
XC430 motor and the tube, with a configurable travel range in mm.

Everything here is plain Python with no rospy / dynamixel_sdk dependency, so
it can be unit tested without ROS or hardware attached.
"""

from dataclasses import dataclass
from typing import Optional

PULSES_PER_REV = 4096
DEG_PER_PULSE = 360.0 / PULSES_PER_REV  # 0.087890625 deg/pulse at the motor shaft

# Dynamixel X-series (XL330/XC430) Profile Velocity / Profile Acceleration LSB
# units, per the ROBOTIS control table (0.229 rev/min and 214.577 rev/min^2
# per unit, respectively — same for both motor models).
PROFILE_VELOCITY_REV_PER_MIN = 0.229
PROFILE_ACCEL_REV_PER_MIN2 = 214.577


@dataclass
class RotationalJointConfig:
    motor_id: int
    gear_ratio: float
    default_speed_deg_s: float = 30.0
    default_accel_deg_s2: Optional[float] = None


@dataclass
class TranslationalJointConfig:
    motor_id: int
    pitch_mm: float
    min_mm: float
    max_mm: float
    default_speed_mm_s: float = 5.0
    default_accel_mm_s2: Optional[float] = None


# --- Rotational (tube-space, accounting for the worm gear ratio) -----------

def rotational_deg_to_pulse(deg: float, gear_ratio: float) -> int:
    """Tube degrees -> motor Goal Position pulses."""
    motor_deg = deg * gear_ratio
    return round(motor_deg / DEG_PER_PULSE)


def rotational_pulse_to_deg(pulse: int, gear_ratio: float) -> float:
    """Motor Present Position pulses -> tube degrees."""
    motor_deg = pulse * DEG_PER_PULSE
    return motor_deg / gear_ratio


def rotational_speed_to_profile_velocity(deg_s: float, gear_ratio: float) -> int:
    """Tube deg/s -> Profile Velocity register value."""
    motor_rev_per_min = (deg_s * gear_ratio) / 6.0
    return max(0, round(motor_rev_per_min / PROFILE_VELOCITY_REV_PER_MIN))


def rotational_profile_velocity_to_speed(value: int, gear_ratio: float) -> float:
    """Profile Velocity register value -> tube deg/s."""
    motor_rev_per_min = value * PROFILE_VELOCITY_REV_PER_MIN
    return (motor_rev_per_min * 6.0) / gear_ratio


def rotational_accel_to_profile_accel(deg_s2: float, gear_ratio: float) -> int:
    """Tube deg/s^2 -> Profile Acceleration register value."""
    motor_rev_per_min2 = (deg_s2 * gear_ratio) * 10.0
    return max(0, round(motor_rev_per_min2 / PROFILE_ACCEL_REV_PER_MIN2))


# --- Translational (lead-screw mm-space) ------------------------------------

def translational_mm_to_pulse(mm: float, pitch_mm: float) -> int:
    """Lead-screw mm -> motor Goal Position pulses."""
    return round((mm / pitch_mm) * PULSES_PER_REV)


def translational_pulse_to_mm(pulse: int, pitch_mm: float) -> float:
    """Motor Present Position pulses -> lead-screw mm."""
    return (pulse / PULSES_PER_REV) * pitch_mm


def translational_speed_to_profile_velocity(mm_s: float, pitch_mm: float) -> int:
    """Lead-screw mm/s -> Profile Velocity register value."""
    motor_rev_per_min = (mm_s / pitch_mm) * 60.0
    return max(0, round(motor_rev_per_min / PROFILE_VELOCITY_REV_PER_MIN))


def translational_profile_velocity_to_speed(value: int, pitch_mm: float) -> float:
    """Profile Velocity register value -> lead-screw mm/s."""
    motor_rev_per_min = value * PROFILE_VELOCITY_REV_PER_MIN
    return (motor_rev_per_min / 60.0) * pitch_mm


def translational_accel_to_profile_accel(mm_s2: float, pitch_mm: float) -> int:
    """Lead-screw mm/s^2 -> Profile Acceleration register value."""
    motor_rev_per_min2 = (mm_s2 / pitch_mm) * 3600.0
    return max(0, round(motor_rev_per_min2 / PROFILE_ACCEL_REV_PER_MIN2))


def validate_translational_range(mm: float, cfg: TranslationalJointConfig) -> None:
    """Raise ValueError if `mm` falls outside the joint's configured travel range.

    Extended Position Control Mode (needed here for multi-turn lead-screw
    travel) does not enforce the firmware's Min/Max Position Limit, so this
    check is the only thing standing between a bad command and the hardstop.
    """
    if not (cfg.min_mm <= mm <= cfg.max_mm):
        raise ValueError(
            f'Target position {mm:.2f} mm is outside the allowed range '
            f'[{cfg.min_mm:.2f}, {cfg.max_mm:.2f}] mm'
        )
