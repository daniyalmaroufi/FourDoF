"""RotationalController: the OOP control module for ITR and OTR.

Each joint is an XL330-M288 motor driving its tube through a 1:40 worm gear,
run in Extended Position Control Mode so continuous multi-turn tube rotation
is tracked without wraparound. No rospy dependency — usable directly from a
script, a test, or a thin ROS1 node wrapper.
"""

from fourdof.base_controller import BaseJointController
from fourdof.joint_math import (
    rotational_accel_to_profile_accel,
    rotational_deg_to_pulse,
    rotational_profile_velocity_to_speed,
    rotational_pulse_to_deg,
    rotational_speed_to_profile_velocity,
)


class RotationalController(BaseJointController):

    def get_position_deg(self, joint_name: str) -> float:
        cfg = self._config(joint_name)
        pulse = self._motor(joint_name).getPresentPosition()
        return rotational_pulse_to_deg(pulse, cfg.gear_ratio)

    def get_state(self, joint_name: str) -> dict:
        """Single-pass feedback read: position (deg), velocity (deg/s), raw pulse, torque."""
        cfg = self._config(joint_name)
        motor = self._motor(joint_name)
        pulse = motor.getPresentPosition()
        velocity_pulse = motor.getPresentVelocity()
        return {
            'position_deg': rotational_pulse_to_deg(pulse, cfg.gear_ratio),
            'velocity_deg_s': rotational_profile_velocity_to_speed(velocity_pulse, cfg.gear_ratio),
            'present_position_pulse': pulse,
            'torque_enabled': bool(motor.torque_status),
        }

    def move_to(self, joint_name: str, position_deg: float,
                speed_deg_s: float = None, accel_deg_s2: float = None) -> None:
        cfg = self._config(joint_name)
        motor = self._motor(joint_name)

        speed = speed_deg_s if speed_deg_s is not None else cfg.default_speed_deg_s
        accel = accel_deg_s2 if accel_deg_s2 is not None else cfg.default_accel_deg_s2

        velocity_pulse = rotational_speed_to_profile_velocity(speed, cfg.gear_ratio) if speed else None
        accel_pulse = rotational_accel_to_profile_accel(accel, cfg.gear_ratio) if accel else None
        self._write_profile(motor, velocity_pulse, accel_pulse)

        motor.setGoalPosition(rotational_deg_to_pulse(position_deg, cfg.gear_ratio))

    def move_by(self, joint_name: str, delta_deg: float,
                speed_deg_s: float = None, accel_deg_s2: float = None) -> None:
        current_deg = self.get_position_deg(joint_name)
        self.move_to(joint_name, current_deg + delta_deg, speed_deg_s, accel_deg_s2)

    def stop(self, joint_name: str) -> None:
        """Hold at the current position (soft stop, mirrors the 1DOF CLI's stop())."""
        self.move_to(joint_name, self.get_position_deg(joint_name))
