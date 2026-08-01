"""TranslationalController: the OOP control module for ITT and OTT.

Each joint is an XL430-W250 motor driving its tube through a lead screw, run
in Extended Position Control Mode since a real travel range spans many motor
turns. Extended Position Mode does not enforce the firmware's position
limits, so each joint's [min_mm, max_mm] travel range is checked in software
before every move. No rospy dependency — usable directly from a script, a
test, or a thin ROS1 node wrapper.
"""

from fourdof.base_controller import BaseJointController
from fourdof.joint_math import (
    translational_accel_to_profile_accel,
    translational_mm_to_pulse,
    translational_profile_velocity_to_speed,
    translational_pulse_to_mm,
    translational_speed_to_profile_velocity,
    validate_translational_range,
)


class TranslationalController(BaseJointController):

    def get_position_mm(self, joint_name: str) -> float:
        cfg = self._config(joint_name)
        pulse = self._motor(joint_name).getPresentPosition()
        return translational_pulse_to_mm(pulse, cfg.pitch_mm)

    def get_state(self, joint_name: str) -> dict:
        """Single-pass feedback read: position (mm), velocity (mm/s), raw pulse, torque."""
        cfg = self._config(joint_name)
        motor = self._motor(joint_name)
        pulse = motor.getPresentPosition()
        velocity_pulse = motor.getPresentVelocity()
        return {
            'position_mm': translational_pulse_to_mm(pulse, cfg.pitch_mm),
            'velocity_mm_s': translational_profile_velocity_to_speed(velocity_pulse, cfg.pitch_mm),
            'present_position_pulse': pulse,
            'torque_enabled': bool(motor.torque_status),
        }

    def move_to(self, joint_name: str, position_mm: float,
                speed_mm_s: float = None, accel_mm_s2: float = None) -> None:
        cfg = self._config(joint_name)
        validate_translational_range(position_mm, cfg)
        motor = self._motor(joint_name)

        speed = speed_mm_s if speed_mm_s is not None else cfg.default_speed_mm_s
        accel = accel_mm_s2 if accel_mm_s2 is not None else cfg.default_accel_mm_s2

        velocity_pulse = translational_speed_to_profile_velocity(speed, cfg.pitch_mm) if speed else None
        accel_pulse = translational_accel_to_profile_accel(accel, cfg.pitch_mm) if accel else None
        self._write_profile(motor, velocity_pulse, accel_pulse)

        motor.setGoalPosition(translational_mm_to_pulse(position_mm, cfg.pitch_mm))

    def move_by(self, joint_name: str, delta_mm: float,
                speed_mm_s: float = None, accel_mm_s2: float = None) -> None:
        current_mm = self.get_position_mm(joint_name)
        self.move_to(joint_name, current_mm + delta_mm, speed_mm_s, accel_mm_s2)

    def stop(self, joint_name: str) -> None:
        """Hold at the current position (soft stop, mirrors the 1DOF CLI's stop())."""
        self.move_to(joint_name, self.get_position_mm(joint_name))
