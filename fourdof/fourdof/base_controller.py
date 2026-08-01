"""Shared, rospy-independent base class for the joint controllers.

Handles Connector/Motor setup (Extended Position Control Mode + torque
enable) and the generic torque/home/shutdown operations that are identical
between RotationalController and TranslationalController. Subclasses add the
physical-unit <-> pulse conversion and public move_to/move_by API for their
own unit (degrees or mm).
"""

from dynamixel_easy_sdk import Connector, OperatingMode
from dynamixel_easy_sdk.data_types import toSignedInt


class BaseJointController:

    def __init__(self, device_name, baud_rate, joints):
        """joints: dict mapping joint name (e.g. "ITR") to a joint config
        dataclass (RotationalJointConfig / TranslationalJointConfig)."""
        self.connector = Connector(device_name, baud_rate)
        self._joint_configs = dict(joints)
        self._motors = {}
        for name, cfg in self._joint_configs.items():
            motor = self.connector.createMotor(cfg.motor_id)
            motor.disableTorque()
            motor.setOperatingMode(OperatingMode.EXTENDED_POSITION)
            motor.enableTorque()
            self._motors[name] = motor

    @property
    def joint_names(self):
        return list(self._joint_configs.keys())

    def _motor(self, joint_name):
        try:
            return self._motors[joint_name]
        except KeyError:
            raise ValueError(f"Unknown joint '{joint_name}'. Expected one of {self.joint_names}")

    def _config(self, joint_name):
        try:
            return self._joint_configs[joint_name]
        except KeyError:
            raise ValueError(f"Unknown joint '{joint_name}'. Expected one of {self.joint_names}")

    def is_torque_enabled(self, joint_name) -> bool:
        return bool(self._motor(joint_name).torque_status)

    def set_torque(self, joint_name, enabled: bool) -> None:
        motor = self._motor(joint_name)
        if enabled:
            motor.enableTorque()
        else:
            motor.disableTorque()

    def set_home(self, joint_name) -> None:
        """Rebase the joint's Homing Offset so the current position reads as 0.

        Present Position = Actual Position + Homing Offset, so shifting the
        offset by -current_present_position zeroes the reading without
        touching the physical position.
        """
        motor = self._motor(joint_name)
        item = motor._getControlTableItem('Homing Offset')
        current_offset = toSignedInt(motor._readData(motor.id, item.address, item.size), item.size)
        current_present_position = motor.getPresentPosition()
        new_offset = current_offset - current_present_position

        motor.disableTorque()
        try:
            motor.setHomingOffset(new_offset)
        finally:
            motor.enableTorque()

    def shutdown(self) -> None:
        for motor in self._motors.values():
            try:
                motor.disableTorque()
            except Exception:
                pass
        self.connector.closePort()

    def _write_profile(self, motor, velocity_pulse, acceleration_pulse=None) -> None:
        if velocity_pulse is not None:
            item = motor._getControlTableItem('Profile Velocity')
            motor._writeData(motor.id, item.address, item.size, velocity_pulse)
        if acceleration_pulse is not None:
            item = motor._getControlTableItem('Profile Acceleration')
            motor._writeData(motor.id, item.address, item.size, acceleration_pulse)
