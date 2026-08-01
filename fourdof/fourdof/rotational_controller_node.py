#!/usr/bin/env python3
"""ROS1 node wrapping RotationalController: exposes ITR/OTR over topics/services.

Params:
  ~config_file: path to a YAML file (see config/rotational_config.yaml)
"""

import os

import rospy
import yaml

from fourdof.joint_math import RotationalJointConfig
from fourdof.msg import JointCommand, JointState
from fourdof.rotational_controller import RotationalController
from fourdof.srv import SetHome, SetHomeResponse, SetTorque, SetTorqueResponse


def _load_joints(joints_cfg: dict) -> dict:
    joints = {}
    for name, cfg in joints_cfg.items():
        joints[name] = RotationalJointConfig(
            motor_id=cfg['motor_id'],
            gear_ratio=cfg.get('gear_ratio', 40.0),
            default_speed_deg_s=cfg.get('default_speed_deg_s', 30.0),
            default_accel_deg_s2=cfg.get('default_accel_deg_s2'),
        )
    return joints


class RotationalControllerNode:

    def __init__(self):
        rospy.init_node('rotational_controller_node')

        config_file = rospy.get_param('~config_file', '')
        if not config_file or not os.path.exists(config_file):
            raise RuntimeError(f"'~config_file' param not set or file not found: '{config_file}'")

        with open(config_file, 'r') as f:
            config = yaml.safe_load(f)

        device_name = config.get('device_name', '/dev/ttyUSB0')
        baud_rate = config.get('baud_rate', 57600)
        joints = _load_joints(config.get('joints', {}))
        if not joints:
            raise RuntimeError(f"No joints configured in '{config_file}'")

        rospy.loginfo(f"Rotational Controller connecting to {device_name} @ {baud_rate} baud "
                      f"(joints: {list(joints)})")
        self.controller = RotationalController(device_name, baud_rate, joints)
        rospy.on_shutdown(self.controller.shutdown)

        rospy.Subscriber('joint_command', JointCommand, self._on_command, queue_size=10)
        self.state_pub = rospy.Publisher('joint_state', JointState, queue_size=10)
        rospy.Service('set_home', SetHome, self._on_set_home)
        rospy.Service('set_torque', SetTorque, self._on_set_torque)

        rospy.loginfo('rotational_controller_node ready.')

    def _on_command(self, msg: JointCommand):
        try:
            velocity = msg.velocity or None
            acceleration = msg.acceleration or None
            if msg.relative:
                self.controller.move_by(msg.joint_name, msg.position, velocity, acceleration)
            else:
                self.controller.move_to(msg.joint_name, msg.position, velocity, acceleration)
        except ValueError as e:
            rospy.logwarn(f'Rejected joint command: {e}')
        except Exception as e:
            rospy.logerr(f'Failed to execute joint command: {e}')

    def _on_set_home(self, req):
        try:
            self.controller.set_home(req.joint_name)
            return SetHomeResponse(success=True, message=f"Zeroed '{req.joint_name}'")
        except Exception as e:
            return SetHomeResponse(success=False, message=str(e))

    def _on_set_torque(self, req):
        try:
            self.controller.set_torque(req.joint_name, req.enable)
            state = 'enabled' if req.enable else 'disabled'
            return SetTorqueResponse(success=True, message=f"Torque {state} for '{req.joint_name}'")
        except Exception as e:
            return SetTorqueResponse(success=False, message=str(e))

    def spin(self):
        rate = rospy.Rate(10)
        while not rospy.is_shutdown():
            for joint_name in self.controller.joint_names:
                try:
                    state = self.controller.get_state(joint_name)
                    msg = JointState()
                    msg.joint_name = joint_name
                    msg.position = state['position_deg']
                    msg.velocity = state['velocity_deg_s']
                    msg.present_position_pulse = state['present_position_pulse']
                    msg.torque_enabled = state['torque_enabled']
                    self.state_pub.publish(msg)
                except Exception as e:
                    rospy.logdebug(f"Failed to read state for '{joint_name}': {e}")
            rate.sleep()


def main():
    try:
        node = RotationalControllerNode()
        node.spin()
    except rospy.ROSInterruptException:
        pass


if __name__ == '__main__':
    main()
