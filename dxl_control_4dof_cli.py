#!/usr/bin/env python
"""
CLI tool for controlling the 4-DOF robot through ROS.

Talks to two separate ROS1 nodes, each on its own serial bus:
  - rotational_controller:    ITR / OTR tube rotation, in degrees (continuous)
  - translational_controller: ITT / OTT tube translation, in mm (bounded range)
"""

import readline  # noqa: F401  (enables arrow-key command history in input())
import sys

import rospy

from fourdof.msg import JointCommand, JointState
from fourdof.srv import SetHome, SetTorque

# Joint name -> (controller namespace, physical unit label)
JOINT_INFO = {
    'ITR': ('rotational_controller', 'deg'),
    'OTR': ('rotational_controller', 'deg'),
    'ITT': ('translational_controller', 'mm'),
    'OTT': ('translational_controller', 'mm'),
}

DEFAULT_SPEEDS = {
    'ITR': 30.0,
    'OTR': 30.0,
    'ITT': 5.0,
    'OTT': 5.0,
}


class FourDoFControlCLI(object):
    """Interactive CLI for controlling all 4 joints (ITR/OTR/ITT/OTT) via ROS."""

    def __init__(self):
        rospy.init_node('four_dof_control_cli', anonymous=True)

        self.states = {}  # joint_name -> latest JointState
        # Each joint tracks its OWN default speed, independently settable
        # (mirrors dxl_control_length_1dof.py's self.default_speed, per joint).
        self.default_speeds = dict(DEFAULT_SPEEDS)

        namespaces = sorted(set(ns for ns, _ in JOINT_INFO.values()))
        self.command_pubs = {
            ns: rospy.Publisher(f'/{ns}/joint_command', JointCommand, queue_size=10)
            for ns in namespaces
        }
        self.state_subs = [
            rospy.Subscriber(f'/{ns}/joint_state', JointState, self._state_callback)
            for ns in namespaces
        ]
        self.set_home_proxies = {ns: rospy.ServiceProxy(f'/{ns}/set_home', SetHome) for ns in namespaces}
        self.set_torque_proxies = {ns: rospy.ServiceProxy(f'/{ns}/set_torque', SetTorque) for ns in namespaces}

        # Wait a bit for ROS connections to establish
        rospy.sleep(0.5)

    def _state_callback(self, msg):
        self.states[msg.joint_name] = msg

    @staticmethod
    def _parse_joint(token):
        joint = token.upper()
        if joint not in JOINT_INFO:
            print(f"Error: Unknown joint '{token}'. Expected one of {list(JOINT_INFO)}")
            return None
        return joint

    # --- command senders -----------------------------------------------

    def _send(self, joint, position, relative, speed):
        ns, unit = JOINT_INFO[joint]
        msg = JointCommand()
        msg.joint_name = joint
        msg.relative = relative
        msg.position = position
        msg.velocity = speed
        msg.acceleration = 0.0
        self.command_pubs[ns].publish(msg)
        verb = 'Relative' if relative else 'Absolute'
        print(f'{verb} move: {joint} -> {position:+.2f} {unit} @ {speed:.2f} {unit}/s')

    def move_absolute(self, joint, value, speed=None):
        self._send(joint, value, relative=False, speed=speed if speed is not None else self.default_speeds[joint])

    def move_relative(self, joint, delta, speed=None):
        self._send(joint, delta, relative=True, speed=speed if speed is not None else self.default_speeds[joint])

    def print_position(self, joint=None):
        for j in ([joint] if joint else list(JOINT_INFO)):
            state = self.states.get(j)
            unit = JOINT_INFO[j][1]
            if state is None:
                print(f'  {j}: position unknown (no feedback received yet)')
            else:
                torque = 'ON' if state.torque_enabled else 'OFF'
                print(f'  {j}: {state.position:+.2f} {unit}  (vel {state.velocity:+.2f} {unit}/s, torque {torque})')

    def set_default_speed(self, joint, speed):
        self.default_speeds[joint] = speed
        unit = JOINT_INFO[joint][1]
        print(f'{joint} default speed set to {speed:.2f} {unit}/s')

    def stop(self, joint):
        state = self.states.get(joint)
        if state is None:
            print(f'Error: No position feedback for {joint} yet, cannot stop.')
            return
        self._send(joint, state.position, relative=False, speed=self.default_speeds[joint])
        print(f'Stop command sent for {joint} (holding at {state.position:.2f} {JOINT_INFO[joint][1]})')

    def home(self, joint):
        ns = JOINT_INFO[joint][0]
        try:
            rospy.wait_for_service(f'/{ns}/set_home', timeout=2.0)
            resp = self.set_home_proxies[ns](joint_name=joint)
            print(resp.message or ('OK' if resp.success else 'Failed'))
        except (rospy.ROSException, rospy.ServiceException) as e:
            print(f'Error calling set_home for {joint}: {e}')

    def torque(self, joint, enable):
        ns = JOINT_INFO[joint][0]
        try:
            rospy.wait_for_service(f'/{ns}/set_torque', timeout=2.0)
            resp = self.set_torque_proxies[ns](joint_name=joint, enable=enable)
            print(resp.message or ('OK' if resp.success else 'Failed'))
        except (rospy.ROSException, rospy.ServiceException) as e:
            print(f'Error calling set_torque for {joint}: {e}')

    # --- CLI plumbing -----------------------------------------------------

    def print_help(self):
        help_text = """
╔════════════════════════════════════════════════════════════════════════════╗
║                    4-DOF ROBOT CONTROL CLI                                  ║
╠════════════════════════════════════════════════════════════════════════════╣
║                                                                            ║
║  JOINTS:                                                                   ║
║  ───────                                                                   ║
║   itr, otr   - Inner/Outer Tube Rotation   (degrees, continuous)          ║
║   itt, ott   - Inner/Outer Tube Translation (mm, bounded range)           ║
║                                                                            ║
║  MOVEMENT COMMANDS:                                                        ║
║  ─────────────────                                                         ║
║   <joint> a <value> [speed]  - Absolute move to <value>                   ║
║   <joint> r <delta> [speed]  - Relative move by <delta>                   ║
║                                                                            ║
║      Example:  itr a 90 20    - Rotate ITR to 90 deg at 20 deg/s          ║
║      Example:  ott r -5       - Move OTT back 5 mm at its default speed  ║
║                                                                            ║
╠════════════════════════════════════════════════════════════════════════════╣
║  UTILITY COMMANDS:                                                         ║
║  ──────────────────                                                        ║
║   p [joint]           - Print position (all 4 joints if omitted)          ║
║   speed <joint> <val>  - Set that joint's own default speed               ║
║   home <joint>         - Zero the joint's current position                ║
║   torque <joint> on|off - Enable/disable torque for that joint            ║
║   s / stop <joint>|all - Hold at current position                        ║
║   h or help            - Show this help                                   ║
║   q or quit or exit    - Exit CLI                                         ║
║                                                                            ║
╠════════════════════════════════════════════════════════════════════════════╣
║  HARDWARE:                                                                 ║
║  ─────────                                                                 ║
║   ITR/OTR: XL330-M288, 1:40 worm gear (40 motor rev = 1 tube rev)         ║
║   ITT/OTT: XC430-W240, 2.54 mm/rev lead screw, joint-specific mm range    ║
╚════════════════════════════════════════════════════════════════════════════╝

⚠  Make sure rotational_controller_node and translational_controller_node are running!
"""
        print(help_text)

    def run(self):
        """Main interactive loop."""
        print('\n' + '=' * 70)
        print('  4-DOF ROBOT CONTROL CLI')
        print('=' * 70)
        print("\nType 'h' for help, 'q' to quit\n")

        print('Waiting for joint state feedback...', end='', flush=True)
        timeout = rospy.Time.now() + rospy.Duration(3.0)
        while len(self.states) < len(JOINT_INFO) and rospy.Time.now() < timeout and not rospy.is_shutdown():
            rospy.sleep(0.1)
            print('.', end='', flush=True)
        print()

        missing = [j for j in JOINT_INFO if j not in self.states]
        if not missing:
            print('✓ Connected! All 4 joints reporting.\n')
        else:
            print(f'⚠ Warning: no feedback yet from {missing}. '
                  'Make sure both controller nodes are running!\n')

        while not rospy.is_shutdown():
            try:
                user_input = input('4dof> ').strip()
                if not user_input:
                    continue

                parts = user_input.split()
                cmd = parts[0].lower()

                if cmd in ('q', 'quit', 'exit'):
                    print('Exiting...')
                    break

                elif cmd in ('h', 'help'):
                    self.print_help()

                elif cmd == 'p':
                    if len(parts) > 1:
                        joint = self._parse_joint(parts[1])
                        if joint is None:
                            continue
                        self.print_position(joint)
                    else:
                        self.print_position()

                elif cmd == 'speed':
                    if len(parts) < 3:
                        print('Error: Missing arguments. Usage: speed <joint> <value>')
                        continue
                    joint = self._parse_joint(parts[1])
                    if joint is None:
                        continue
                    try:
                        value = float(parts[2])
                        if value <= 0:
                            print('Error: Speed must be positive')
                        else:
                            self.set_default_speed(joint, value)
                    except ValueError:
                        print('Error: Speed must be a number')

                elif cmd == 'home':
                    if len(parts) < 2:
                        print('Error: Missing joint. Usage: home <joint>')
                        continue
                    joint = self._parse_joint(parts[1])
                    if joint is not None:
                        self.home(joint)

                elif cmd == 'torque':
                    if len(parts) < 3 or parts[2].lower() not in ('on', 'off'):
                        print('Error: Usage: torque <joint> on|off')
                        continue
                    joint = self._parse_joint(parts[1])
                    if joint is not None:
                        self.torque(joint, parts[2].lower() == 'on')

                elif cmd in ('s', 'stop'):
                    if len(parts) < 2:
                        print('Error: Missing joint. Usage: stop <joint>|all')
                        continue
                    if parts[1].lower() == 'all':
                        for j in JOINT_INFO:
                            self.stop(j)
                    else:
                        joint = self._parse_joint(parts[1])
                        if joint is not None:
                            self.stop(joint)

                elif cmd.upper() in JOINT_INFO:
                    joint = cmd.upper()
                    if len(parts) < 3:
                        print(f'Error: Missing arguments. Usage: {cmd} a|r <value> [speed]')
                        continue
                    move_cmd = parts[1].lower()
                    try:
                        value = float(parts[2])
                        speed = float(parts[3]) if len(parts) > 3 else None
                    except ValueError:
                        print('Error: Invalid number format')
                        continue

                    if move_cmd == 'a':
                        self.move_absolute(joint, value, speed)
                    elif move_cmd == 'r':
                        self.move_relative(joint, value, speed)
                    else:
                        print(f"Error: Unknown move type '{move_cmd}'. Use 'a' (absolute) or 'r' (relative).")

                else:
                    print(f"Unknown command: '{cmd}'. Type 'h' for help.")

            except KeyboardInterrupt:
                print("\n\nInterrupted. Type 'q' to quit or continue entering commands.")
            except EOFError:
                print('\nExiting...')
                break
            except Exception as e:
                print(f'Error: {e}')

        print('CLI terminated.')


def main():
    try:
        cli = FourDoFControlCLI()
        cli.run()
    except rospy.ROSInterruptException:
        print('ROS interrupted')
    except Exception as e:
        print(f'Error: {e}')
        sys.exit(1)


if __name__ == '__main__':
    main()
