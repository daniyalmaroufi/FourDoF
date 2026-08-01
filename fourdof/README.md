# fourdof

ROS1 (Noetic) package for the 4-DOF robot: 2x XL330-M288 rotational motors
(ITR, OTR) and 2x XC430-W240 translational motors (ITT, OTT), driven over the
Dynamixel SDK via the vendored `dynamixel_easy_sdk` wrapper.

## Architecture

The two joint pairs each have their own controller, on their own serial bus,
run as independent ROS1 nodes:

- **Rotational Controller** (`rotational_controller_node.py`) — ITR / OTR tube
  rotation, in degrees, through a 1:40 worm gear. Continuous, no travel limit.
- **Translational Controller** (`translational_controller_node.py`) — ITT / OTT
  tube translation, in mm, through a 2.54 mm/rev lead screw. Each joint has its
  own configured `[min_mm, max_mm]` travel range (see
  `config/translational_config.yaml` — these **must** be set to the real
  hardware limits before running; the node refuses to start otherwise).

Both nodes wrap a plain-Python, rospy-independent control module
(`fourdof/rotational_controller.py`, `fourdof/translational_controller.py`,
built on `fourdof/base_controller.py` + `fourdof/joint_math.py`) on top of the
`dynamixel_easy_sdk` layer in `libs/`.

## Dependencies not covered by rosdep

`libs/dynamixel_easy_sdk` imports `dynamixel_sdk` and `pyserial`. Install them
with pip if `rosdep install` doesn't resolve them on your system:
```bash
pip install dynamixel-sdk pyserial
```

## Running

```bash
roslaunch fourdof four_dof.launch
```
This starts both controller nodes (see `config/rotational_config.yaml` and
`config/translational_config.yaml` for device ports, motor IDs, gear
ratio/pitch, and per-joint default speeds). Run `tests/scan_motors.py` first
to confirm the actual motor IDs on each bus before trusting the config
defaults.

Then, from the project root, run the interactive CLI:
```bash
python3 dxl_control_4dof_cli.py
```
It drives all 4 joints with commands like `itr a 90 20`, `ott r -5`,
`p`, `home itt`, `torque otr off`. Run `h` inside the CLI for the full list.

## Notes
- Both nodes run their motors in **Extended Position Control Mode** (needed
  for continuous rotation and multi-turn lead-screw travel).
- Ensure your user has permissions to access the serial port (`sudo usermod -aG dialout $USER`).
