# FourDoF

Control stack for a 4-DOF robot built from two nested tubes:

| Joint | Motor | Motion | Mechanism |
|---|---|---|---|
| ITR (Inner Tube Rotation) | XL330-M288 | degrees, continuous | 1:40 worm gear |
| OTR (Outer Tube Rotation) | XL330-M288 | degrees, continuous | 1:40 worm gear |
| ITT (Inner Tube Translation) | XC430-W240 | mm, bounded range | 2.54 mm/rev lead screw |
| OTT (Outer Tube Translation) | XC430-W240 | mm, bounded range | 2.54 mm/rev lead screw |

ITR/OTR and ITT/OTT are wired to **two separate serial buses** and run as two
independent ROS1 nodes — the **Rotational Controller** and the
**Translational Controller**.

## Layout

```
FourDoF/
├── dxl_control_4dof_cli.py       # interactive CLI for all 4 joints (run this)
├── dxl_control_length_1dof.py    # older single-DOF reference example (unrelated to the 4-DOF stack)
└── fourdof/                      # the ROS1 catkin package
    ├── fourdof/                  # control module (rospy-independent classes + thin ROS1 node wrappers)
    ├── libs/dynamixel_easy_sdk/  # vendored ROBOTIS Dynamixel SDK wrapper (Connector/Motor/OperatingMode)
    ├── msg/, srv/                # JointCommand/JointState messages, SetHome/SetTorque services
    ├── config/                  # rotational_config.yaml, translational_config.yaml
    ├── launch/four_dof.launch   # starts both controller nodes
    └── tests/                   # test_joint_math.py (unit tests), scan_motors.py (bus discovery)
```

See [fourdof/README.md](fourdof/README.md) for package-level details.

## Prerequisites

- Ubuntu with **ROS1 Noetic** installed and a catkin workspace.
- Python 3, `pip3`.
- The `dynamixel_easy_sdk` layer needs `dynamixel_sdk` and `pyserial`, which
  aren't always resolvable via `rosdep`:
  ```bash
  pip3 install dynamixel-sdk pyserial
  ```
- If you launch through the repo's virtual environment, that interpreter also
   needs `rospkg` because `rospy` imports it at startup:
   ```bash
   pip install rospkg
   ```
- Serial port access:
  ```bash
  sudo usermod -aG dialout $USER   # log out/in (or reboot) for this to take effect
  ```

## Setup

1. **Add the package to a catkin workspace.** `fourdof/` is the catkin
   package itself, so it needs to live under a workspace's `src/`. Either
   symlink it in or clone directly there:
   ```bash
   mkdir -p ~/fourdof_catkin_ws/src
   ln -s /path/to/FourDoF/fourdof ~/fourdof_catkin_ws/src/fourdof
   ```

2. **Build.**
   ```bash
   cd ~/fourdof_catkin_ws
   rosdep install --from-paths src --ignore-src -r -y
   catkin_make        # or: catkin build
   source devel/setup.zsh   # use setup.bash only from bash
   ```

3. **Identify the serial ports.** Plug in both controllers and check which
   `/dev/ttyUSB*` (or `/dev/ttyACM*`) each one enumerates as, e.g.:
   ```bash
   dmesg | grep -i tty
   ```
   Edit `device_name` in `fourdof/config/rotational_config.yaml` and
   `fourdof/config/translational_config.yaml` to match.

4. **Confirm motor IDs.** The configs assume ID 1/2 on each bus. Verify with
   the bus-scan utility (edit the `port`/`baud` constants at the top of the
   script for whichever bus you're checking):
   ```bash
   source ~/fourdof_catkin_ws/devel/setup.zsh 
   python3 ~/fourdof_catkin_ws/src/fourdof/tests/scan_motors.py
   ```
   Update `motor_id` in the config files if they don't match.

5. **Set the translational travel limits.** `fourdof/config/translational_config.yaml`
   ships with `min_mm`/`max_mm` set to `null` for ITT and OTT — the node
   refuses to start until these are set to the real hardware travel limits
   (Extended Position Mode doesn't enforce firmware position limits, so this
   is the only safety guard against a lead screw running into its hardstop).

## Running

Start both controller nodes:
```bash
source ~/fourdof_catkin_ws/devel/setup.zsh 
roslaunch fourdof four_dof.launch
```

Expected Output:
```
➜  fourdof git:(main) ✗ roslaunch fourdof four_dof.launch
... logging to /home/daniyal/.ros/log/39147f32-6e68-11f1-906a-4b41e9977cdc/roslaunch-tina-1156137.log
Checking log directory for disk usage. This may take a while.
Press Ctrl-C to interrupt

started roslaunch server http://172.31.1.146:46321/

SUMMARY
========

PARAMETERS
 * /rosdistro: noetic
 * /rosversion: 1.17.0
 * /rotational_controller/rotational_controller/config_file: /home/daniyal/fou...
 * /translational_controller/translational_controller/config_file: /home/daniyal/fou...

NODES
  /rotational_controller/
    rotational_controller (fourdof/rotational_controller_node.py)
  /translational_controller/
    translational_controller (fourdof/translational_controller_node.py)

ROS_MASTER_URI=http://localhost:11311

process[rotational_controller/rotational_controller-1]: started with pid [1156151]
process[translational_controller/translational_controller-2]: started with pid [1156152]
[INFO] [1786559981.611970]: Rotational Controller connecting to /dev/ttyUSB2 @ 57600 baud (joints: ['ITR', 'OTR'])
[INFO] [1786559981.617089]: Translational Controller connecting to /dev/ttyUSB3 @ 57600 baud (joints: ['ITT', 'OTT'])
[INFO] [1786559981.681574]: rotational_controller_node ready.
[INFO] [1786559981.693663]: translational_controller_node ready.
```

In another terminal (with the same workspace sourced), run the CLI from the
project root:
```bash
source ~/fourdof_catkin_ws/devel/setup.zsh 
python3 dxl_control_4dof_cli.py
```

Example commands inside the CLI (`h` for the full list):
```
itr a 90 20      # rotate ITR to 90 deg at 20 deg/s
ott r -5         # move OTT back 5 mm at its own default speed
p                # print all 4 joint positions
speed itt 8      # change ITT's own default speed to 8 mm/s
home itr         # zero ITR's current position
torque otr off   # disable torque on OTR
stop all         # hold all 4 joints at their current position
```

## Verifying without hardware

The conversion math (gear ratio / lead-screw pitch, degree/mm ⇄ pulse, and
translational range clamping) is unit tested with no ROS/hardware
dependency:
```bash
cd ~/fourdof_catkin_ws/src/fourdof/tests
python3 -m unittest test_joint_math -v
```
# FourDoF
