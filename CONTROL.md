# Control of the 4-DoF Concentric Tube Steerable Drilling Robot

## Abstract

This document describes the control architecture of a four-degree-of-freedom (4-DoF) concentric tube steerable drilling robot (CT-SDR) composed of two telescoping, pre-curved nitinol tubes. Each tube contributes one rotational and one translational degree of freedom, giving a joint space of dimension four. The document covers the definition of the joint (DoF) space, the distinction between absolute and relative motion commands, the actuation and sensing hardware, the low-level and supervisory control layers, and the kinematic (Cosserat-rod) model used to relate joint-space commands to the backbone shape and tip pose of the robot in task space.

---

## 1. Introduction

Concentric tube robots (CTRs) are a class of continuum manipulators formed from nested, pre-curved elastic tubes that are individually actuated in rotation and translation at their proximal (base) end. Relative rotation and axial translation between tubes elastically deform the assembly, and the superposition of each tube's intrinsic curvature produces a continuously curving backbone whose shape can be varied without any distal actuation, wiring, or joints along the tube length. This property makes CTRs attractive for minimally invasive and steerable-drilling applications, where a slender, dexterous instrument must navigate through confined or curved channels.

The robot considered here consists of two tubes — an **outer tube** and an **inner tube** — arranged concentrically, with the inner tube passing through the lumen of the outer tube. Each tube is independently driven at its base in two directions: rotation about its own axis, and axial translation along that axis. This yields **four independently commandable degrees of freedom**, which constitute the joint space of the system.

The control system is organized in three conceptual layers:

1. A **kinematic/model layer** that relates the four joint variables to the physical backbone shape and tip pose of the robot, using a geometrically exact (Cosserat) rod formulation.
2. A **supervisory command layer** that expresses desired joint motion as either absolute or relative set-points and converts physical (tube-space) units into actuator (motor-space) units.
3. A **low-level actuation layer**, implemented in the servo firmware, that closes a position/velocity feedback loop around each individual joint using integrated magnetic encoders.

---

## 2. Degrees of Freedom

### 2.1 Definition

The four degrees of freedom are conventionally labeled by a two-letter tube identifier and a one-letter motion type:

| Symbol | Name | Motion type | Unit | Range character |
|:------:|------|:-----------:|:----:|:----------------:|
| **ITR** | Inner Tube Rotation | Rotation about tube axis | degrees (°) | Continuous (unbounded) |
| **OTR** | Outer Tube Rotation | Rotation about tube axis | degrees (°) | Continuous (unbounded) |
| **ITT** | Inner Tube Translation | Axial translation | millimeters (mm) | Bounded (finite stroke) |
| **OTT** | Outer Tube Translation | Axial translation | millimeters (mm) | Bounded (finite stroke) |

The full joint-space configuration of the robot is therefore the vector

$$
\mathbf{q} = \begin{bmatrix} q_{\mathrm{OTT}} & q_{\mathrm{ITT}} & q_{\mathrm{OTR}} & q_{\mathrm{ITR}} \end{bmatrix}^{\top} \in \mathbb{R}^2 \times \mathbb{S}^1 \times \mathbb{S}^1,
$$

reflecting that the two translational variables live in a bounded real interval, while the two rotational variables are angles on the circle (unbounded in the sense that the actuator itself supports continuous multi-turn rotation, even though the physical tube motion is periodic in $2\pi$).

### 2.2 Kinematic coupling between joints

Because the inner tube is mechanically nested inside — and carried by — the outer tube, a translation of the outer tube rigidly transports the inner tube's carriage along with it. Consequently the *tube-space* insertion depth of the inner tube relative to the outer tube's distal end depends on **both** $q_{\mathrm{OTT}}$ and $q_{\mathrm{ITT}}$; a pure outer-tube translation does not, by itself, change the relative overlap of the two tubes, but it does change the inner tube's translation in the fixed (world/base) frame. This coupling must be accounted for whenever joint commands are converted into the relative insertion and rotation variables used by the kinematic model (Section 7).

### 2.3 Sign convention

A consistent right-handed sign convention is required to relate a positive joint command to a physical direction of motion. For this robot, the validated convention (established by comparison of the kinematic model against external tip-tracking measurements, Section 8) is:

- **Translation** ($q_{\mathrm{OTT}}, q_{\mathrm{ITT}}$): positive values correspond to the tube advancing distally, i.e., increasing insertion out of the guiding channel.
- **Rotation** ($q_{\mathrm{OTR}}, q_{\mathrm{ITR}}$): positive commanded rotation produces a **left-handed** rotation of the tube about its own distal axis (i.e., the opposite sign from the naive right-hand-rule expectation). This inversion is a compound effect of (a) the handedness introduced by the rotary transmission between actuator and tube, and (b) the handedness of the reference frame in which tip pose is measured.

Both sign conventions are treated as calibrated constants of the physical system rather than free parameters — they must be re-verified whenever the transmission, actuator wiring, or tracking-frame definition is changed.

---

## 3. Absolute and Relative Motion

Every joint accepts motion commands in one of two modes. This distinction is orthogonal to the joint identity (ITR/OTR/ITT/OTT) — any of the four DoFs can be driven absolutely or relatively.

### 3.1 Absolute motion

An **absolute** command specifies a target position directly in the joint's physical unit (degrees for rotation, millimeters for translation), referenced to the joint's calibrated zero (Section 3.3):

$$
q_j(t^{+}) = q_j^{\text{target}}
$$

where $q_j(t^{+})$ denotes the commanded destination for joint $j \in \{\mathrm{ITR, OTR, ITT, OTT}\}$. Absolute motion is the natural mode for tasks specified in terms of a desired configuration (e.g., "set the outer tube rotation to $45°$"), and is the mode required whenever the destination must be reproducible independent of the joint's current state.

### 3.2 Relative motion

A **relative** (incremental) command specifies a signed displacement $\Delta q_j$ to be applied on top of the joint's instantaneous measured position:

$$
q_j(t^{+}) = q_j(t) + \Delta q_j
$$

where $q_j(t)$ is read back from the joint's position feedback at the time the command is issued. Relative motion is the natural mode for exploratory or corrective motions (e.g., jogging a joint by a small increment) where the operator or planner reasons in terms of a change from the present state rather than a fixed target.

Because a relative command is resolved against the instantaneously measured position rather than a previously commanded target, relative motions are **not commutative with concurrent motion** in the strict sense — a relative command issued while the joint is still moving toward a prior target will compound with whatever position has actually been reached at the instant the new command is evaluated, not with the prior target itself.

### 3.3 Zero-reference (homing) and its relation to absolute motion

Absolute motion is only meaningful relative to a defined zero reference. Zeroing (homing) redefines the joint's reported position to read zero at its current physical pose, without commanding any motion:

$$
q_j^{\text{offset,new}} = q_j^{\text{offset,old}} - q_j^{\text{present}}
$$

so that the joint's reported position is subsequently given by

$$
q_j^{\text{present}} = q_j^{\text{actual}} + q_j^{\text{offset}}.
$$

This operation is used to establish a repeatable, operator-defined reference configuration (typically a fully retracted, de-rotated state) after physical setup or re-mounting of the tubes, and all subsequent absolute commands are interpreted relative to that reference.

---

## 4. Actuation and Sensing Hardware

### 4.1 Actuators

Each of the four joints is driven by an independent rotary servo actuator with an integrated magnetic encoder, brushed DC motor, and embedded closed-loop position/velocity controller:

- **Rotational joints (ITR, OTR):** compact rotary servo actuators, driving the tube through a **1:40 worm-gear reduction**. The worm-gear reduction provides both the large speed reduction needed to resolve sub-degree tube motion from a small, fast actuator and a self-locking transmission that resists back-driving of the tube under external load.
- **Translational joints (ITT, OTT):** rotary servo actuators of a larger frame size, driving a **lead-screw carriage** with a **2.54 mm lead-screw pitch**, converting the actuator's rotary output into linear translation of the tube carriage.

### 4.2 Feedback sensing

Each actuator's internal magnetic encoder provides absolute angular position feedback at the motor shaft with a resolution of **4096 counts per revolution**, i.e., an angular resolution of

$$
\Delta\theta_{\text{enc}} = \frac{360^{\circ}}{4096} \approx 0.0879^{\circ}/\text{count}
$$

at the motor shaft. This resolution is subsequently amplified by the transmission ratio (worm gear or lead screw) when referred back to the tube, as described in Section 5.

### 4.3 External tip-pose sensing (validation only)

For model validation and offline performance characterization, the distal tip of the inner tube can additionally be instrumented with a passive retro-reflective marker and tracked with an external optical/electromagnetic tracking system (a commercial surgical navigation-grade tracker) at update rates on the order of tens of hertz. This tip-tracking measurement is **not** part of the real-time control loop; it is used exclusively to validate the kinematic model against ground-truth tip pose (Section 8) and is not fed back into any controller.

### 4.4 Communication and bus architecture

The four actuators are partitioned across **two independent serial buses**, one dedicated to the rotational pair (ITR, OTR) and one to the translational pair (ITT, OTT). Each bus uses a half-duplex, packet-based serial protocol with individually addressable device IDs, allowing each joint to be commanded and read back independently while sharing bus bandwidth with its paired joint. Partitioning rotational and translational joints onto separate buses decouples the update-rate and traffic characteristics of the two motion types (rotation is commanded far more frequently and at higher bandwidth than translation in typical steering sequences).

---

## 5. Joint-Space to Actuator-Space Mapping

Because the physically meaningful quantities (tube rotation angle in degrees, tube translation in millimeters) are related to the actuator's native units (encoder counts, internal velocity/acceleration profile units) through fixed mechanical ratios, every commanded motion is passed through a deterministic unit conversion before being written to the actuator, and every read-back position is passed through the inverse conversion before being reported as tube-space state.

### 5.1 Rotational conversion

Let $N$ denote the signed transmission (gear) ratio between the actuator shaft and the tube ($|N| = 40$, with sign fixed by the convention of Section 2.3). A commanded tube rotation $\theta_{\text{tube}}$ (degrees) corresponds to a motor-shaft rotation

$$
\theta_{\text{motor}} = N \cdot \theta_{\text{tube}},
$$

which is converted to an encoder count target via

$$
c = \frac{\theta_{\text{motor}}}{\Delta\theta_{\text{enc}}} = \frac{N \cdot \theta_{\text{tube}}}{360^{\circ}/4096}.
$$

### 5.2 Translational conversion

Let $p$ denote the signed lead-screw pitch (mm per actuator revolution, $|p| = 2.54\ \text{mm}$). A commanded tube translation $x_{\text{tube}}$ (mm) corresponds to an actuator revolution count $x_{\text{tube}}/p$, and thus to an encoder count target

$$
c = \frac{x_{\text{tube}}}{p}\cdot 4096.
$$

### 5.3 Velocity and acceleration profiles

Commanded motions are executed under a trapezoidal (or triangular, for short moves) velocity profile, parameterized by a maximum profile velocity and profile acceleration expressed in the actuator's native units — increments of $0.229\ \text{rev/min}$ for velocity and $214.577\ \text{rev/min}^2$ for acceleration. A desired tube-space rate $\dot{q}$ is converted to actuator revolutions per minute using the same gear ratio or lead-screw pitch as above (e.g., for rotation, $\omega_{\text{motor}} = N\dot\theta_{\text{tube}}/6$ rev/min), and then divided by the velocity LSB to obtain the profile-velocity register value. This ensures that a desired tube-space speed and acceleration limit is enforced consistently regardless of the transmission ratio of the joint being commanded.

---

## 6. Control System Architecture

### 6.1 Low-level closed-loop control

Position (and, transiently, velocity) control is closed at the actuator level: each servo runs an internal feedback loop, using its own magnetic encoder as the sole feedback sensor, that drives the motor current/PWM so as to track a commanded goal position under a bounded velocity/acceleration profile. This inner loop executes entirely within the actuator's own firmware and is not exposed to, or overridden by, the supervisory layer; the supervisory layer only writes goal positions and profile limits, and reads back present position, present velocity, and present load.

Because the physical stroke of the translational joints and the essentially unbounded rotation of the rotational joints both require multi-turn operation of the actuator (i.e., encoder counts that exceed a single mechanical revolution of the actuator shaft), all four joints operate in an **extended, multi-turn position-control mode** rather than a single-turn position mode. This is a necessary consequence of the transmission ratios in Section 5: even a moderate tube rotation or translation corresponds to many actuator revolutions once passed through the worm gear or lead screw.

### 6.2 Supervisory command layer

Above the per-joint closed loop, a supervisory layer is responsible for:

- Converting task-relevant joint commands (absolute or relative, Section 3) into actuator-native goal positions and profile limits (Section 5);
- Maintaining a software-enforced record of each joint's calibrated zero reference (Section 3.3);
- Aggregating individual joint states into a single robot joint-state vector $\mathbf{q}(t)$ published at a fixed update rate; and
- Enforcing safety interlocks (Section 6.3) before any command is dispatched to the actuator.

This layer does **not** currently close a loop around the *task-space* tip pose: the mapping from joint command to physical tip position is realized open-loop, relying on the accuracy of the kinematic model (Section 7) and the fidelity of the low-level position control, rather than on real-time tip-pose feedback. Closing this outer loop — using, e.g., the external tip-tracking measurement of Section 4.3 as a feedback signal — is identified as a natural extension of the present architecture (Section 9).

### 6.3 Safety interlocks

Because extended (multi-turn) position-control mode does not itself enforce a mechanical travel limit, the translational joints (ITT, OTT), which have a genuine finite physical stroke bounded by hard mechanical stops, are protected by a **software-enforced travel-range check**: every commanded translational target is validated against a calibrated $[\,q_{\min}, q_{\max}\,]$ interval before being dispatched, and out-of-range commands are rejected rather than being sent to the actuator. This software check is the only safeguard against driving a translational joint into its hard stop, and the supervisory layer refuses to initialize a translational joint for which these limits have not been explicitly calibrated. The rotational joints, being continuous, carry no analogous travel limit.

### 6.4 Torque enable/disable

Each joint's torque (i.e., whether the actuator is actively holding/driving position, versus free to be back-driven by hand) can be independently enabled or disabled. Disabling torque is used during manual setup, homing, and inspection, when the tubes must be repositioned by hand without commanding motion through the control loop.

---

## 7. Kinematic Model: From Joint Space to Backbone Shape

### 7.1 Purpose

While the low-level control loop only ever needs to know the four joint values, understanding and predicting the robot's *physical* behavior — the three-dimensional curve traced by the backbone, and the pose of the distal tip — requires a mechanics-based model of how the two pre-curved, elastically interacting tubes deform under a given combination of relative rotation and translation. This model is used offline, for planning, validation, and post-hoc analysis; it is not evaluated inside the real-time control loop.

### 7.2 Geometrically exact (Cosserat) rod formulation

Each tube is modeled as a geometrically exact (Cosserat/Kirchhoff) rod: an inextensible, unshearable elastic curve parameterized by arc length $s$, described by a position $\mathbf{p}(s) \in \mathbb{R}^3$ and an orientation frame $\mathbf{R}(s) \in SO(3)$. The frame evolves along the backbone according to

$$
\mathbf{R}'(s) = \mathbf{R}(s)\,\widehat{\mathbf{u}}(s), \qquad \mathbf{p}'(s) = \mathbf{R}(s)\,\mathbf{e}_3,
$$

where $\mathbf{u}(s) \in \mathbb{R}^3$ is the local curvature/torsion vector expressed in the material frame, $\widehat{(\cdot)}$ denotes the skew-symmetric cross-product matrix, and $\mathbf{e}_3$ is the local tangent direction. Each tube possesses an intrinsic (stress-free) precurvature $\mathbf{u}^{*}(s)$, and the internal bending/torsional moment is related to the *deviation* from that precurvature by a linear constitutive law,

$$
\mathbf{m}(s) = \mathbf{K}(s)\,\big(\mathbf{u}(s) - \mathbf{u}^{*}(s)\big), \qquad \mathbf{K}(s) = \operatorname{diag}\!\big(EI,\ EI,\ GJ\big),
$$

where $E$ is the tube material's elastic modulus, $I$ the second moment of area of the (typically annular) cross-section, $G$ the shear modulus, and $J$ the polar moment of area. Static equilibrium of internal force $\mathbf{n}(s)$ and moment $\mathbf{m}(s)$, in the absence of distributed external load, requires

$$
\mathbf{n}'(s) = \mathbf{0}, \qquad \mathbf{m}'(s) + \mathbf{p}'(s)\times \mathbf{n}(s) = \mathbf{0}.
$$

### 7.3 Coupling between concentric tubes

Where the two tubes overlap, they are constrained to share a common backbone position and tangent (they are threaded one inside the other), while each tube contributes its own curvature and precurvature to the pair. This produces a superposition of curvature along the overlapping region and a **torsional interaction**: rotating one tube at its base does not rotate the overlapping segment of that tube rigidly, because the tube itself twists elastically along its torsionally-compliant length before that rotation is fully transmitted to the curved, engaged region distally. This distributed torsional compliance — "torsional windup" — means that a commanded base rotation of, say, $180°$ is generally **not** realized as a full $180°$ reorientation of the curved distal section; a portion of the commanded rotation is absorbed as elastic twist along the tube's straight, torsionally-free transmission length before it reaches the curved, load-bearing segment.

### 7.4 Joint variables as boundary conditions

The four joint variables enter the model as boundary conditions at the proximal (base) end of each tube: $q_{\mathrm{OTT}}$ and $q_{\mathrm{ITT}}$ set the axial insertion of each tube's curved segment relative to the fixed base, and $q_{\mathrm{OTR}}$ and $q_{\mathrm{ITR}}$ set the base rotation angle applied to each tube. As noted in Section 2.2, because the inner tube's carriage is physically transported by the outer tube's translation stage, the *relative* insertion depth of the inner tube with respect to the outer tube's distal end used by the model is a function of both $q_{\mathrm{OTT}}$ and $q_{\mathrm{ITT}}$, not $q_{\mathrm{ITT}}$ alone. Solving the coupled boundary-value problem above for a given $\mathbf{q}$ yields the full backbone shape $\mathbf{p}(s)$, from which the distal tip position and orientation are read off directly.

---

## 8. Model–Hardware Correspondence

The kinematic model of Section 7 was validated by comparing its predicted tip trajectory against externally tracked tip positions (Section 4.3) over a series of hardware trials spanning the reachable joint space. This comparison served two purposes: (i) fixing the sign convention of Section 2.3, which cannot be derived from the mechanical drawings alone because it depends on the as-built handedness of the worm-gear transmission and the handedness of the external tracking frame; and (ii) identifying and attributing the dominant sources of model–hardware discrepancy, chiefly torsional windup in the inner tube's torsionally-free transmission length (Section 7.3), with curvature superposition effects between the two tubes contributing a smaller secondary component. The validated model achieves tip-position agreement with hardware on the order of a few millimeters root-mean-square error across the trial set, which is treated as the present accuracy bound of the open-loop joint-to-tip-pose mapping.

---

## 9. Limitations and Future Work

The present control architecture closes a position/velocity loop only at the level of the individual actuator; there is no closed loop around the mechanical state of the tube (encoder feedback cannot observe elastic deflection, backlash in the worm gear, or torsional windup) nor around the task-space tip pose. Consequently, the accuracy of any commanded tip pose is bounded by the fidelity of the open-loop kinematic model of Section 7, and by any unmodeled elastic or friction effects in the transmission. Two directions follow naturally from this limitation:

1. **Improved open-loop kinematic modeling** — extending the model to account more completely for curvature superposition and distributed torsional compliance between the two tubes, so that a commanded joint vector $\mathbf{q}$ more accurately predicts the realized tip pose without requiring any additional sensing.
2. **Closed-loop tip feedback** — using the external tip-tracking measurement (or an equivalent onboard sensing modality) as a real-time feedback signal to a task-space controller, closing the loop that is presently open between commanded joint values and realized tip pose.

Either direction would move the system from its current open-loop, model-mediated joint-to-tip mapping toward a controller whose task-space accuracy is bounded by sensing and feedback bandwidth rather than by model fidelity alone.
