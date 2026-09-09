# Validating the Cosserat model against the ICRA2027 NDI trials

**What was done:** the geometrically exact Cosserat model in [`modeling/ctr/`](ctr/)
was run through the exact commanded joint sequence of every free-space
experiment in [`ICRA2027/exp_info.md`](../ICRA2027/exp_info.md) and compared
against the NDI tip recordings — 6 sets × 3 trials = **18 trials, 47 782
samples**. Exp 4 is excluded (logged as "Random Exp — don't count").

Per-experiment figures are in [`figures/validation/`](figures/validation/), one
per experiment number (`exp1.png` … `exp19.png`), plus four summary figures and
[`metrics.csv`](figures/validation/metrics.csv) with every number below.

```bash
cd modeling
python3 validation/run_validation.py            # all 18 trials + summaries
python3 validation/run_validation.py --fit      # re-identify geometry first
python3 validation/windup_vs_superposition.py   # the 180° attribution
python3 validation/friction_study.py            # friction, calibrated on set2
```

---

## Headline results

1. **The model reproduces the measured tip paths to ~1.5 mm RMS** across five of
   six sets, and to **0.70 mm** on the cleanest set — over a 34 mm advance *and*
   a full 360° revolution. Trial-to-trial repeatability in the existing
   evaluation is 0.08–0.95 mm, so the model is running within about 2× the
   measurement's own noise floor.
2. **The 180° rotation shortfall is torsional windup, not curvature
   superposition.** This *reverses* the conclusion in
   [`EVALUATION_REPORT.md`](../ICRA2027/EVALUATION_REPORT.md) §3.5. See §3.
3. **The commanded rotation sense is inverted** relative to a right-handed
   convention about the guide axis. Wrong-way round, set2 scores 16.1 mm; the
   right way round, 1.5 mm. See §4.
4. **~2.9 mm of inner tube is already deployed at ITT = 0.** Identifying it took
   set2 from 4.00 mm to 0.70 mm RMS. See §5.
5. **The tubes are not the 50 mm arcs the design description specifies.** They
   were measured at **57 mm** (2026-09-08), which corroborates the model for the
   outer tube but conflicts sharply with what `set2` says about the inner one
   (~39.5 mm, pinned by four independent observables). See §10. The two-tube
   assembly's curvature also *tapers* along an advance in a way no
   constant-curvature pair can produce. See §6.
6. **Translation undershoots by 0.66 mm on every one of 8 segments**, and the
   model cannot explain it — it predicts the commanded length exactly. See §7.
7. **Friction is now modelled**, calibrated independently on `set2` (595 N/m of
   guide contact reproduces its 8.2° loss exactly). It is a large effect but
   **not one-signed**, so it does not simply close §3's gap — and it predicts an
   11 mm hysteresis loop that would be easy to measure. See §9.

---

## 1. Method

Model and measurement live in unrelated frames — the NDI tracker was never
registered to the robot's guide — so shapes are compared after a **rigid
(rotation + translation, no scaling) Kabsch alignment of the whole trial**.
What survives that alignment is real shape disagreement. Alongside it we report
**frame-invariant scalars** (path length, fitted arc radius, swept angle) that
need no registration and so cannot be flattered by it.

Correspondence is by *fraction of each commanded step completed*, not global arc
length: the steps are separated by operator pauses of arbitrary length, and one
global parameter would smear a disagreement in one step across the others.

Data conditioning (spike rejection, smoothing, motion segmentation) is
**imported from `ICRA2027/ndi_tip_analysis.py`** rather than reimplemented, so
the model is scored against exactly the traces the existing evaluation used.

**Geometry was identified from translation segments only.** Everything about
rotation — which is where the interesting physics is — is therefore an
out-of-sample prediction, not a fit.

### What the tracker measures

The marker is on the **inner tube's tip**: in set2 the outer tube never moves and
the tracker still records a 34 mm advance. That matters, because while the outer
tube is advanced ahead of the inner one, the inner tip lies *part way along* the
backbone, not at its distal end. Model paths are sampled with
`Solution.position_at(inner_deployed_length)` accordingly.

---

## 2. Agreement per experiment set

Whole-trial registration RMS, mean of the 3 trials (`summary_error.png`):

"Measured tubes" is the config as shipped — both tubes at the measured 57 mm
radius, no home offset. "Identified" keeps that 57 mm for the outer tube and
adds the two things the data pins independently: `R_inner` = 39.5 mm and the
2.9 mm inner-tube home offset (§5, §10).

| Set | Commanded sequence | Measured tubes | Identified | Rotation step |
| :-- | :-- | --: | --: | :-- |
| `set1`  | T 20 mm                        | 0.19 mm | **0.18 mm** | — |
| `set2`  | T 35 mm → ITR 360°             | 5.05 mm | **0.70 mm** | inner tube alone |
| `set3a` | T 35 mm → OTR 360°             | 13.44 mm | **15.12 mm** | outer tube only |
| `set3b` | T 35 mm → OTR+ITR 360°         | 5.35 mm | **3.13 mm** | co-rotation |
| `set4a` | T 35 mm → ITR 180° → T 35 mm   | 1.65 mm | **2.31 mm** | inner vs overlapped outer |
| `set4b` | T 17.5 mm → ITR 180° → T 17.5 mm | 1.23 mm | **1.21 mm** | same, half overlap |
| **All** | | 4.49 mm | 3.77 mm | |
| **Excluding `set3a`** | | 2.69 mm | **1.51 mm** | |

`set2` (`exp5.png`) is the strongest single result: the model tracks a 34 mm
advance followed by a complete revolution with **0.70 mm** RMS, matching the
rotation circle radius to 16.8 vs 16.9 mm and the advance arc radius to 39.5 vs
40.1 mm.

**`set3a` is the one clear failure**, and it fails in one specific way: the model
predicts 302° of tip sweep where 356° was measured. The deviation trace in
`exp8.png` oscillates twice per revolution — the signature of a *phase* error on
a circle of about the right size, not a shape error. The model is over-predicting
how much of the outer tube's rotation is lost.

---

## 3. The 180° shortfall is windup, not superposition

`EVALUATION_REPORT.md` §3.5 reports that inner-tube rotation reaches the tip
almost perfectly at OTT = 0 (97.7 %) but only 33–42 % once the outer tube is
advanced, and attributes this to curvature superposition, explicitly ruling out
torsional windup on the grounds that the 360° trials bound windup at 2–8°.

The model contains both mechanisms, so it can separate them rather than argue
from bounds. Running it twice — once as built, once with `GJ → ∞` so that every
commanded degree arrives at the deployed section and only superposition can
act (`windup_vs_superposition.png`):

| Set | Rotates | OTT | Cmd | Measured | Model (full) | Model (torsionally rigid) |
| :-- | :-- | --: | --: | --: | --: | --: |
| `set2`  | ITR | 0 mm    | 360° | 351.8° (97.7 %) | 360.0° (100 %) | 360.0° (100 %) |
| `set3a` | OTR | 35 mm   | 360° | 356.3° (99.0 %) | 300.8° (83.6 %) | 360.0° (100 %) |
| `set3b` | both | 35 mm  | 360° | 358.3° (99.5 %) | 360.0° (100 %) | 360.0° (100 %) |
| `set4b` | ITR | 17.5 mm | 180° | 76.2° (42.3 %) | 55.1° (30.6 %) | 179.6° (**99.8 %**) |
| `set4a` | ITR | 35 mm   | 180° | 60.0° (33.3 %) | 36.9° (**20.5 %**) | 178.3° (**99.1 %**) |

**With the tubes made torsionally rigid, `set4a` delivers 99.1 % of the
commanded rotation.** Curvature superposition accounts for roughly **1
percentage point** of the loss; torsional compliance accounts for the rest. The
full model predicts 21 % where 33 % was measured — the right regime and the
right mechanism, from geometry and mechanics alone with nothing fitted to
rotation data. The rigid model's 99 % is not close to anything measured.

The mechanism is explicit in the model's state: at `set4a`'s final pose only
**−26.9° of the −180° commanded relative roll** has reached the deployed
overlap. The other 153° is wound into the inner tube's ~273 mm of torsionally
free transmission behind the guide.

**Why the earlier bound failed.** Windup is not a constant angular offset — it
is driven by the torque the two tubes exert on each other, which is *zero* when
they are aligned or co-rotating, and zero when only one tube is deployed. So the
360° trials genuinely have near-zero windup (`set3b`: 0.0° of relative roll
between the tubes throughout; `set2`: only one tube deployed, no coupling at
all), and that tells you nothing about the 180° trials, where the tubes are
driven into strong opposition. The 2–8° bound was measured on exactly the
configurations where the mechanism switches off.

**The model over-predicts windup** in every case where it acts (21 % vs 33 % on
`set4a`, 31 % vs 42 % on `set4b`, 84 % vs 99 % on `set3a`). That is the expected
direction of error: the model treats the whole retracted length as
torsionally free and frictionless, whereas the real tubes are supported by, and
rub against, each other, the guide and the carriages.

---

## 4. The rotation sense is inverted

Modelled with positive OTR/ITR as right-handed about +z, the model walks set2's
rotation circle **backwards**: the whole-trial residual is 16.1 mm, roughly
twice the circle's 16.8 mm radius, which is what you get from a circle traversed
the wrong way. Reversing it gives 1.5 mm. A rigid registration uses proper
rotations only, so it *cannot* absorb a handedness error — this is a real
chirality mismatch.

All four sign combinations, whole-trial RMS (mm):

| OTR, ITR | `set2` | `set3a` | `set3b` | `set4a` |
| :-- | --: | --: | --: | --: |
| +1, +1 | 16.10 | **3.71** | 15.51 | **2.08** |
| +1, −1 | **1.56** | **3.71** | 15.68 | 2.21 |
| −1, +1 | 16.10 | 15.36 | 4.99 | **2.08** |
| −1, −1 | **1.56** | 15.36 | **2.78** | 2.21 |

`set2` (a single tube — unambiguous) and `set3b` (a pure rigid roll — kinematically
trivial) independently both demand a **negative** sense, by 10× and 5×. That is
now set in [`config/ct_sdr.yaml`](config/ct_sdr.yaml) and locked by a test.

`set3a` prefers the opposite sense for OTR, which is the same failure as §2 —
but `set3a` is also the one configuration where the model's physics already
disagrees, so it is the contaminated test, and reversing OTR alone would
contradict `set3b`.

⚠️ **Two causes are indistinguishable from this data:** the 1:40 worm gear
reversing the commanded sense, or the NDI recordings being in a left-handed
frame. Check the worm-gear direction against a physical rotation before relying
on the absolute sense.

---

## 5. The inner tube is ~2.9 mm deployed at ITT = 0

The advances constrain the tubes' *curvature* but are blind to how much tube is
already out — advancing along a uniform arc traces the same circle wherever you
start. The rotation circles are the opposite: their radius is the tip's distance
from the guide axis, which depends on total deployed length.

Solving `set2`'s measured 16.83 mm rotation circle for the offset gives
**2.90 mm**, and that same number then reproduces, without further adjustment:

| Quantity | Measured | Model |
| :-- | --: | --: |
| Advance path length | 34.40 mm | 34.39 mm |
| Advance chord | 33.40 mm | 33.33 mm |
| Advance arc radius | 39.5 mm | 39.5 mm |
| Rotation circle radius | 16.83 mm | 16.83 mm |

The outer tube's offset is **not** identifiable here: once it covers the inner
tube's tip, advancing it further changes nothing the tracker can see. It is
pinned to zero because `set2` runs at OTT = 0 and is the calibration case — a
non-zero value leaves a stub of outer tube stiffening the first few millimetres
there and pulls that circle to 15.06 mm.

---

## 6. What the tubes trace when deployed

Fitting circles to nested sub-arcs *within* each advance (so the trend is
within-trial, free of any cross-session confound) gives the tip path's mean
radius as a function of how far the tubes have advanced — `curvature_profile.png`.
A tube of uniform curvature gives a flat line.

**Inner tube alone** (`set2` T1, `set4a` T3, `set4b` T3 — 9 segments): roughly
flat at **39–42 mm** — against 50 mm in the design description and **57 mm
measured on the free tubes** (§10). The identified 39.5 mm reproduces set2 to
0.70 mm RMS across four independent observables, so it is well determined; the
disagreement with the direct measurement is the subject of §10.

**Both tubes advancing** (15 segments): the measured radius **rises from ~36 mm
at 12 mm of advance to ~61 mm at 30 mm** — the assembly starts as curved as the
inner tube alone and progressively straightens. This is highly repeatable across
9 trials in three different sets, and **no constant-curvature pair of tubes
reproduces it**: the best fit over (R_outer, tip lead-in, R_inner) improves the
turn-angle residual only from 5.25° to 4.78°, and the fitted curves in
`curvature_profile.png` are flat where the data slopes.

Candidate explanations, none of which the current data can separate:

* the outer tube's pre-curved section is not uniform, or does not start at its tip;
* the 0.25 mm-per-side clearance between the tubes lets them lag each other,
  which the zero-clearance model cannot represent;
* OTT and ITT did not advance the way the log implies (see §8).

Two further internal tensions worth knowing about:

* `set1` (OTT 20, ITT 20) has a curvature profile essentially **indistinguishable
  from inner-tube-only deployment**, unlike `set3a`/`3b`/`4a`, which all show
  clear outer-tube straightening at the same advance.
* For `set3b`, the advance path's own shape implies the tip ends ~11 mm off the
  guide axis, but its rotation circle measures 17.1 mm. Those two measurements
  of the same pose are not consistent with any rigid configuration. An
  un-calibrated marker offset would explain it — but §10 rules that out for the
  *inner* tube on a path-length argument, so for `set3b` it remains open.

**A tube of `tip_straight_length` (a straight lead-in between the curved section
and the distal tip) was added to the model** while chasing this, since the
measured curvature sits *proximal* — near the guide — with the distal end nearly
straight. It fits `set3a`'s profile well (3.1 mm vs 8–25 mm for the
alternatives) but contradicts `set1`, so it is **not** claimed as the
explanation. The capability is in `ctr/tube.py` with tests; the shipped config
leaves it at zero.

---

## 7. What the model does *not* explain

* **Translation undershoot.** All 8 translation segments undershoot the command,
  by 0.66 mm on average (`set1` 19.45 vs 20, `set2` 34.48 vs 35, `set4b` T3
  16.41 vs 17.5). The model predicts the commanded length *exactly* — it is
  inextensible with a rigid transmission. This is a lead-screw scale factor,
  backlash or axial compliance, and it is calibratable.
* **The 8.2° lost in `set2`'s revolution.** Only one tube is deployed, so the
  model predicts exactly 360° with zero coupling. The 2.3 % shortfall is
  transmission backlash or tube-to-guide friction — outside a frictionless
  elastostatic model.
* **`set3a`** (§2, §3).

---

## 8. What should be done better

Ordered by how much they would improve the next round.

1. **Log joint positions alongside the NDI stream.** `ndi_tip_recorder.py`
   records `stamp, elapsed, x, y, z, frame_id` and nothing else. Almost every
   ambiguity in this report — whether OTT and ITT moved together or in sequence,
   whether the commanded rotation actually happened, where the home positions
   were — would be settled by a synchronised joint channel. The ROS nodes already
   publish `JointState`; the recorder just doesn't subscribe. **Highest payoff
   per line of code in this entire analysis.**

2. **Pivot-calibrate the marker and register the NDI frame to the guide.**
   Everything here needed a rigid Procrustes fit, which absorbs real error, and
   §6 shows evidence of an uncalibrated marker offset. A one-time registration
   turns "shape agreement" into *absolute accuracy*, which is what a drilling
   robot actually needs to report.

3. **Measure the tubes' free shape directly** — photograph them against a grid,
   or CT them. §6 is the largest unresolved discrepancy and one afternoon with a
   camera would close it. Record the curved length and where the curve starts,
   not just the radius.

4. **Verify the rotation sense physically** (§4) and record the tubes' home
   offsets (§5) as part of the setup procedure.

5. **Resolve `set3a` vs `set3b`.** They demand opposite outer-tube rotation
   senses. Re-run both with joint logging.

6. **Record the drilling experiments.** `exp_info.md` lists three drill runs with
   no CSVs. The model's externally-loaded branch — the whole point of using this
   paper rather than a constant-curvature model — is completely untested against
   hardware. Even a single instrumented run with a force reading would exercise it.

7. **Sweep rotation in small steps with dwells, and sweep back.** The model
   predicts snap-through and hysteresis (see [README](README.md) and
   `figures/workspace_itr_sweep.png`). Continuous 360° sweeps at one speed
   cannot see either. Step ITR in ~15° increments with a 2 s dwell, then reverse:
   if the return path differs from the outbound one, that is the hysteresis, and
   it bounds the friction the model currently omits.

8. **~~Add friction to the transmission model.~~** *Done — see §9. The result
   did not go the way this recommendation predicted, and the recommendation as
   originally written was wrong: friction is not a one-signed correction and
   fitting it to `set4a`/`set4b` would not have been legitimate anyway, since
   `set2` calibrates it independently. What replaces it is §9's last paragraph:
   run the hysteresis experiment.*

9. **Calibrate the translation scale factor** (§7) — 24/24 one-sided undershoot
   is free accuracy.

---

## 9. Friction

Added after the first pass of this report, in [`ctr/friction.py`](ctr/friction.py),
with the study in [`validation/friction_study.py`](validation/friction_study.py)
and `figures/validation/friction_study.png`.

### What is modelled

Two contacts, with different materials and different coefficients:

| Contact | Where | Dry μ | Wetted μ |
| :-- | :-- | --: | --: |
| nitinol on nitinol | inner tube inside outer, over the overlap | **0.35** | 0.15 |
| nitinol on stainless steel | tubes inside the rigid guide | **0.25** | 0.12 |

Dry values are appropriate for a robot drilling in air or bone. NiTi-on-NiTi is
the higher of the two because it is adhesive and galls readily — that ordering
is not incidental and is locked by a test. These are order-of-magnitude
literature values, **not measurements of your tubes**; the study treats the
contact force as the free parameter and calibrates it.

Friction enters each tube's torsion ODE as a distributed axial moment, equal
and opposite between contacting tubes, so the axial moment balance
`Σ GJ_i u_i,z = (R e₃)·m` is preserved exactly — friction moves torque between
tubes, it does not create any. The sliding sense comes from the change in
commanded *relative* roll, so a co-rotation or a pure translation correctly
produces no sliding at all.

### Calibration: `set2` is the free lunch

`set2` deploys **only** the inner tube, so there is no second tube to couple to.
The frictionless model therefore predicts a perfect 360° with exactly zero
windup and has *no mechanism whatsoever* for the 8.2° the tracker recorded.
Everything in that shortfall is friction against the guide bore, which makes it
a clean calibration — the same role `set2` played for the tube geometry in §5.

The response is linear, and **595 N/m of guide contact force reproduces the
measured 8.2° exactly**. That number is physically sensible rather than a
fudge: it corresponds to the guide reacting the inner tube's 1.66 N·m
straightening moment over a ≈53 mm effective grip.

### Prediction: large, and *not* one-signed

Applying that calibrated friction to the two 180° rotations — out of sample,
nothing here was fitted to them:

| Set | Measured | Frictionless | With calibrated friction |
| :-- | --: | --: | --: |
| `set4a` (OTT 35) | 60.0° | 36.9° | **70.5°** |
| `set4b` (OTT 17.5) | 76.2° | 55.2° | **38.7°** |

Friction closes most of `set4a`'s gap (error 23.1° → 10.5°) and opens `set4b`'s
(21.0° → 37.5°). It is a large effect — tens of degrees — but its **sign is
configuration-dependent**, so it is not the correction §3's gap was waiting for.

The reason is that friction resists *both* legs of a sweep. Over a commanded
0 → 180°, the relative roll at the deployed section is not monotonic: it rises,
peaks, then relaxes back toward alignment. On the way up friction is a pure
torque sink between actuator and tip and costs delivered rotation; on the way
down it resists the relaxation and *holds* the assembly wound. Which effect
dominates depends on the configuration. Held at one fixed pose part way up
(ITR = 90°, a-priori coefficients) the first effect is all you see, and it is
small: 14.26° → 14.00° delivered, a 12 µm tip shift.

⚠️ The model uses **one sliding sense per commanded step**, taken from the
change in commanded relative roll. That is right on a monotonic leg and wrong
where the local sliding reverses part way through — which is exactly what
happens here. A proper stick–slip treatment (a complementarity solve, with the
local relative rotation rate setting each contact's state) is what would make
the sign trustworthy. Until then, treat §9's `set4a`/`set4b` numbers as showing
*scale*, not as a validated prediction.

### The prediction worth testing

Friction makes the model **path-dependent**, which the frictionless model
cannot be. Sweeping ITR out to 180° and back with the sliding sense reversed on
the return leg gives two different tip paths:

| Set | Mean gap | Max gap |
| :-- | --: | --: |
| `set4a` | 9.02 mm | **11.08 mm** |
| `set4b` | 2.15 mm | 4.06 mm |

Against a measured trial-to-trial repeatability of 0.08–0.95 mm, an 11 mm loop
is enormous — this is a cheap, unambiguous experiment. **Sweep ITR out and
back, in steps with dwells, and record the return path.** If the loop is there,
it calibrates the friction properly and settles the sign question; if the return
retraces the outbound path, friction is far smaller than these coefficients say
and §3's gap has to be explained by the transmission's free length instead.
Either answer is worth more than any amount of further modelling.

## 10. The measured 57 mm radius, and where it disagrees

The tubes were **measured directly on 2026-09-08: both 57 mm radius of
curvature**, superseding the 50 mm in the design description.
[`config/ct_sdr.yaml`](config/ct_sdr.yaml) now carries 57 mm.

That measurement lands differently on the two tubes.

### The outer tube: corroborated

§6 identified `R_outer` ≈ 53.6 mm and flagged it as **poorly determined** — the
both-tube advances taper in a way no constant-curvature pair reproduces, so the
fit had little to grip. 57 mm is close to that and is a direct measurement, so
it simply wins. It also *improves* the two sets that lean hardest on the outer
tube (`set3a` 15.3 → 13.4 mm, `set4a` 2.3 → 1.7 mm).

### The inner tube: a real conflict

`set2` deploys the inner tube alone and pins its radius through four
observables at once. 39.5 mm matches all four; 57 mm misses all four:

| `set2` observable | Measured | Model at 39.5 mm | Model at 57 mm |
| :-- | --: | --: | --: |
| advance path length | 34.40 mm | 34.39 | 35.00 |
| advance chord | 33.40 mm | 33.87 | 34.45 |
| advance arc radius | 39.5 mm | **39.5** | **57.0** |
| rotation circle radius | 16.83 mm | **16.83** | **10.41** |
| whole-trial RMS | — | **0.70 mm** | **5.04 mm** |

**The obvious escape is ruled out.** If the tracker marker sat off the tube
centreline by 17.5 mm toward the centre of curvature, a 57 mm tube would trace a
39.5 mm arc — which would reconcile everything. But the same offset would make
the marker travel only `35 × 39.5/57 = 24.3 mm` while the tube advances 35 mm.
The marker actually traces **34.40 mm, 98 % of the advance**, so it is on the
centreline and the traced arc radius *is* the deployed radius.

### What it costs

Mean whole-trial RMS over all 18 trials (and over the five sets excluding
`set3a`, which fails for unrelated reasons — §2):

| Tube parameters | set1 | set2 | set3a | set3b | set4a | set4b | mean | excl. set3a |
| :-- | --: | --: | --: | --: | --: | --: | --: | --: |
| measured 57 / 57 | 0.19 | 5.05 | 13.44 | 5.35 | 1.65 | 1.23 | 4.49 | 2.69 |
| 57 / 57 + 2.9 mm home | 0.19 | 3.72 | 13.92 | 4.00 | 1.59 | 1.12 | 4.09 | 2.12 |
| **57 outer / 39.5 inner + home** | 0.18 | **0.70** | 15.12 | 3.13 | 2.31 | 1.21 | **3.77** | **1.51** |
| fitted 53.6 / 39.5 + home | 0.18 | 0.70 | 15.34 | 2.76 | 2.28 | 1.19 | 3.74 | 1.42 |

**Recommendation: use your measured 57 mm for the outer tube and keep 39.5 mm
for the inner**, with the 2.9 mm home offset. That hybrid is as good as the
fully fitted parameter set (3.77 vs 3.74 mm mean) while replacing the one
badly-determined fitted number with a measurement. Note also that the 2.9 mm
home offset is worth 0.4 mm of mean RMS *regardless* of which radius is right —
it is a separable finding.

### How to settle it

The two numbers describe different things and could both be right:

* **How was the 57 mm measured?** If a circle was fitted to the whole curved
  section including the straight transitions at each end, it reads high, while
  the distal ~35 mm that `set2` actually deploys is tighter. That would make
  both numbers correct and fits the non-uniform curvature §6 already found.
  Measuring centreline-vs-outer-surface is only worth 1.3 mm, so it is not that.
* **Is the inner tube's curvature uniform?** §6 says the inner-only traces are
  flat at 39–42 mm over deployments spanning 0–70 mm from its tip. If the free
  tube really is a uniform 57 mm over that same span, the two cannot both hold.

**The decisive experiment is one photograph.** Retract the outer tube fully,
advance the inner tube ~35 mm, and photograph the deployed shape against a grid
or ruler. That measures the deployed radius directly, in the configuration that
matters, and takes minutes.

---

## Files

| Path | What |
| :-- | :-- |
| `validation/ndi_experiments.py` | Commanded sequences, trial loading, model trajectories, registration |
| `validation/run_validation.py` | Per-trial comparison, geometry identification, figures, `metrics.csv` |
| `validation/windup_vs_superposition.py` | The §3 attribution |
| `ctr/friction.py` | Coulomb friction, tube-tube and tube-guide (§9) |
| `validation/friction_study.py` | §9: calibrate on set2, predict on set4a/b |
| `figures/validation/friction_study.png` | §9 |
| `figures/validation/exp*.png` | One figure per experiment number (18) |
| `figures/validation/summary_error.png` | §2 |
| `figures/validation/windup_vs_superposition.png` | §3 |
| `figures/validation/curvature_profile.png` | §6 |
| `figures/validation/rotation_delivery.png` | Commanded vs delivered rotation |
| `figures/validation/metrics.csv` | Every number above |
