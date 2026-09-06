"""
ndi_experiments.py -- the ICRA2027 free-space trials, as model inputs.

Bridges three things:

* the commanded joint sequences in ``ICRA2027/exp_info.md``,
* the NDI tip recordings and their conditioning pipeline, which is imported
  from ``ICRA2027/ndi_tip_analysis.py`` rather than reimplemented -- that
  module's spike rejection, smoothing and motion segmentation are already
  validated against this dataset and its segment count matches the commanded
  step count for all 18 trials,
* the Cosserat model in ``modeling/ctr``.

What the tracker actually measures
----------------------------------
The marker is on the **inner tube's tip**: in `set2` the outer tube never
moves and the tracker still records a 34 mm advance, so it cannot be on the
outer tube.  That matters for the model, because while the outer tube is
advanced ahead of the inner one the inner tip is *part way along* the
backbone, not at its distal end.  Model paths are therefore sampled with
``Solution.position_at(inner_deployed_length)``, never at the distal tip.

Translation steps: which tube moves when
----------------------------------------
`dxl_control_4dof_cli.py` drives one joint at a time, so a step written
"OTT = 35, ITT = 35" is two commands, not one simultaneous move.  Under that
reading the outer tube goes out first, during which the inner tip does not
move at all and the tracker sees nothing; the recorded motion is the inner
tube advancing inside an outer tube already at its final extension.  Every
translation step therefore reduces to the same thing -- an ITT advance at a
fixed OTT -- which is what ``SEQUENTIAL`` generates.

``SIMULTANEOUS`` generates the alternative (both tubes advancing together) so
the two can be scored against the data instead of assumed; see
``run_validation.py --advance``.
"""

from __future__ import annotations

import os
import sys
from dataclasses import dataclass, field
from typing import List, Sequence, Tuple

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_MODELING = os.path.dirname(_HERE)
_REPO = os.path.dirname(_MODELING)
ICRA_DIR = os.path.join(_REPO, "ICRA2027")

sys.path.insert(0, _MODELING)
sys.path.insert(0, ICRA_DIR)

import ndi_tip_analysis as ndi  # noqa: E402  (the existing, validated pipeline)

from ctr import CTSDR, Joints  # noqa: E402
from ctr.robot import INNER, OUTER  # noqa: E402

SEQUENTIAL = "sequential"
SIMULTANEOUS = "simultaneous"


# --- commanded sequences ----------------------------------------------------

@dataclass(frozen=True)
class Step:
    """One commanded motion.

    kind
        ``"T"`` translation or ``"R"`` rotation, matching ``exp_info.md``.
    ott, itt, otr, itr
        Joint values at the **end** of the step.  The start is the previous
        step's end (or all zeros for the first).
    moves
        Which joints this step actually drives, for labelling.
    """

    kind: str
    ott: float
    itt: float
    otr: float
    itr: float
    moves: str

    def joints(self) -> Joints:
        return Joints(ott=self.ott, itt=self.itt, otr=self.otr, itr=self.itr)


@dataclass(frozen=True)
class ExperimentSet:
    key: str
    name: str
    config: str
    steps: Tuple[Step, ...]
    trials: Tuple[Tuple[str, str], ...]
    note: str = ""


def _seq(*rows) -> Tuple[Step, ...]:
    return tuple(Step(*r) for r in rows)


#: The six sets, transcribed from exp_info.md.  Trial file choices (which
#: duplicate capture to use, which run to drop) follow ndi_tip_analysis.py so
#: the model is compared against exactly the traces the existing evaluation
#: report used.  Exp 4 is excluded there as "Random Exp - don't count".
EXPERIMENTS: List[ExperimentSet] = [
    ExperimentSet(
        key="set1", name="1st Free Space Test",
        config="OTT 20, ITT 20, OTR 0, ITR 0",
        steps=_seq(("T", 20.0, 20.0, 0.0, 0.0, "OTT+ITT")),
        trials=(("Exp 1", "exp1_20260811_215059.csv"),
                ("Exp 2", "exp2_20260811_215200.csv"),
                ("Exp 3", "exp3_20260811_215253.csv")),
    ),
    ExperimentSet(
        key="set2", name="2nd Free Space Test",
        config="OTT 0, ITT 35, ITR 360 (outer tube never moves)",
        steps=_seq(("T", 0.0, 35.0, 0.0, 0.0, "ITT"),
                   ("R", 0.0, 35.0, 0.0, 360.0, "ITR")),
        trials=(("Exp 5", "exp5_20260811_222950.csv"),
                ("Exp 6", "exp6_20260811_223714.csv"),
                ("Exp 7", "exp7_20260811_224018.csv")),
        note="Only the inner tube is deployed, so there is no tube-to-tube "
             "torsional coupling at all: the model predicts exactly zero windup.",
    ),
    ExperimentSet(
        key="set3a", name="3rd Free Space Test (Top Right)",
        config="OTT 35, ITT 35, OTR 360, ITR 0",
        steps=_seq(("T", 35.0, 35.0, 0.0, 0.0, "OTT+ITT"),
                   ("R", 35.0, 35.0, 360.0, 0.0, "OTR")),
        trials=(("Exp 8", "exp8_20260811_230128.csv"),
                ("Exp 9", "exp9_20260811_230810.csv"),
                ("Exp 10", "exp10_20260811_231038.csv")),
        note="The stiff outer tube rotates; the resultant curvature follows it.",
    ),
    ExperimentSet(
        key="set3b", name="3rd Free Space Test (Bottom Left)",
        config="OTT 35, ITT 35, OTR 360, ITR 360",
        steps=_seq(("T", 35.0, 35.0, 0.0, 0.0, "OTT+ITT"),
                   ("R", 35.0, 35.0, 360.0, 360.0, "OTR+ITR")),
        trials=(("Exp 11", "exp11_20260811_232209.csv"),
                ("Exp 12", "exp12_20260811_232456.csv"),
                ("Exp 13", "exp13_20260811_232803.csv")),
        note="Both tubes co-rotate, so the whole shape turns rigidly about +z.",
    ),
    ExperimentSet(
        key="set4a", name="4th Free Space Test (Bottom Mid)",
        config="OTT 35, ITT 35, then ITR 180, then ITT +35",
        steps=_seq(("T", 35.0, 35.0, 0.0, 0.0, "OTT+ITT"),
                   ("R", 35.0, 35.0, 0.0, 180.0, "ITR"),
                   ("T", 35.0, 70.0, 0.0, 180.0, "ITT")),
        trials=(("Exp 14", "exp14_20260811_233834.csv"),
                ("Exp 15", "exp15_20260811_234353.csv"),
                ("Exp 16", "exp16_20260811_235022.csv")),
        note="The weaker inner tube rotates against a fully overlapped outer tube.",
    ),
    ExperimentSet(
        key="set4b", name="4th Free Space Test (Bottom Right)",
        config="OTT 17.5, ITT 17.5, then ITR 180, then ITT +17.5",
        steps=_seq(("T", 17.5, 17.5, 0.0, 0.0, "OTT+ITT"),
                   ("R", 17.5, 17.5, 0.0, 180.0, "ITR"),
                   ("T", 17.5, 35.0, 0.0, 180.0, "ITT")),
        trials=(("Exp 17", "exp17_20260812_001514.csv"),
                ("Exp 18", "exp18_20260812_010331.csv"),
                ("Exp 19", "exp19_20260812_010534.csv")),
        note="Same as set4a at half the overlap.",
    ),
]

#: Trial label -> (set, csv filename), for looking an experiment up by number.
TRIAL_INDEX = {label: (es, fn) for es in EXPERIMENTS for label, fn in es.trials}


# --- measured trajectories --------------------------------------------------

@dataclass
class MeasuredTrial:
    label: str
    set_key: str
    path: str
    t: np.ndarray
    p: np.ndarray                  # (M,3) mm, smoothed, referenced to sample 0
    p_raw: np.ndarray
    segments: List[Tuple[int, int]]
    n_spikes: int = 0

    def segment(self, k: int) -> np.ndarray:
        a, b = self.segments[k]
        return self.p[a:b + 1]


def load_trial(label: str, filename: str, set_key: str,
               smooth_s: float = 0.5, jump_mm: float = 2.0,
               speed_thresh: float = 0.3, min_move_s: float = 1.0,
               min_pause_s: float = 0.8) -> MeasuredTrial:
    """Load and condition one NDI recording using the existing pipeline."""
    full = os.path.join(ICRA_DIR, filename)
    t, p, n_spikes = ndi.drop_spikes(*ndi.load_csv(full), jump_mm=jump_mm)
    t = t - t[0]
    p = p - p[0]
    ps = ndi.smooth(t, p, smooth_s)
    mov, _, _ = ndi.segment_motion(t, ps, speed_thresh, min_move_s, min_pause_s)
    return MeasuredTrial(label=label, set_key=set_key, path=full, t=t, p=ps,
                         p_raw=p, segments=mov, n_spikes=n_spikes)


# --- modelled trajectories --------------------------------------------------

def _waypoints(es: ExperimentSet, k: int, n: int, advance: str) -> List[Joints]:
    """Joint waypoints spanning step ``k``, at the resolution the tracker sees.

    A translation step emits only the part of the motion the tracker can
    observe -- see the module docstring.
    """
    step = es.steps[k]
    prev = es.steps[k - 1] if k > 0 else Step("T", 0.0, 0.0, 0.0, 0.0, "")

    if step.kind == "R":
        f = np.linspace(0.0, 1.0, n)
        return [
            Joints(ott=step.ott, itt=step.itt,
                   otr=prev.otr + a * (step.otr - prev.otr),
                   itr=prev.itr + a * (step.itr - prev.itr))
            for a in f
        ]

    both = step.ott != prev.ott and step.itt != prev.itt
    if both and advance == SIMULTANEOUS:
        f = np.linspace(0.0, 1.0, n)
        return [
            Joints(ott=prev.ott + a * (step.ott - prev.ott),
                   itt=prev.itt + a * (step.itt - prev.itt),
                   otr=step.otr, itr=step.itr)
            for a in f
        ]
    # Sequential: the outer tube is already at its final extension and only
    # the inner tube's advance is visible to the tracker.
    return [
        Joints(ott=step.ott, itt=itt, otr=step.otr, itr=step.itr)
        for itt in np.linspace(prev.itt, step.itt, n)
    ]


def model_trajectory(robot: CTSDR, es: ExperimentSet, n_per_step: int = 41,
                     advance: str = SEQUENTIAL, steps: Sequence[int] = None):
    """Modelled inner-tube-tip path for a whole experiment set.

    Returns ``(p_mm, seg_slices, solutions)`` with ``p_mm`` an ``(M,3)`` array
    referenced to the first sample, matching how the measured traces are
    referenced, and ``seg_slices`` one ``slice`` per requested step.

    ``steps`` restricts which commanded steps are simulated.  Rotation sweeps
    dominate the cost, so geometry fitting -- which only scores advances --
    passes ``steps=[0]`` and skips them.  Each step's joint waypoints are built
    from the commanded sequence, not from the previous step's simulation, so
    skipping a step does not corrupt the ones after it; only the continuation
    warm start is lost.
    """
    pts: List[np.ndarray] = []
    slices: List[slice] = []
    sols = []
    guess = None
    wanted = range(len(es.steps)) if steps is None else steps
    for k in wanted:
        wps = _waypoints(es, k, n_per_step, advance)
        start = len(pts)
        for j in wps:
            inner_tip_s = robot.deployed_lengths(j)[INNER]
            if inner_tip_s <= 1e-9:
                # Nothing of the inner tube is out, so its tip is exactly at
                # the guide exit.  No boundary value problem to solve, and at
                # zero deployment there would not be one to solve anyway.
                pts.append(np.zeros(3))
                sols.append(None)
                continue
            sol = robot.solve(j, guess=guess, n_points=60)
            guess = sol.u_z0
            pts.append(sol.position_at(inner_tip_s))
            sols.append(sol)
        slices.append(slice(start, len(pts)))
    p = np.asarray(pts) * 1e3
    return p - p[0], slices, sols


# --- registration and metrics ----------------------------------------------

def kabsch(a: np.ndarray, b: np.ndarray):
    """Rigid transform (no scaling) taking ``a`` onto ``b``, least squares.

    Returns ``(R, t, rms)`` with ``b ~ a @ R.T + t``.  Rotation is *not*
    fitted away as a nuisance for its own sake -- the tracker frame and the
    robot's guide frame are genuinely unrelated and were never registered
    during the experiments, so a rigid alignment is the only way to compare
    shapes at all.  What survives it is real shape disagreement.
    """
    ca, cb = a.mean(axis=0), b.mean(axis=0)
    h = (a - ca).T @ (b - cb)
    u, _, vt = np.linalg.svd(h)
    d = np.sign(np.linalg.det(vt.T @ u.T))
    r = vt.T @ np.diag([1.0, 1.0, d]) @ u.T
    t = cb - r @ ca
    resid = b - (a @ r.T + t)
    return r, t, float(np.sqrt(np.mean(np.sum(resid ** 2, axis=1))))


def resample_by_arclength(p: np.ndarray, n: int) -> np.ndarray:
    """Resample a polyline to ``n`` points equally spaced in arc length."""
    d = np.r_[0.0, np.cumsum(np.linalg.norm(np.diff(p, axis=0), axis=1))]
    if d[-1] <= 0:
        return np.repeat(p[:1], n, axis=0)
    q = np.linspace(0.0, d[-1], n)
    return np.column_stack([np.interp(q, d, p[:, k]) for k in range(3)])


def correspond(measured: MeasuredTrial, model_p: np.ndarray,
               model_slices: Sequence[slice], n_per_seg: int = 120):
    """Put measured and modelled paths in correspondence, segment by segment.

    Pairing is by *fraction of each commanded step completed*, not by time and
    not by fraction of the whole path: the steps are separated by operator
    pauses of arbitrary length, and a single global arc-length parameter would
    smear a disagreement in one step across all the others.
    """
    n = min(len(measured.segments), len(model_slices))
    meas, mod, bounds = [], [], []
    for k in range(n):
        a = resample_by_arclength(measured.segment(k), n_per_seg)
        b = resample_by_arclength(model_p[model_slices[k]], n_per_seg)
        bounds.append((len(meas) * 1, len(meas) + n_per_seg))
        meas.append(a)
        mod.append(b)
    return np.vstack(meas), np.vstack(mod), n


def arc_metrics(p: np.ndarray):
    """Path length, chord, best-fit circle radius/residual and swept angle."""
    length = float(np.linalg.norm(np.diff(p, axis=0), axis=1).sum())
    chord = float(np.linalg.norm(p[-1] - p[0]))
    c, u_ax, v_ax, normal, plane_rms = ndi.fit_plane(p)
    u, v = ndi.project_plane(p, c, u_ax, v_ax)
    cu, cv, r, resid = ndi.fit_circle_2d(u, v)
    ang = np.unwrap(np.arctan2(v - cv, u - cu))
    swept = float(np.degrees(abs(ang[-1] - ang[0])))
    return dict(length=length, chord=chord, radius=r, circle_resid=resid,
                swept=swept, plane_rms=plane_rms)
