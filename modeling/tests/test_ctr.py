#!/usr/bin/env python3
"""
test_ctr.py -- verification of the Cosserat concentric-tube model.

Run from anywhere:

    python3 -m unittest discover -s modeling/tests -v
    python3 modeling/tests/test_ctr.py

The tests check the model against cases with closed-form answers rather than
against itself: a single pre-curved tube must trace an exact circular arc, two
identical tubes at 180 deg must cancel to a straight line, dissimilar tubes
must combine with the stiffness-weighted average curvature, and a straight
tube under a small tip load must reproduce the Euler-Bernoulli cantilever.
Together those exercise every term of the model except the distributed-moment
one, which is checked for consistency against an equivalent tip moment.
"""

import os
import sys
import unittest

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from ctr import (  # noqa: E402
    DEFAULT_CONFIG, CTSDR, CosseratModel, ExternalLoad, Joints, ShootingError, Tube,
)
from ctr.robot import INNER, OUTER  # noqa: E402

E_NITINOL = 60e9
NU = 0.33


def make_tube(name, od, wall, length, curved_length, roc, E=E_NITINOL):
    return Tube(
        name=name,
        outer_diameter=od,
        inner_diameter=od - 2.0 * wall,
        length=length,
        curved_length=curved_length,
        curvature=(0.0 if roc is None else 1.0 / roc),
        youngs_modulus=E,
        poisson_ratio=NU,
    )


def analytic_arc(kappa, s):
    """Tip of a tube with pre-curvature about +x, released from the origin."""
    if kappa == 0.0:
        return np.array([0.0, 0.0, s])
    return np.array([0.0, (np.cos(kappa * s) - 1.0) / kappa, np.sin(kappa * s) / kappa])


class TestTubeGeometry(unittest.TestCase):
    """Cross-section properties for the two CT-SDR nitinol tubes."""

    def test_outer_tube_bore(self):
        t = make_tube("outer", 3.6e-3, 0.25e-3, 0.178, 0.07854, 0.050)
        self.assertAlmostEqual(t.inner_diameter * 1e3, 3.1, places=9)
        self.assertAlmostEqual(t.wall_thickness * 1e3, 0.25, places=9)

    def test_inner_tube_bore(self):
        t = make_tube("inner", 2.6e-3, 0.2e-3, 0.308, 0.07854, 0.050)
        self.assertAlmostEqual(t.inner_diameter * 1e3, 2.2, places=9)

    def test_second_moment_matches_annulus_formula(self):
        t = make_tube("outer", 3.6e-3, 0.25e-3, 0.178, 0.07854, 0.050)
        expected = (np.pi / 64.0) * ((3.6e-3) ** 4 - (3.1e-3) ** 4)
        self.assertAlmostEqual(t.second_moment, expected, places=18)
        self.assertAlmostEqual(t.polar_moment, 2.0 * expected, places=18)

    def test_shear_modulus(self):
        t = make_tube("inner", 2.6e-3, 0.2e-3, 0.308, 0.07854, 0.050)
        self.assertAlmostEqual(t.shear_modulus, E_NITINOL / (2.0 * (1.0 + NU)), places=3)

    def test_precurvature_is_distal_only(self):
        t = make_tube("outer", 3.6e-3, 0.25e-3, 0.178, 0.07854, 0.050)
        self.assertEqual(t.precurvature(0.0), (0.0, 0.0))
        self.assertEqual(t.precurvature(0.05), (0.0, 0.0))
        self.assertAlmostEqual(t.precurvature(0.178)[0], 20.0, places=9)
        self.assertAlmostEqual(t.precurvature(0.15)[0], 20.0, places=9)
        self.assertEqual(t.precurvature(0.2), (0.0, 0.0))  # past the tip

    def test_rejects_impossible_geometry(self):
        with self.assertRaises(ValueError):
            make_tube("bad", 2.0e-3, 1.5e-3, 0.1, 0.05, 0.05)  # wall > radius
        with self.assertRaises(ValueError):
            make_tube("bad", 3.0e-3, 0.2e-3, 0.1, 0.2, 0.05)  # curved > length


class TestStraightTipLeadIn(unittest.TestCase):
    """A straight run between the curved section and the tube's distal tip.

    Identified from the NDI translation segments for this robot's outer tube
    (see modeling/VALIDATION_REPORT.md).  Because the emerged part of a tube is
    its distal part, a lead-in puts the curvature *proximal* -- near the guide
    exit -- with the straight run beyond it, which is the opposite of a
    distally-curved tube and is what the measured data shows.
    """

    def setUp(self):
        self.kappa = 20.0
        self.S = 0.010
        self.tube = Tube(
            name="leadin", outer_diameter=2.6e-3, inner_diameter=2.2e-3,
            length=0.308, curved_length=0.07854, curvature=self.kappa,
            youngs_modulus=E_NITINOL, poisson_ratio=NU, tip_straight_length=self.S,
        )
        self.model = CosseratModel([self.tube])

    def _deploy(self, d):
        return self.model.solve([-self.tube.length + d], [0.0])

    def test_zero_lead_in_matches_the_plain_tube(self):
        plain = make_tube("plain", 2.6e-3, 0.2e-3, 0.308, 0.07854, 0.050)
        self.assertEqual(plain.tip_straight_length, 0.0)
        self.assertAlmostEqual(plain.curve_end, plain.length, places=12)
        self.assertAlmostEqual(plain.deployable_length, plain.curved_length, places=12)

    def test_deploying_less_than_the_lead_in_stays_straight(self):
        for d in (0.003, 0.006, 0.0099):
            sol = self._deploy(d)
            np.testing.assert_allclose(sol.tip_position, [0.0, 0.0, d], atol=1e-9)

    def test_beyond_the_lead_in_it_is_an_arc_then_a_tangent(self):
        d = 0.020
        arc = d - self.S
        sol = self._deploy(d)
        end = analytic_arc(self.kappa, arc)
        tangent = np.array([0.0, -np.sin(self.kappa * arc), np.cos(self.kappa * arc)])
        np.testing.assert_allclose(sol.tip_position, end + self.S * tangent, atol=1e-9)

    def test_the_curvature_sits_proximal_not_distal(self):
        sol = self._deploy(0.030)
        mag = np.linalg.norm(sol.u[:, :2], axis=1)
        near = mag[sol.s < 0.030 - self.S - 1e-6]
        far = mag[sol.s > 0.030 - self.S + 1e-6]
        np.testing.assert_allclose(near, self.kappa, atol=1e-9)
        np.testing.assert_allclose(far, 0.0, atol=1e-9)

    def test_deployable_length_includes_the_lead_in(self):
        self.assertAlmostEqual(self.tube.deployable_length, 0.07854 + self.S, places=12)

    def test_lead_in_longer_than_the_tube_is_rejected(self):
        with self.assertRaises(ValueError):
            Tube(name="bad", outer_diameter=2.6e-3, inner_diameter=2.2e-3,
                 length=0.05, curved_length=0.04, curvature=20.0,
                 youngs_modulus=E_NITINOL, poisson_ratio=NU,
                 tip_straight_length=0.02)


class TestSingleTubeArc(unittest.TestCase):
    """A lone pre-curved tube must trace an exact circular arc."""

    def setUp(self):
        self.roc = 0.050
        self.tube = make_tube("solo", 2.6e-3, 0.2e-3, 0.308, 0.07854, self.roc)
        self.model = CosseratModel([self.tube])

    def _deploy(self, d, alpha=0.0):
        return self.model.solve([-self.tube.length + d], [alpha])

    def test_tip_matches_circle(self):
        for d in (0.005, 0.020, 0.035, 0.070):
            sol = self._deploy(d)
            want = analytic_arc(1.0 / self.roc, d)
            np.testing.assert_allclose(sol.tip_position, want, atol=1e-9)

    def test_whole_backbone_lies_on_the_circle(self):
        sol = self._deploy(0.070)
        radial = np.linalg.norm(sol.p - np.array([0.0, -self.roc, 0.0]), axis=1)
        np.testing.assert_allclose(radial, self.roc, atol=1e-9)

    def test_backbone_is_inextensible(self):
        """|p'| = 1 exactly; the sampled polyline may only lose O(h^2) chord."""
        sol = self._deploy(0.070)
        np.testing.assert_allclose(np.linalg.norm(sol.R[:, :, 2], axis=1), 1.0, atol=1e-9)
        path = np.sum(np.linalg.norm(np.diff(sol.p, axis=0), axis=1))
        deficit = (sol.arc_length - path) / sol.arc_length
        self.assertGreaterEqual(deficit, 0.0)
        # Chord vs arc for a circle: 1 - sinc(dphi/2) ~ dphi^2 / 24.
        dphi = sol.arc_length / self.roc / (sol.s.size - 1)
        self.assertLess(deficit, 2.0 * dphi ** 2 / 24.0)

    def test_no_torsion_without_a_second_tube(self):
        sol = self._deploy(0.070, alpha=np.radians(137.0))
        np.testing.assert_allclose(sol.u_z[:, 0], 0.0, atol=1e-12)

    def test_axial_rotation_rotates_the_arc(self):
        d, phi = 0.050, np.radians(137.0)
        base = self._deploy(d).tip_position
        rolled = self._deploy(d, alpha=phi).tip_position
        c, s = np.cos(phi), np.sin(phi)
        want = np.array([c * base[0] - s * base[1], s * base[0] + c * base[1], base[2]])
        np.testing.assert_allclose(rolled, want, atol=1e-9)

    def test_tip_tangent_matches_the_arc_angle(self):
        d = 0.060
        sol = self._deploy(d)
        theta = d / self.roc
        np.testing.assert_allclose(
            sol.tip_tangent, [0.0, -np.sin(theta), np.cos(theta)], atol=1e-9
        )


class TestTwoTubeSuperposition(unittest.TestCase):
    """Bending compatibility: tubes share one curvature, stiffness-weighted."""

    def setUp(self):
        self.roc = 0.050
        self.a = make_tube("a", 3.6e-3, 0.25e-3, 0.178, 0.07854, self.roc)
        self.b = make_tube("b", 2.6e-3, 0.2e-3, 0.308, 0.07854, self.roc)

    def _solve(self, tubes, d, alphas):
        model = CosseratModel(tubes)
        betas = [-t.length + d for t in tubes]
        return model.solve(betas, alphas)

    def test_identical_tubes_aligned_behave_as_one(self):
        twin = [self.b, self.b]
        d = 0.050
        both = self._solve(twin, d, [0.0, 0.0])
        np.testing.assert_allclose(both.tip_position, analytic_arc(1.0 / self.roc, d), atol=1e-9)

    def test_identical_tubes_opposed_cancel_to_a_straight_line(self):
        twin = [self.b, self.b]
        d = 0.050
        sol = self._solve(twin, d, [0.0, np.pi])
        np.testing.assert_allclose(sol.tip_position, [0.0, 0.0, d], atol=1e-9)
        np.testing.assert_allclose(sol.u[:, :2], 0.0, atol=1e-9)

    def test_dissimilar_tubes_opposed_give_stiffness_weighted_curvature(self):
        d = 0.035
        sol = self._solve([self.a, self.b], d, [0.0, np.pi])
        ka, kb = self.a.bending_stiffness, self.b.bending_stiffness
        kappa = (ka - kb) * (1.0 / self.roc) / (ka + kb)
        np.testing.assert_allclose(sol.tip_position, analytic_arc(kappa, d), atol=1e-8)

    def test_dissimilar_tubes_aligned_give_the_full_curvature(self):
        d = 0.035
        sol = self._solve([self.a, self.b], d, [0.0, 0.0])
        np.testing.assert_allclose(sol.tip_position, analytic_arc(1.0 / self.roc, d), atol=1e-9)

    def test_curvature_magnitude_is_bounded_by_the_precurvature(self):
        model = CosseratModel([self.a, self.b])
        for itr in np.radians([0, 30, 60, 90, 120, 150]):
            sol = model.solve([-self.a.length + 0.035, -self.b.length + 0.035], [0.0, itr])
            mag = np.linalg.norm(sol.u[:, :2], axis=1)
            self.assertLessEqual(mag.max(), 1.0 / self.roc + 1e-9)


class TestTorsion(unittest.TestCase):
    """Torsional compliance: the paper's u_z ODE and its boundary conditions."""

    def setUp(self):
        self.robot = CTSDR.from_yaml()

    def test_free_tube_ends_carry_no_torsional_moment(self):
        sol = self.robot.solve(Joints(ott=35, itt=35, itr=90))
        # Outer tube ends partway along; inner tube ends at the distal tip.
        outer_end = np.nanargmax(np.where(np.isnan(sol.u_z[:, OUTER]), -np.inf, sol.s))
        self.assertLess(abs(sol.u_z[outer_end, OUTER]), 1e-6)
        self.assertLess(abs(sol.u_z[-1, INNER]), 1e-6)

    def test_aligned_tubes_do_not_wind_up(self):
        sol = self.robot.solve(Joints(ott=35, itt=35, itr=0))
        np.testing.assert_allclose(np.nan_to_num(sol.u_z), 0.0, atol=1e-9)

    def test_windup_lags_the_commanded_rotation(self):
        """Elastic windup means the tip twists less than the actuator does.

        Compared as magnitudes: the shipped config carries a negative rotation
        sense (identified from the NDI trials), so a positive ITR command
        produces a negative roll.  What the physics fixes is that the roll
        delivered is smaller than the roll commanded, not its handedness.
        """
        for itr in (30.0, 60.0, 90.0):
            sol = self.robot.solve(Joints(ott=35, itt=35, itr=itr))
            tip_rel = np.degrees(sol.theta[-1, INNER] - sol.theta[-1, OUTER])
            commanded = np.degrees(sol.alphas[INNER] - sol.alphas[OUTER])
            self.assertGreater(abs(tip_rel), 0.0)
            self.assertLess(abs(tip_rel), abs(commanded))
            self.assertEqual(np.sign(tip_rel), np.sign(commanded))

    def test_torsionally_rigid_limit_recovers_the_commanded_angle(self):
        """As GJ -> infinity the tubes stop winding and theta tracks alpha."""
        stiff = [
            Tube(
                name=t.name,
                outer_diameter=t.outer_diameter,
                inner_diameter=t.inner_diameter,
                length=t.length,
                curved_length=t.curved_length,
                curvature=t.curvature,
                youngs_modulus=t.youngs_modulus,
                poisson_ratio=-0.999999,  # G = E / (2(1+nu)) -> huge
            )
            for t in self.robot.tubes
        ]
        model = CosseratModel(stiff)
        d, itr = 0.035, np.radians(90.0)
        sol = model.solve([-stiff[0].length + d, -stiff[1].length + d], [0.0, itr])
        tip_rel = sol.theta[-1, INNER] - sol.theta[-1, OUTER]
        self.assertAlmostEqual(tip_rel, itr, delta=np.radians(0.05))

        # ...and the shape is then the constant-curvature superposition.
        ka, kb = model.k_b
        kappa = np.linalg.norm(
            [ka * 20.0 + kb * 20.0 * np.cos(itr), kb * 20.0 * np.sin(itr)]
        ) / (ka + kb)
        self.assertAlmostEqual(
            np.linalg.norm(sol.tip_position), np.linalg.norm(analytic_arc(kappa, d)),
            delta=2e-5,
        )

    def test_transmission_windup_dominates_for_the_long_inner_tube(self):
        """Most of the lost rotation is in the retracted length, not the arc.

        The inner tube has ~273 mm of torsionally free transmission behind the
        guide exit at ITT = 35 mm against only 35 mm of deployed length, so the
        model must attribute the bulk of the windup to s < 0.
        """
        sol = self.robot.solve(Joints(ott=35, itt=35, itr=90))
        betas, alphas = self.robot.joints_to_model(Joints(ott=35, itt=35, itr=90))
        theta0_rel = sol.theta[0, INNER] - sol.theta[0, OUTER]
        transmission = np.radians(90.0) - theta0_rel
        deployed = theta0_rel - (sol.theta[-1, INNER] - sol.theta[-1, OUTER])
        self.assertGreater(transmission, deployed)
        # theta_i(0) = alpha_i - beta_i u_i,z(0), the rigid-guide condition.
        np.testing.assert_allclose(
            sol.theta[0], alphas - betas * sol.u_z0, atol=1e-9
        )


class TestAxialMomentBalance(unittest.TestCase):
    """sum_i GJ_i u_i,z = 0 wherever no external axial torque is applied.

    This is an exact identity, not an approximation: the z-component of
    ``m = sum_i R_i K_i (u_i - u_i*)`` in the reference frame is exactly
    ``sum_i k_ti u_i,z``, and an unloaded robot has ``m = 0`` everywhere.  It
    exercises the torsion ODE and the tube-end boundary conditions together,
    so a sign error in either shows up here.
    """

    def setUp(self):
        self.robot = CTSDR.from_yaml()
        self.k_t = np.array([t.torsional_stiffness for t in self.robot.tubes])

    def test_axial_torque_balances_everywhere_when_unloaded(self):
        sol = self.robot.solve(Joints(ott=35, itt=35, itr=90))
        torque = np.nansum(np.nan_to_num(sol.u_z) * self.k_t, axis=1)
        np.testing.assert_allclose(torque, 0.0, atol=1e-9)

    def test_the_two_tubes_carry_equal_and_opposite_torque(self):
        sol = self.robot.solve(Joints(ott=35, itt=35, itr=120))
        t_out = self.k_t[OUTER] * sol.u_z[0, OUTER]
        t_in = self.k_t[INNER] * sol.u_z[0, INNER]
        self.assertGreater(abs(t_out), 1e-3)  # non-trivial
        self.assertAlmostEqual(t_out, -t_in, places=9)

    def test_applied_tip_torque_shows_up_in_the_balance(self):
        M = 0.01  # N.m about the tip tangent
        j = Joints(ott=35, itt=35, itr=90)
        free = self.robot.solve(j)
        load = self.robot.drilling_load(0.0, torque=M, solution=free)
        sol = self.robot.solve(j, load=load, guess=free.u_z0)
        torque = np.nansum(np.nan_to_num(sol.u_z) * self.k_t, axis=1)
        axial = np.einsum("kj,kj->k", sol.R[:, :, 2], sol.m)
        np.testing.assert_allclose(torque, axial, atol=1e-8)

    def test_drill_torque_goes_into_the_tube_it_is_attached_to(self):
        """The other tube's distal end stays free even under drill torque."""
        j = Joints(ott=35, itt=70, itr=90)  # inner tube reaches the tip alone
        free = self.robot.solve(j)
        load = self.robot.drilling_load(0.0, torque=0.01, solution=free)
        sol = self.robot.solve(j, load=load, guess=free.u_z0)
        tangent = sol.R[-1][:, 2]
        expected = float(tangent @ load.tip_moment) / self.k_t[INNER]
        self.assertAlmostEqual(sol.u_z[-1, INNER], expected, places=8)
        self.assertNotAlmostEqual(expected, 0.0, places=4)
        # The outer tube ends earlier and is unloaded there.
        last_outer = np.flatnonzero(~np.isnan(sol.u_z[:, OUTER]))[-1]
        self.assertLess(abs(sol.u_z[last_outer, OUTER]), 1e-6)

    def test_the_carrying_tube_can_be_chosen(self):
        j = Joints(ott=35, itt=35, itr=45)
        free = self.robot.solve(j)
        base = self.robot.drilling_load(0.0, torque=0.01, solution=free)
        outer_load = ExternalLoad(tip_force=base.tip_force, tip_moment=base.tip_moment,
                                  axial_torque_tube=OUTER)
        sol = self.robot.solve(j, load=outer_load, guess=free.u_z0)
        self.assertLess(abs(sol.u_z[-1, INNER]), 1e-6)
        self.assertGreater(abs(sol.u_z[-1, OUTER]), 1e-4)

    def test_choosing_a_retracted_tube_is_rejected(self):
        j = Joints(ott=0, itt=35)
        load = ExternalLoad(tip_moment=[0, 0, 0.01], axial_torque_tube=OUTER)
        with self.assertRaises(ValueError):
            self.robot.solve(j, load=load)


class TestMaterialLimits(unittest.TestCase):
    """The linear-elastic model has to say when it is being extrapolated."""

    def setUp(self):
        self.robot = CTSDR.from_yaml()

    def test_bending_strain_is_curvature_times_outer_radius(self):
        t = self.robot.tubes[OUTER]
        self.assertAlmostEqual(t.bending_strain(20.0), 20.0 * 1.8e-3, places=12)

    def test_shear_stress_is_g_times_twist_times_radius(self):
        t = self.robot.tubes[INNER]
        self.assertAlmostEqual(t.shear_stress(4.0), t.shear_modulus * 4.0 * 1.3e-3, places=6)

    def test_precurvature_alone_already_leaves_the_linear_branch(self):
        """3.6 % surface strain: nitinol is on its superelastic plateau here."""
        sol = self.robot.solve(Joints(ott=35, itt=35, itr=0))
        rows = self.robot.material_limits(sol)
        self.assertAlmostEqual(rows[OUTER]["bending_strain"], 0.036, delta=0.001)

    def test_torque_signs_are_preserved_and_opposite(self):
        sol = self.robot.solve(Joints(ott=35, itt=35, itr=90))
        rows = self.robot.material_limits(sol)
        self.assertAlmostEqual(rows[OUTER]["torque"], -rows[INNER]["torque"], places=9)

    def test_heavily_wound_configurations_are_flagged(self):
        sol = self.robot.solve_path(
            [Joints(ott=50, itt=50, itr=a) for a in np.linspace(0, 340, 35)]
        )[-1]
        rows = self.robot.material_limits(sol, shear_limit=250e6)
        self.assertFalse(rows[INNER]["shear_ok"])

    def test_absent_tubes_report_none(self):
        sol = self.robot.solve(Joints(ott=0, itt=35))
        rows = self.robot.material_limits(sol)
        self.assertIsNone(rows[OUTER])
        self.assertIsNotNone(rows[INNER])


class TestSnapThrough(unittest.TestCase):
    """Branch folds are real behaviour; the solver must survive them."""

    def setUp(self):
        self.robot = CTSDR.from_yaml()

    def test_the_folded_configuration_still_solves(self):
        """At 10 mm deployment the wound branch dies between ITR 190 and 200."""
        sol = self.robot.solve(Joints(ott=10, itt=10, itr=200))
        self.assertLess(sol.residual_norm, 1e-6)

    def test_restart_is_reported_when_the_first_branch_fails(self):
        prev = self.robot.solve(Joints(ott=10, itt=10, itr=190))
        sol = self.robot.solve(Joints(ott=10, itt=10, itr=200), guess=prev.u_z0)
        self.assertGreater(sol.restarts, 0)
        self.assertLess(sol.residual_norm, 1e-6)

    def test_a_full_sweep_never_fails(self):
        for d in (10.0, 30.0):
            sols = self.robot.solve_path(
                [Joints(ott=d, itt=d, itr=a) for a in np.linspace(0, 360, 25)]
            )
            for s in sols:
                self.assertLess(s.residual_norm, 1e-6)

    def test_solutions_found_after_a_restart_are_genuine_roots(self):
        sol = self.robot.solve(Joints(ott=10, itt=10, itr=200))
        k_t = np.array([t.torsional_stiffness for t in self.robot.tubes])
        # A stalled least-squares point would not satisfy the moment balance.
        np.testing.assert_allclose(
            np.nansum(np.nan_to_num(sol.u_z) * k_t, axis=1), 0.0, atol=1e-9
        )


class TestExternalLoads(unittest.TestCase):
    """The externally loaded half of the paper: n, m and the tip wrench."""

    def setUp(self):
        self.L = 0.060
        self.tube = make_tube("straight", 2.6e-3, 0.2e-3, 0.308, 0.0, None)
        self.model = CosseratModel([self.tube])
        self.betas = [-self.tube.length + self.L]

    def test_unloaded_straight_tube_stays_straight(self):
        sol = self.model.solve(self.betas, [0.0])
        np.testing.assert_allclose(sol.tip_position, [0.0, 0.0, self.L], atol=1e-12)
        np.testing.assert_allclose(sol.n, 0.0, atol=1e-14)
        np.testing.assert_allclose(sol.m, 0.0, atol=1e-14)

    def test_small_tip_load_matches_euler_bernoulli_cantilever(self):
        EI = self.tube.bending_stiffness
        for F in (0.05, 0.1, 0.2):
            sol = self.model.solve(self.betas, [0.0], load=ExternalLoad(tip_force=[F, 0, 0]))
            linear = F * self.L ** 3 / (3.0 * EI)
            self.assertAlmostEqual(sol.tip_position[0] / linear, 1.0, delta=0.02)

    def test_tip_moment_bends_the_tube_into_a_circle(self):
        """A pure end moment gives constant curvature M / EI, exactly."""
        EI = self.tube.bending_stiffness
        M = 0.05
        sol = self.model.solve(self.betas, [0.0], load=ExternalLoad(tip_moment=[M, 0, 0]))
        np.testing.assert_allclose(sol.u[:, 0], M / EI, rtol=1e-6)
        np.testing.assert_allclose(sol.tip_position, analytic_arc(M / EI, self.L), atol=1e-8)

    def test_internal_wrench_satisfies_the_tip_boundary_condition(self):
        load = ExternalLoad(tip_force=[0.3, -0.1, 0.05], tip_moment=[0.0, 0.002, 0.0])
        sol = self.model.solve(self.betas, [0.0], load=load)
        np.testing.assert_allclose(sol.n[-1], load.tip_force, atol=1e-8)
        np.testing.assert_allclose(sol.m[-1], load.tip_moment, atol=1e-8)
        # No distributed load, so n is constant and m' = -p' x n.
        np.testing.assert_allclose(
            sol.n, np.broadcast_to(load.tip_force, sol.n.shape), atol=1e-8
        )

    def test_uniform_distributed_force_equals_its_resultant_at_the_base(self):
        w = np.array([0.4, 0.0, 0.0])  # N/m
        sol = self.model.solve(self.betas, [0.0], load=ExternalLoad(force_density=w))
        np.testing.assert_allclose(sol.n[0], w * self.L, atol=1e-7)
        np.testing.assert_allclose(sol.n[-1], 0.0, atol=1e-8)

    def test_distributed_moment_integrates_to_the_equivalent_tip_moment(self):
        """A uniform l over the length must load the base like l*L applied at the tip."""
        l_dens = np.array([0.02, 0.0, 0.0])  # N.m/m
        dist = self.model.solve(self.betas, [0.0], load=ExternalLoad(moment_density=l_dens))
        base_m = dist.m[0]
        np.testing.assert_allclose(base_m, l_dens * self.L, atol=1e-7)

    def test_gravity_deflects_the_ct_sdr_downward(self):
        robot = CTSDR.from_yaml()
        j = Joints(ott=35, itt=35, itr=0)
        free = robot.solve(j)
        sag = robot.solve(j, load=robot.gravity_load(j, gravity=(0, 0, -9.81)))
        # Guide axis is +z, so gravity along -z pulls the tip back toward -z.
        self.assertLess(sag.tip_position[2], free.tip_position[2])
        self.assertLess(abs(sag.tip_position[2] - free.tip_position[2]), 1e-5)

    def test_load_stiffens_predictably_with_deployment(self):
        """A cantilever gets much softer as it grows: delta ~ L^3."""
        F = ExternalLoad(tip_force=[0.05, 0, 0])
        short = self.model.solve([-self.tube.length + 0.030], [0.0], load=F)
        long_ = self.model.solve([-self.tube.length + 0.060], [0.0], load=F)
        self.assertAlmostEqual(long_.tip_position[0] / short.tip_position[0], 8.0, delta=0.25)


class TestCTSDRConfig(unittest.TestCase):
    """The shipped configuration reproduces the CT-SDR design description."""

    def setUp(self):
        self.robot = CTSDR.from_yaml()

    def test_tube_dimensions(self):
        outer, inner = self.robot.tubes
        self.assertAlmostEqual(outer.outer_diameter * 1e3, 3.6, places=6)
        self.assertAlmostEqual(outer.wall_thickness * 1e3, 0.25, places=6)
        self.assertAlmostEqual(outer.length * 1e3, 178.0, places=6)
        self.assertAlmostEqual(inner.outer_diameter * 1e3, 2.6, places=6)
        self.assertAlmostEqual(inner.wall_thickness * 1e3, 0.2, places=6)
        self.assertAlmostEqual(inner.length * 1e3, 308.0, places=6)

    def test_shipped_rotation_sense_is_negative(self):
        """Identified from the NDI trials, not a convention chosen for taste.

        With a right-handed sense the model walks set2's rotation circle
        backwards and scores 16.1 mm against 1.5 mm; a rigid registration uses
        proper rotations only, so it cannot absorb that.  Locked in a test so a
        future config edit has to confront the evidence.
        """
        np.testing.assert_allclose(self.robot.signs, [1.0, 1.0, -1.0, -1.0])

    def test_radius_of_curvature_is_50_mm(self):
        for t in self.robot.tubes:
            self.assertAlmostEqual(1.0 / t.curvature * 1e3, 50.0, places=6)

    def test_inner_tube_fits_inside_the_outer_bore(self):
        outer, inner = self.robot.tubes
        self.assertLess(inner.outer_diameter, outer.inner_diameter)

    def test_outer_tube_is_the_stiffer_one(self):
        outer, inner = self.robot.tubes
        self.assertGreater(outer.bending_stiffness, inner.bending_stiffness)

    def test_zero_joints_deploy_nothing(self):
        np.testing.assert_allclose(self.robot.deployed_lengths(Joints()), 0.0)
        with self.assertRaises(ValueError):
            self.robot.solve(Joints())

    def test_translation_maps_one_to_one_to_deployed_length(self):
        d = self.robot.deployed_lengths(Joints(ott=20, itt=35))
        np.testing.assert_allclose(d * 1e3, [20.0, 35.0], atol=1e-9)

    def test_home_offset_shifts_deployment_without_changing_the_traced_arc(self):
        """A non-zero home offset adds deployed length but keeps the same circle.

        The NDI validation identifies ~2.9 mm of inner tube already protruding
        at ITT = 0.  It has to move the tip along the tube's arc without
        changing that arc's radius, which is exactly why the advances are blind
        to it and the rotation circles are not.
        """
        outer, inner = self.robot.tubes
        shifted = CTSDR(outer, inner,
                        base_offsets=[-outer.length + 2.9e-3,
                                      -inner.length + 2.9e-3])
        np.testing.assert_allclose(
            shifted.deployed_lengths(Joints(ott=20, itt=20)) * 1e3, [22.9, 22.9],
            atol=1e-9)

        d = 0.020
        base = self.robot.solve(Joints(ott=20, itt=20))
        moved = shifted.solve(Joints(ott=20, itt=20))
        # Same circle: both tips sit one radius from the same centre.
        roc = 1.0 / inner.curvature
        for sol in (base, moved):
            self.assertAlmostEqual(
                np.linalg.norm(sol.tip_position - np.array([0.0, -roc, 0.0])),
                roc, delta=1e-9)
        # But further along it, by the offset.
        self.assertAlmostEqual(moved.arc_length - base.arc_length, 2.9e-3, places=9)

    def test_config_base_offset_is_read(self):
        import tempfile

        import yaml as _yaml
        with open(DEFAULT_CONFIG) as fh:
            cfg = _yaml.safe_load(fh)
        cfg["tubes"]["inner"]["base_offset_mm"] = -300.0
        with tempfile.NamedTemporaryFile("w", suffix=".yaml", delete=False) as fh:
            _yaml.safe_dump(cfg, fh)
            path = fh.name
        try:
            robot = CTSDR.from_yaml(path)
            self.assertAlmostEqual(robot.base_offsets[INNER], -0.300, places=9)
            np.testing.assert_allclose(
                robot.deployed_lengths(Joints(itt=10))[INNER] * 1e3, 18.0, atol=1e-9)
        finally:
            os.unlink(path)

    def test_rotation_signs_are_applied(self):
        """config/ct_sdr.yaml carries a negative rotation sense; honour it."""
        outer, inner = self.robot.tubes
        flipped = CTSDR(outer, inner, base_offsets=self.robot.base_offsets,
                        signs=(1.0, 1.0, -1.0, -1.0))
        _, alphas = flipped.joints_to_model(Joints(ott=35, itt=35, otr=30, itr=90))
        np.testing.assert_allclose(np.degrees(alphas), [-30.0, -90.0], atol=1e-9)

    def test_check_rejects_over_deployment(self):
        with self.assertRaises(ValueError):
            self.robot.check(Joints(ott=100, itt=100))
        self.robot.check(Joints(ott=35, itt=70))  # inside the 78.54 mm arc


class TestSegmentation(unittest.TestCase):
    """Tube starts, tube ends and curvature transitions all become breakpoints."""

    def setUp(self):
        self.robot = CTSDR.from_yaml()
        self.model = self.robot.model

    def test_only_the_inner_tube_is_active_beyond_the_outer_tip(self):
        betas, _ = self.robot.joints_to_model(Joints(ott=20, itt=35))
        segs, L = self.model.segments(betas)
        self.assertAlmostEqual(L, 0.035, places=9)
        self.assertEqual(segs[0][2], (OUTER, INNER))
        self.assertEqual(segs[-1][2], (INNER,))
        self.assertAlmostEqual(segs[-1][0], 0.020, places=9)

    def test_segments_tile_the_domain(self):
        betas, _ = self.robot.joints_to_model(Joints(ott=20, itt=35))
        segs, L = self.model.segments(betas)
        self.assertAlmostEqual(segs[0][0], 0.0, places=12)
        self.assertAlmostEqual(segs[-1][1], L, places=12)
        for (_, b), (a, _) in zip([(s[0], s[1]) for s in segs[:-1]],
                                  [(s[0], s[1]) for s in segs[1:]]):
            self.assertAlmostEqual(a, b, places=12)

    def test_retracted_robot_is_rejected(self):
        with self.assertRaises(ValueError):
            self.model.segments([-0.178, -0.308])

    def test_solution_marks_absent_tubes_as_nan(self):
        sol = self.robot.solve(Joints(ott=20, itt=35))
        beyond = sol.s > 0.020 + 1e-9
        self.assertTrue(np.all(np.isnan(sol.theta[beyond, OUTER])))
        self.assertFalse(np.any(np.isnan(sol.theta[:, INNER])))


class TestContinuation(unittest.TestCase):
    """Warm starting must not change the answer, only how it is found."""

    def setUp(self):
        self.robot = CTSDR.from_yaml()

    def test_warm_start_reproduces_the_cold_solution(self):
        j = Joints(ott=35, itt=35, itr=60)
        cold = self.robot.solve(j)
        seed = self.robot.solve(Joints(ott=35, itt=35, itr=55))
        warm = self.robot.solve(j, guess=seed.u_z0)
        np.testing.assert_allclose(warm.tip_position, cold.tip_position, atol=1e-9)

    def test_solve_path_is_continuous(self):
        path = [Joints(ott=35, itt=35, itr=a) for a in np.linspace(0, 120, 13)]
        sols = self.robot.solve_path(path)
        tips = np.array([s.tip_position for s in sols])
        steps = np.linalg.norm(np.diff(tips, axis=0), axis=1)
        self.assertLess(steps.max(), 2e-3)  # no branch jumps
        for s in sols:
            self.assertLess(s.residual_norm, 1e-6)

    def test_full_turn_is_the_same_configuration(self):
        a = self.robot.solve(Joints(ott=35, itt=35, itr=45))
        b = self.robot.solve(Joints(ott=35, itt=35, itr=45 + 360))
        np.testing.assert_allclose(a.tip_position, b.tip_position, atol=1e-9)


class TestExperimentConfigurations(unittest.TestCase):
    """The joint values actually commanded in the ICRA2027 free-space tests."""

    def setUp(self):
        self.robot = CTSDR.from_yaml()

    def test_all_recorded_configurations_solve(self):
        cases = [
            Joints(ott=20, itt=20, otr=0, itr=0),
            Joints(ott=0, itt=35, otr=0, itr=360),
            Joints(ott=35, itt=35, otr=360, itr=0),
            Joints(ott=35, itt=35, otr=360, itr=360),
            Joints(ott=35, itt=70, otr=0, itr=180),
            Joints(ott=17.5, itt=35, otr=0, itr=180),
        ]
        for j in cases:
            self.robot.check(j)
            sol = self.robot.solve(j)
            self.assertLess(sol.residual_norm, 1e-6, msg=str(j))
            self.assertTrue(np.all(np.isfinite(sol.p)), msg=str(j))
            # Deployed backbone length equals the further of the two tubes.
            self.assertAlmostEqual(
                sol.arc_length, self.robot.deployed_lengths(j).max(), delta=1e-9
            )

    def test_common_rotation_of_both_tubes_rotates_the_whole_shape(self):
        """Co-rotating both tubes is a rigid roll about the guide axis.

        This is the kinematics `set3b` exercises, and the one that pins the
        rotation sign in the NDI validation, so the expected roll is taken from
        the configured sense rather than assumed right-handed.
        """
        base = self.robot.solve(Joints(ott=35, itt=35, itr=90))
        rolled = self.robot.solve(Joints(ott=35, itt=35, otr=70, itr=90 + 70))
        phi = float(rolled.alphas[OUTER] - base.alphas[OUTER])
        c, s = np.cos(phi), np.sin(phi)
        Rz = np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])
        np.testing.assert_allclose(rolled.tip_position, Rz @ base.tip_position, atol=1e-8)

    def test_shape_is_planar_when_the_tubes_are_aligned(self):
        sol = self.robot.solve(Joints(ott=35, itt=35, otr=0, itr=0))
        np.testing.assert_allclose(sol.p[:, 0], 0.0, atol=1e-12)

    def test_advancing_the_inner_tube_extends_past_the_outer_one(self):
        sol = self.robot.solve(Joints(ott=35, itt=70, itr=180))
        self.assertAlmostEqual(sol.arc_length, 0.070, places=9)
        self.assertGreater(np.linalg.norm(sol.tip_position), 0.030)


class TestFailureModes(unittest.TestCase):

    def test_shooting_error_is_raised_not_silently_returned(self):
        robot = CTSDR.from_yaml()
        with self.assertRaises(ShootingError):
            robot.solve(Joints(ott=35, itt=35, itr=90), max_nfev=1, tol=1e-12)

    def test_non_strict_returns_the_unconverged_solution(self):
        robot = CTSDR.from_yaml()
        sol = robot.solve(
            Joints(ott=35, itt=35, itr=90), max_nfev=1, tol=1e-12, strict=False
        )
        self.assertTrue(np.all(np.isfinite(sol.p)))

    def test_alphas_length_is_validated(self):
        robot = CTSDR.from_yaml()
        with self.assertRaises(ValueError):
            robot.model.solve([-0.14, -0.27], [0.0])


if __name__ == "__main__":
    unittest.main(verbosity=2)
