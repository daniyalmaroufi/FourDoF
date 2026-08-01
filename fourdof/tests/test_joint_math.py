import os
import sys
import unittest

# Current file is in fourdof/tests/; joint_math.py lives in fourdof/fourdof/.
# joint_math has no rospy/dynamixel_sdk dependency, so this test needs no
# mocking and no ROS/hardware.
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../fourdof')))

from joint_math import (  # noqa: E402
    TranslationalJointConfig,
    rotational_accel_to_profile_accel,
    rotational_deg_to_pulse,
    rotational_profile_velocity_to_speed,
    rotational_pulse_to_deg,
    rotational_speed_to_profile_velocity,
    translational_mm_to_pulse,
    translational_profile_velocity_to_speed,
    translational_pulse_to_mm,
    translational_speed_to_profile_velocity,
    validate_translational_range,
)

GEAR_RATIO = 40.0
PITCH_MM = 2.54


class TestRotationalConversions(unittest.TestCase):

    def test_one_motor_revolution_per_gear_ratio_tube_degrees(self):
        # 9 tube deg * gear ratio 40 = 360 motor deg = exactly one full motor revolution.
        self.assertEqual(rotational_deg_to_pulse(9.0, GEAR_RATIO), 4096)

    def test_pulse_to_deg_round_trip(self):
        for deg in (-720.0, -90.0, 0.0, 45.5, 360.0, 2304.0):
            pulse = rotational_deg_to_pulse(deg, GEAR_RATIO)
            self.assertAlmostEqual(rotational_pulse_to_deg(pulse, GEAR_RATIO), deg, places=2)

    def test_speed_round_trip(self):
        for deg_s in (5.0, 30.0, 90.0):
            pulse = rotational_speed_to_profile_velocity(deg_s, GEAR_RATIO)
            self.assertGreater(pulse, 0)
            self.assertAlmostEqual(
                rotational_profile_velocity_to_speed(pulse, GEAR_RATIO), deg_s, delta=0.5
            )

    def test_accel_never_negative(self):
        self.assertGreaterEqual(rotational_accel_to_profile_accel(0.0, GEAR_RATIO), 0)
        self.assertGreater(rotational_accel_to_profile_accel(50.0, GEAR_RATIO), 0)


class TestTranslationalConversions(unittest.TestCase):

    def test_one_motor_revolution_per_pitch_travel(self):
        # Moving exactly one lead-screw pitch = one full motor revolution.
        self.assertEqual(translational_mm_to_pulse(PITCH_MM, PITCH_MM), 4096)

    def test_pulse_to_mm_round_trip(self):
        for mm in (-50.0, 0.0, 12.7, 100.0):
            pulse = translational_mm_to_pulse(mm, PITCH_MM)
            self.assertAlmostEqual(translational_pulse_to_mm(pulse, PITCH_MM), mm, places=2)

    def test_speed_round_trip(self):
        for mm_s in (1.0, 5.0, 20.0):
            pulse = translational_speed_to_profile_velocity(mm_s, PITCH_MM)
            self.assertGreater(pulse, 0)
            self.assertAlmostEqual(
                translational_profile_velocity_to_speed(pulse, PITCH_MM), mm_s, delta=0.1
            )

    def test_range_validation_accepts_in_range(self):
        cfg = TranslationalJointConfig(motor_id=1, pitch_mm=PITCH_MM, min_mm=0.0, max_mm=100.0)
        validate_translational_range(0.0, cfg)
        validate_translational_range(100.0, cfg)
        validate_translational_range(50.0, cfg)

    def test_range_validation_rejects_out_of_range(self):
        cfg = TranslationalJointConfig(motor_id=1, pitch_mm=PITCH_MM, min_mm=0.0, max_mm=100.0)
        with self.assertRaises(ValueError):
            validate_translational_range(-0.1, cfg)
        with self.assertRaises(ValueError):
            validate_translational_range(100.1, cfg)


if __name__ == '__main__':
    unittest.main()
