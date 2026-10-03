"""Small executable checks for the public code skeleton."""

from __future__ import annotations

import unittest

import numpy as np
import pandas as pd

from cpet_skeleton import (
    PhaseMarkers,
    SupportedTrajectory,
    compare_trajectories,
    transform_record,
)


class RepresentationTests(unittest.TestCase):
    def test_boundaries_are_not_extrapolated(self) -> None:
        record = pd.DataFrame(
            {
                "time_s": [0, 2, 4, 6, 8, 10, 12, 17, 22],
                "workload_w": [20, 40, 60, 80, 100, 100, 0, 0, 0],
                "VO2": [400, 500, 600, 700, 800, 810, 700, 600, 520],
            }
        )
        represented = transform_record(
            record,
            PhaseMarkers(0.0, 10.0, 100.0, recovery_anchor_s=12.0),
            channels=["VO2"],
            rolling_window=1,
        )

        ramp = represented.ramp["VO2"]
        recovery = represented.recovery["VO2"]
        self.assertFalse(ramp.defined[0])
        self.assertTrue(ramp.defined[20])
        self.assertTrue(ramp.defined[100])
        self.assertTrue(recovery.defined[0])
        self.assertTrue(recovery.defined[10 // 5])
        self.assertFalse(recovery.defined[15 // 5])

    def test_pairwise_rms_uses_only_joint_support(self) -> None:
        grid = np.array([0.0, 1.0, 2.0])
        first = SupportedTrajectory(
            grid,
            np.array([1.0, 2.0, np.nan]),
            np.array([True, True, False]),
            np.array([True, True, False]),
        )
        second = SupportedTrajectory(
            grid,
            np.array([2.0, 4.0, 99.0]),
            np.array([True, True, True]),
            np.array([True, True, True]),
        )
        result = compare_trajectories(
            first, second, support_fraction=0.0, minimum_locations=1
        )

        self.assertEqual(result.n_supported, 2)
        self.assertTrue(np.isnan(result.pointwise_difference[2]))
        self.assertAlmostEqual(result.rms, np.sqrt((1.0**2 + 2.0**2) / 2.0))


if __name__ == "__main__":
    unittest.main()
