"""Analytic and missing-phase cases independent of study data."""

import unittest

import numpy as np
import pandas as pd

from cpet_skeleton import PhaseMarkers, compare_records, compare_trajectories, infer_recovery_anchor, transform_record


def ramp(values, workloads=None, **kwargs):
    n = len(values)
    record = pd.DataFrame({"time_s": np.arange(n),
                           "workload_w": np.arange(n) if workloads is None else workloads,
                           "signal": values})
    return transform_record(record, PhaseMarkers(0, n - 1, 100), channels=["signal"], **kwargs)


class MethodContractTests(unittest.TestCase):
    def test_single_arithmetic_mean_pass_with_truncated_edges(self):
        # A median would return zeros; a second mean pass would spread the impulse further.
        result = ramp([0, 0, 30, 0, 0, 0, 0]).ramp["signal"]
        np.testing.assert_allclose(result.values[:7], [10, 7.5, 6, 6, 6, 0, 0])
        self.assertFalse(result.defined[7:].any())

    def test_missing_centre_is_restored_then_internal_bin_is_interpolated(self):
        result = ramp([0, 0, np.nan, 0, 50]).ramp["signal"]
        # Smoothed observed values: 0, 0, NA, 50/3, 25. Interpolated centre: 25/3.
        self.assertFalse(result.observed_bin[2])
        self.assertTrue(result.interpolated[2])
        self.assertAlmostEqual(result.values[2], 25 / 3)
        edges = ramp([np.nan, 10, 10, 10, np.nan]).ramp["signal"]
        self.assertTrue(np.isnan(edges.values[[0, 4]]).all())

    def test_mean_is_confined_to_phase_before_windowing(self):
        record = pd.DataFrame({"time_s": [0, 1, 2, 3, 4, 5, 6, 7],
                               "workload_w": [20, 60, 100, 50, 10, 10, 10, 10],
                               "signal": [1, 1, 1, 9999, 7, 7, 7, 9999]})
        result = transform_record(record, PhaseMarkers(0, 2, 100), channels=["signal"],
                                  recovery_horizon_s=2, recovery_grid=np.arange(3.))
        np.testing.assert_allclose(result.ramp["signal"].values[20:], 1)
        np.testing.assert_allclose(result.recovery["signal"].values, 7)
        self.assertEqual(result.markers.recovery_anchor_s, 4)

    def test_median_per_bin_remains_distinct(self):
        result = ramp([0, 0, 90, 10], workloads=[20, 20, 20, 100], rolling_window=1)
        self.assertEqual(result.ramp["signal"].values[20], 0)
        self.assertEqual(result.ramp["signal"].values[60], 5)

    def test_workload_downsteps_use_cumulative_envelope(self):
        result = ramp([1, 2, 8, 10], workloads=[20, 60, 40, 100], rolling_window=1)
        self.assertEqual(result.ramp["signal"].values[60], 5)
        self.assertFalse(result.ramp["signal"].observed_bin[40])

    def test_original_units_are_retained(self):
        result = ramp([700, 700, 700], workloads=[20, 60, 100])
        np.testing.assert_allclose(result.ramp["signal"].values[20:], 700)

    def test_linear_truth_on_internal_grid(self):
        result = ramp([45, 125, 205], workloads=[20, 60, 100], rolling_window=1)
        trajectory = result.ramp["signal"]
        np.testing.assert_allclose(trajectory.values[20:], 2 * trajectory.grid[20:] + 5)
        self.assertFalse(trajectory.defined[:20].any())

    def test_one_occupied_bin_is_unavailable(self):
        result = ramp([1, 2, 3], workloads=[40, 40, 40])
        self.assertFalse(result.ramp["signal"].defined.any())
        self.assertFalse(result.ramp["signal"].observed_bin.any())

    def test_all_missing_channel_is_unavailable(self):
        result = ramp([np.nan] * 3, workloads=[20, 60, 100])
        self.assertTrue(np.isnan(result.ramp["signal"].values).all())

    def test_rounding_ties_to_even_and_native_domain_exclusion(self):
        result = ramp([10, 30, 999], workloads=[.5, 2.5, 100.4], rolling_window=1)
        trajectory = result.ramp["signal"]
        np.testing.assert_allclose(trajectory.values[:3], [10, 20, 30])
        self.assertTrue(trajectory.observed_bin[[0, 2]].all())
        self.assertFalse(trajectory.defined[3:].any())

    def test_anchor_is_strictly_post_peak_without_fallback(self):
        record = pd.DataFrame({"time_s": [0, 1, 2], "workload_w": [10, 50, 10]})
        self.assertEqual(infer_recovery_anchor(record, 0), 2)
        self.assertIsNone(infer_recovery_anchor(record, 2))

    def test_missing_phase_does_not_erase_other_records_recovery(self):
        full = pd.DataFrame({"time_s": [0, 5, 10, 12, 17, 22],
                             "workload_w": [20, 60, 100, 10, 10, 10],
                             "signal": [1, 2, 3, 3, 2, 1]})
        markers = PhaseMarkers(0, 10, 100)
        available = transform_record(full, markers, channels=["signal"])
        before = available.recovery["signal"].values.copy()
        missing = transform_record(full.iloc[:3], markers, channels=["signal"])
        self.assertIsNone(missing.markers.recovery_anchor_s)
        self.assertTrue(missing.ramp["signal"].defined.any())
        self.assertFalse(missing.recovery["signal"].defined.any())
        pair = compare_records(missing, available, phase="recovery", channel="signal")
        self.assertFalse(pair.eligible)
        self.assertEqual(pair.n_supported, 0)
        self.assertTrue(np.isnan(pair.rms))
        np.testing.assert_equal(available.recovery["signal"].values, before)

    def test_explicit_protocol_anchor_and_explicit_absence(self):
        record = pd.DataFrame({"time_s": [0, 5, 10, 12, 17],
                               "workload_w": [20, 60, 100, 20, 20], "signal": [1] * 5})
        explicit = transform_record(record, PhaseMarkers(0, 10, 100, 12), channels=["signal"])
        self.assertEqual(explicit.recovery["signal"].defined.sum(), 2)
        inferred = transform_record(record, PhaseMarkers(0, 10, 100), channels=["signal"], recovery_threshold_w=20)
        self.assertEqual(inferred.markers.recovery_anchor_s, 12)
        absent = transform_record(record, PhaseMarkers(0, 10, 100), channels=["signal"],
                                  recovery_threshold_w=20, infer_recovery_if_missing=False)
        self.assertFalse(absent.recovery["signal"].defined.any())

    def test_default_pair_support_thresholds_and_known_offset(self):
        first = ramp([7, 7, 7], workloads=[0, 50, 100])
        second = ramp([10, 10, 10], workloads=[0, 50, 100])
        pair = compare_records(first, second, phase="ramp", channel="signal")
        self.assertEqual(pair.n_minimum, 21)
        self.assertEqual(pair.rms, 3)
        recovery = compare_records(first, second, phase="recovery", channel="signal")
        self.assertEqual(recovery.n_minimum, 8)

    def test_input_not_modified_and_rows_sorted(self):
        record = pd.DataFrame({"time_s": [2, 0, 1], "workload_w": [100, 20, 60], "signal": [3, 1, 2]})
        original = record.copy(deep=True)
        result = transform_record(record, PhaseMarkers(0, 2, 100), channels=["signal"], rolling_window=1)
        pd.testing.assert_frame_equal(record, original)
        np.testing.assert_allclose(result.ramp["signal"].values[[20, 60, 100]], [1, 2, 3])

    def test_invalid_grids_markers_and_windows_are_rejected(self):
        for grid in ([0], [0, 2, 3], [1, 2], [0, np.nan], [0, -1], [0, 101]):
            with self.subTest(grid=grid), self.assertRaises(ValueError):
                ramp([1, 2, 3], ramp_grid=np.array(grid))
        for window in (0, 1.5, True):
            with self.subTest(window=window), self.assertRaises(ValueError):
                ramp([1, 2, 3], rolling_window=window)
        with self.assertRaises(ValueError):
            ramp([1, 2, 3], recovery_horizon_s=-1)

    def test_duplicate_and_nonfinite_structure_are_not_silently_deleted(self):
        for times, loads in [([0, 0, 1], [20, 60, 100]), ([0, 1, 2], [20, np.inf, 100])]:
            record = pd.DataFrame({"time_s": times, "workload_w": loads, "signal": [1, 2, 3]})
            with self.subTest(times=times, loads=loads), self.assertRaises(ValueError):
                transform_record(record, PhaseMarkers(0, 2, 100), channels=["signal"])

    def test_empty_pair_cannot_be_made_eligible_by_zero_minimum(self):
        trajectory = ramp([np.nan] * 3).ramp["signal"]
        with self.assertRaises(ValueError):
            compare_trajectories(trajectory, trajectory, minimum_locations=0, support_fraction=0)


if __name__ == "__main__":
    unittest.main()
