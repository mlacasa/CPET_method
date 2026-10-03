"""Decimal endpoints must retain valid bins without admitting outside events."""
import unittest
import numpy as np
import pandas as pd
from cpet_skeleton import PhaseMarkers, transform_record


class DecimalGridTests(unittest.TestCase):
    def recovery(self, times, values):
        record = pd.DataFrame({
            'time_s': [0., 1., 2., *times],
            'workload_w': [0., 50., 100., *([10.] * len(times))],
            'signal': [1., 2., 3., *values],
        })
        return transform_record(
            record, PhaseMarkers(0, 2, 100, 3), channels=['signal'],
            rolling_window=1, recovery_horizon_s=.3,
            recovery_grid=np.array([0., .1, .2, .3]),
        ).recovery['signal']

    def test_reported_decimal_endpoint_keeps_all_four_bins(self):
        curve = self.recovery([3., 3.1, 3.2, 3.3], [10., 11., 12., 13.])
        np.testing.assert_allclose(curve.values, [10., 11., 12., 13.])
        self.assertTrue(curve.defined.all())
        self.assertTrue(curve.observed_bin.all())
        self.assertFalse(curve.interpolated.any())

    def test_two_valid_bins_including_decimal_endpoint_are_sufficient(self):
        curve = self.recovery([3.1, 3.3], [11., 13.])
        np.testing.assert_allclose(curve.values, [np.nan, 11., 12., 13.])
        np.testing.assert_array_equal(curve.defined, [False, True, True, True])
        np.testing.assert_array_equal(curve.observed_bin, [False, True, False, True])
        np.testing.assert_array_equal(curve.interpolated, [False, False, True, False])

    def test_outside_decimal_domain_cannot_round_into_endpoint(self):
        # Within the selected phase, but the ramp coordinate exceeds the grid.
        record = pd.DataFrame({'time_s': [0, 1, 2],
                               'workload_w': [0., .1, .301], 'signal': [10., 11., 999.]})
        curve = transform_record(
            record, PhaseMarkers(0, 2, 100), channels=['signal'], rolling_window=1,
            ramp_grid=np.array([0., .1, .2, .3]),
        ).ramp['signal']
        np.testing.assert_allclose(curve.values, [10., 11., np.nan, np.nan])
        np.testing.assert_array_equal(curve.observed_bin, [True, True, False, False])

    def test_decimal_ramp_endpoint_and_bin_median(self):
        record = pd.DataFrame({'time_s': [0, 1, 2, 3, 4],
                               'workload_w': [0., .1, .3, .3, .3],
                               'signal': [10., 11., 13., 13., 999.]})
        curve = transform_record(
            record, PhaseMarkers(0, 4, 100), channels=['signal'], rolling_window=1,
            ramp_grid=np.array([0., .1, .2, .3]),
        ).ramp['signal']
        np.testing.assert_allclose(curve.values, [10., 11., 12., 13.])
        np.testing.assert_array_equal(curve.observed_bin, [True, True, False, True])


if __name__ == '__main__':
    unittest.main()
