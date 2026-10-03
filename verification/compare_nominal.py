"""Compare an installed candidate with a source checkout on synthetic nominal grids.

Example: python verification/compare_nominal.py --baseline-src ../baseline/src
This software regression check is separate from manuscript reconstruction.
"""
from pathlib import Path
import argparse
import json
import os
import subprocess
import sys
import tempfile
import numpy as np


def snapshot(path):
    import pandas as pd
    from cpet_skeleton import PhaseMarkers, compare_records, transform_record
    rng = np.random.default_rng(20261003)
    output, results = {}, []
    for i in range(40):
        n = int(rng.integers(40, 220))
        duration = float(rng.uniform(180, 540))
        t = np.r_[0., np.sort(rng.uniform(0, duration, n - 2)), duration]
        w = np.linspace(0, 100, n)
        w[5::13] -= 5  # Non-monotone measurements exercise the workload envelope.
        tau = np.r_[0., np.sort(rng.uniform(0, 180, 45)), 180.]
        times = np.r_[t, duration + 2 + tau]
        record = pd.DataFrame({'time_s':times, 'workload_w':np.r_[w, np.full(len(tau),10.)]})
        for j in range(3):
            y = rng.normal(300 + 100*j, 30, len(times))
            y[rng.random(len(y)) < .2] = np.nan
            if j == 2 and i % 4 == 0:
                y[:] = np.nan
            record[f'channel{j}'] = y
        if i % 5 == 0:
            record = record.iloc[:n].copy()
        r = transform_record(record, PhaseMarkers(0, duration, 100), channels=['channel0','channel1','channel2'])
        results.append(r)
        for phase in ('ramp','recovery'):
            for channel, c in getattr(r, phase).items():
                for field in ('grid','values','defined','observed_bin','interpolated'):
                    output[f'{i}_{phase}_{channel}_{field}'] = getattr(c, field)
    for i in range(0, len(results), 2):
        for phase in ('ramp','recovery'):
            for channel in ('channel0','channel1','channel2'):
                pair = compare_records(results[i], results[i+1], phase=phase, channel=channel)
                for field in ('pointwise_difference','jointly_defined','rms','eligible','n_supported','n_minimum'):
                    output[f'pair{i}_{phase}_{channel}_{field}'] = np.asarray(getattr(pair, field))
    np.savez(path, **output)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--baseline-src', type=Path)
    p.add_argument('--worker', type=Path, help=argparse.SUPPRESS)
    args = p.parse_args()
    if args.worker:
        snapshot(args.worker)
        return
    if args.baseline_src is None:
        p.error('--baseline-src is required')
    assert (args.baseline_src / 'cpet_skeleton/representation.py').is_file()
    with tempfile.TemporaryDirectory(prefix='cpet-nominal-') as temp:
        outputs = []
        for name in ('baseline', 'candidate'):
            dest = Path(temp) / (name + '.npz')
            env = os.environ.copy()
            env.pop('PYTHONPATH', None)
            if name == 'baseline':
                env['PYTHONPATH'] = str(args.baseline_src.resolve())
            subprocess.run([sys.executable, str(Path(__file__).resolve()), '--worker', str(dest)], env=env, check=True)
            outputs.append(np.load(dest, allow_pickle=False))
        first, second = outputs
        assert set(first.files) == set(second.files)
        for key in first.files:
            np.testing.assert_array_equal(first[key], second[key], err_msg=key)
        count = len(first.files)
        first.close()
        second.close()
    print(json.dumps(dict(records=40, channels=3, phases=2, pairs=20, seed=20261003,
                          arrays_and_scalars_compared=count, exact_match=True,
                          grids={'ramp':'0:1:100','recovery':'0:5:180'}), indent=2))


if __name__ == '__main__':
    main()
