"""Seeded synthetic CPET example; writes inspectable CSV and JSON outputs."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
from cpet_skeleton import PhaseMarkers, compare_records, transform_record


def synthetic_record(seed: int, response_shift: float = 0.0):
    rng = np.random.default_rng(seed)
    time_s = np.cumsum(rng.uniform(2.0, 4.5, size=240))
    workload = np.where(time_s < 120, 10, 10 + (time_s - 120) * .25)
    workload[time_s > 520] = 10
    peak_index = int(np.argmax(workload))
    peak_time = float(time_s[peak_index])
    peak_workload = float(workload[peak_index])
    vo2 = 350 + 9 * workload
    vco2 = 300 + 8.5 * workload
    recovery = time_s > peak_time
    vo2[recovery] = 350 + 9 * peak_workload * np.exp(-(time_s[recovery] - peak_time) / 60)
    vco2[recovery] = 300 + 8.5 * peak_workload * np.exp(-(time_s[recovery] - peak_time) / 70)
    vo2 += response_shift + rng.normal(0, 25, len(time_s))
    vco2 += response_shift + rng.normal(0, 30, len(time_s))
    vo2[20::23] = np.nan
    vco2[25::29] = np.nan
    frame = pd.DataFrame({'time_s': time_s, 'workload_w': workload, 'VO2': vo2, 'VCO2': vco2})
    return frame, PhaseMarkers(120, peak_time, peak_workload)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, default=Path('example_output'))
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    represented, rows, marker_rows = [], [], []
    for index, (seed, shift) in enumerate([(101, 0), (202, 30)], start=1):
        record, markers = synthetic_record(seed, shift)
        record.to_csv(args.output_dir / f'synthetic_cpet{index}.csv', index=False)
        result = transform_record(record, markers, channels=['VO2', 'VCO2'])
        represented.append(result)
        marker_rows.append(vars(result.markers))
        for phase in ('ramp', 'recovery'):
            for channel, curve in getattr(result, phase).items():
                rows.extend(dict(record=f'synthetic_cpet{index}', phase=phase, channel=channel,
                                 coordinate=float(g), value=float(v), defined=bool(m),
                                 observed_bin=bool(o), interpolated=bool(i))
                            for g, v, m, o, i in zip(curve.grid, curve.values, curve.defined,
                                                     curve.observed_bin, curve.interpolated))
    pd.DataFrame(rows).to_csv(args.output_dir / 'trajectories.csv', index=False)
    summary = {'data': 'synthetic only', 'seed': [101, 202], 'markers': marker_rows, 'paired_VCO2': {}}
    differences = []
    for phase in ('ramp', 'recovery'):
        pair = compare_records(*represented, phase=phase, channel='VCO2')
        summary['paired_VCO2'][phase] = dict(n_supported=pair.n_supported,
                                           n_minimum=pair.n_minimum, rms=pair.rms)
        differences.extend(dict(phase=phase, coordinate=float(g), difference=float(v), jointly_defined=bool(m))
                           for g, v, m in zip(pair.grid, pair.pointwise_difference, pair.jointly_defined))
    pd.DataFrame(differences).to_csv(args.output_dir / 'paired_VCO2.csv', index=False)
    full, markers = synthetic_record(101)
    missing = transform_record(full.loc[full.time_s <= markers.peak_time_s], markers, channels=['VCO2'])
    missing_pair = compare_records(missing, represented[1], phase='recovery', channel='VCO2')
    assert missing.markers.recovery_anchor_s is None
    assert missing_pair.n_supported == 0 and np.isnan(missing_pair.rms)
    assert represented[1].recovery['VCO2'].defined.any()
    summary['missing_recovery_case'] = dict(anchor=None, joint_support=0, rms=None,
                                           other_record_recovery_retained=True)
    text = json.dumps(summary, indent=2, allow_nan=False)
    (args.output_dir / 'summary.json').write_text(text + chr(10), encoding='utf-8')
    print(text)
    print(f'Synthetic inputs, curves, masks and differences saved to {args.output_dir}')


if __name__ == '__main__':
    main()
