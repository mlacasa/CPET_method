"""Four/eight-minute synthetic ramps, equally explicit recovery, and CSV exports.

Run after installing .[tutorial]. No clinical data or study paths are used.
"""
from __future__ import annotations
import argparse
from dataclasses import asdict
from importlib.metadata import version
import json
from pathlib import Path
import shutil
import tempfile
import zipfile

import numpy as np
import pandas as pd
from cpet_skeleton import PhaseMarkers, compare_records, transform_record


def analytic_response(coordinate, phase):
    """Invented VCO2 in mL/min; not fitted physiological ground truth."""
    x = np.asarray(coordinate, dtype=float)
    return 300 + 800 * (x / 100 if phase == 'ramp' else np.exp(-x / 60))


def synthetic_record(duration_s, peak_w, n_ramp, n_recovery, seed):
    rng = np.random.default_rng(seed)

    def irregular_times(duration, count):
        times = np.linspace(0, duration, count)
        times[1:-1] += rng.uniform(-.35, .35, count - 2) * duration / (count - 1)
        return times

    ramp_t = irregular_times(duration_s, n_ramp)
    tau = irregular_times(180, n_recovery)
    anchor = duration_s + 2
    # Start at 20% of the supplied Wpeak, reach Wpeak at the audited peak time.
    workload = peak_w * (.2 + .8 * ramp_t / duration_s)
    workload[-1] = peak_w
    ramp_y = analytic_response(100 * workload / peak_w, 'ramp')
    recovery_y = analytic_response(tau, 'recovery')
    # Deliberate internal and terminal channel gaps; retain the native rows.
    ramp_y[8::17] = np.nan
    recovery_y[8::17] = np.nan
    recovery_y[-2:] = np.nan
    frame = pd.DataFrame({
        'time_s': np.r_[ramp_t, anchor + tau],
        'workload_w': np.r_[workload, np.full(n_recovery, 10.)],
        'VCO2': np.r_[ramp_y, recovery_y],
    })
    return frame, PhaseMarkers(0., float(duration_s), float(peak_w), float(anchor))


def make_cases():
    short = synthetic_record(240, 100, 81, 61, 401)
    long = synthetic_record(480, 160, 161, 91, 801)
    shifted = long[0].copy()
    shifted['VCO2'] += 30  # Preserve coordinates and NaNs; isolate an amplitude change.
    absent = long[0].loc[long[0].time_s <= long[1].peak_time_s].copy()
    absent_markers = PhaseMarkers(0., 480., 160.)
    return {'ramp_4min': short, 'ramp_8min': long,
            'ramp_8min_plus30': (shifted, long[1]),
            'ramp_8min_no_recovery': (absent, absent_markers)}


def transform_cases(cases):
    # There is no paired input or reference record in this loop.
    return {name: transform_record(record, markers, channels=['VCO2'])
            for name, (record, markers) in cases.items()}


def export_results(cases, results, output_dir):
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    rows, support, metadata = [], [], []
    for name, (record, markers) in cases.items():
        record.to_csv(out / f'native_{name}.csv', index=False)
        metadata.append(dict(record=name, **asdict(results[name].markers),
                             ramp_duration_s=markers.peak_time_s - markers.ramp_start_s,
                             ramp_slope_w_per_min=.8 * markers.peak_workload_w /
                             (markers.peak_time_s - markers.ramp_start_s) * 60,
                             ramp_observations=int((record.time_s <= markers.peak_time_s).sum()),
                             recovery_observations=int((record.time_s > markers.peak_time_s).sum())))
        for phase in ('ramp', 'recovery'):
            c = getattr(results[name], phase)['VCO2']
            support.append(dict(record=name, phase=phase, grid_size=len(c.grid),
                                defined=int(c.defined.sum()), observed_bin=int(c.observed_bin.sum()),
                                interpolated=int(c.interpolated.sum())))
            rows.extend(dict(record=name, phase=phase, channel='VCO2', coordinate=float(g),
                             value=float(v), defined=bool(d), observed_bin=bool(o), interpolated=bool(i))
                        for g, v, d, o, i in zip(c.grid, c.values, c.defined, c.observed_bin, c.interpolated))
    pd.DataFrame(rows).to_csv(out / 'trajectories.csv', index=False)
    pd.DataFrame(metadata).to_csv(out / 'markers_and_counts.csv', index=False)
    pd.DataFrame(support).to_csv(out / 'support.csv', index=False)
    # Every matrix uses the same named row order and explicit coordinate columns.
    for phase in ('ramp', 'recovery'):
        grid = getattr(next(iter(results.values())), phase)['VCO2'].grid
        for field in ('values', 'defined', 'observed_bin', 'interpolated'):
            matrix = pd.DataFrame([getattr(getattr(r, phase)['VCO2'], field) for r in results.values()],
                                  index=list(results), columns=grid)
            matrix.index.name = 'record'
            matrix.to_csv(out / f'{phase}_{field}.csv')
    paired, pair_rows = {}, []
    for phase in ('ramp', 'recovery'):
        paired[phase] = {}
        for label, first, second in [('alignment', 'ramp_4min', 'ramp_8min'),
                                     ('known_offset', 'ramp_8min', 'ramp_8min_plus30'),
                                     ('missing_recovery', 'ramp_8min', 'ramp_8min_no_recovery')]:
            pair = compare_records(results[first], results[second], phase=phase, channel='VCO2')
            paired[phase][label] = dict(n_supported=pair.n_supported, n_minimum=pair.n_minimum,
                                       rms=float(pair.rms) if np.isfinite(pair.rms) else None)
            pair_rows.extend(dict(case=label, phase=phase, coordinate=float(g), difference=float(v),
                                  jointly_defined=bool(d))
                             for g, v, d in zip(pair.grid, pair.pointwise_difference, pair.jointly_defined))
    pd.DataFrame(pair_rows).to_csv(out / 'paired_differences.csv', index=False)
    summary = dict(package_version=version('cpet-functional-skeleton'), data='synthetic only',
                   seeds=[401, 801], channel_units={'VCO2':'mL/min'},
                   coordinate_units={'ramp':'%Wpeak', 'recovery':'seconds from anchor'},
                   rolling_window=5, paired=paired)
    (out / 'summary.json').write_text(json.dumps(summary, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    return pd.DataFrame(metadata), pd.DataFrame(support), summary


def plot_results(cases, results, output_dir):
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
    from matplotlib.patches import Patch
    out = Path(output_dir)
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
    for name in ('ramp_4min', 'ramp_8min'):
        record, markers = cases[name]
        for row, phase in enumerate(('ramp', 'recovery')):
            part = record.loc[record.time_s <= markers.peak_time_s] if phase == 'ramp' else record.loc[record.time_s >= markers.recovery_anchor_s]
            native_time = part.time_s - (markers.ramp_start_s if phase == 'ramp' else markers.recovery_anchor_s)
            axes[row, 0].plot(native_time, part.VCO2, '.-', label=name, markersize=3)
            curve = getattr(results[name], phase)['VCO2']
            axes[row, 1].plot(curve.grid, curve.values, label=name)
    for row, phase in enumerate(('ramp', 'recovery')):
        grid = getattr(results['ramp_4min'], phase)['VCO2'].grid
        axes[row, 1].plot(grid, analytic_response(grid, phase), 'k--', label='analytic target')
        for col in (0, 1):
            axes[row, col].set_ylabel('Synthetic VCO2 (mL/min)')
            axes[row, col].set_title(f'{phase.capitalize()}: ' + ('native observations' if col == 0 else 'represented'))
            axes[row, col].legend()
        axes[row, 0].set_xlabel('Seconds from ' + ('ramp onset' if phase == 'ramp' else 'recovery anchor'))
        axes[row, 1].set_xlabel('%Wpeak' if phase == 'ramp' else 'Seconds from recovery anchor')
    fig.savefig(out / 'native_and_represented.png', dpi=140)
    mask_fig, mask_axes = plt.subplots(2, 1, figsize=(12, 5), constrained_layout=True)
    colors = ['#dedede', '#2166ac', '#f4a340']
    for ax, phase in zip(mask_axes, ('ramp', 'recovery')):
        curves = [getattr(r, phase)['VCO2'] for r in results.values()]
        states = np.array([c.observed_bin.astype(int) + 2 * c.interpolated for c in curves])
        ax.imshow(states, interpolation='nearest', aspect='auto', cmap=ListedColormap(colors), vmin=0, vmax=2)
        ax.set_yticks(range(len(results)), labels=list(results))
        ticks = np.linspace(0, len(curves[0].grid)-1, 6, dtype=int)
        ax.set_xticks(ticks, labels=curves[0].grid[ticks])
        ax.set_xlabel('Ramp %Wpeak' if phase == 'ramp' else 'Recovery seconds from anchor')
    mask_axes[0].legend(handles=[Patch(color=c, label=l) for c, l in zip(colors, ['unavailable', 'occupied bin', 'interpolated'])],
                        loc='lower center', bbox_to_anchor=(.5, 1.02), ncol=3)
    mask_fig.savefig(out / 'support_masks.png', dpi=140)
    return fig, mask_fig


def archive_results(output_dir):
    out = Path(output_dir)
    # Build and verify on a local temporary filesystem before copying to cloud drives.
    with tempfile.TemporaryDirectory() as temp:
        archive = Path(shutil.make_archive(str(Path(temp) / 'cpet_synthetic_results'), 'zip', out))
        with zipfile.ZipFile(archive) as handle:
            assert handle.testzip() is None
        target = out.parent / (out.name + '.zip')
        shutil.copyfile(archive, target)
    return target


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, default=Path('unequal_duration_output'))
    args = parser.parse_args()
    cases = make_cases()
    results = transform_cases(cases)
    metadata, support, summary = export_results(cases, results, args.output_dir)
    import matplotlib
    matplotlib.use('Agg')
    plot_results(cases, results, args.output_dir)
    print(metadata.to_string(index=False))
    print(support.to_string(index=False))
    print(json.dumps(summary, indent=2))
    print('Downloadable bundle:', archive_results(args.output_dir))


if __name__ == '__main__':
    main()
