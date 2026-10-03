# Synthetic examples and interpretation

## First example: deliberately unequal durations

Use [the tutorial notebook](../notebooks/CPET_method_tutorial.ipynb) or
`python examples/unequal_duration_example.py --output-dir unequal_duration_output`.
Both use the same formulas, seeds and default transformation settings.

| Property | Four-minute record | Eight-minute record |
|---|---:|---:|
| Ramp duration | 240 s | 480 s |
| Native ramp observations | 81 | 161 |
| Ramp start / peak workload | 20 / 100 W | 32 / 160 W |
| Linear workload slope | 20 W/min | 16 W/min |
| Recovery anchor | 242 s | 482 s |
| Recorded recovery duration | 180 s | 180 s |
| Native recovery observations | 61 | 91 |
| Recovery workload | 10 W | 10 W |

Seeds 401 and 801 generate jittered interior timestamps with fixed phase
endpoints. Ramp and recovery have deliberate missing channel positions, with
the final two recovery values missing. Native counts include those rows.
No measurement noise is added. NaNs stay in the native record.

The analytic synthetic VCO2 target (mL/min) is `300 + 800*p/100` on ramp and
`300 + 800*exp(-tau/60)` on recovery. Here `p` is %Wpeak and `tau` is seconds
from recovery anchor. It is a teaching construction, not physiological truth.
At 50% of peak the two records are at 50 and 80 W. Their matching analytic
responses were chosen deliberately and do not establish physiological equivalence.

Three comparisons are kept separate:

1. **Alignment:** four versus eight minutes, with the same known target functions
   and different native counts in both phases. Smoothing and binning of unequal
   samples can produce nonzero differences; identical output curves are not assumed.
2. **Known amplitude difference:** a copy of the eight-minute record gains exactly
   30 mL/min at every finite channel observation; timestamps, workloads, markers
   and NaNs are unchanged. Both phase differences equal 30 on joint support.
3. **Missing recovery:** another copy ends at peak. Its ramp remains available,
   its recovery remains unavailable and other records' recovery is retained.

The nominal five-position mean is used in every case. There is no smoothing
across the phase boundary, scaling of channel amplitudes or imputation of an
absent phase. The default 101/37 candidate grids do not imply full support.

## Inspect and download

| Files | Contents |
|---|---|
| `native_*.csv` | Four invented records, including the two deliberately modified copies |
| `markers_and_counts.csv` | Audited synthetic markers, durations, load slopes and native counts |
| `trajectories.csv` | Phase/channel/coordinate/value and three masks for every record |
| `support.csv` | Available, occupied and interpolated counts per record/phase |
| `ramp_*.csv`, `recovery_*.csv` | `values`, `defined`, `observed_bin`, `interpolated` matrices; named record rows and explicit coordinate columns |
| `paired_differences.csv` | Second-minus-first differences and joint masks for the three illustrations |
| `summary.json` | Package version, seeds, units, settings and paired support/RMS |
| `native_and_represented.png`, `support_masks.png` | Both phases before/after mapping and mask provenance |

The sibling ZIP contains all these outputs. Matrix masks are Boolean, and blank
numerical cells are unavailable values. The notebook enables a Colab download
or supplies a local link. Re-running overwrites the generated files.

## Original seeded example, retained unchanged

`python examples/synthetic_example.py --output-dir example_output` retains seeds
101 and 202, two channels, random amplitude noise and its original near-equal
ramp durations (397.666 and 398.259 s). It is useful as a software regression
case, not as the demonstration of unequal durations. Its five outputs are
`synthetic_cpet1.csv`, `synthetic_cpet2.csv`, `trajectories.csv`, `paired_VCO2.csv`
and `summary.json`.

VCO2 RMS is 32.95382346 on 91 jointly supported ramp positions (minimum 21),
and 33.62799838 on 36 recovery positions (minimum 8). The missing-recovery
case has zero joint support and unavailable RMS while retaining its partner's
recovery. These synthetic values are illustrative, not physiological validation.

Neither example reproduces the [historical reconstruction experiment](../validation/README.md).
