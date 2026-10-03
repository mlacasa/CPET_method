# Input, processing and support contract

## Input contract

Pass one pandas `DataFrame` per record:

| Column | Required meaning |
|---|---|
| `time_s` | Finite elapsed time in seconds, one unique time per native observation |
| `workload_w` | Finite delivered workload in watts |
| Named channels, e.g. `VO2`, `VCO2` | Numerical measurements in their original units; missing values allowed |

Alternative time/workload column names are configurable. Rows are sorted by time;
the input is not modified. Duplicate times and nonfinite structural fields raise
an error: resolve device-specific harmonisation upstream instead of silently
collapsing or deleting breaths. Non-numeric or infinite channel values become
missing. Apply physiological quality-control masks upstream and keep affected
breath rows, so that the five-position windows retain their intended meaning.

Supply audited ramp onset, peak time and positive peak workload. The library does
not estimate these markers or impose a physiological QC policy.

```python
import pandas as pd
from cpet_skeleton import PhaseMarkers, transform_record

record = pd.DataFrame({
    "time_s": [0, 5, 10, 12, 17, 22],
    "workload_w": [20, 60, 100, 10, 10, 10],
    "VO2": [400, 600, 800, 750, 680, 600],
})
represented = transform_record(
    record,
    PhaseMarkers(ramp_start_s=0, peak_time_s=10, peak_workload_w=100),
    channels=["VO2"],
)
curve = represented.ramp["VO2"]
print(curve.grid, curve.values, curve.defined)
print(represented.markers.recovery_anchor_s)  # 12.0
```

## Processing and defaults

1. Select the ramp observations from onset through peak, inclusive. Infer recovery
   as the first observation **strictly after peak** with workload <=10 W, or use
   an explicit audited anchor. Select `[anchor, anchor + 180 s]` before smoothing.
2. Apply one centred, five-observation **arithmetic moving mean**, with stride one
   and at least one finite value per window, separately within each phase/channel.
   Truncate windows at phase edges and restore the original channel missing mask.
3. Map ramp using `100 * cumulative_max(workload) / peak_workload`. Map recovery
   using `time - anchor`. Channel amplitudes retain their original units.
4. Retain native coordinates inside the grid domain before assigning nearest
   centres. Rounding uses ties-to-even, as in NumPy. Summarise each bin with the
   **median** of its available smoothed values; this is distinct from smoothing.
5. Require at least **two occupied bins per phase/channel**. Interpolate linearly
   between occupied centres, with no extrapolation outside their bounds.

| Argument | Default | Meaning |
|---|---|---|
| `rolling_window` | `5` | Positions per mean window; `1` disables smoothing |
| `ramp_grid` | `0, 1, ..., 100` | Percentage points of attained peak workload |
| `recovery_grid` | `0, 5, ..., 180` | Seconds from anchor |
| `recovery_horizon_s` | `180` | Native recovery selection duration |
| `recovery_threshold_w` | `10` | Threshold for automatic transition detection |
| `infer_recovery_if_missing` | `True` | Infer an anchor when its supplied value is `None` |

Custom grids must be finite, uniform, strictly increasing, start at zero and have
at least two points. Ramp cannot exceed 100; recovery cannot exceed its horizon.
If only the horizon changes, the default recovery grid ends at its last 5-second
multiple. A horizon shorter than 5 seconds requires an explicit finer grid.
Changing defaults is a researcher-defined adaptation, not an evaluated alternative
protocol or an automatic reproduction of the manuscript results.

### Recovery absence and protocol adaptation

When no transition is recorded, `recovery_anchor_s` stays `None`; all recovery
values are NaN and all recovery masks are false. There is **no peak-time fallback**
and no paired-record rule. Ramp remains available when supported. An observed
recovery can also lack sufficient channel support; phase availability and
trajectory availability are different concepts.

For another protocol, supply an audited `recovery_anchor_s` strictly after peak
or adapt the threshold explicitly. To mark a phase unavailable without inference,
use `PhaseMarkers(..., recovery_anchor_s=None)` with
`infer_recovery_if_missing=False`. Cross-protocol validity requires separate
evaluation.

## Output and support

`RecordRepresentation` has `ramp` and `recovery` dictionaries keyed by channel,
plus the record's resolved `markers`. Every `SupportedTrajectory` contains:

| Field | Meaning |
|---|---|
| `grid` | Fixed candidate coordinate array |
| `values` | Values in input units; NaN outside defined support |
| `defined` | Boolean availability mask |
| `observed_bin` | Defined positions supplied by occupied-bin summaries |
| `interpolated` | Defined positions estimated between occupied bins |

With fewer than two occupied bins the entire trajectory, including its
`observed_bin` mask, is unavailable. Bin summaries already include smoothing and
aggregation; they are not raw measurements. Internal interpolation may bridge
missing bins, but does not create independent observations. No gap-length limit
is imposed by the manuscript method.

To construct a matrix for a channel and phase, stack the `.values` arrays from
independently processed records, and stack `.defined` in the same record order.
Keep record identifiers and markers alongside the matrices.

## Optional paired comparison

```python
from cpet_skeleton import compare_records

# first and second are independently transformed RecordRepresentation objects.
pair = compare_records(first, second, phase="recovery", channel="VCO2")
print(pair.pointwise_difference)  # second minus first; NaN outside joint support
print(pair.jointly_defined, pair.n_supported, pair.eligible, pair.rms)
```

Both inputs must use identical grids. RMS is computed on joint support only when
`n_supported >= max(minimum_locations, ceil(support_fraction * grid_size))`.
Defaults are 5 locations and 0.20, giving 21/101 for ramp and 8/37 for recovery.
An ineligible RMS is NaN. This helper uses the full supplied grid; it does not
derive or apply the manuscript's cohort-level display restriction. Differences,
RMS, ratios, clustering and other downstream analyses are optional, separate
research choices.


## Interpretation

Ordered correspondence does not establish physiological equivalence, superiority
to alternative coordinates, clinical subtypes or validated classifiers. Workload
normalisation changes the coordinate, not the physiological channel amplitudes.
The output can support differences, ratios, functional analysis or other
researcher-chosen analyses; those choices require their own validation.
