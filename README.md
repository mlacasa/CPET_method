# Protocol-aware functional representation of breath-by-breath CPET

This folder contains a repository-ready code skeleton accompanying the manuscript
*From Irregular Breath-by-Breath Records to Comparable Functional CPET Trajectories: A
Protocol-Aware Phase-Specific Framework*.

The code demonstrates the general methodological contribution: transforming one
irregular CPET record at a time into two supported, comparably indexed trajectories:

- an incremental-ramp trajectory indexed by relative workload progression; and
- a recovery trajectory indexed by time from a record-specific post-exercise anchor.

The paired CPET2--CPET1 difference and RMS calculation are implemented in a separate
module because they are examples of downstream use, not required components of the
record-level transformation.

## What is included

```text
github/
|-- README.md
|-- pyproject.toml
|-- src/cpet_skeleton/
|   |-- __init__.py
|   |-- representation.py   # General record-level transformation
|   `-- paired.py           # Optional paired difference and RMS example
|-- examples/
|   `-- synthetic_example.py
`-- tests/
    `-- test_core.py
```

No participant-level data, clinical identifiers, institutional paths, fitted reference
population, rarity cut-points, or diagnostic rules are included.

## Installation and checks

Python 3.10 or later is recommended.

```bash
python -m venv .venv
python -m pip install -e .
python examples/synthetic_example.py
python -m unittest discover -s tests -v
```

The only runtime dependencies are NumPy and pandas.

## Expected input

Each record is supplied as one pandas `DataFrame`. At minimum it must contain:

| Field | Meaning |
|---|---|
| `time_s` | Irregular elapsed time in seconds |
| `workload_w` | Externally delivered workload in watts |
| one or more channel columns | Physiological signals such as `VO2` or `VCO2` |

The caller also supplies audited ramp onset, peak time, and peak workload. A recovery
anchor may be supplied; otherwise the helper selects the first post-peak observation at
or below 10 W and falls back to that record's peak time.

```python
from cpet_skeleton import PhaseMarkers, transform_record

markers = PhaseMarkers(
    ramp_start_s=300.0,
    peak_time_s=720.0,
    peak_workload_w=110.0,
)

represented = transform_record(
    record,
    markers,
    channels=["VO2", "VCO2"],
)

vo2_ramp = represented.ramp["VO2"]
print(vo2_ramp.grid)
print(vo2_ramp.values)
print(vo2_ramp.defined)       # Explicit support mask
print(vo2_ramp.observed_bin)  # Direct bin summaries versus internal interpolation
```

The primary defaults use a 0--100% ramp grid in 1-percentage-point steps and a 0--180 s
recovery grid in 5-s steps. Median binning is followed by bounded internal linear
interpolation. Values outside the first and last defined bins remain unavailable; no
extrapolation is performed.

## Optional paired demonstration

`paired.py` illustrates one operation that can be performed after two records have been
represented independently:

```python
from cpet_skeleton import compare_records

result = compare_records(
    represented_cpet1,
    represented_cpet2,
    phase="ramp",
    channel="VCO2",
)

print(result.pointwise_difference)  # CPET2 minus CPET1
print(result.jointly_defined)
print(result.rms)
```

RMS is returned only when the jointly defined grid meets the configured support rule.
It remains in the channel's native units and is not a cross-channel importance score or
a universal inter-individual distance.

## Scope and reproducibility boundary

This is a readable methodological skeleton for external inspection and adaptation. It
does not replace the audited study pipeline. Before use with real data, investigators
must define and validate device-specific import rules, physiological QC masks, phase
markers, protocol compatibility, and missing-data handling. The supplied functions do
not create physiological equivalence, normative reference limits, clinical classes, or
validated classifiers.

The synthetic example is generated locally and contains no study data. Numerical
results reported in the manuscript require the restricted source records and the frozen
study analysis pipeline.

## Public-release checklist

Before publishing this folder as a standalone GitHub repository:

1. choose an institutionally approved open-source license;
2. add the final journal citation and repository DOI;
3. create a versioned release matching the submitted or accepted manuscript; and
4. replace the staged code-availability wording in the manuscript with the public URL.

