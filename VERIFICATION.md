# Verification of version 0.2.0

Checked locally on 1 October 2026. This records implementation checks, not
physiological validation or a published GitHub release.

## Behaviour tests and installation

All **20 tests** passed in both environments below:

| Environment | Python | NumPy | pandas |
|---|---|---|---|
| Study verification environment | 3.13.9 | 2.3.5 | 2.3.3 |
| Fresh virtual environment, regular package installation | 3.13.9 | 2.5.3 | 3.0.6 |

The fresh installation built a wheel with `python -m pip install .` from a
separate source copy. Imports resolved to the environment's `site-packages`,
not the study workspace. `python -m pip check` found no broken requirements.
The synthetic example completed and wrote all five documented outputs.
Testing was on Windows; other operating systems and Python versions were not
exercised in this check. Python >=3.10 is the declared package requirement.

The optional [constraints file](constraints-verified.txt) records exact dependency
versions for the fresh Python 3.13 environment. To repeat that dependency setup:

```bash
python -m pip install -c constraints-verified.txt .
```

The seeded example produced the following VCO2 summaries (illustrative synthetic
values only; RMS is in the input channel's units):

| Phase | Jointly defined positions | Required positions | RMS |
|---|---:|---:|---:|
| Ramp | 91 | 21 | 32.95382346 |
| Recovery | 36 | 8 | 33.62799838 |

The missing-recovery case produced zero jointly defined positions and an
unavailable RMS, while retaining the other record's recovery.

## Private study-reference concordance

The public implementation was applied locally to all **120 records, 17 channels,
two phases and both QC configurations** using audited ramp/peak markers and the
already harmonised, QC-masked inputs. Inferred anchors matched all 120 audited
anchor/absence entries in each QC configuration. No participant data are included
in this repository.

For both record curves and their paired differences, the masks matched the
corrected nominal study reference exactly. The maximum absolute numerical
difference was **0** in both phases and both QC configurations. Comparisons used
`rtol=1e-10`, `atol=1e-8`, with missing values compared as missing. The default
101-point ramp and 37-point recovery grids were checked.

This is a concordance check of the nominal representation. It does not rerun the
historical controlled experiments, all anchor perturbations, clinical analyses,
cohort display selection, or the entire manuscript pipeline. Those private
checks require restricted inputs and are recorded separately by the authors.
