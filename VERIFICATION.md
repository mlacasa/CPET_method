# Verification of version 0.2.1

Checked locally on **3 October 2026**, on Windows 11, before publication. These
are software checks and teaching examples. The manuscript reconstruction study
and restricted clinical concordance were **not rerun**. Their dates and scope
remain in the [historical 0.2.0 record](docs/VERIFICATION_0.2.0.md) and
[reconstruction index](validation/README.md).

Machine-readable results and test logs are in [verification/2026-10-03](verification/2026-10-03).

## Baseline and decimal-bin correction

Before modifying the example, the 20 tests at reviewed commit
`4b424221b0d5e47d7d7eb05e8028263a12a8999d` passed. The unchanged seeded example
matched the documented VCO2 RMS values 32.95382346049664 (91 ramp locations)
and 33.627998383469894 (36 recovery locations).

The reported decimal recovery grid `[0, 0.1, 0.2, 0.3]` reproduced the defect:
`[10, 11, 12, NaN]`. Version 0.2.1 returns `[10, 11, 12, 13]`, with all four
occupied/defined positions and no interpolated positions. Three new endpoint
regressions fail on 0.2.0 and pass on 0.2.1; a fourth guards native-domain
exclusion. The two-bin endpoint case retains `[NaN, 11, 12, 13]`, with the
interior value marked interpolated. Existing tests still cover ties-to-even,
median binning, the minimum support rule and phase-local smoothing.

## Current environments and checks

| Environment | Python | NumPy | pandas | Test result |
|---|---|---|---|---|
| Existing study environment, source import | 3.13.9 | 2.3.5 | 2.3.3 | 24 passed |
| New isolated virtual environment, regular installed wheel | 3.13.9 | 2.5.3 | 3.0.6 | 24 passed |

The clean environment used Matplotlib 3.11.2 and ipykernel 7.4.0. Imports resolved
to its `site-packages`, not the source tree; `pip check` found no broken
requirements. The original five example files matched the baseline exactly
(numerical CSV contents and JSON values). Windows was tested; declared support
for Python >=3.10 is not a claim that all versions/platforms were exercised.

The additional nominal regression compares 0.2.0 with the installed candidate
on 40 synthetic records, three channels, both phases and 20 pairs, seed 20261003.
It exercises irregular times, channel gaps, absent recoveries and workload
downsteps. All **1,920 arrays/scalars matched exactly**, including values, grids,
support masks, differences, support counts, eligibility and RMS. This supports
nominal-grid regression coverage; it is not a new clinical concordance check.

Commands, from the repository root after activating a new virtual environment:

```bash
python -m pip install ".[tutorial]" nbclient nbformat ipykernel
python -m pip check
python -m unittest discover -s tests -v
python examples/synthetic_example.py --output-dir example_output
python examples/unequal_duration_example.py --output-dir unequal_duration_output

git clone https://github.com/mlacasa/CPET_method.git ../cpet-baseline
git -C ../cpet-baseline checkout 4b424221b0d5e47d7d7eb05e8028263a12a8999d
python verification/compare_nominal.py --baseline-src ../cpet-baseline/src
```

[constraints-verified.txt](constraints-verified.txt) retains the exact NumPy/pandas
versions from the historical regular-install check; the new clean installation
also resolved those versions. It is optional, not a lockfile for all notebook
dependencies.

## Unequal-duration example and notebook

The four/eight-minute example completed with native ramp counts 81/161 and
recovery counts 61/91. Its outputs preserve 81 ramp and 36 recovery positions
for each complete record. Missing terminal recovery values do not become
extrapolated values. The absent-recovery copy has zero recovery support.

| Illustration | Ramp RMS (mL/min) | Recovery RMS (mL/min) |
|---|---:|---:|
| Same analytic targets, different sampling | 2.73246683 | 8.20703253 |
| Known +30 mL/min copy | 30 | 30 |
| Truncated copy compared with its complete source | 0 | Unavailable |

The alignment illustration is not required to have zero RMS: it uses unequal
native samples and the default smoothing. The known-offset comparison isolates
the deliberate amplitude change while keeping coordinates/missingness fixed.

All **seven code cells** of the notebook executed in order, without cell errors,
using the regular-install environment as kernel. Before the tag existed, the
documented `CPET_TUTORIAL_PACKAGE` override selected the local 0.2.1 source build.
The default notebook installation is the public Git tag `v0.2.1`. The local
execution controller used nbclient 0.10.2 / nbformat 5.10.4; the isolated runtime
also has nbclient 0.11.0 / nbformat 5.11.1 installed. Published notebook outputs
are cleared; the authors retain the executed notebook in the local audit.

**Verificado localmente; ejecución en Colab pendiente.** A Colab link is provided,
but this release does not claim an execution in an actual Google Colab runtime.
To execute locally with a registered Python kernel:

```bash
python -m ipykernel install --user --name cpet-review
python -c "import nbformat; from nbclient import NotebookClient; n=nbformat.read('notebooks/CPET_method_tutorial.ipynb', as_version=4); NotebookClient(n, kernel_name='cpet-review', timeout=240).execute(); nbformat.write(n, 'CPET_method_tutorial_executed.ipynb')"
```

The example and notebook each exported 19 files. CSV/JSON contents agreed
exactly; matrix shapes, coordinate labels, record order, Boolean masks, NaN
support and known differences were checked. Both ZIP archives passed integrity
checks and contained exactly their respective exported files. Ramp/recovery
plots and mask figures were inspected visually. Download handling was exercised
through the local fallback; the Colab `files.download` branch remains pending.

## Publication scope and historical evidence

The validation selection includes 17 existing design/code/aggregate files and
an explanatory index. SHA-256 provenance is recorded for every selected file;
only an absolute source path and an internal task-description field were
redacted from two configuration files. Scientific parameters and stored
estimates were preserved. Frozen source snapshots are for inspection and retain
their documented restricted dependencies.

No clinical observations, per-record sampling templates, per-template error
tables, author filesystem paths or reviewer correspondence are included in this
selection. The exact remaining dependencies are listed in the validation index.
The new tutorial is separate from that historical experiment.
