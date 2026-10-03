# Protocol-aware CPET functional representation

Transform irregular breath-by-breath CPET records of different lengths into
comparable **ramp and recovery matrices with explicit availability masks**.
Each record is transformed independently: ramp uses % of its attained peak
workload; recovery uses seconds from its protocol transition. Channel amplitudes
retain their original units.

**Version 0.2.1 — manuscript and software under peer review.** The authors make
the code publicly available now so reviewers and researchers can inspect and run
it. Updates may follow review. The final journal citation, archival DOI and reuse
license are pending; public access does not establish an open reuse license.
See [release status](RELEASE.md).

## Start with the notebook

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/mlacasa/CPET_method/blob/v0.2.1/notebooks/CPET_method_tutorial.ipynb)

[Open the notebook on GitHub](notebooks/CPET_method_tutorial.ipynb). It installs
the tagged version, generates four/eight-minute synthetic ramps, gives recovery
the same treatment, plots native and represented signals, explains support masks
and downloads the results. No clinical records, credentials or author-specific
paths are required. **Verificado localmente; ejecución en Colab pendiente.**

## What you supply

- One harmonised pandas `DataFrame` per record: unique finite `time_s` (seconds),
  finite `workload_w` (watts), and named numerical channels in their original units.
- Audited ramp onset, peak time and positive peak workload; an audited recovery
  anchor or a protocol-specific transition rule. Missing channel values are NaNs.
- Upstream device harmonisation and physiological quality control. Keep channel
  gaps as missing values without silently deleting their native breath positions.

The package provides the transformation; it does not automatically perform the
entire ingestion, QC and marker-auditing workflow.

## Install and run locally

Python >=3.10; core dependencies are NumPy >=1.24 and pandas >=2.0. The tutorial
extra adds Matplotlib. From a terminal:

```bash
git clone --branch v0.2.1 https://github.com/mlacasa/CPET_method.git
cd CPET_method
python -m venv .venv
```

Activate with `source .venv/bin/activate` on Linux/macOS, or
`.\.venv\Scripts\Activate.ps1` on Windows PowerShell, then:

```bash
python -m pip install ".[tutorial]"
python -m unittest discover -s tests -v
python examples/unequal_duration_example.py --output-dir unequal_duration_output
```

If activation is unavailable, use `.venv/bin/python` (Linux/macOS) or
`.\.venv\Scripts\python.exe` (Windows). For development use `pip install -e ".[tutorial]"`.
The notebook also supports local Jupyter execution; see [verification](VERIFICATION.md).

## What you receive

The default grids have **101 ramp positions** (0–100 %Wpeak) and **37 recovery
positions** (0–180 seconds). The grid is a common set of candidate coordinates;
availability remains specific to each record and channel.

| Output | Meaning |
|---|---|
| `values` | Channel values; NaN outside supported bounds |
| `defined` | Positions with available output |
| `observed_bin` | Available positions supplied by occupied-bin summaries |
| `interpolated` | Available positions estimated between occupied bins |

The example exports native observations, markers/counts, coordinate-labelled
matrices and all masks for both phases, optional paired differences, figures and
a ZIP bundle. Empty numerical CSV cells mean unavailable, not zero. Output files
are replaced when the same output directory is reused.

The four/eight-minute comparison illustrates alignment, a separate +30 mL/min
copy isolates a known amplitude change, and a truncated copy illustrates absent
recovery. Equal %Wpeak need not mean equal watts or equal physiological state.
See [example interpretation and outputs](docs/EXAMPLES.md).

## Inspect the method and evidence

| Material | Scope |
|---|---|
| [Method contract and API examples](docs/METHOD.md) | Smoothing, binning, interpolation, markers, masks and optional pairing |
| [Current verification](VERIFICATION.md) | Installation, 24 software tests, nominal regression and notebook execution |
| [Historical reconstruction materials](validation/README.md) | Analytic functions, scenario parameters, metric definitions, frozen code and aggregate results; restricted dependencies identified |
| [Historical 0.2.0 verification](docs/VERIFICATION_0.2.0.md) | Checks dated 1 October 2026, including restricted study-reference concordance; not rerun here |
| [Changelog](CHANGELOG.md) | Scientific-policy changes in 0.2.0 and decimal-bin correction in 0.2.1 |

The teaching example is separate from the manuscript reconstruction experiment.
The public materials do not reproduce every manuscript result from public inputs.
Interpolation supplies comparable positions, not independent measurements or
clinical validation.

## Citation and feedback

[CITATION.cff](CITATION.cff) identifies this software version. Report both the
version and Git commit; no final journal citation or DOI is claimed. Report bugs
with a small synthetic example, Python/NumPy/pandas versions and expected versus
observed behaviour. Do not post participant records or clinical identifiers.
