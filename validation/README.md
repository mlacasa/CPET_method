# Historical reconstruction study: inspectable public materials

This directory publishes a selected, documented part of the existing manuscript
experiment. **No simulations were rerun for version 0.2.1.** It contains code,
analytic definitions, aggregate scales, configuration and aggregate results;
it contains no individual clinical observations or empirical timestamp/mask
arrays. [PROVENANCE.json](PROVENANCE.json) records original and public SHA-256
hashes, historical run identifiers and the two metadata redactions.

## Four different kinds of evidence

| Evidence | What it establishes | Access and timing |
|---|---|---|
| [Software tests](../tests) | Defined numerical behaviour, boundary handling and masks | Public synthetic cases; 24 tests run for 0.2.1 |
| [Teaching example](../docs/EXAMPLES.md) | How to prepare inputs, map records and inspect outputs | Public synthetic four/eight-minute example; new in 0.2.1 |
| Reconstruction experiment, indexed below | Error against analytic targets on historical empirical observation/missingness patterns | Frozen September/October 2026 code and aggregates; exact rerun needs restricted inputs |
| [Clinical reference concordance](../docs/VERIFICATION_0.2.0.md) | Agreement with corrected nominal study objects | Historical check on 1 October 2026, using restricted inputs; not rerun for 0.2.1 |

## Index

| File | Purpose |
|---|---|
| [Analytic dictionary](design/truth_model_dictionary.json) | Channel-specific mathematical targets, features and known session changes |
| [Scenario manifest](design/scenario_manifest.csv) and [design configuration](design/analysis_config.json) | Exact 31 scenarios, 10 families, five replicates, five filter candidates, seed 20260930 |
| [Channel scales](design/channel_scales.csv) | One aggregate offset/amplitude per channel; no per-record scale inputs |
| [Data/metric dictionary](design/data_dictionary.json) | Units, normalization, keys and feature-error definitions |
| [Frozen study code](source/study_code_frozen.py) | `baseline`, `make_truth`, `scenarios`, `apply_filter`, `grid_matrix`, `metric_columns` and simulation/evaluation logic |
| [Frozen summary code](source/summarize_validation_final.py) | Family-balanced summaries, paired contrasts and conditional bootstrap intervals |
| [Summary configuration](aggregates/analysis_config.json) | Aggregation, finite-value handling, seed 20261001 and 2000 resamples |
| [All-channel aggregates](aggregates/validation_by_channel.csv) | Curve and paired-difference RMSE, bias, coverage and nominal intervals |
| [Family aggregates](aggregates/validation_by_family.csv) | Same endpoints by scenario family, including adverse feature cases |
| [Contrasts](aggregates/validation_contrasts.csv) | Mean-one-pass minus no-presmoothing contrasts; no p-values |
| [Main table values](aggregates/main_table_values.csv), [main figure values](aggregates/main_figure_values.csv) | Stored VCO2 selections for the manuscript |
| [Historical environment](aggregates/environment.json), [QA](aggregates/quality_assurance.json), [completion record](design/simulation_complete.json), [contract checks](design/contract_tests.csv) | Original execution evidence, not current reruns |

## Definitions and interpretation

Analytic truth is evaluated at native coordinates and at target grid centres;
it is not defined as a smoothed or bin-averaged truth. With finite prediction and
truth, RMSE is `sqrt(mean((prediction - truth)**2))`, bias is
`mean(prediction - truth)`, and coverage is the finite comparison count divided
by the full candidate count. Paired error compares the reconstructed difference
with the known analytic difference. Unsupported values stay missing.

Errors are normalized by the channel's aggregate amplitude scale, **not** by the
true paired change or a percentage. Physical-unit RMSE/bias multiply normalized
errors by that scale; VCO2 uses 862.9 mL/min. `*_low` and `*_high` are nominal
95% percentile intervals; `*_n` gives the finite template count. Coverage is
stored as a proportion. The supplied dictionaries and code retain the original
feature definitions and the separate historical composite selection score.

The reconstruction uses 60 paired empirical templates (120 records), 17 channels,
two phases, 31 scenarios, five replicates and five candidate filters. Replicates
are averaged before equal session/scenario weights within families, equal family
weights within templates, and equal template weights. The secondary bootstrap
uses 2000 paired-template resamples, keeping dependent observations together.
Its uncertainty is conditional on the fixed synthetic design and stored
replicate averages; it is not clinical or Monte Carlo uncertainty.

The heterogeneous templates supply observed coordinates and missingness patterns,
not clinical ground truth. The historical design includes one pair's recovery
segmentation from before the later anchor correction. It was not silently
replaced with the corrected 119-recovery/59-pair configuration. Published results
must be interpreted with that historical limitation. These materials do not show
that reconstruction error is independent of duration.

The operational manuscript method uses one five-position arithmetic mean.
All five candidates and adverse scenario families are retained here; neither
universal optimality nor preservation of every narrow feature is claimed. The
historical composite ranking and separate reconstruction endpoints should not
be conflated with a new filter-selection exercise.

## Exact limits of public reproduction

The frozen `.py` files support **inspection**, not a standalone public rerun.
They preserve the historical implementation, including its nominal-grid binning;
the package's 0.2.1 decimal-grid fix has not been backported into this evidence.
They are outside the installed package and are never imported by the notebook.

The original study driver imports an unavailable project module
`audit_reproduction` (`HERE`, `ROOT`, `M`, `D`, `S`, `UNITS`, `sha`, `dump`). Its
input hash audit belongs to reference run `20260930T101337+0200`. It reads:

- `cpet_long_all_patients_v3_mainHardNA.csv` and
  `cpet_long_all_patients_v3_sensHardPlusSusNA.csv`: individual harmonised native
  observations, times, workloads, QC masks and record identifiers.
- `step3_master_by_test_v3.csv`: record-specific audited phase/peak markers and
  the historical recovery segmentation.

Recreating synthetic observation patterns requires the MAIN empirical timestamps,
workloads, per-channel finite masks and historical markers. Recreating empirical
scale estimates and clinical impact additionally needs the observed channel
values (and SENS inputs for that QC comparison). Those files and individual
template manifests are not published here.

The secondary summary driver requires the original per-template
`synthetic/metrics_template_*.parquet` and `synthetic/paired_template_*.parquet`
files, along with its source dictionaries/configuration and stored
`tables/T-SR4_known_paired_change.csv`, `T-SR4_by_phase_variable.csv` and
`T-SR4_global.csv` for internal regression checks. Feature reanalysis also needs
`features_template_*.parquet`. These per-template inputs are not included.
The aggregate CSVs permit inspection of reported estimates, but not reconstruction
of the bootstrap draws or all figures/tables from scratch. No missing result or
parameter has been invented to fill those gaps. Restricted access follows the
manuscript's Data availability / corresponding-author route.

The synthetic notebook has its own invented coordinates and is therefore an
educational example, not a substitute reproduction of this experiment.
