# Changelog

## 0.2.1 — 2026-10-03

- Assign bins by integer index and use validated grid centres, preventing a
  decimal endpoint such as 0.3 from being discarded through multiplication
  roundoff. Preserve native-domain exclusion, ties-to-even, medians, the
  two-occupied-bin rule, internal interpolation and support meanings.
- Add four decimal-grid regressions; retain the 20 existing tests and unchanged
  original seeded example. Record nominal-grid comparisons with 0.2.0.
- Add a self-contained, tagged-version notebook and an unequal-duration example
  with explicit ramp/recovery plots, masks, known differences, absent recovery
  and downloadable CSV/matrix exports.
- Publish selected historical reconstruction definitions, code and aggregate
  results with provenance and exact restricted-input limitations; no rerun of
  the manuscript experiment or clinical concordance is claimed.
- Shorten the README path to first execution; retain the detailed method
  contract and historical verification in secondary documentation. Keep
  under-review status, pending license and actual-Colab limitation explicit.


## 0.2.0 — manuscript revision

- Use a single centred arithmetic moving mean (five positions by default),
  restoring the original channel mask, instead of a rolling median.
- Infer recovery strictly after the audited peak. Return an unavailable phase
  when no transition exists; remove the peak fallback.
- Require two occupied bins and exclude out-of-domain native coordinates before
  rounding. Retain bin medians and bounded linear interpolation.
- Reject ambiguous duplicate times, invalid structural data and unsupported grid
  shapes. Provide explicit protocol anchors, threshold configuration and an
  option to retain declared recovery absence without inference.
- Preserve independent record processing and native channel units. Clarify
  optional full-grid paired summaries and their support requirements.
- Add contract tests, inspectable synthetic CSV/JSON outputs and installation,
  input/output, citation and release documentation.

This changes numerical behaviour relative to 0.1.0; it is not a cosmetic update.

## 0.1.0

Initial local methodological skeleton; not a published release.
