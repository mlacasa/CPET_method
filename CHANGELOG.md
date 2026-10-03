# Changelog

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
