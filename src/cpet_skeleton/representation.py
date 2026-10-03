"""General record-level CPET functional representation.

The module ends with supported ramp and recovery functions. It deliberately does not
define a distance, reference population, rarity rule, cluster, or classifier.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Iterable

import numpy as np
import pandas as pd


@dataclass(frozen=True)
class PhaseMarkers:
    """Protocol markers determined independently for one CPET record."""

    ramp_start_s: float
    peak_time_s: float
    peak_workload_w: float
    recovery_anchor_s: float | None = None


@dataclass(frozen=True)
class SupportedTrajectory:
    """One channel on a fixed grid with explicit support provenance."""

    grid: np.ndarray
    values: np.ndarray
    defined: np.ndarray
    observed_bin: np.ndarray

    @property
    def interpolated(self) -> np.ndarray:
        """Locations defined by bounded interpolation rather than a bin summary."""

        return self.defined & ~self.observed_bin


@dataclass(frozen=True)
class RecordRepresentation:
    """Phase-specific output for one independently processed CPET record."""

    ramp: dict[str, SupportedTrajectory]
    recovery: dict[str, SupportedTrajectory]
    markers: PhaseMarkers


def _numeric_record(
    record: pd.DataFrame,
    *,
    time_col: str,
    workload_col: str,
    channels: Iterable[str],
) -> pd.DataFrame:
    """Validate harmonised structure without deleting native breath positions."""

    channels = list(channels)
    required = [time_col, workload_col, *channels]
    if len(set(required)) != len(required) or not record.columns.is_unique:
        raise ValueError("Time, workload and channel columns must be distinct and unique.")
    missing = [column for column in required if column not in record.columns]
    if missing:
        raise KeyError(f"Missing required columns: {missing}")

    frame = record.loc[:, required].copy()
    for column in required:
        frame[column] = pd.to_numeric(frame[column], errors="coerce")
    frame = frame.replace([np.inf, -np.inf], np.nan)
    if frame[[time_col, workload_col]].isna().any().any():
        raise ValueError("Time and workload must be finite; harmonise the record first.")
    frame = frame.sort_values(time_col, kind="mergesort")
    if frame[time_col].duplicated().any():
        raise ValueError("Duplicate times require an explicit upstream harmonisation policy.")
    if frame.empty:
        raise ValueError("No structurally valid observations remain in the record.")
    if not frame[time_col].is_monotonic_increasing:
        raise ValueError("Elapsed time must be monotonic after duplicate handling.")
    return frame


def infer_recovery_anchor(
    record: pd.DataFrame,
    peak_time_s: float,
    *,
    time_col: str = "time_s",
    workload_col: str = "workload_w",
    threshold_w: float = 10.0,
) -> float | None:
    """Return the first post-peak time at or below the workload threshold.

    Post-peak means strictly later than the audited peak time. Return None when
    no transition is recorded. No peak fallback or paired record is consulted.
    """

    if not np.isfinite(peak_time_s) or not np.isfinite(threshold_w):
        raise ValueError("Peak time and recovery threshold must be finite.")
    frame = _numeric_record(record, time_col=time_col, workload_col=workload_col, channels=[])
    candidates = frame.loc[
        (frame[time_col] > peak_time_s) & (frame[workload_col] <= threshold_w),
        time_col,
    ]
    return float(candidates.iloc[0]) if not candidates.empty else None


def _smooth_preserving_missing(values: pd.Series, window: int) -> np.ndarray:
    """One centred arithmetic-mean pass, restoring the original missing mask."""

    smoothed = values.rolling(window, center=True, min_periods=1).mean()
    smoothed = smoothed.where(values.notna())
    return smoothed.to_numpy(dtype=float)


def _supported_grid_curve(
    coordinate: np.ndarray,
    values: np.ndarray,
    grid: np.ndarray,
    step: float,
) -> SupportedTrajectory:
    """Median binning plus linear interpolation within observed-bin bounds only."""

    # Restrict native coordinates BEFORE rounding; outside events cannot round in.
    valid = (np.isfinite(coordinate) & np.isfinite(values)
             & (coordinate >= grid[0]) & (coordinate <= grid[-1]))
    output = np.full(grid.shape, np.nan, dtype=float)
    observed = np.zeros(grid.shape, dtype=bool)
    if not valid.any():
        return SupportedTrajectory(grid.copy(), output, np.isfinite(output), observed)

    binned = pd.DataFrame(
        {
            "bin": np.round(coordinate[valid] / step) * step,
            "value": values[valid],
        }
    )
    summary = binned.groupby("bin", as_index=False)["value"].median().sort_values("bin")
    summary = summary.loc[
        (summary["bin"] >= grid.min()) & (summary["bin"] <= grid.max())
    ]
    if len(summary) < 2:
        return SupportedTrajectory(grid.copy(), output, np.isfinite(output), observed)

    x = summary["bin"].to_numpy(dtype=float)
    y = summary["value"].to_numpy(dtype=float)
    inside = (grid >= x.min()) & (grid <= x.max())
    output[inside] = np.interp(grid[inside], x, y)
    for bin_value in x:
        observed |= np.isclose(grid, bin_value, rtol=0, atol=1e-10)

    return SupportedTrajectory(
        grid=grid.copy(),
        values=output,
        defined=np.isfinite(output),
        observed_bin=observed & np.isfinite(output),
    )


def _validate_grid(grid: np.ndarray, name: str, upper: float) -> np.ndarray:
    """Require a finite, increasing uniform grid with the implemented zero origin."""
    grid = np.asarray(grid, dtype=float)
    if (grid.ndim != 1 or len(grid) < 2 or not np.isfinite(grid).all()
            or grid[0] != 0 or grid[-1] > upper):
        raise ValueError(f"{name} must start at zero and contain >=2 finite points within [0, {upper}].")
    steps = np.diff(grid)
    if steps[0] <= 0 or not np.allclose(steps, steps[0], rtol=1e-10, atol=1e-12):
        raise ValueError(f"{name} must be strictly increasing and uniformly spaced.")
    return grid


def transform_record(
    record: pd.DataFrame,
    markers: PhaseMarkers,
    *,
    channels: Iterable[str],
    time_col: str = "time_s",
    workload_col: str = "workload_w",
    rolling_window: int = 5,
    ramp_grid: np.ndarray | None = None,
    recovery_grid: np.ndarray | None = None,
    recovery_horizon_s: float = 180.0,
    recovery_threshold_w: float = 10.0,
    infer_recovery_if_missing: bool = True,
) -> RecordRepresentation:
    """Map one CPET record to supported ramp and recovery trajectories.

    Every argument is record-specific. The function never inspects a comparison partner.
    QC masking should be applied by the caller before this transformation.
    None triggers anchor inference by default. Set infer_recovery_if_missing=False
    to retain an explicitly unavailable phase. A missing phase returns all-NaN
    values and all-False masks, while the available ramp is retained.
    """

    channels = list(channels)
    if not channels:
        raise ValueError("At least one physiological channel is required.")
    if isinstance(rolling_window, bool) or not isinstance(rolling_window, (int, np.integer)) or rolling_window < 1:
        raise ValueError("rolling_window must be a positive integer.")
    if not np.isfinite([markers.ramp_start_s, markers.peak_time_s, markers.peak_workload_w]).all():
        raise ValueError("Ramp and peak markers must be finite.")
    if markers.peak_workload_w <= 0:
        raise ValueError("peak_workload_w must be positive.")
    if markers.peak_time_s < markers.ramp_start_s:
        raise ValueError("peak_time_s must not precede ramp_start_s.")
    if not np.isfinite(recovery_horizon_s) or recovery_horizon_s <= 0:
        raise ValueError("recovery_horizon_s must be positive and finite.")
    if not np.isfinite(recovery_threshold_w):
        raise ValueError("recovery_threshold_w must be finite.")
    if markers.recovery_anchor_s is not None and (
        not np.isfinite(markers.recovery_anchor_s)
        or markers.recovery_anchor_s <= markers.peak_time_s
    ):
        raise ValueError("An explicit recovery anchor must be finite and strictly after peak time.")

    frame = _numeric_record(
        record,
        time_col=time_col,
        workload_col=workload_col,
        channels=channels,
    )
    if markers.recovery_anchor_s is None and infer_recovery_if_missing:
        anchor = infer_recovery_anchor(
            frame,
            markers.peak_time_s,
            time_col=time_col,
            workload_col=workload_col,
            threshold_w=recovery_threshold_w,
        )
        markers = replace(markers, recovery_anchor_s=anchor)

    ramp_grid = (
        np.arange(0.0, 101.0, 1.0) if ramp_grid is None else np.asarray(ramp_grid, float)
    )
    recovery_grid = (
        np.arange(int(np.floor(recovery_horizon_s / 5.0)) + 1, dtype=float) * 5.0
        if recovery_grid is None
        else np.asarray(recovery_grid, float)
    )
    ramp_grid = _validate_grid(ramp_grid, "ramp_grid", 100.0)
    recovery_grid = _validate_grid(recovery_grid, "recovery_grid", recovery_horizon_s)

    ramp_frame = frame.loc[
        (frame[time_col] >= markers.ramp_start_s)
        & (frame[time_col] <= markers.peak_time_s)
    ].copy()
    anchor = markers.recovery_anchor_s
    recovery_frame = frame.iloc[:0].copy() if anchor is None else frame.loc[
        (frame[time_col] >= anchor)
        & (frame[time_col] <= anchor + recovery_horizon_s)
    ].copy()
    if ramp_frame.empty:
        raise ValueError("No observations fall within the defined ramp phase.")

    workload_envelope = np.maximum.accumulate(ramp_frame[workload_col].to_numpy(float))
    ramp_coordinate = 100.0 * workload_envelope / markers.peak_workload_w
    recovery_coordinate = (
        np.empty(0) if anchor is None else recovery_frame[time_col].to_numpy(float) - anchor
    )

    ramp: dict[str, SupportedTrajectory] = {}
    recovery: dict[str, SupportedTrajectory] = {}
    ramp_step = float(np.median(np.diff(ramp_grid)))
    recovery_step = float(np.median(np.diff(recovery_grid)))
    for channel in channels:
        ramp_values = _smooth_preserving_missing(ramp_frame[channel], rolling_window)
        recovery_values = _smooth_preserving_missing(
            recovery_frame[channel], rolling_window
        )
        ramp[channel] = _supported_grid_curve(
            ramp_coordinate, ramp_values, ramp_grid, ramp_step
        )
        recovery[channel] = _supported_grid_curve(
            recovery_coordinate, recovery_values, recovery_grid, recovery_step
        )

    return RecordRepresentation(ramp=ramp, recovery=recovery, markers=markers)
