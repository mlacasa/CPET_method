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
    """Select, coerce, sort, and deterministically collapse duplicate times."""

    channels = list(channels)
    required = [time_col, workload_col, *channels]
    missing = [column for column in required if column not in record.columns]
    if missing:
        raise KeyError(f"Missing required columns: {missing}")

    frame = record.loc[:, required].copy()
    for column in required:
        frame[column] = pd.to_numeric(frame[column], errors="coerce")
    frame = frame.dropna(subset=[time_col, workload_col])
    frame = frame.sort_values(time_col, kind="mergesort")

    # A median collapse is transparent and deterministic. Production pipelines should
    # replace it if their device-specific duplicate policy differs.
    frame = frame.groupby(time_col, as_index=False, sort=True).median(numeric_only=True)
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
) -> float:
    """Return the first post-peak time at or below the workload threshold.

    The fallback is the peak time from the same record. No paired record is consulted.
    """

    frame = record.loc[:, [time_col, workload_col]].copy()
    frame[time_col] = pd.to_numeric(frame[time_col], errors="coerce")
    frame[workload_col] = pd.to_numeric(frame[workload_col], errors="coerce")
    frame = frame.dropna().sort_values(time_col, kind="mergesort")
    candidates = frame.loc[
        (frame[time_col] >= peak_time_s) & (frame[workload_col] <= threshold_w),
        time_col,
    ]
    return float(candidates.iloc[0]) if not candidates.empty else float(peak_time_s)


def _smooth_preserving_missing(values: pd.Series, window: int) -> np.ndarray:
    """Centred rolling median that leaves originally missing locations missing."""

    smoothed = values.rolling(window, center=True, min_periods=1).median()
    smoothed = smoothed.where(values.notna())
    return smoothed.to_numpy(dtype=float)


def _supported_grid_curve(
    coordinate: np.ndarray,
    values: np.ndarray,
    grid: np.ndarray,
    step: float,
) -> SupportedTrajectory:
    """Median binning plus linear interpolation within observed-bin bounds only."""

    valid = np.isfinite(coordinate) & np.isfinite(values)
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
    if summary.empty:
        return SupportedTrajectory(grid.copy(), output, np.isfinite(output), observed)

    x = summary["bin"].to_numpy(dtype=float)
    y = summary["value"].to_numpy(dtype=float)
    if len(summary) == 1:
        matches = np.isclose(grid, x[0])
        output[matches] = y[0]
        observed[matches] = True
    else:
        inside = (grid >= x.min()) & (grid <= x.max())
        output[inside] = np.interp(grid[inside], x, y)
        for bin_value in x:
            observed |= np.isclose(grid, bin_value)

    return SupportedTrajectory(
        grid=grid.copy(),
        values=output,
        defined=np.isfinite(output),
        observed_bin=observed & np.isfinite(output),
    )


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
) -> RecordRepresentation:
    """Map one CPET record to supported ramp and recovery trajectories.

    Every argument is record-specific. The function never inspects a comparison partner.
    QC masking should be applied by the caller before this transformation.
    """

    channels = list(channels)
    if not channels:
        raise ValueError("At least one physiological channel is required.")
    if rolling_window < 1:
        raise ValueError("rolling_window must be at least 1.")
    if markers.peak_workload_w <= 0:
        raise ValueError("peak_workload_w must be positive.")
    if markers.peak_time_s < markers.ramp_start_s:
        raise ValueError("peak_time_s must not precede ramp_start_s.")

    frame = _numeric_record(
        record,
        time_col=time_col,
        workload_col=workload_col,
        channels=channels,
    )
    if markers.recovery_anchor_s is None:
        anchor = infer_recovery_anchor(
            frame,
            markers.peak_time_s,
            time_col=time_col,
            workload_col=workload_col,
        )
        markers = replace(markers, recovery_anchor_s=anchor)

    ramp_grid = (
        np.arange(0.0, 101.0, 1.0) if ramp_grid is None else np.asarray(ramp_grid, float)
    )
    recovery_grid = (
        np.arange(0.0, 181.0, 5.0)
        if recovery_grid is None
        else np.asarray(recovery_grid, float)
    )

    ramp_frame = frame.loc[
        (frame[time_col] >= markers.ramp_start_s)
        & (frame[time_col] <= markers.peak_time_s)
    ].copy()
    recovery_frame = frame.loc[
        (frame[time_col] >= markers.recovery_anchor_s)
        & (frame[time_col] <= markers.recovery_anchor_s + recovery_horizon_s)
    ].copy()
    if ramp_frame.empty:
        raise ValueError("No observations fall within the defined ramp phase.")

    workload_envelope = np.maximum.accumulate(ramp_frame[workload_col].to_numpy(float))
    ramp_coordinate = 100.0 * workload_envelope / markers.peak_workload_w
    recovery_coordinate = (
        recovery_frame[time_col].to_numpy(float) - markers.recovery_anchor_s
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

