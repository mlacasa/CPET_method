"""Optional paired operations applied after independent record representation."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np

from .representation import RecordRepresentation, SupportedTrajectory


@dataclass(frozen=True)
class PairedDiscrepancy:
    """Pointwise CPET2-minus-CPET1 difference and an optional RMS summary."""

    grid: np.ndarray
    pointwise_difference: np.ndarray
    jointly_defined: np.ndarray
    rms: float
    n_supported: int
    n_minimum: int
    eligible: bool


def compare_trajectories(
    cpet1: SupportedTrajectory,
    cpet2: SupportedTrajectory,
    *,
    support_fraction: float = 0.20,
    minimum_locations: int = 5,
) -> PairedDiscrepancy:
    """Compare two already represented trajectories on their jointly defined grid."""

    if cpet1.grid.shape != cpet2.grid.shape or not np.allclose(
        cpet1.grid, cpet2.grid
    ):
        raise ValueError("Trajectories must use the same phase-specific grid.")
    if not 0 <= support_fraction <= 1:
        raise ValueError("support_fraction must lie between 0 and 1.")

    joint = (
        cpet1.defined
        & cpet2.defined
        & np.isfinite(cpet1.values)
        & np.isfinite(cpet2.values)
    )
    difference = np.full(cpet1.grid.shape, np.nan, dtype=float)
    difference[joint] = cpet2.values[joint] - cpet1.values[joint]
    n_supported = int(joint.sum())
    n_minimum = max(minimum_locations, int(np.ceil(support_fraction * len(joint))))
    eligible = n_supported >= n_minimum
    rms = float(np.sqrt(np.mean(difference[joint] ** 2))) if eligible else float("nan")

    return PairedDiscrepancy(
        grid=cpet1.grid.copy(),
        pointwise_difference=difference,
        jointly_defined=joint,
        rms=rms,
        n_supported=n_supported,
        n_minimum=n_minimum,
        eligible=eligible,
    )


def compare_records(
    cpet1: RecordRepresentation,
    cpet2: RecordRepresentation,
    *,
    phase: Literal["ramp", "recovery"],
    channel: str,
    support_fraction: float = 0.20,
    minimum_locations: int = 5,
) -> PairedDiscrepancy:
    """Select one phase and channel, then call :func:`compare_trajectories`."""

    first = getattr(cpet1, phase)
    second = getattr(cpet2, phase)
    if channel not in first or channel not in second:
        raise KeyError(f"Channel {channel!r} is unavailable in one or both records.")
    return compare_trajectories(
        first[channel],
        second[channel],
        support_fraction=support_fraction,
        minimum_locations=minimum_locations,
    )

