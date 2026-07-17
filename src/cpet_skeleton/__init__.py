"""Public interface for the CPET functional-representation skeleton."""

from .paired import PairedDiscrepancy, compare_records, compare_trajectories
from .representation import (
    PhaseMarkers,
    RecordRepresentation,
    SupportedTrajectory,
    infer_recovery_anchor,
    transform_record,
)

__all__ = [
    "PairedDiscrepancy",
    "PhaseMarkers",
    "RecordRepresentation",
    "SupportedTrajectory",
    "compare_records",
    "compare_trajectories",
    "infer_recovery_anchor",
    "transform_record",
]

