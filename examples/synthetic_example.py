"""Run the framework on two synthetic, irregular breath-by-breath CPET records."""

from __future__ import annotations

import numpy as np
import pandas as pd

from cpet_skeleton import PhaseMarkers, compare_records, transform_record


def synthetic_record(seed: int, response_shift: float = 0.0) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    time_s = np.cumsum(rng.uniform(2.0, 4.5, size=220))
    ramp_start = 120.0
    peak_time = 520.0

    workload = np.where(
        time_s < ramp_start,
        10.0,
        np.minimum(110.0, 10.0 + (time_s - ramp_start) * 0.25),
    )
    workload = np.where(time_s > peak_time, 0.0, workload)
    vo2 = 350.0 + 9.0 * workload + response_shift + rng.normal(0.0, 35.0, len(time_s))
    vco2 = 300.0 + 8.5 * workload + response_shift + rng.normal(0.0, 40.0, len(time_s))
    recovery = np.maximum(time_s - peak_time, 0.0)
    vo2 -= np.where(recovery > 0, np.minimum(650.0, 3.5 * recovery), 0.0)
    vco2 -= np.where(recovery > 0, np.minimum(600.0, 3.2 * recovery), 0.0)

    return pd.DataFrame(
        {"time_s": time_s, "workload_w": workload, "VO2": vo2, "VCO2": vco2}
    )


def main() -> None:
    markers = PhaseMarkers(
        ramp_start_s=120.0,
        peak_time_s=520.0,
        peak_workload_w=110.0,
    )
    first = transform_record(
        synthetic_record(101), markers, channels=["VO2", "VCO2"]
    )
    second = transform_record(
        synthetic_record(202, response_shift=30.0),
        markers,
        channels=["VO2", "VCO2"],
    )

    ramp = compare_records(first, second, phase="ramp", channel="VCO2")
    recovery = compare_records(first, second, phase="recovery", channel="VCO2")

    print("Synthetic example; no study data are used.")
    print(f"Ramp defined locations, CPET1: {first.ramp['VCO2'].defined.sum()}")
    print(f"Recovery defined locations, CPET1: {first.recovery['VCO2'].defined.sum()}")
    print(f"Ramp joint support / minimum: {ramp.n_supported} / {ramp.n_minimum}")
    print(f"Ramp RMS: {ramp.rms:.2f}")
    print(
        f"Recovery joint support / minimum: "
        f"{recovery.n_supported} / {recovery.n_minimum}"
    )
    print(f"Recovery RMS: {recovery.rms:.2f}")


if __name__ == "__main__":
    main()

