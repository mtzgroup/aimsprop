from pathlib import Path

import numpy as np
import pytest

import aimsprop as ap


@pytest.fixture(scope="module")
def trajectory():
    # 1.  Parse a series of FMS90 trajectories that Hayley has run for ethylene
    trajs = [
        ap.parse_fms90(Path(__file__).parent / "test_data" / f"000{x}") for x in [2, 3]
    ]
    # 2. Merge the trajectories into one super-big Trajectory with uniform weights
    traj = ap.Trajectory.merge(trajs, ws=[1.0 / len(trajs)] * len(trajs), labels=[2, 3])

    # 3. Interpolate trajectory with ~1 fs intervals, removing adaptive timesteps from AIMS
    ts = np.arange(0.0, max(traj.ts), 40.0)
    traj = traj.interpolate_nearest(ts)
    yield traj
