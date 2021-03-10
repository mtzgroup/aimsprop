import numpy as np

import aimsprop as ai

# => Parse/Align/Weight/Clean <= #

# Parse a series of FMS90 trajectories that Hayley has run for Stilbene
ICs = list(range(1, 16 + 1))
trajs = [
    ai.parse_fms90("/home/hweir/stilbene/5-aims/s0_extension/aims_%04d/job1" % x)
    for x in ICs
]
# Merge the trajectories into one super-big Trajectory with uniform weights
traj = ai.Trajectory.merge(trajs, ws=[1.0 / len(trajs)] * len(trajs), labels=ICs)
# Compute properties at ~1 fs intervals, removing nonsense due to adaptive timesteps
ts = np.arange(0.0, max(traj.ts), 40.0)
trajs = [traj2.interpolate_nearest(ts) for traj2 in trajs]
traj = traj.interpolate_nearest(ts)

# => Population Plot <= #

# If you want the detailed population curves (corresponding to traj.ts), e.g., for decay constant fitting:
# populations = ai.compute_population(traj)
# If you just want a plot of the populations
plt = ai.plot_population(
    "P.pdf",
    traj,
    trajs,
    time_units="fs",
    state_colors=["r", "b"],
)
