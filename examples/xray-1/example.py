import numpy as np

import aimsprop as ai

# => Parse/Align/Weight/Clean <= #

# Parse a series of FMS90 trajectories that Hayley has run for Stilbene
# Troubles: 11, 12 (not finished), 15
# trajs = [ai.parse_fms90('/home/hweir/stilbene/5-aims/s0_extension/aims_%04d/job1' % x) for x in range(1,16+1) if x not in [11, 12, 15]]
trajs = [ai.parse_fms90("/home/parrish/chem/stil/6-aims/%04d" % x) for x in [1, 2]]
# TODO: Align trajectories to IC transition dipole moment (on z) and weight by oscillator strength at IC
# Merge the trajectories into one super-big Trajectory with uniform weights
traj = ai.Trajectory.merge(trajs, ws=[1.0 / len(trajs)] * len(trajs), labels=[1, 2])
# Compute properties at ~1 fs intervals, removing nonsense due to adaptive timesteps
ts = np.arange(0.0, max(traj.ts), 400.0)  # TODO: Cleaner edges
traj = traj.interpolate_nearest(ts)

# => Tag with X-Ray Scattering Signal <= #

# Grab form factors for all atoms in this trajectory
factors = ai.iam.AtomicFormFactor.build_factors(traj.frames[0], mode="xray")
# The q values to compute scattering cross section at (in A^-1)
q = np.linspace(0.5, 3.0, 100)
# Compute the diffraction pattern moments (l=0,2,4)
ai.iam.compute_diffraction_moments(
    traj=traj,
    key="xray",
    q=q,
    factors=factors,
    nlebedev=74,
    nlebedev2=74,
    nomega2=12,
    nlegendre=4,
    print_level=True,
)
# Plots of the result
ai.plot_vector(
    "I0.pdf",
    traj,
    "xray-0",
    y=q,
    ylabel=r"$q [\AA{}^{-1}]$",
    time_units="fs",
    diff=True,
)
ai.plot_vector(
    "I2.pdf",
    traj,
    "xray-2",
    y=q,
    ylabel=r"$q [\AA{}^{-1}]$",
    time_units="fs",
    diff=True,
)
ai.plot_vector(
    "I4.pdf",
    traj,
    "xray-4",
    y=q,
    ylabel=r"$q [\AA{}^{-1}]$",
    time_units="fs",
    diff=True,
)  # should be zero


# Compute the diffraction pattern moments (l=0 only)
ai.iam.compute_diffraction_moment0(
    traj=traj,
    key="xray0",
    q=q,
    factors=factors,
    nlebedev=74,
    print_level=True,
)
# Plots of the result
ai.plot_vector(
    "I00.pdf",
    traj,
    "xray0-0",
    y=q,
    ylabel=r"$q [\AA{}^{-1}]$",
    time_units="fs",
    diff=True,
)
