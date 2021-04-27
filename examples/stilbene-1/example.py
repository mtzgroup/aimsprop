import numpy as np

import aimsprop as ai

# Parse a series of FMS90 trajectories that Hayley has run for Stilbene
trajs = [
    ai.parse_fms90("/home/hweir/stilbene/5-aims/s0_extension/aims_%04d/job1" % x)
    for x in [1, 2, 3]
]
# Merge the trajectories into one super-big Bundle with uniform weights
traj = ai.Bundle.merge(trajs, ws=[1.0 / len(trajs)] * len(trajs), labels=[1, 2, 3])
print((traj.labels))
# Compute properties at ~1 fs intervals, removing nonsense due to adaptive timesteps
ts = np.arange(0.0, max(traj.ts), 40.0)  # TODO: Cleaner edges
traj = traj.interpolate_nearest(ts)

# => Showoff Plot of Bond Distances (Spaghetti + Blur) <= #

# Compute the a bond distance property for all Frames in traj
# This particular bond is the one that leads to a three-ring complex
ai.compute_bond(traj, "R01", 7, 8)

# Blur the bond distance
R = np.linspace(1.0, 6.0, 50)
ai.blur_property(traj, "R01", "Rblur", R, alpha=8.0)

# Plot the heat map of blurred bond distance
# ai.plot_vector('R2.pdf', traj, 'Rblur', y=R, ylabel=r'$R [\AA{}]$', time_units='fs', twosided=False, cmap=plt.get_cmap('viridis'))
ai.plot_vector(
    "R.pdf", traj, "Rblur", y=R, ylabel=r"$R [\AA{}]$", time_units="fs", nlevel=64
)
ai.plot_scalar(
    "R.pdf",
    traj,
    "R01",
    ylabel=r"$R_{CC} [\AA{}]$",
    time_units="fs",
    state_colors=["r", "b"],
    clf=False,
    plot_average=False,
)

# => Showoff Plot of Torsion Angle (Spaghetti + Blur) <= #

# Compute the a torsion angle property for all Frames in traj
# This particular torsion is the one the leads to cis-trans isomerization
ai.compute_torsion(traj, "T4019", 4, 0, 1, 9)
ai.unwrap_property(traj, "T4019", 360.0)

# Blur the torsion
T = np.linspace(-180.0, +180.0, 100)
ai.blur_property(traj, "T4019", "Tblur", T, alpha=0.02)

# Plot the heat map of blurred torison
ai.plot_vector(
    "T.pdf",
    traj,
    "Tblur",
    y=T,
    ylabel=r"$Theta [^{\circ{}}]$",
    time_units="fs",
    nlevel=64,
)
ai.plot_scalar(
    "T.pdf",
    traj,
    "T4019",
    ylabel=r"$\Theta [^{\circ{}}]$",
    time_units="fs",
    state_colors=["r", "b"],
    clf=False,
    plot_average=False,
)

# => UED Cross Section <= #

# Compute the "simple" form of the UED cross section in R
R = np.linspace(1.0, 6.0, 50)
ai.compute_ued_simple(traj, "UED", R=R, alpha=8.0)

# Plot the heat map of the UED cross section detailed above
ai.plot_vector(
    "U.pdf", traj, "UED", y=R, ylabel=r"$R [\AA{}]$", time_units="fs", diff=True
)
