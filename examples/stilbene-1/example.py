import aimsprop as ai
import numpy as np

# Parse a series of FMS90 trajectories that Hayley has run for Stilbene
trajs = [ai.parse_fms90('/home/hweir/stilbene/5-aims/aims_%04d/job2' % x) for x in [0,1,2,3]]
# Merge the trajectories into one super-big Trajectory with uniform weights
traj = ai.Trajectory.merge(trajs, [1.0 / len(trajs)] * len(trajs))
# Compute properties at ~1 fs intervals, removing nonsense due to adaptive timesteps
ts = np.arange(0.0, max(traj.ts), 40.0) # TODO: Cleaner edges
traj = traj.interpolate_nearest(ts)

# => Spaghetti Plot of Bond Distances <= #

# Compute the a bond distance property for all Frames in traj
# This particular bond is the one that leads to a three-ring complex
ai.compute_bond(traj, 'R01', 7, 8)

# Plot the "spaghetti plot" of the bond distance computed above
ai.plot_scalar('R.pdf', traj, 'R01', ylabel=r'$R_{CC} [\AA{}]$', time_units='fs', state_colors=['r', 'b'])

# => Showoff Plot of Bond Distances (Spaghetti + Blur) <= #    

# Blur the bond distance
R = np.linspace(1.0,6.0,50)
ai.blur_property(traj, 'R01', 'Rblur', R, alpha=8.0)

# Plot the heat map of blurred bond distance
# ai.plot_vector('R2.pdf', traj, 'Rblur', y=R, ylabel=r'$R [\AA{}]$', time_units='fs', twosided=False, cmap=plt.get_cmap('viridis'))
ai.plot_vector('R2.pdf', traj, 'Rblur', y=R, ylabel=r'$R [\AA{}]$', time_units='fs', nlevel=64)
ai.plot_scalar('R.pdf', traj, 'R01', ylabel=r'$R_{CC} [\AA{}]$', time_units='fs', state_colors=['r', 'b'], clf=False)

# => UED Cross Section <= #

# Compute the "simple" form of the UED cross section in R
R = np.linspace(1.0,6.0,50)
ai.compute_ued_simple(traj, 'UED', R=R, alpha=8.0)

# Plot the heat map of the UED cross section detailed above
ai.plot_vector('U.pdf', traj, 'UED', y=R, ylabel=r'$R [\AA{}]$', time_units='fs', diff=True)
