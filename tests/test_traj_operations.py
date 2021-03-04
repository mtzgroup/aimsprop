import os
import tempfile

import aimsprop as ai
import numpy as np


def test_major_trajectory_operations():
    # 1.  Parse a series of FMS90 trajectories that Hayley has run for ethylene
    trajs = [ai.parse_fms90(
        '%s/test_data/%04d' %
        (os.path.dirname(__file__), x)) for x in [2, 3]]
    assert len(trajs) == 2

    # 2. Merge the trajectories into one super-big Trajectory with uniform weights
    traj = ai.Trajectory.merge(
        trajs, ws=[1.0 / len(trajs)] * len(trajs), labels=[2, 3])
    assert traj.labels == [(2, 1), (2, 2), (2, 3), (3, 1), (3, 2)]

    # 3. Interpolate trajectory with ~1 fs intervals, removing adaptive timesteps from AIMS
    ts = np.arange(0.0, max(traj.ts), 40.0)
    traj = traj.interpolate_nearest(ts)

    # 4. Plot Bond Distances (Spaghetti + Blur)
    # Compute the a bond distance property for all Frames in traj

    tmp_dir = tempfile.gettempdir()
    ai.compute_bond(traj, 'R01', 0, 1)
    ai.plot_scalar(tmp_dir + 'R.pdf', traj, 'R01', ylabel=r'$R_{CC} [\AA{}]$', time_units='fs', state_colors=[
                   'r', 'b'], plot_average=True)

    # Blur the bond distance (convolve)
    R = np.linspace(0.5, 3.0, 50)
    ai.blur_property(traj, 'R01', 'Rblur', R, alpha=8.0)
    #    Plot the heat map of blurred bond distance
    ai.plot_vector(tmp_dir + 'Rblur.pdf', traj, 'Rblur', y=R,
                   ylabel=r'$R [\AA{}]$', time_units='fs', nlevel=64)

    # 5. Plot of Torsion Angle (Spaghetti + Blur)

    # Compute the a torsion angle property for all Frames in traj
    ai.compute_torsion(traj, 'T0123', 0, 1, 2, 3)
    ai.unwrap_property(traj, 'T0123', 360.0)
    ai.plot_scalar(tmp_dir + 'T.pdf', traj, 'T0123', ylabel=r'$\Theta [^{\circ{}}]$', time_units='fs', state_colors=[
                   'r', 'b'], plot_average=True)

    # Blur the torsion
    T = np.linspace(-180.0, +180.0, 100)
    ai.blur_property(traj, 'T0123', 'Tblur', T, alpha=0.02)
    #    Plot the heat map of blurred torison
    ai.plot_vector(tmp_dir + 'Tblur.pdf', traj, 'Tblur', y=T,
                   ylabel=r'$Theta [^{\circ{}}]$', time_units='fs', nlevel=64)

    # 6. UED Cross Section

    # Compute the "simple" form of the UED cross section in R
    R = np.linspace(1.0, 6.0, 50)
    ai.compute_ued_simple(traj, 'UED', R=R, alpha=8.0)

    # Plot the heat map of the UED cross section detailed above
    ai.plot_vector(tmp_dir + 'UED.pdf', traj, 'UED', y=R,
                   ylabel=r'$R [\AA{}]$', time_units='fs', diff=True)
