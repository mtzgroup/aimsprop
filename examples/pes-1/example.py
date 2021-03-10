import get_dyson
import get_molpro
import numpy as np
import prep_dyson
import prep_molpro

import aimsprop as ai

# => Parse/Align/Weight/Clean <= #

# Parse a series of FMS90 trajectories that Hayley has run for Stilbene
# Troubles: 11, 12 (not finished), 15
# trajs = [ai.parse_fms90('/home/hweir/stilbene/5-aims/s0_extension/aims_%04d/job1' % x) for x in range(1,16+1) if x not in [11, 12, 15]]
trajs = [
    ai.parse_fms90("/home/hweir/stilbene/5-aims/s0_extension/aims_%04d/job1" % x)
    for x in [1]
]
# TODO: Align trajectories to IC transition dipole moment (on z) and weight by oscillator strength at IC
# Merge the trajectories into one super-big Trajectory with uniform weights
traj = ai.Trajectory.merge(trajs, [1.0 / len(trajs)] * len(trajs))
# Compute properties at ~10 fs intervals, removing nonsense due to adaptive timesteps
ts = np.arange(0.0, max(traj.ts), 400.0)  # TODO: Cleaner edges
traj = traj.interpolate_nearest(ts)

# => Obtain Ionization Potential with Molpro (user defined process) <= #

# divide molpro calculations into reasonable number of jobs
njobs = 25

# molpro reference input files
neutral_input = "/home/monikaw/examples/pes_stilbene/ref/neutral.dat"
cation_input = "/home/monikaw/examples/pes_stilbene/ref/cation.dat"

# cation states to run, depending on frame state (1 = groundstate for both molpro and fms90)
# the active space for this system only permits 1 excited state
cation_states = {
    1: 1,
    2: 2,
}

# reference submission scripts
submit = "/home/monikaw/examples/pes_stilbene/ref/sbatch.sh"

# output working directory
out_dir = "/home/monikaw/examples/pes_stilbene/IP/"
opt_in = {
    "neutral_input": neutral_input,
    "cation_input": cation_input,
    "submit_input": submit,
    "out_dir": out_dir,
    "traj": traj,
    "cation_states": cation_states,
    "submit": False,
    "njobs": njobs,
}

# prep and / or run molpro jobs
# prep_molpro.prep_molpro(opt_in)

# assign ionization potentials as frame properties
# output working directory
out_dir = "/home/monikaw/examples/pes_stilbene/"
traj = get_molpro.get_IP_molpro(
    out_dir=out_dir,
    cation_out="cation.out",
    neutral_out="neutral.out",
    traj=traj,
    njobs=njobs,
)

# => Prepare Jobs for Calculation of Dyson Orbitals in Lightspeed <= #

out_dir_dyson = "/home/monikaw/examples/pes_stilbene/dyson/"
submit = "/home/monikaw/examples/pes_stilbene/ref/dyson/sbatch.sh"
opt_dyson = {
    "submit_input": submit,
    "out_dir": out_dir_dyson,
    "traj": traj,
    "submit": False,
    "njobs": njobs,
}
# prep and / or run dyson calculations
# prep_dyson.prep_dyson(opt_dyson)
traj = get_dyson.get_dyson(out_dir_dyson, traj, njobs)

# => Compute PES <= #

# kinetic energy of electron (in eV) (y axis of plot)
eKT = np.linspace(0.0, 15.0, 1000)
# experimental carrier frequency (eV)
carrier_frequency = 3.0
# gaussian blur, broadening parameter
alpha = 1.0
traj = ai.pes.compute_pes(
    traj=traj,
    carrier_frequency=carrier_frequency,
    alpha=alpha,
    eKT=eKT,
)

# => Plot Results <= #

import os

import matplotlib.pyplot as plt

# change to plotting directory
filename = "/home/monikaw/examples/pes_stilbene/pes2D.pdf"
ai.plot.plot_vector(
    filename=filename,
    traj=traj,
    key="pes",
    y=eKT,
    ylabel="$\epsilon [eV]$",
    twosided=False,
    time_units="fs",
    nlevel=65,
    cmap=plt.cm.jet,
)
