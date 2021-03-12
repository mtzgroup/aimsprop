import copy
import glob
import re

import matplotlib
import numpy as np

import aimsprop as ai

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def build_traj(
    ICs,
):

    # => Parse/Align/Weight/Clean <= #

    # Parse a series of FMS90 trajectories
    main_dir = "/home/monikaw/chem/stilbene/4-aims"

    broken_frame = {}
    trajs = []
    labels = []
    for IC in ICs:
        input_dirs = [name for name in glob.glob("%s/%04d/*-run" % (main_dir, IC))]
        trajs_IC = []
        for run in input_dirs:
            mobj = re.match("\S+(\d+)-run", run)

        restart_index = np.array(
            [int(re.match("\S+(\d+)-run", run).group(1)) for run in input_dirs]
        )
        restart_inds = np.argsort(
            np.array(
                [int(re.match("\S+(\d+)-run", run).group(1)) for run in input_dirs]
            )
        )
        for ind in restart_inds:
            run = input_dirs[ind]
            ind1 = restart_index[ind]
            print(("Reading:", run))
            if (IC, ind1) in list(broken_frame.keys()):
                traj = ai.parse_fms90_dumpfile(
                    "%s" % run, scheme="mulliken", cutoff_time=broken_frame[(IC, ind1)]
                )
            else:
                traj = ai.parse_fms90_dumpfile("%s" % run, scheme="mulliken")
            trajs_IC.append(copy.copy(traj))
        print("merging restarts")
        traj_IC = ai.Trajectory.merge(
            trajs_IC, [1.0] * len(trajs_IC), [IC] * len(trajs_IC)
        )

        print("removing duplicate frames")
        traj_IC = traj_IC.remove_duplicates()

        print("interpolating across IC")
        ts = np.arange(0.0, max(traj_IC.ts), 400.0)
        traj_IC = traj_IC.interpolate_linear(ts)

        # Each dat contains np.array([[dE (eV), Osc, Tx, Ty, Tz]])
        dat = np.load("/home/monikaw/chem/stilbene/4-spectra/%04d/align.npz" % IC)[
            "dat"
        ]
        # Extract the transition dipole moments for alignment purposes
        tdm = dat[0, 2:]
        # Rotate each AIMS trajectory so the transition dipole moment vector is in aligned in z at t=0
        traj_IC = ai.rotate_frames_to_z(traj_IC, tdm)

        labels += [IC] * len(traj_IC.frames)
        trajs.append(copy.copy(traj_IC))

    # Parse the oscillator strength weights and transition dipole alignments
    # from the 4-spectra dir (computed by the align.py script)
    # Each dat contains np.array([[dE (eV), Osc, Tx, Ty, Tz]])
    dats = [
        np.load("/home/monikaw/chem/stilbene/4-spectra/%04d/align.npz" % x)["dat"]
        for x in ICs
    ]
    # Compute the normalized oscillator strength weights
    ws = np.array([dat[0, 1] for dat in dats])
    ws /= np.sum(ws)
    # Merge the trajectories into one super-big Trajectory with uniform weights
    print("merging ICs")
    traj = ai.Trajectory.merge(trajs, ws)
    print(("nframes", len(traj.frames)))

    return traj, trajs


if __name__ == "__main__":

    time_scale = 1.0 / 41.3413745758  # au to fs

    ICs = np.arange(100, 132, 1)
    traj, trajs = build_traj(ICs)
    ts = traj.ts

    # => UED Signal <= #

    print("calculating raw UED signal")
    R = np.linspace(0.0, 6.0, 200)
    traj = ai.compute_ued_simple(traj, "UED", R=R, alpha=8.0)

    # => Bootstrap Results <= #

    print("resampling trajectories")
    resampled_trajs = ai.bootstrap(traj, 1000, ICs)

    # => Plotting Population with Standard Deviation Error Bars <= #

    # Getting population statistics
    I1 = []
    I2 = []
    for ind, retraj in enumerate(resampled_trajs):
        pop1 = ai.pop.compute_population(retraj)
        I1.append(pop1[1])
        I2.append(pop1[2])
    std1 = np.std(I1, axis=0)
    std2 = np.std(I2, axis=0)
    avg_pop1 = np.average(I1, axis=0)
    avg_pop2 = np.average(I2, axis=0)

    # Plotting population with error bars
    plt.clf()
    plt.errorbar(ts * time_scale, avg_pop1, yerr=std1, fmt="none", ecolor="r")
    plt.errorbar(ts * time_scale, avg_pop2, yerr=std2, fmt="none", ecolor="b")
    plt.plot(ts * time_scale, avg_pop1, "r", label="$S_{0}$")
    plt.plot(ts * time_scale, avg_pop2, "b", label="$S_{1}$")

    plt.xlabel("t [fs]")
    plt.ylabel("Population")
    plt.axis([time_scale * min(ts), max(ts) * time_scale, -0.1, 1.1])
    plt.legend()
    plt.tight_layout()
    plt.savefig("Population.pdf")

    # => Plotting UED stats Heatmap <= #

    print("extracting UED statistics")
    avg, std = ai.extract_stats(resampled_trajs, "UED", diff=True)

    nlevel = 65

    plt.clf()
    vmax = np.max(np.abs(std))
    levels = np.linspace(0.0, vmax, nlevel)
    cticks = [0, +int(vmax)]
    plt.contourf(ts * time_scale, R, std.T, levels=levels)
    plt.colorbar(ticks=cticks)
    plt.xlabel("t [fs]")
    plt.ylabel("R[$\AA$]")
    plt.tight_layout()
    plt.savefig("ued_std.pdf")

    plt.clf()
    vmax = np.max(np.abs(avg))
    levels = np.linspace(-vmax, +vmax, nlevel)
    cticks = [-int(vmax), 0, +int(vmax)]
    plt.contourf(ts * time_scale, R, avg.T, cmap=plt.cm.bwr, levels=levels)
    plt.colorbar(ticks=cticks)
    plt.xlabel("t [fs]")
    plt.ylabel("R[$\AA$]")
    plt.tight_layout()
    plt.savefig("ued_avg.pdf")
