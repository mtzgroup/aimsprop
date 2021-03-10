import os

import numpy as np

import aimsprop as ai
from aimsprop import atom_data
from util import manage_xyz


def vis_aims(
    traj,
):

    # => Writing vmd scripts <= #
    os.system("mkdir xyzs")

    # => Write Trajectory Files (for vmd to read) <= #
    filenames = []
    opacities = {}
    states = {}
    for index, label in enumerate(traj.labels):

        # all timesteps for trajectories with label x
        trajI = traj.subset_by_label(label)
        xyzs = np.array([frame.xyz for frame in trajI.frames])
        # Ns = np.array([frame.N for frame in trajI.frames])
        ts = np.array([frame.t for frame in trajI.frames])
        ws = np.array([frame.w for frame in trajI.frames])

        # creating and concatenating dummy frames
        # adding previous steps
        new_ts1 = np.arange(ts[1] - ts[0], ts[0], ts[1] - ts[0])
        new_ws1 = np.zeros_like(new_ts1)
        nts1 = len(new_ts1)
        natoms = np.shape(xyzs)[1]
        new_xyzs1 = np.zeros((nts1, natoms, 3))
        tot_xyzs = np.vstack((new_xyzs1, xyzs))
        tot_ws = np.hstack((new_ws1, ws))

        # output coordinates as xyz file for vmd
        symbols = [atom_data.atom_symbol_table[N] for N in frame.N]
        geoms = manage_xyz.xyzs_to_geom(symbols, tot_xyzs)
        manage_xyz.write_xyzs(
            "xyzs/x%02d-%02d.xyz" % (label[0], label[1]), geoms, scale=1.0
        )
        filenames.append("x%02d-%02d.xyz" % (label[0], label[1]))

        opacities[index] = tot_ws
        states[index] = frame.I - 1
        nframes = len(tot_ws)

    colors = [
        [1.0, 0.0, 0.0],
        [0.0, 0.0, 1.0],
        [0.0, 0.5, 0.5],
    ]
    ai.vis_aims.write_vmd(colors, filenames, states)
    options = {
        "opacities": opacities,
        "end_frame": nframes,
        "increment_frame": 1,
        "sleep": 0.5,
        "ntrajs": len(filenames),
    }
    ai.vis_aims.write_render(options)
    ai.vis_aims.run_vmd_render()
    ai.vis_aims.get_gif()


def test():

    # => Make a trajectory object <= #

    ICs = [3]
    trajs = [
        ai.parse_fms90("/home/monikaw/chem/o-np/6-aims/%04d/1-run" % x) for x in ICs
    ]

    # Merge the trajectories into one super-big Trajectory with uniform weights
    traj = ai.Trajectory.merge(trajs, [1.0 / len(trajs)] * len(trajs), labels=ICs)

    # Compute properties at ~1 fs intervals, removing nonsense due to adaptive timesteps
    ts = np.arange(np.min(traj.ts), np.max(traj.ts), 800.0)
    traj = traj.interpolate_linear(ts)
    print(f"nframes: {len(traj.ts)}")

    vis_aims(
        traj,
    )


if __name__ == "__main__":

    test()
