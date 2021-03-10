import os
import re
import sys

import numpy as np


def get_dyson(out_dir, traj, njobs):

    # TODO: verify frame corresponds to molpro run based on ID

    frames_per_job = len(traj.frames) / njobs
    nframes = len(traj.frames)
    for job in range(njobs + 1):
        ind1 = frames_per_job * job
        ind2 = frames_per_job * (job + 1)
        if job == njobs:
            ind2 = nframes
        for ind, frame in enumerate(traj.frames[ind1:ind2]):
            data = np.loadtxt("%s%02djob/%04d/norms.out" % (out_dir, job, ind + ind1))
            dyson_norms = {}
            if len(data) == 3:
                dyson_norms = np.array(
                    [
                        [int(data[1]) + 1, data[2]],
                    ]
                )
                frame.properties["dyson_norms"] = dyson_norms
            else:
                dyson_norms = np.array([(int(d[1]) + 1, d[2]) for d in data])
            frame.properties["dyson_norms"] = dyson_norms

    return traj
