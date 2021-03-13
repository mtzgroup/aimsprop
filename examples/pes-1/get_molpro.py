import os
import re
import sys

import numpy as np


# TODO: make energies frame property
def get_caspt2_energies(
    fn_in,
):

    ev_per_h = 27.211385
    nm_per_bohr = 0.052917721092

    # TODO: check whether computation finished
    all_lines = open(fn_in).readlines()
    energies = {}
    for n, line in enumerate(all_lines):
        mobj = re.match(r"^\s*\!RSPT2 STATE (\d)\.\d Energy\s*(\S+)", line)
        if mobj:
            state = int(mobj.group(1))
            energy = float(mobj.group(2))
            energies[state] = energy * ev_per_h

    return energies


def get_IP_molpro(
    out_dir,
    cation_out,
    neutral_out,
    traj,
    njobs,
):

    # TODO: verify frame corresponds to molpro run based on ID

    frames_per_job = len(traj.frames) / njobs
    nframes = len(traj.frames)
    for job in range(njobs + 1):
        os.chdir(out_dir)
        os.chdir("IP/%02djob" % job)
        ind1 = frames_per_job * job
        ind2 = frames_per_job * (job + 1)
        if job == njobs:
            ind2 = nframes
        for ind, frame in enumerate(traj.frames[ind1:ind2]):
            cEs = get_caspt2_energies(
                "%sIP/%02djob/%04d/%s" % (out_dir, job, ind + ind1, cation_out)
            )
            nEs = get_caspt2_energies(
                "%sIP/%02djob/%04d/%s" % (out_dir, job, ind + ind1, neutral_out)
            )
            nE = nEs[frame.I]
            IP = np.array([(key, nE - cE) for key, cE in list(cEs.items())])
            frame.properties["IP"] = IP
        os.chdir("../")
    os.chdir("../")

    return traj
