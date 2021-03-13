import os
import re
import sys

import numpy as np

import aimsprop.atom_data as atom_data
import aimsprop.manage_xyz as manage_xyz

# TODO: make more invariant to filename (sbatch.sh)

# run molpro calculations to obtain ionization potential information
def prep_dyson(
    opt_in,
):

    """Runs lightspeed calculations from template files and traj object

    Params:
        opt_in - [dict], options to override defaults

    """

    options = {
        # all paths must be absolute!
        "submit_input": None,  # sbatch file to edit
        "submit_input2": "sbatch.sh",  # modified sbatch file to execute
        "out_dir": None,  # dir for output
        "submit": False,  # whether to submit or just prep the job
        "traj": None,  # traj object to run (all frames)
        "njobs": 25,  # number of jobs
    }

    for key, val in list(opt_in.items()):
        options[key] = val

    # TODO: throw error messages for invalid combination of options
    # override options
    submit_input = options["submit_input"]
    submit_input2 = options["submit_input2"]

    out_dir = options["out_dir"]
    submit = options["submit"]
    traj = options["traj"]
    njobs = options["njobs"]

    # loop over jobs
    frames_per_job = len(traj.frames) / njobs
    nframes = len(traj.frames)
    for job in range(njobs + 1):

        # setup directories
        os.chdir(out_dir)
        os.system("mkdir %02djob" % job)
        os.system("cp %s %02djob" % (submit_input, job))
        os.system("cp ../ref/dyson_fomo.py %02djob" % job)
        os.chdir("%02djob" % job)

        # divide frames per job (nframes % njobs goes into last job)
        ind1 = frames_per_job * job
        ind2 = frames_per_job * (job + 1)
        if job == njobs:
            ind2 = nframes

        states = [frame.I - 1 for frame in traj.frames[ind1:ind2]]
        nstates = "[ "
        for state in states:
            nstates += str(state) + ", "
        nstates += "]"
        # update submission script
        update_sbatch("%d" % job, submit_input2, ind1=ind1, ind2=ind2)
        update_dyson(nstates=nstates)

        # loop over frames
        for ind, frame in enumerate(traj.frames[ind1:ind2]):
            # set up directory
            os.system("mkdir %04d" % (ind + ind1))
            os.chdir("%04d" % (ind + ind1))

            # modify neutral molpro params
            symbols = [atom_data.atom_symbol_table[N] for N in frame.N]
            write_xyz(frame, symbols, "geom.xyz")

            # create file to identify frame with job
            os.system(
                "echo 'labels: %4d %4d\nt: %12.6f' > %s"
                % (frame.label[0], frame.label[1], frame.t, "%04d_ID.txt" % ind)
            )

            # return to working directory
            os.chdir("../")

        # if requested, submit the job
        if submit == True:
            os.system("sbatch %s" % submit_input2)
        os.chdir("../")


def update_sbatch(
    new_txt,
    sbatch_in,
    ind1,
    ind2,
):

    # jobname
    os.system('sed -i "s/ID/%s/g" %s' % (new_txt, sbatch_in))

    # frame indices
    os.system('sed -i "s/ind1/%s/g" %s' % (ind1, sbatch_in))
    os.system('sed -i "s/ind2/%s/g" %s' % (ind2, sbatch_in))


def update_dyson(
    nstates,
    filename="dyson_fomo.py",
):

    os.system('sed -i "s/neutral_state_values/%s/g" %s' % (nstates, filename))


def write_xyz(
    frame,
    symbols,
    filename,
):

    geom = []
    for ind, atom in enumerate(frame.xyz):
        geom.append((symbols[ind], atom[0], atom[1], atom[2]))
    manage_xyz.write_xyz(filename, geom, scale=1.0)
