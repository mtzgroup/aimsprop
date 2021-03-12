import os
import re
import sys

import numpy as np

import aimsprop.atom_data as atom_data
import aimsprop.manage_xyz as manage_xyz

# TODO: make more invariant to filename (sbatch.sh)

# run molpro calculations to obtain ionization potential information
def prep_molpro(
    opt_in,
):

    """Runs molpro calculations from template files and traj object

    Params:
        opt_in - [dict], options to override defaults

    """

    options = {
        # all paths must be absolute!
        "neutral_input": None,  # molpro input file for neutral
        "cation_input": None,  # molpro input file for cation
        "submit_input": None,  # sbatch file (runs neutral & cation in same job)
        "out_dir": None,  # dir for output
        "submit": False,  # whether to submit or just prep the job
        "cation_states": None,  # dict cation states corresponding to frame states
        "traj": None,  # traj object to run (all frames)
        "neutral_input2": "neutral.dat",  # modified molpro input file for neutral (must match sbatch.sh
        "cation_input2": "cation.dat",  # modified molpro input file for cation
        "submit_input2": "sbatch.sh",  # modified input script
        "njobs": 25,  # modified input script
    }

    for key, val in list(opt_in.items()):
        options[key] = val

    # TODO: throw error messages for invalid combination of options
    # override options
    neutral_input = options["neutral_input"]
    cation_input = options["cation_input"]
    submit_input = options["submit_input"]

    neutral_input2 = options["neutral_input2"]
    cation_input2 = options["cation_input2"]
    submit_input2 = options["submit_input2"]

    out_dir = options["out_dir"]
    cation_states = options["cation_states"]
    submit = options["submit"]
    traj = options["traj"]
    njobs = options["njobs"]

    frames_per_job = len(traj.frames) / njobs
    nframes = len(traj.frames)
    for job in range(njobs + 1):
        os.chdir(out_dir)
        os.system("mkdir %02djob" % job)
        os.system("cp %s %02djob" % (submit_input, job))
        os.system("cp run_molpro.py %02djob" % job)
        os.chdir("%02djob" % job)
        ind1 = frames_per_job * job
        ind2 = frames_per_job * (job + 1)
        if job == njobs:
            ind2 = nframes
        # update submission script
        update_sbatch("%d" % job, submit_input2, ind1=ind1, ind2=ind2)
        for ind, frame in enumerate(traj.frames[ind1:ind2]):
            # set up directory
            os.system("mkdir %04d" % (ind + ind1))
            os.chdir("%04d" % (ind + ind1))

            # modify neutral molpro params
            symbols = [atom_data.atom_symbol_table[N] for N in frame.N]
            insert_geom(frame.xyz, symbols, neutral_input, neutral_input2)
            update_state(frame.I, neutral_input2)
            write_xyz(frame, symbols, "geom.xyz")

            # modify cation molpro params
            symbols = [atom_data.atom_symbol_table[N] for N in frame.N]
            insert_geom(frame.xyz, symbols, cation_input, cation_input2)
            cI = cation_states[frame.I]
            update_state(cI, cation_input2)

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


def insert_geom(
    xyz,
    symbols,
    molpro_in,
    molpro_out,
):

    lines = open(molpro_in).readlines()
    fh = open(molpro_out, "w")
    for line in lines:
        mobj = re.match("\s*geometry=\{\s*", line)
        fh.write(line)
        if mobj:
            for ind, atom in enumerate(xyz):
                fh.write(
                    "%6s %12.6f %12.6f %12.6f\n"
                    % (symbols[ind], atom[0], atom[1], atom[2])
                )


def update_state(
    I,
    molpro_in,
):

    # Note: molpro is 1 based in state indexing (1 = groundstate, 2 = first excited state)
    # Note: fms is also 1 based in state indexing
    os.system('sed -i "s/xstate/%s/g" %s' % (I, molpro_in))


def update_sbatch(
    new_txt,
    sbatch_in,
    ind1,
    ind2,
):

    os.system('sed -i "s/ID/%s/g" %s' % (new_txt, sbatch_in))
    os.system('sed -i "s/ind1/%s/g" %s' % (ind1, sbatch_in))
    os.system('sed -i "s/ind2/%s/g" %s' % (ind2, sbatch_in))


def write_xyz(
    frame,
    symbols,
    filename,
):

    geom = []
    for ind, atom in enumerate(frame.xyz):
        geom.append((symbols[ind], atom[0], atom[1], atom[2]))

    manage_xyz.write_xyz(filename, geom)
