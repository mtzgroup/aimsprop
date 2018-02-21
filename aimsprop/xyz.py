import numpy as np
import traj
import re
import atom_data

# TODO: Maybe should be in atom data
_N_table = { val: key for key, val in atom_data.atom_symbol_table.iteritems() }

def parse_xyz(
    filename,
    label=1,
    w=1.0,
    I=1,
    t0=0.0,
    dt=20.0,
    ts=None,
    N_table=None,
    ):

    """ Parse an XYZ adiabatic trajectory file directly into a Trajectory.

        filename (str) - the absolute or relative path to the xyz file.
        label (hashable) - the label of this trajectory
        w (float) - the weight of this trajectory
        I (int) - electronic state label
        t0 (float) - the initial time in au
        dt (float) - the timestep in au 
        ts (list of float) - an explicit list of times in au, overrides t0 and dt
        N_table (dict of str : int) - an optional dictionary mapping atomic
            symbol to atomic number, used for non-standard atom names.        

    Returns:
        trajectory (Trajectory) - the Trajectory object.
    """

    lines = open(filename).readlines()

    natom = int(lines[0]) # This should always work
    if len(lines) % (natom+2): raise ValueError('Invalid number of lines in xyz file')
    nframe = len(lines) / (natom+2)
    
    xyzs = []
    Zs = []
    for frame in range(nframe):
        lines2 = lines[frame*(natom+2)+2:(frame+1)*(natom+2)]
        Z = []
        xyz = []
        for line in lines2:
            mobj = re.match(r'^\s*(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s*$',line)
            Z.append(mobj.group(1))
            xyz.append([float(mobj.group(x)) for x in (2,3,4)])
        xyz = np.array(xyz)
        xyzs.append(xyz)
        Zs.append(Z)

    # User symbol table or default?
    N_table2 = N_table if N_table else _N_table

    frames2 = []
    for ind, xyz in enumerate(xyzs):
        Z = Zs[ind]
        Ns = [N_table2[key] for key in Z]
        frame2 = traj.Frame(
            label=label,
            t=dt*ind + t0 if ts is None else ts[ind],
            w=w,
            I=I,
            N=Ns,
            xyz=xyz,
            )
        frames2.append(frame2)

    trajectory = traj.Trajectory(frames2)

    return trajectory
