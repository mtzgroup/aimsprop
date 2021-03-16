import glob
import os
import re

import numpy as np

from . import atom_data, traj

_N_table = {val: key for key, val in list(atom_data.atom_symbol_table.items())}


def parse_positions(filepath, cutoff_time=None):
    """Parse information in positions.*.xyz.

    Arguments:
        filepath (str): the path to the FMS run directory
        cutoff_time (float): cutoff time to stop reading trajectory info after
            (None reads all times).
    Returns:
        (N2s, xyz2s): Tuple of atomic indices and list of xyzs for each timme and state (dict of dict)
    """

    posfiles = glob.glob(f"{filepath}/positions.*.xyz")

    posfiles = {
        int(re.match(r"%s/positions\.(\d+)\.xyz" % filepath, path).group(1)): path
        for path in posfiles
    }

    # Read the Amp and positions files
    N2s = {}  # list of Ns
    xyz2s = {}  # list of xyzs
    for I, posfile in list(posfiles.items()):

        Ns = {}
        xyzs = {}
        lines = open(posfile).readlines()
        natom = int(lines[0])
        if len(lines) % (natom + 2) != 0:
            raise RuntimeError("Position file does not have correct number of lines.")
        nframe = len(lines) // (natom + 2)
        for A in range(nframe):
            lines2 = lines[A * (natom + 2) : (A + 1) * (natom + 2)]
            mobj = re.match(r"^\s*Time:\s+(\S+),\s+Trajectory:\d+\s*$", lines2[1])
            t = float(mobj.group(1))
            if cutoff_time and t > cutoff_time:
                continue
            N = []
            xyz = np.zeros((natom, 3))
            for A, line in enumerate(lines2[2:]):
                mobj = re.match(r"^\s*(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s*$", line)
                N.append(_N_table[mobj.group(1).upper()])
                xyz[A, 0] = float(mobj.group(2))
                xyz[A, 1] = float(mobj.group(3))
                xyz[A, 2] = float(mobj.group(4))
            Ns[t] = N
            xyzs[t] = xyz
        N2s[I] = Ns
        xyz2s[I] = xyzs

        # N2s have the list of atoms repeated - why???

    return N2s, xyz2s


def parse_Cs(filepath, cutoff_time=None):
    """Parse information in Amps.*.

    Arguments:
        filepath (str): the path to the FMS run directory
        cutoff_time (float): cutoff time to stop reading trajectory info after
            (None reads all times).
    Returns:
        C2s: Amps (dict of dicts)

    """

    Cfiles = glob.glob(f"{filepath}/Amp.*")

    Cfiles = {
        int(re.match(r"%s/Amp\.(\d+)" % filepath, path).group(1)): path
        for path in Cfiles
    }

    # Read the Amp files
    C2s = {}

    for I, Cfile in list(Cfiles.items()):

        lines = open(Cfile).readlines()[1:]
        Cs = {}
        for line in lines:
            # (t) Norm (Real) (Imag)
            mobj = re.match(r"^\s*(\S+)\s+\S+\s+(\S+)\s+(\S+)\s*$", line)
            if mobj:
                if cutoff_time and float(mobj.group(1)) > cutoff_time:
                    continue
                Cs[float(mobj.group(1))] = complex(
                    float(mobj.group(2)), float(mobj.group(3))
                )
            else:
                raise RuntimeError(f"Invalid Amp line: {line}")

        C2s[I] = Cs

    return C2s


def parse_Ss(filepath, cutoff_time=None):
    """Parse information in S.dat.

    Arguments:
        filepath (str): the path to the FMS run directory
        cutoff_time (float): cutoff time to stop reading trajectory info after
            (None reads all times).
    Returns:
        S2s: Overlap matrix (dict of dicts)
    """

    Sfile = os.path.join(filepath, "S.dat")

    # Read the overlap matrix
    Ss = {}
    lines = open(Sfile).readlines()
    tinds = []
    ts = []
    for lind, line in enumerate(lines):
        mobj = re.match(r"^\s*(\S+)\s*$", line)
        if mobj:
            tinds.append(lind)
            ts.append(float(mobj.group(1)))
    tinds.append(len(lines))
    for A, t in enumerate(ts):
        if cutoff_time and t > cutoff_time:
            continue
        lines2 = lines[tinds[A] + 1 : tinds[A + 1]]
        # This thing wraps at a certain point.
        # Just read in a 1D array, and then reshape to 2D
        Spart = []
        for line in lines2:
            toks = line.split()
            for B in range(len(toks)):
                if B % 2 != 0:
                    continue
                Sval = complex(float(toks[B]), float(toks[B + 1]))
                Spart.append(Sval)
        Smat = np.array(Spart)
        # Reshape to square array
        n = int(np.sqrt(len(Smat)))
        if n ** 2 != len(Smat):
            raise RuntimeError("S matrix is not square")
        Smat = np.reshape(Smat, (n, n))
        Ss[t] = Smat

    return Ss


def parse_spawnlog(filepath, cutoff_time=None):
    """Parse information in Spawn.log.

    Arguments:
        filepath (str): the path to the FMS run directory
        cutoff_time (float): cutoff time to stop reading trajectory info after
            (None reads all times).
    Returns:
        states: electronic states
    """

    Spawnfile = os.path.join(filepath, "Spawn.log")

    # Read the Spawn.log to figure out electronic states (not read
    states = {}
    if not os.path.isfile(Spawnfile):
        if initial_I is None:
            raise RuntimeError("No Spawn.log and no initial_I")
        states = {k: initial_I for k in list(Cfiles.keys())}
    else:
        lines = open(Spawnfile).readlines()[1:]
        for lind, line in enumerate(lines):
            # match groups are TBF ID, state, parent TBF ID, parent TBF state
            mobj = re.match(
                r"^\s*\S+\s+\S+\s+\S+\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)", line
            )
            states[int(mobj.group(1))] = int(mobj.group(2))
            states[int(mobj.group(3))] = int(
                mobj.group(4)
            )  # Redundant, but captures IC TBFs

    return states


def swap_dic_axis(x2s):
    """Swap axis of dict to make it faster"""
    x3s = {}
    for I, x2 in list(x2s.items()):
        for t, x in list(x2.items()):
            x3s.setdefault(t, {})[I] = x

    return x3s


def create_frames_mulliken(Ss, C3s, xyz3s, N3s, states):
    """Build list of Frames from parsed data using mulliken scheme

    Arguments:
        Ss (dic): Overlap
        C3s (dic): Amps
        xyz3s (dic): positions
        N3s (dic): atomic numbers
        states: electronic states
    Returns:
        Frames: list of frame objects
    """
    frames = []
    for t, S in list(Ss.items()):
        if t not in C3s:
            # Sometimes timestamps do not match because Amp.* only holds 2 decimal digits
            # E.g., 1000.875 (in S) vs. 1000.88 (in Amp)
            print(
                f"Warning: Time {t} not in amplitudes (OK if very small adaptive timestep)"
            )
            continue

        Cs = C3s[t]
        xyzs = xyz3s[t]
        Ns = N3s[t]

        if len(Cs) != len(S):
            raise RuntimeError("C and S are different sizes")

        Isort = list(sorted(Cs.keys()))

        for I2, I in enumerate(Isort):
            q = 0.0
            for J2, J in enumerate(Isort):
                if states[I] != states[J]:
                    continue  # Electronic orthogonality
                q += np.real(
                    0.5 * np.conj(Cs[I]) * S[I2, J2] * Cs[J]
                    + 0.5 * np.conj(Cs[J]) * S[J2, I2] * Cs[I]
                )
            frame = traj.Frame(
                label=I,
                t=t,
                w=q,
                I=states[I],
                N=Ns[I],
                xyz=xyzs[I],
            )
            frames.append(frame)

    return frames


def create_frames_saddle(Ss, C3s, xyz3s, N3s, states, cutoff_saddle):
    """Build list of Frames from parsed data using saddle scheme

    Arguments:
        Ss (dic): Overlap
        C3s (dic): Amps
        xyz3s (dic): positions
        N3s (dic): atomic numbers
        states: electronic states
        cutoff_saddle (float): cutoff for centroid TBF pair in the saddle point approach.
    Returns:
        Frames: list of frame objects
    """

    frames = []
    for t, S in list(Ss.items()):
        if t not in C3s:
            # Sometimes timestamps do not match because Amp.* only holds 2 decimal digits
            # E.g., 1000.875 (in S) vs. 1000.88 (in Amp)
            print(
                f"Warning: Time {t} not in amplitudes (OK if very small adaptive timestep)"
            )
            continue

        Cs = C3s[t]
        xyzs = xyz3s[t]
        Ns = N3s[t]

        if len(Cs) != len(S):
            raise RuntimeError("C and S are different sizes")

        Isort = list(sorted(Cs.keys()))

        for I2, I in enumerate(Isort):
            for J2, J in enumerate(Isort):
                if states[I] != states[J]:
                    continue  # Electronic orthogonality
                if J > I:
                    continue
                q = np.real(
                    0.5 * np.conj(Cs[I]) * S[I2, J2] * Cs[J]
                    + 0.5 * np.conj(Cs[J]) * S[J2, I2] * Cs[I]
                )
                if J < I:
                    q *= 2.0  # Upper/lower triangle
                if abs(q) < cutoff_saddle:
                    continue  # Vanishing weight
                frame = traj.Frame(
                    label=(I, J),
                    t=t,
                    w=q,
                    I=states[I],
                    N=Ns[I],
                    xyz=0.5 * (xyzs[I] + xyzs[J]),  # centroid
                )
                frames.append(frame)

    return frames


def parse_fms90(
    filepath,
    scheme="mulliken",
    cutoff_time=None,
    cutoff_saddle=1.0e-4,
    initial_I=None,
):

    """Parse an FMS90 results directory into a Trajectory.

    This uses information in positions.*.xyz, Amps.*, S.dat, and Spawn.log to
    populate a density matrix by mulliken or saddle point rules.

    Arguments:
        filepath (str): the path to the FMS run directory
        scheme (str): 'mulliken' or 'saddle' to indicate the approximation
            used for property evaluation.
        cutoff_time (float): cutoff time to stop reading trajectory info after
            (None reads all times).
        cutoff_saddle (float): cutoff for centroid TBF pair in the saddle
            point approach.
        initial_I (int): initial electronic state, used only if there is
            no Spawn.log (e.g., if no spawning has happened yet) to place
            electronic label.
    Returns:
        Trajectory: The Trajectory object.
    """

    # Read in FMS output files: positions*, Amp.*, S.dat and Spawn.log
    N2s, xyz2s = parse_positions(filepath, cutoff_time)

    C2s = parse_Cs(filepath, cutoff_time)

    Ss = parse_Ss(filepath, cutoff_time)

    states = parse_spawnlog(filepath, cutoff_time)

    if len(C2s) != len(N2s):
        raise RuntimeError("xyz and C files not same number of TBF")

    # Swap to put time on slow axis (C3s[t][I] instead of C2s[I][t])
    C3s, xyz3s, N3s = [swap_dic_axis(x) for x in [C2s, xyz2s, N2s]]

    if scheme == "mulliken":
        frames = create_frames_mulliken(Ss, C3s, xyz3s, N3s, states)
    elif scheme == "saddle":
        frames = create_frames_saddle(Ss, C3s, xyz3s, N3s, states, cutoff_saddle)
    else:
        raise RuntimeError(f"Invalid scheme: {scheme}")

    trajectory = traj.Trajectory(frames)

    return trajectory
