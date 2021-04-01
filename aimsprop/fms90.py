import re
from pathlib import Path
from typing import Any, Dict, List

import numpy as np

from . import atom_data, traj

_N_table = {val: key for key, val in list(atom_data.atom_symbol_table.items())}


def _get_file_index(filename: Path) -> int:
    """
    extracts the number X from filename fanculo.X.ext
    """
    return int(filename.name.split(".")[1])


def _parse_positions(filepath: Path, cutoff_time: float = None) -> Dict[Any, Any]:
    """Parse information in positions.*.xyz.

    Arguments:
        filepath: the path to the FMS run directory
        cutoff_time: cutoff time to stop reading trajectory info after
            (None reads all times).
    Returns:
        (N2s, xyz2s): Tuple of atomic indices and list of xyzs for each timme and state (dict of dict)
    """

    posfiles = (
        (_get_file_index(path), path)
        for path in filepath.iterdir()
        if path.match("positions.*.xyz")
    )

    # Read the Amp and positions files
    N2s = {}  # dict of Ns
    xyz2s = {}  # dict of xyzs
    for index, posfile in posfiles:
        Ns = {}
        xyzs = {}
        with open(posfile) as f:
            lines = f.readlines()
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
        N2s[index] = Ns
        xyz2s[index] = xyzs

        # N2s have the list of atoms repeated - why???

    return N2s, xyz2s


def _parse_Cs(filepath: Path, cutoff_time: float = None) -> Dict[Any, Any]:
    """Parse information in Amps.*.

    Arguments:
        filepath: the path to the FMS run directory
        cutoff_time: cutoff time to stop reading trajectory info after
            (None reads all times).
    Returns:
        C2s: Amps (dict of dicts)

    """

    Cfiles = (
        (_get_file_index(path), path)
        for path in filepath.iterdir()
        if path.match("Amp.*")
    )
    C2s = {}
    for index, Cfile in Cfiles:

        # AV 2021, I do not want to further modify this parser, but of course it is stupid
        # to read the file in lines and THEN loop over them.
        with open(Cfile) as f:
            lines = f.readlines()[1:]

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

        C2s[index] = Cs

    return C2s


def _parse_Ss(filepath: Path, cutoff_time: float = None) -> Dict[Any, Any]:
    """Parse information in S.dat.

    Arguments:
        filepath: the path to the FMS run directory
        cutoff_time: cutoff time to stop reading trajectory info after
            (None reads all times).
    Returns:
        S2s: Overlap matrix (dict of dicts)
    """

    Sfile = filepath / "S.dat"

    # Read the overlap matrix
    Ss = {}
    with open(Sfile) as f:
        lines = f.readlines()
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


def _parse_spawnlog(filepath: Path, initial_I: int = None) -> Dict[Any, Any]:
    """Parse information in Spawn.log.

    Arguments:
        filepath: the path to the FMS run directory
    Returns:
        states: electronic states
    """

    Spawnfile = filepath / "Spawn.log"

    # Read the Spawn.log to figure out electronic states (not read
    states = {}
    # AV: why this if? Why don't we ALWAYS use the Amp files to get the states variable?
    if not Spawnfile.is_file():
        if initial_I is None:
            raise RuntimeError("No Spawn.log and no initial_I")
        Cfiles_indexes = (
            _get_file_index(path) for path in filepath.iterdir() if path.match("Amp.*")
        )
        states = {k: initial_I for k in Cfiles_indexes}
    else:
        with open(Spawnfile) as f:
            lines = f.readlines()[1:]

    for lind, line in enumerate(lines):
        # match groups are TBF ID, state, parent TBF ID, parent TBF state
        mobj = re.match(r"^\s*\S+\s+\S+\s+\S+\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)", line)
        states[int(mobj.group(1))] = int(mobj.group(2))
        # Redundant, but captures IC TBFs
        states[int(mobj.group(3))] = int(mobj.group(4))

    return states


def _swap_dic_axis(x2s: Dict[Any, Any]) -> Dict[Any, Any]:
    """Swap axis of dict to make it faster"""
    x3s = {}
    for index, x2 in list(x2s.items()):
        for t, x in list(x2.items()):
            x3s.setdefault(t, {})[index] = x

    return x3s


def _create_frames_mulliken(
    Ss: Dict[Any, Any],
    C3s: Dict[Any, Any],
    xyz3s: Dict[Any, Any],
    N3s: Dict[Any, Any],
    states: Dict[Any, Any],
) -> List[traj.Frame]:
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


def _create_frames_saddle(
    Ss: Dict[Any, Any],
    C3s: Dict[Any, Any],
    xyz3s: Dict[Any, Any],
    N3s: Dict[Any, Any],
    states: Dict[Any, Any],
    cutoff_saddle: float,
) -> List[traj.Frame]:
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
    filepath: Path,
    scheme: str = "mulliken",
    cutoff_time: float = None,
    cutoff_saddle: float = 1.0e-4,
    initial_I: int = None,
) -> traj.Trajectory:
    """Parse an FMS90 results directory into a Trajectory.

    This uses information in positions.*.xyz, Amps.*, S.dat, and Spawn.log to
    populate a density matrix by mulliken or saddle point rules.

    Arguments:
        filepath: the path to the FMS run directory
        scheme: 'mulliken' or 'saddle' to indicate the approximation
            used for property evaluation.
        cutoff_time: cutoff time to stop reading trajectory info after
            (None reads all times).
        cutoff_saddle: cutoff for centroid TBF pair in the saddle
            point approach.
        initial_I: initial electronic state, used only if there is
            no Spawn.log (e.g., if no spawning has happened yet) to place
            electronic label.
    Returns:
        Trajectory: The Trajectory object.
    """

    # Read in FMS output files: positions*, Amp.*, S.dat and Spawn.log
    N2s, xyz2s = _parse_positions(filepath, cutoff_time)

    C2s = _parse_Cs(filepath, cutoff_time)

    Ss = _parse_Ss(filepath, cutoff_time)

    states = _parse_spawnlog(filepath, initial_I)

    if len(C2s) != len(N2s):
        raise RuntimeError("xyz and C files not same number of TBF")

    # Swap to put time on slow axis (C3s[t][I] instead of C2s[I][t])
    C3s, xyz3s, N3s = [_swap_dic_axis(x) for x in [C2s, xyz2s, N2s]]

    if scheme == "mulliken":
        frames = _create_frames_mulliken(Ss, C3s, xyz3s, N3s, states)
    elif scheme == "saddle":
        frames = _create_frames_saddle(Ss, C3s, xyz3s, N3s, states, cutoff_saddle)
    else:
        raise RuntimeError(f"Invalid scheme: {scheme}")

    trajectory = traj.Trajectory(frames)

    return trajectory
