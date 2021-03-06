import numpy as np
from . import traj
import glob
import re
import os
from . import atom_data

_N_table = { val: key for key, val in atom_data.atom_symbol_table.items() }

def parse_fms90(
    filepath,
    scheme='mulliken',
    cutoff_time=None,
    cutoff_saddle=1.0E-4,
    initial_I=None,
    ):

    """ Parse an FMS90 results directory into a Trajectory.

    This uses information in positions.*.xyz, Amps.*, S.dat, and Spawn.log to
    populate a density matrix by mulliken or saddle point rules.

    Params:
        filepath (str) - the path to the FMS run directory
        scheme (str) - 'mulliken' or 'saddle' to indicate the approximation
            used for property evaluation.
        cutoff_time (float) - cutoff time to stop reading trajectory info after
            (None reads all times).
        cutoff_saddle (float) - cutoff for centroid TBF pair in the saddle
            point approach.
        initial_I (int) - initial electronic state, used only if there is
            no Spawn.log (e.g., if no spawning has happened yet) to place
            electronic label.
    Returns:
        trajectory (Trajectory) - the Trajectory object.
    """

    # The files (these are standard textual output from FMS90)
    posfiles = glob.glob('%s/positions.*.xyz' % filepath)
    Cfiles = glob.glob('%s/Amp.*' % filepath)
    Sfile = '%s/S.dat' % filepath
    Spawnfile = '%s/Spawn.log' % filepath
    if len(posfiles) != len(Cfiles):
        raise RuntimeError('xyz and C files not same number of TBF')

    # Make a dict to index files
    posfiles = { int(re.match(r'%s/positions\.(\d+)\.xyz' % filepath, path).group(1)) : path for path in posfiles }
    Cfiles = { int(re.match(r'%s/Amp\.(\d+)' % filepath, path).group(1)) : path for path in Cfiles }

    # Read the Amp and positions files
    C2s = {}
    N2s = {}
    xyz2s = {}
    for I, Cfile in Cfiles.items():
        posfile = posfiles[I]
        lines = open(Cfile).readlines()[1:]
        Cs = {}
        for line in lines:
            # (t) Norm (Real) (Imag)
            mobj = re.match(r'^\s*(\S+)\s+\S+\s+(\S+)\s+(\S+)\s*$', line)
            if mobj:
                if cutoff_time and float(mobj.group(1)) > cutoff_time: continue
                Cs[float(mobj.group(1))] = complex(float(mobj.group(2)), float(mobj.group(3)))
            else:
                raise RuntimeError('Invalid Amp line: %s' % line)
        
        Ns = {}
        xyzs = {}
        lines = open(posfile).readlines()
        natom = int(lines[0])
        if len(lines) % (natom + 2) != 0:
            raise RuntimeError('Position file does not have correct number of lines.')
        nframe = len(lines) // (natom + 2)
        for A in range(nframe):
            lines2 = lines[A * (natom + 2) : (A+1) * (natom+2)]
            mobj = re.match(r'^\s*Time:\s+(\S+),\s+Trajectory:\d+\s*$', lines2[1])
            t = float(mobj.group(1))
            if cutoff_time and t > cutoff_time: continue
            N = []
            xyz = np.zeros((natom, 3))
            for A, line in enumerate(lines2[2:]):
                mobj = re.match(r'^\s*(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s*$', line)
                N.append(_N_table[mobj.group(1).upper()])
                xyz[A,0] = float(mobj.group(2))
                xyz[A,1] = float(mobj.group(3))
                xyz[A,2] = float(mobj.group(4))
            Ns[t] = N
            xyzs[t] = xyz
        C2s[I] = Cs
        N2s[I] = Ns
        xyz2s[I] = xyzs

    # Read the overlap matrix
    Ss = {}
    lines = open(Sfile).readlines()
    tinds = []
    ts = []
    for lind, line in enumerate(lines):
        mobj = re.match(r'^\s*(\S+)\s*$', line)
        if mobj:
            tinds.append(lind)
            ts.append(float(mobj.group(1)))
    tinds.append(len(lines))
    for A, t in enumerate(ts):
        if cutoff_time and t > cutoff_time: continue
        lines2 = lines[tinds[A]+1:tinds[A+1]]
        # This thing wraps at a certain point. 
        # Just read in a 1D array, and then reshape to 2D
        Spart = []
        for line in lines2:
            toks = line.split()
            for B in range(len(toks)):
                if B % 2 != 0: continue
                Sval = complex(float(toks[B]), float(toks[B+1]))
                Spart.append(Sval)
        Smat = np.array(Spart)
        # Reshape to square array
        n = int(np.sqrt(len(Smat)))
        if n**2 != len(Smat): raise RuntimeError('S matrix is not square')
        Smat = np.reshape(Smat, (n,n))
        Ss[t] = Smat

    # Read the Spawn.log to figure out electronic states (not read
    states = {}
    if not os.path.isfile(Spawnfile):
        if initial_I is None: raise RuntimeError("No Spawn.log and no initial_I")
        states = { k : initial_I for k in list(Cfiles.keys()) }
    else:
        lines = open(Spawnfile).readlines()[1:]
        for lind, line in enumerate(lines):
            # match groups are TBF ID, state, parent TBF ID, parent TBF state
            mobj = re.match(r'^\s*\S+\s+\S+\s+\S+\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)', line)
            states[int(mobj.group(1))] = int(mobj.group(2))
            states[int(mobj.group(3))] = int(mobj.group(4)) # Redundant, but captures IC TBFs

    # Swap to put time on slow axis (C3s[t][I] instead of C2s[I][t])
    C3s = {}
    for I, C2 in C2s.items():
        for t, C in C2.items():
            C3s.setdefault(t, {})[I] = C
    xyz3s = {}
    for I, xyz2 in xyz2s.items():
        for t, xyz in xyz2.items():
            xyz3s.setdefault(t, {})[I] = xyz
    N3s = {}
    for I, N2 in N2s.items():
        for t, N in N2.items():
            N3s.setdefault(t, {})[I] = N

    # Build Frames from parsed data
    frames = []
    for t, S in Ss.items():
        if t not in C3s:
            # Sometimes timestamps do not match because Amp.* only holds 2 decimal digits
            # E.g., 1000.875 (in S) vs. 1000.88 (in Amp)
            print('Warning: Time %r not in amplitudes (OK if very small adaptive timestep)' % t)
            continue
        Cs = C3s[t]
        xyzs = xyz3s[t]
        Ns = N3s[t]
        Isort = list(sorted(Cs.keys()))
        if scheme == 'mulliken':
            for I2, I in enumerate(Isort):
                q = 0.0
                for J2, J in enumerate(Isort):
                    if states[I] != states[J]: continue # Electronic orthogonality
                    q += np.real(0.5 * np.conj(Cs[I]) * S[I2, J2] * Cs[J] + \
                                 0.5 * np.conj(Cs[J]) * S[J2, I2] * Cs[I])
                frame = traj.Frame(
                    label=I,
                    t=t,
                    w=q,
                    I=states[I],
                    N=Ns[I],
                    xyz=xyzs[I],
                    )
                frames.append(frame)
        elif scheme == 'saddle':
            for I2, I in enumerate(Isort):
                for J2, J in enumerate(Isort):
                    if states[I] != states[J]: continue # Electronic orthogonality
                    if J > I: continue
                    q = np.real(0.5 * np.conj(Cs[I]) * S[I2, J2] * Cs[J] + \
                                0.5 * np.conj(Cs[J]) * S[J2, I2] * Cs[I])
                    if J < I: q *= 2.0 # Upper/lower triangle
                    if abs(q) < cutoff_saddle: continue # Vanishing weight
                    frame = traj.Frame(
                        label=(I,J),
                        t=t,
                        w=q,
                        I=states[I],
                        N=Ns[I],
                        xyz=0.5*(xyzs[I]+xyzs[J]), # centroid
                        )
                    frames.append(frame)
        else:
            raise RuntimeError('Invalid scheme: %s' % scheme)
    
    trajectory = traj.Trajectory(frames)
    return trajectory

def parse_fms90_dumpfile(
    filepath,
    scheme='mulliken',
    cutoff_time=None,
    cutoff_saddle=1.0E-4,
    ):

    """ Parse an FMS90 results directory into a Trajectory.

    This uses information in positions.*.xyz, S.dat, and TrajDump to
    populate a density matrix by mulliken or saddle point rules.

    Note: TrajDump contains a lower level of precision than other output files

    Params:
        filepath (str) - the path to the FMS run directory
        scheme (str) - 'mulliken' or 'saddle' to indicate the approximation
            used for property evaluation.
        cutoff_time (float) - cutoff time to stop reading trajectory info after
            (None reads all times).
        cutoff_saddle (float) - cutoff for centroid TBF pair in the saddle
            point approach.
    Returns:
        trajectory (Trajectory) - the Trajectory object.
    """

    # The files (these are standard textual output from FMS90)
    posfiles = glob.glob('%s/positions.*.xyz' % filepath)
    dumpfiles = glob.glob('%s/TrajDump.*' % filepath)
    Sfile = '%s/S.dat' % filepath
    if len(posfiles) != len(dumpfiles):
        raise RuntimeError('xyz and dump files not same number of TBF')

    # Make a dict to index files
    posfiles = { int(re.match(r'%s/positions\.(\d+)\.xyz' % filepath, path).group(1)) : path for path in posfiles }
    dumpfiles = { int(re.match(r'%s/TrajDump\.(\d+)' % filepath, path).group(1)) : path for path in dumpfiles }

    # Read the Amp and positions files
    C2s = {}
    N2s = {}
    xyz2s = {}
    states = {}
    for I, dumpfile in dumpfiles.items():
        posfile = posfiles[I]
        data = np.loadtxt(dumpfile, skiprows=1, ndmin=2)
        ts = data[:,0]
        if cutoff_time is not None:
            cut_tind = np.array(np.where(ts >= cutoff_time)).flatten()[0]
            ts = data[:cut_tind, 0]
        C2s[I] = {t: complex(data[tind,-4], data[tind,-3]) for tind, t in enumerate(ts)}
        states[I] = data[0,-1]

        Ns = {}
        xyzs = {}
        lines = open(posfile).readlines()
        natom = int(lines[0])
        if len(lines) % (natom + 2) != 0:
            raise RuntimeError('Position file does not have correct number of lines.')
        nframe = len(lines) / (natom + 2)
        for A in range(nframe):
            lines2 = lines[A * (natom + 2) : (A+1) * (natom+2)]
            mobj = re.match(r'^\s*Time:\s+(\S+),\s+Trajectory:\d+\s*$', lines2[1])
            t = float(mobj.group(1))
            if cutoff_time and t > cutoff_time: continue
            N = []
            xyz = np.zeros((natom, 3))
            for A, line in enumerate(lines2[2:]):
                mobj = re.match(r'^\s*(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s*$', line)
                N.append(_N_table[mobj.group(1).upper()])
                xyz[A,0] = float(mobj.group(2))
                xyz[A,1] = float(mobj.group(3))
                xyz[A,2] = float(mobj.group(4))
            Ns[t] = N
            xyzs[t] = xyz
        N2s[I] = Ns
        xyz2s[I] = xyzs

    # Read the overlap matrix
    Ss = {}
    lines = open(Sfile).readlines()
    tinds = []
    ts = []
    for lind, line in enumerate(lines):
        mobj = re.match(r'^\s*(\S+)\s*$', line)
        if mobj:
            tinds.append(lind)
            ts.append(float(mobj.group(1)))
    tinds.append(len(lines))
    for A, t in enumerate(ts):
        if cutoff_time and t > cutoff_time: continue
        lines2 = lines[tinds[A]+1:tinds[A+1]]
        # This thing wraps at a certain point. 
        # Just read in a 1D array, and then reshape to 2D
        Spart = []
        for line in lines2:
            toks = line.split()
            for B in range(len(toks)):
                if B % 2 != 0: continue
                Sval = complex(float(toks[B]), float(toks[B+1]))
                Spart.append(Sval)
        Smat = np.array(Spart)
        # Reshape to square array
        n = int(np.sqrt(len(Smat)))
        if n**2 != len(Smat): raise RuntimeError('S matrix is not square')
        Smat = np.reshape(Smat, (n,n))
        Ss[t] = Smat

    # Swap to put time on slow axis (C3s[t][I] instead of C2s[I][t])
    C3s = {}
    for I, C2 in C2s.items():
        for t, C in C2.items():
            C3s.setdefault(t, {})[I] = C
    xyz3s = {}
    for I, xyz2 in xyz2s.items():
        for t, xyz in xyz2.items():
            xyz3s.setdefault(t, {})[I] = xyz
    N3s = {}
    for I, N2 in N2s.items():
        for t, N in N2.items():
            N3s.setdefault(t, {})[I] = N
    # Build Frames from parsed data
    frames = []
    # for t, S in Ss.iteritems():
    ts = np.sort(list(Ss.keys()))
    for t in ts:
        S = Ss[t]
        # TODO: find a better way around this?
        t = float('%.2f' % t)
        if t not in C3s:
            # Sometimes timestamps do not match because Amp.* only holds 2 decimal digits
            # E.g., 1000.875 (in S) vs. 1000.88 (in Amp)
            print('Warning: Time %r not in amplitudes (OK if very small adaptive timestep)' % t)
            continue
        Cs = C3s[t]
        xyzs = xyz3s[t]
        Ns = N3s[t]
        Isort = list(sorted(Cs.keys()))
        if scheme == 'mulliken':
            for I2, I in enumerate(Isort):
                q = 0.0
                for J2, J in enumerate(Isort):
                    if states[I] != states[J]: continue # Electronic orthogonality
                    q += np.real(0.5 * np.conj(Cs[I]) * S[I2, J2] * Cs[J] + \
                                 0.5 * np.conj(Cs[J]) * S[J2, I2] * Cs[I])
                # print 'I', I,'t', t, 'Ns', Ns.keys()
                frame = traj.Frame(
                    label=I,
                    t=t,
                    w=q,
                    I=states[I],
                    N=Ns[I],
                    xyz=xyzs[I],
                    )
                frames.append(frame)
        elif scheme == 'saddle':
            for I2, I in enumerate(Isort):
                for J2, J in enumerate(Isort):
                    if states[I] != states[J]: continue # Electronic orthogonality
                    if J > I: continue
                    q = np.real(0.5 * np.conj(Cs[I]) * S[I2, J2] * Cs[J] + \
                                0.5 * np.conj(Cs[J]) * S[J2, I2] * Cs[I])
                    if J < I: q *= 2.0 # Upper/lower triangle
                    if abs(q) < cutoff_saddle: continue # Vanishing weight
                    frame = traj.Frame(
                        label=(I,J),
                        t=t,
                        w=q,
                        I=states[I],
                        N=Ns[I],
                        xyz=0.5*(xyzs[I]+xyzs[J]), # centroid
                        )
                    frames.append(frame)
        else:
            raise RuntimeError('Invalid scheme: %s' % scheme)

    trajectory = traj.Trajectory(frames)

    return trajectory

