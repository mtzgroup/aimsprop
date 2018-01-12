import numpy as np
import traj
import glob
import re

# TODO: Fill atom numbers
_N_table = {
    'H' : 1,
    'C' : 6,
}

def parse_fms90(
    filepath,
    scheme='mulliken',
    ):

    """ Parse an FMS90 results directory into a Trajectory.

    Params:
        filepath (str) - the path to the FMS run directory
        scheme (str) - 'mulliken' or 'saddle' to indicate the approximation
            used for property evaluation.
    Returns:
        trajectory (Trajectory) - the Trajectory object.
    """

    # The files
    posfiles = glob.glob('%s/positions.*.xyz' % filepath)
    Cfiles = glob.glob('%s/Amp.*' % filepath)
    Sfile = '%s/S.dat' % filepath
    if len(posfiles) != len(Cfiles):
        raise RuntimeError('xyz and C files not same number of TBF')

    # Make a dict to index files
    posfiles = { int(re.match(r'%s/positions\.(\d+)\.xyz' % filepath, path).group(1)) : path for path in posfiles }
    Cfiles = { int(re.match(r'%s/Amp\.(\d+)' % filepath, path).group(1)) : path for path in Cfiles }

    C2s = {}
    N2s = {}
    xyz2s = {}
    for I, Cfile in Cfiles.iteritems():
        posfile = posfiles[I]

        lines = open(Cfile).readlines()[1:]
        Cs = {}
        for line in lines:
            mobj = re.match(r'^\s*(\S+)\s+\S+\s+(\S+)\s+(\S+)\s*$', line)
            if mobj:
                Cs[float(mobj.group(1))] = complex(float(mobj.group(2)), float(mobj.group(3)))
            else:
                raise RuntimeError('Invalid Amp line: %s' % line)
        
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
        lines2 = lines[tinds[A]+1:tinds[A+1]]
        Spart = []
        for line in lines2:
            toks = line.split()
            Spart2 = []
            for B in range(len(toks)):
                if B % 2 != 0: continue
                Sval = complex(float(toks[B]), float(toks[B+1]))
                Spart2.append(Sval)
            Spart.append(Spart2)
        Smat = np.array(Spart)
        Ss[t] = Smat

    # Swap to put time on slow axis
    C3s = {}
    for I, C2 in C2s.iteritems():
        for t, C in C2.iteritems():
            C3s.setdefault(t, {})[I] = C
    xyz3s = {}
    for I, xyz2 in xyz2s.iteritems():
        for t, xyz in xyz2.iteritems():
            xyz3s.setdefault(t, {})[I] = xyz
    N3s = {}
    for I, N2 in N2s.iteritems():
        for t, N in N2.iteritems():
            N3s.setdefault(t, {})[I] = N
     
    # Build frames
    frames = []
    for t, S in Ss.iteritems():
        Cs = C3s[t]
        xyzs = xyz3s[t]
        Ns = N3s[t]

        Isort = list(sorted(Cs.keys()))

        if scheme == 'mulliken':
            for I2, I in enumerate(Isort):
                q = 0.0
                for J2, J in enumerate(Isort):
                    q += np.real(0.5 * np.conj(Cs[I]) * S[I2, J2] * Cs[J] + \
                                 0.5 * np.conj(Cs[J]) * S[J2, I2] * Cs[I])
                frame = traj.Frame(
                    label=I,
                    t=t,
                    w=q,
                    I=0,
                    N=Ns,
                    xyz=xyzs,
                    )
                frames.append(frame)
        elif scheme == 'saddle':
            raise NotImplementedError('This')
        else:
            raise RuntimeError('Invalid scheme: %s' % scheme)
    
    trajectory = traj.Trajectory(frames)
    return trajectory

if __name__ == '__main__':

    traj = parse_fms90('/home/hweir/stilbene/5-aims/aims_0000/job2')
    
    for t in traj.ts: 
        print sum([x.w for x in traj.subset_by_t(t).frames])
