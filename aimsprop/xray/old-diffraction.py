import lebedev
import factors
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

class Diffraction(object):

    def __init__(
        self,
        traj,
        ):

        self.traj = traj

        # Build atomic form factors
        self.factors = [factors.AtomicFormFactor.build(x.element.symbol.upper()) for x in traj.topology.atoms]

    def compute_diffraction(
        self,
        q,
        nlebedev,
        ):

        leb = lebedev.Lebedev.build(nlebedev)
        qx = np.outer(q, leb.x)
        qy = np.outer(q, leb.y)
        qz = np.outer(q, leb.z)

        R = np.zeros((self.traj.n_frames, len(q)))
    
        for T in range(self.traj.n_frames):
            f = np.zeros(qx.shape, dtype=complex) 
            for A in range(self.traj.n_atoms):
                xA = self.traj.xyz[T,A,0] * 10.0 # In Angstrom
                yA = self.traj.xyz[T,A,1] * 10.0 # In Angstrom
                zA = self.traj.xyz[T,A,2] * 10.0 # In Angstrom
                factorA = self.factors[A]
                f += factorA.evaluate(qx, qy, qz) * np.exp(-1.j * (qx * xA + qy * yA + qz * zA))
            I = (np.abs(f)**2).real
            R[T,:] = np.sum(np.outer(np.ones(q.shape),leb.w / (4.0 * np.pi)) * I, 1)

        return R

    def compute_diffraction_full(
        self,
        q,
        leb,
        ):

        qx = np.outer(q, leb.x)
        qy = np.outer(q, leb.y)
        qz = np.outer(q, leb.z)

        R = np.zeros((self.traj.n_frames, len(q), len(leb.x)))
    
        for T in range(self.traj.n_frames):
            f = np.zeros(qx.shape, dtype=complex) 
            for A in range(self.traj.n_atoms):
                xA = self.traj.xyz[T,A,0] * 10.0 # In Angstrom
                yA = self.traj.xyz[T,A,1] * 10.0 # In Angstrom
                zA = self.traj.xyz[T,A,2] * 10.0 # In Angstrom
                factorA = self.factors[A]
                f += factorA.evaluate(qx, qy, qz) * np.exp(-1.j * (qx * xA + qy * yA + qz * zA))
            I = (np.abs(f)**2).real
            R[T,:,:] = I

        return R
