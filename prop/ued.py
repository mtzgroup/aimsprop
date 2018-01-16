import numpy as np
import math

# UED cross sections (TODO: Fix - more atoms, where does it come from)
_ued_cross_sections = {
     1 : 0.0,
     6 : 8.75e-3,
     9 : 1.03e-2,
    53 : 2.07e-1,
}

def compute_ued_simple(
    traj,
    key,
    R,
    alpha,
    ABpairs=None,
    ):

    """ Compute the simple pairwise-distance form of the UED cross section,
        with Gaussian blurring in R.

    Params:
        traj (Trajectory) - the Trajectory object to compute the property for (modified in
            place)
        key (str) - the name of the property
        R (np.ndarray of distances) - the distances to collocate the
            UED cross section to.
        alpha (float) - the Guassian blurring exponent
        ABpairs (list of (int, int)) - a restricted list of atom pair indices
            to include in the computation, or None for all atom pairs.
    Result/Return:
        traj - reference to the input Trajectory object. The property
            key is set to computed UED property.
    """
    
    for frame in traj.frames: 
        N = frame.N
        xyz = frame.xyz
        # Which pair indices?
        if ABpairs is None:
            ABpairs2 = []
            for A in range(len(N)):
                for B in range(A):
                    ABpairs2.append((A,B))
        else:
            ABpairs2 = ABpairs
        # Compute UED cross section
        V = np.zeros_like(R)
        for A,B in ABpairs2:
            rAB = xyz[A,:] - xyz[B,:]
            RAB = math.sqrt(sum(rAB**2))
            SAB = math.sqrt(_ued_cross_sections[N[A]] * _ued_cross_sections[N[B]]) / RAB
            V += SAB * math.sqrt(alpha / math.pi) * np.exp(-alpha * (R - RAB)**2)
        frame.properties[key] = V
    return traj
            
        
        

    
    
