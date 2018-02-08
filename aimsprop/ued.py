import numpy as np
import math

# UED cross sections (computed by ELSEPA for a 3.7 MeV e- beam with default settings)
_ued_cross_sections = {
  1 : 3.92943E-04,
  2 : 5.96348E-04,
  3 : 3.89833E-03,
  4 : 6.17327E-03,
  5 : 7.76737E-03,
  6 : 8.74560E-03,
  7 : 9.42320E-03,
  8 : 9.92602E-03,
  9 : 1.03156E-02,
 10 : 1.06265E-02,
 11 : 1.69220E-02,
 12 : 2.20789E-02,
 13 : 2.80195E-02,
 14 : 3.22411E-02,
 15 : 3.54220E-02,
 16 : 3.79121E-02,
 17 : 3.99156E-02,
 18 : 4.15608E-02,
 19 : 5.48441E-02,
 20 : 6.56685E-02,
 21 : 6.76687E-02,
 22 : 6.88909E-02,
 23 : 6.97234E-02,
 24 : 6.49267E-02,
 25 : 7.07457E-02,
 26 : 7.10577E-02,
 27 : 7.12812E-02,
 28 : 7.14359E-02,
 29 : 6.67471E-02,
 30 : 7.15914E-02,
 31 : 7.91437E-02,
 32 : 8.51549E-02,
 33 : 9.02497E-02,
 34 : 9.46627E-02,
 35 : 9.85306E-02,
 36 : 1.01948E-01,
 37 : 1.20694E-01,
 38 : 1.36372E-01,
 39 : 1.42990E-01,
 40 : 1.47529E-01,
 41 : 1.44643E-01,
 42 : 1.47227E-01,
 43 : 1.49397E-01,
 44 : 1.51232E-01,
 45 : 1.52791E-01,
 46 : 1.47081E-01,
 47 : 1.55245E-01,
 48 : 1.63144E-01,
 49 : 1.74926E-01,
 50 : 1.84575E-01,
 51 : 1.92955E-01,
 52 : 2.00383E-01,
 53 : 2.07039E-01,
 54 : 2.13039E-01,
 55 : 2.40272E-01,
 56 : 2.62970E-01,
 57 : 2.73268E-01,
 58 : 2.64265E-01,
 59 : 2.64055E-01,
 60 : 2.63588E-01,
 61 : 2.62944E-01,
 62 : 2.62170E-01,
 63 : 2.61295E-01,
 64 : 2.68502E-01,
 65 : 2.59327E-01,
 66 : 2.58262E-01,
 67 : 2.57156E-01,
 68 : 2.56016E-01,
 69 : 2.54849E-01,
 70 : 2.53659E-01,
 71 : 2.60687E-01,
 72 : 2.65547E-01,
 73 : 2.69569E-01,
 74 : 2.73027E-01,
 75 : 2.76042E-01,
 76 : 2.78691E-01,
 77 : 2.81022E-01,
 78 : 2.76923E-01,
 79 : 2.78661E-01,
 80 : 2.86460E-01,
 81 : 3.00666E-01,
 82 : 3.12359E-01,
 83 : 3.22665E-01,
 84 : 3.31940E-01,
 85 : 3.40371E-01,
 86 : 3.48076E-01,
 87 : 3.78187E-01,
 88 : 4.03532E-01,
 89 : 4.18951E-01,
 90 : 4.30842E-01,
 91 : 4.24330E-01,
 92 : 4.25599E-01,
 93 : 4.26351E-01,
 94 : 4.17340E-01,
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
            
        
        

    
    
