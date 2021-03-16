import math

import numpy as np

from .traj import Trajectory

# UED cross sections (computed by ELSEPA for a 3.7 MeV e- beam with default settings)
_ued_cross_sections = {
    1: 3.92943e-04,
    2: 5.96348e-04,
    3: 3.89833e-03,
    4: 6.17327e-03,
    5: 7.76737e-03,
    6: 8.74560e-03,
    7: 9.42320e-03,
    8: 9.92602e-03,
    9: 1.03156e-02,
    10: 1.06265e-02,
    11: 1.69220e-02,
    12: 2.20789e-02,
    13: 2.80195e-02,
    14: 3.22411e-02,
    15: 3.54220e-02,
    16: 3.79121e-02,
    17: 3.99156e-02,
    18: 4.15608e-02,
    19: 5.48441e-02,
    20: 6.56685e-02,
    21: 6.76687e-02,
    22: 6.88909e-02,
    23: 6.97234e-02,
    24: 6.49267e-02,
    25: 7.07457e-02,
    26: 7.10577e-02,
    27: 7.12812e-02,
    28: 7.14359e-02,
    29: 6.67471e-02,
    30: 7.15914e-02,
    31: 7.91437e-02,
    32: 8.51549e-02,
    33: 9.02497e-02,
    34: 9.46627e-02,
    35: 9.85306e-02,
    36: 1.01948e-01,
    37: 1.20694e-01,
    38: 1.36372e-01,
    39: 1.42990e-01,
    40: 1.47529e-01,
    41: 1.44643e-01,
    42: 1.47227e-01,
    43: 1.49397e-01,
    44: 1.51232e-01,
    45: 1.52791e-01,
    46: 1.47081e-01,
    47: 1.55245e-01,
    48: 1.63144e-01,
    49: 1.74926e-01,
    50: 1.84575e-01,
    51: 1.92955e-01,
    52: 2.00383e-01,
    53: 2.07039e-01,
    54: 2.13039e-01,
    55: 2.40272e-01,
    56: 2.62970e-01,
    57: 2.73268e-01,
    58: 2.64265e-01,
    59: 2.64055e-01,
    60: 2.63588e-01,
    61: 2.62944e-01,
    62: 2.62170e-01,
    63: 2.61295e-01,
    64: 2.68502e-01,
    65: 2.59327e-01,
    66: 2.58262e-01,
    67: 2.57156e-01,
    68: 2.56016e-01,
    69: 2.54849e-01,
    70: 2.53659e-01,
    71: 2.60687e-01,
    72: 2.65547e-01,
    73: 2.69569e-01,
    74: 2.73027e-01,
    75: 2.76042e-01,
    76: 2.78691e-01,
    77: 2.81022e-01,
    78: 2.76923e-01,
    79: 2.78661e-01,
    80: 2.86460e-01,
    81: 3.00666e-01,
    82: 3.12359e-01,
    83: 3.22665e-01,
    84: 3.31940e-01,
    85: 3.40371e-01,
    86: 3.48076e-01,
    87: 3.78187e-01,
    88: 4.03532e-01,
    89: 4.18951e-01,
    90: 4.30842e-01,
    91: 4.24330e-01,
    92: 4.25599e-01,
    93: 4.26351e-01,
    94: 4.17340e-01,
}


def compute_ued_simple(
    traj: Trajectory,
    key: str,
    R: np.ndarray,
    alpha: float,
    ABpairs=None,
) -> Trajectory:
    """Compute the simple pairwise-distance form of the UED cross section,
        with Gaussian blurring in R.

    Params:
        traj: the Trajectory object to compute the property for (modified in
            place)
        key: the name of the property
        R: the distances to collocate the
            UED cross section to.
        alpha: the Guassian blurring exponent
        ABpairs: a restricted list of atom pair indices
            to include in the computation, or None for all atom pairs.
    Return:
        traj: reference to the input Trajectory object. The property
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
                    ABpairs2.append((A, B))
        else:
            ABpairs2 = ABpairs
        # Compute UED cross section
        V = np.zeros_like(R)
        for A, B in ABpairs2:
            rAB = xyz[A, :] - xyz[B, :]
            RAB = math.sqrt(sum(rAB ** 2))
            SAB = math.sqrt(_ued_cross_sections[N[A]] * _ued_cross_sections[N[B]]) / RAB
            V += SAB * math.sqrt(alpha / math.pi) * np.exp(-alpha * (R - RAB) ** 2)
        frame.properties[key] = V
    return traj
