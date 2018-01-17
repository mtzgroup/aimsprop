import math

""" Compute common geometric properties of Trajectory objects, such as bond
    distances, bond angles, torsion angles, and out-of-plane angles. 
"""
    
def compute_bond(
    traj,
    key,
    A,
    B,
    ):

    """ Compute the a bond-length property for a Trajectory.

    Params:
        traj - the Trajectory object to compute the property for (modified in
            place)
        key - the name of the property
        A - the index of the first atom
        B - the index of the second atom
    Result/Return:
        traj - reference to the input Trajectory object. The property
            key is set to the float value of the bond length for the
            indices A and B.
    """

    for frame in traj.frames:
        xyz = frame.xyz
        rAB = xyz[B,:] - xyz[A,:]
        frame.properties[key] = math.sqrt(sum(rAB**2))
    return traj
            
def compute_angle(
    traj,
    key,
    A,
    B,
    C,
    ):

    """ Compute the a bond-angle property for a Trajectory (in degrees).

    Params:
        traj - the Trajectory object to compute the property for (modified in
            place)
        key - the name of the property
        A - the index of the first atom
        B - the index of the second atom
        C - the index of the third atom
    Result/Return:
        traj - reference to the input Trajectory object. The property
            key is set to the float value of the bond angle for the
            indices A, B and C
    """

    for frame in traj.frames:
        xyz = frame.xyz
        rAB = xyz[B,:] - xyz[A,:]
        rCB = xyz[B,:] - xyz[C,:]
        frame.properties[key] = 180.0 / math.pi * math.acos(sum(rAB * rCB) / math.sqrt(sum(rAB**2) * sum(rCB**2)))
    return traj
            
def compute_torsion(
    traj,
    key,
    A,
    B,
    C,
    D,
    ):

    """ Compute the a torsion-angle property for a Trajectory (in degrees).

    Params:
        traj - the Trajectory object to compute the property for (modified in
            place)
        key - the name of the property
        A - the index of the first atom
        B - the index of the second atom
        C - the index of the third atom
        D - the index of the fourth atom
    Result/Return:
        traj - reference to the input Trajectory object. The property
            key is set to the float value of the torsion angle for the
            indices A, B, C, and D
    """

    # TODO: Google Sherrill Geometry Analysis Program for nice definitions
    raise RuntimeError('Not Implemented')

def compute_oop(
    traj,
    key,
    A,
    B,
    C,
    D,
    ):

    """ Compute the a out-of-plane-angle property for a Trajectory (in degrees).
    
    # TODO: Document which angle this corresponds to 

    Params:
        traj - the Trajectory object to compute the property for (modified in
            place)
        key - the name of the property
        A - the index of the first atom
        B - the index of the second atom
        C - the index of the third atom
        D - the index of the fourth atom
    Result/Return:
        traj - reference to the input Trajectory object. The property
            key is set to the float value of the out-of-plane angle for the
            indices A, B, C, and D
    """

    # TODO
    raise RuntimeError('Not Implemented')

