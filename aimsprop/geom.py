import math

import numpy as np

from .traj import Trajectory

# => Utility Math functions <= #


def _normalize(
    vec,
):
    """ Return a normalized version of vec """

    return vec / math.sqrt(sum(vec ** 2))


def _dot(
    vec1,
    vec2,
):
    """ Dot product between vec1 and vec2 """

    return sum(vec1 * vec2)


def _cross(
    vec1,
    vec2,
):
    """ Cross product between vec1 and vec2 in R^3 """

    vec3 = np.zeros((3,))
    vec3[0] = +(vec1[1] * vec2[2] - vec1[2] * vec2[1])
    vec3[1] = -(vec1[0] * vec2[2] - vec1[2] * vec2[0])
    vec3[2] = +(vec1[0] * vec2[1] - vec1[1] * vec2[0])
    return vec3


# => Geometric Properties <= #

""" Compute common geometric properties of Trajectory objects, such as bond
    distances, bond angles, torsion angles, and out-of-plane angles. 
"""


def compute_bond(
    traj: Trajectory,
    key: str,
    A: int,
    B: int,
) -> Trajectory:
    """Compute the a bond-length property for a Trajectory.

    Params:
        traj: the Trajectory object to compute the property for (modified in
            place)
        key: the name of the property
        A: the index of the first atom
        B: the index of the second atom
    Return:
        traj - reference to the input Trajectory object. The property
            key is set to the float value of the bond length for the
            indices A and B.
    """

    for frame in traj.frames:
        xyz = frame.xyz
        rAB = xyz[B, :] - xyz[A, :]
        frame.properties[key] = math.sqrt(sum(rAB ** 2))
    return traj


def compute_angle(
    traj: Trajectory,
    key: str,
    A: int,
    B: int,
    C: int,
) -> Trajectory:
    """Compute the a bond-angle property for a Trajectory (in degrees).

    Params:
        traj: the Trajectory object to compute the property for (modified in
            place)
        key: the name of the property
        A: the index of the first atom
        B: the index of the second atom
        C: the index of the third atom
    Return:
        traj: reference to the input Trajectory object. The property
            key is set to the float value of the bond angle for the
            indices A, B and C
    """

    for frame in traj.frames:
        xyz = frame.xyz
        rAB = xyz[B, :] - xyz[A, :]
        rCB = xyz[B, :] - xyz[C, :]
        frame.properties[key] = (
            180.0
            / math.pi
            * math.acos(sum(rAB * rCB) / math.sqrt(sum(rAB ** 2) * sum(rCB ** 2)))
        )
    return traj


def compute_torsion(
    traj: Trajectory,
    key: str,
    A: int,
    B: int,
    C: int,
    D: int,
) -> Trajectory:
    """Compute the a torsion-angle property for a Trajectory (in degrees).

    Params:
        traj: the Trajectory object to compute the property for (modified in
            place)
        key: the name of the property
        A: the index of the first atom
        B: the index of the second atom
        C: the index of the third atom
        D: the index of the fourth atom
    Return:
        traj: reference to the input Trajectory object. The property
            key is set to the float value of the torsion angle for the
            indices A, B, C, and D
    """

    for frame in traj.frames:
        xyz = frame.xyz
        rAB = xyz[B, :] - xyz[A, :]
        rBC = xyz[C, :] - xyz[B, :]
        rCD = xyz[D, :] - xyz[C, :]
        eAB = _normalize(rAB)
        eBC = _normalize(rBC)
        eCD = _normalize(rCD)
        n1 = _normalize(_cross(rAB, rBC))
        n2 = _normalize(_cross(rBC, rCD))
        m1 = _cross(n1, eBC)
        x = _dot(n1, n2)
        y = _dot(m1, n2)
        theta = 180.0 / math.pi * math.atan2(y, x)
        frame.properties[key] = theta

    return traj


def compute_oop(
    traj: Trajectory,
    key: str,
    A: int,
    B: int,
    C: int,
    D: int,
) -> Trajectory:
    """Compute the a out-of-plane-angle property for a Trajectory (in degrees).

    Params:
        traj: the Trajectory object to compute the property for (modified in
            place)
        key: the name of the property
        A: the index of the first atom (OOP atom)
        B: the index of the second atom
        C: the index of the third atom
        D: the index of the fourth atom (Defines vector to A)
    Return:
        traj: reference to the input Trajectory object. The property
            key is set to the float value of the out-of-plane angle for the
            indices A, B, C, and D
    """

    traj = compute_angle(traj, "AngleBDC", B, D, C)

    for frame in traj.frames:
        xyz = frame.xyz
        rDA = xyz[A, :] - xyz[D, :]
        rDB = xyz[B, :] - xyz[D, :]
        rDC = xyz[C, :] - xyz[D, :]
        eDA = _normalize(rDA)
        eDB = _normalize(rDB)
        eDC = _normalize(rDC)

        thetaBDC = math.radians(frame.properties["AngleBDC"])

        cross = _cross(eDB, eDC)
        OOP = math.degrees(math.asin(_dot(cross / math.sin(thetaBDC), eDA)))

        frame.properties[key] = OOP

    return traj


def compute_transfer_coord(
    traj: Trajectory,
    key: str,
    A: int,
    B: int,
    C: int,
) -> Trajectory:
    """Compute the a proton transfer coordinate property for a Trajectory (au).

    Params:
        traj: the Trajectory object to compute the property for (modified in
            place)
        key: the name of the property
        A: the index of the first atom
        B: the index of the second atom
        C: the index of the transfered atom
    Return:
        traj: reference to the input Trajectory object. The property
            key is set to the proton transfer coordinate for the
            indices A, B, C
    """

    traj = compute_bond(traj, "dAC", A, C)
    traj = compute_bond(traj, "dBC", B, C)
    traj = compute_bond(traj, "dAB", A, B)
    for frame in traj.frames:
        dAC = frame.properties["dAC"]
        dBC = frame.properties["dBC"]
        dAB = frame.properties["dAB"]
        tau = (dBC - dAC) / dAB
        frame.properties[key] = np.array([tau])
    return traj
