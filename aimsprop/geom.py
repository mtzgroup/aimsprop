import math

import numpy as np

from .bundle import Bundle

# => Utility Math functions <= #


def _normalize(
    vec,
):
    """Return a normalized version of vec"""

    return vec / math.sqrt(sum(vec ** 2))


def _dot(
    vec1,
    vec2,
):
    """Dot product between vec1 and vec2"""

    return sum(vec1 * vec2)


def _cross(
    vec1,
    vec2,
):
    """Cross product between vec1 and vec2 in R^3"""

    vec3 = np.zeros((3,))
    vec3[0] = +(vec1[1] * vec2[2] - vec1[2] * vec2[1])
    vec3[1] = -(vec1[0] * vec2[2] - vec1[2] * vec2[0])
    vec3[2] = +(vec1[0] * vec2[1] - vec1[1] * vec2[0])
    return vec3


# => Geometric Properties <= #

""" Compute common geometric properties of Bundle objects, such as bond
    distances, bond angles, torsion angles, and out-of-plane angles. 
"""


def compute_bond(
    bundle: Bundle,
    key: str,
    A: int,
    B: int,
) -> Bundle:
    """Compute the a bond-length property for a Bundle.

    Params:
        bundle: the Bundle object to compute the property for (modified in
            place)
        key: the name of the property
        A: the index of the first atom
        B: the index of the second atom
    Return:
        bundle - reference to the input Bundle object. The property
            key is set to the float value of the bond length for the
            indices A and B.
    """

    for frame in bundle.frames:
        xyz = frame.xyz
        rAB = xyz[B, :] - xyz[A, :]
        frame.properties[key] = math.sqrt(sum(rAB ** 2))
    return bundle


def compute_angle(
    bundle: Bundle,
    key: str,
    A: int,
    B: int,
    C: int,
) -> Bundle:
    """Compute the a bond-angle property for a Bundle (in degrees).

    Params:
        bundle: the Bundle object to compute the property for (modified in
            place)
        key: the name of the property
        A: the index of the first atom
        B: the index of the second atom
        C: the index of the third atom
    Return:
        bundle: reference to the input Bundle object. The property
            key is set to the float value of the bond angle for the
            indices A, B and C
    """

    for frame in bundle.frames:
        xyz = frame.xyz
        rAB = xyz[B, :] - xyz[A, :]
        rCB = xyz[B, :] - xyz[C, :]
        frame.properties[key] = (
            180.0
            / math.pi
            * math.acos(sum(rAB * rCB) / math.sqrt(sum(rAB ** 2) * sum(rCB ** 2)))
        )
    return bundle


def compute_torsion(
    bundle: Bundle,
    key: str,
    A: int,
    B: int,
    C: int,
    D: int,
) -> Bundle:
    """Compute the a torsion-angle property for a Bundle (in degrees).

    Params:
        bundle: the Bundle object to compute the property for (modified in
            place)
        key: the name of the property
        A: the index of the first atom
        B: the index of the second atom
        C: the index of the third atom
        D: the index of the fourth atom
    Return:
        bundle: reference to the input Bundle object. The property
            key is set to the float value of the torsion angle for the
            indices A, B, C, and D
    """

    for frame in bundle.frames:
        xyz = frame.xyz
        rAB = xyz[B, :] - xyz[A, :]
        rBC = xyz[C, :] - xyz[B, :]
        rCD = xyz[D, :] - xyz[C, :]
        _normalize(rAB)
        eBC = _normalize(rBC)
        _normalize(rCD)
        n1 = _normalize(_cross(rAB, rBC))
        n2 = _normalize(_cross(rBC, rCD))
        m1 = _cross(n1, eBC)
        x = _dot(n1, n2)
        y = _dot(m1, n2)
        theta = 180.0 / math.pi * math.atan2(y, x)
        frame.properties[key] = theta

    return bundle


def compute_oop(
    bundle: Bundle,
    key: str,
    A: int,
    B: int,
    C: int,
    D: int,
) -> Bundle:
    """Compute the a out-of-plane-angle property for a Bundle (in degrees).

    Params:
        bundle: the Bundle object to compute the property for (modified in
            place)
        key: the name of the property
        A: the index of the first atom (OOP atom)
        B: the index of the second atom
        C: the index of the third atom
        D: the index of the fourth atom (Defines vector to A)
    Return:
        bundle: reference to the input Bundle object. The property
            key is set to the float value of the out-of-plane angle for the
            indices A, B, C, and D
    """

    bundle = compute_angle(bundle, "AngleBDC", B, D, C)

    for frame in bundle.frames:
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

    return bundle


def compute_transfer_coord(
    bundle: Bundle,
    key: str,
    A: int,
    B: int,
    C: int,
) -> Bundle:
    """Compute the a proton transfer coordinate property for a Bundle (au).

    Params:
        bundle: the Bundle object to compute the property for (modified in
            place)
        key: the name of the property
        A: the index of the first atom
        B: the index of the second atom
        C: the index of the transfered atom
    Return:
        bundle: reference to the input Bundle object. The property
            key is set to the proton transfer coordinate for the
            indices A, B, C
    """

    bundle = compute_bond(bundle, "dAC", A, C)
    bundle = compute_bond(bundle, "dBC", B, C)
    bundle = compute_bond(bundle, "dAB", A, B)
    for frame in bundle.frames:
        dAC = frame.properties["dAC"]
        dBC = frame.properties["dBC"]
        dAB = frame.properties["dAB"]
        tau = (dBC - dAC) / dAB
        frame.properties[key] = np.array([tau])
    return bundle
