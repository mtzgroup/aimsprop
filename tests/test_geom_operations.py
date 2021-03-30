import pytest as pt

import aimsprop as ai


def test_compute_bond(trajectory):

    # 1. Compute the bond distance between atoms 0 and 1 for the trajectory.
    # Test if the first frame's bond length is correct.

    ai.compute_bond(trajectory, "R01", 0, 1)
    assert trajectory.frames[0].properties["R01"] == pt.approx(
        1.282298373, rel=1e-10
    )  # 0767554


def test_compute_angle(trajectory):

    # 2. Compute the angle between atoms 0, 1, 2 for the trajectory.
    #    Test if the first frame's angle is correct.

    ai.compute_angle(trajectory, "A012", 0, 1, 2)
    assert trajectory.frames[0].properties["A012"] == pt.approx(
        122.541809807, rel=1e-10
    )  # 4228


def test_compute_torsion(trajectory):

    # 3. Compute the torsion between atoms 0, 1, 2, 3 for the trajectory
    #    Test if the first frame's dihedral angle is correct.

    ai.compute_torsion(trajectory, "T0123", 0, 1, 2, 3)
    assert trajectory.frames[0].properties["T0123"] == pt.approx(
        160.465152603, rel=1e-10
    )  # 32674


def test_compute_oop(trajectory):

    # 4. Compute the out of plane torsion between atoms 0, 1, 2, 3 for the trajectory
    #    Test if the first frame's out of plane angle is correct.

    ai.compute_oop(trajectory, "OOP0123", 0, 1, 2, 3)
    assert trajectory.frames[0].properties["OOP0123"] == pt.approx(
        9.1198640147, rel=1e-11
    )  # 02102


def test_compute_transfer_coord(trajectory):

    # 5. Compute the proton transfer coordinate between atoms 0, 1, 2, 3 for the trajectory
    #   Test if the first frame's proton transfer coordinate matches the original

    ai.compute_transfer_coord(trajectory, "PT012", 0, 1, 2)
    coordinate_0 = float(trajectory.frames[0].properties["PT012"])
    assert coordinate_0 == pt.approx(-0.777663303, rel=1e-10)  # 0639367
