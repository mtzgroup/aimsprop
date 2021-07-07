import pytest as pt

import aimsprop as ai


def test_rotate_frames(bundle):

    # 1. Rotate the coordinates according to conversion matrix R.
    #    Test if the first rotated coordinate is correct.

    R = [[0, 1, 0], [1, 0, 1], [0, 1, 1]]
    ai.rotate_frames(bundle, R)
    assert bundle.frames[0].xyz[0][0] == pt.approx(0.035608866, rel=1e-9)


def test_rotate_frames_to_z(bundle):

    # 1. Rotate the coordinates according to the z vector for TDP.
    #    Test if the first rotated coordinate is correct.

    z = [0.2, 0.1, 0.1]
    ai.rotate_frames_to_z(bundle, z)
    assert bundle.frames[0].xyz[0][0] == pt.approx(-0.2099543620557265, rel=1e-9)
