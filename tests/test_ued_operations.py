import numpy as np
import pytest as pt

import aimsprop as ai


def test_ued_simple(bundle):

    # 1. Compute the UED spectrum.
    #    Test if the computed UED spectrum at time 0 is correct.
    R = np.linspace(0.0, 6.0, 200)
    bundle = ai.compute_ued_simple(bundle, "UED", R=R, alpha=8.0)
    assert bundle.frames[0].properties["UED"][0] == pt.approx(8.95739037e-07, rel=1e-7)
