import pytest as pt

import aimsprop as ai


def test_compute_population(bundle):

    # 1. Compute the populations for states.
    # Test if the computed population at time 25 is correct.
    pop = ai.compute_population(bundle)
    assert pop[1][24] == pt.approx(0.441434205, rel=1e-7)
