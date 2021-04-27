import numpy as np

from .bundle import Bundle


def compute_population(
    bundle: Bundle,
):
    """Compute the state populations for a bundle

    Params:
        bundle: the Bundle object to compute the property for
    Returns:
        populations: dictionary from I to time population of I. The indices of each vector in populations
        corresponds to the set of bundle.ts
    """

    ts = bundle.ts
    populations = {}
    for I in bundle.Is:
        bundle2 = bundle.subset_by_I(I)
        populations[I] = np.zeros((len(ts),))
        for tind, t in enumerate(ts):
            bundle3 = bundle2.subset_by_t(t)
            populations[I][tind] += sum([frame.w for frame in bundle3.frames])
    return populations
