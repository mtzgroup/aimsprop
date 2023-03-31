import numpy as np
from .bundle import Bundle
from collections import defaultdict


def compute_population(bundle: Bundle):
    """Compute the state populations for a bundle
    Params:
        bundle: the Bundle object to compute the property for
    Returns:
        populations: dictionary from I to time population of I. The indices of each vector in populations
        corresponds to the set of bundle.ts
    """

    ts = bundle.ts
    populations = defaultdict(lambda: np.zeros(len(ts),))

    # Group frames by I and t
    frames_by_I_and_t = defaultdict(list)
    for frame in bundle.frames:
        I, t = frame.I, frame.t
        frames_by_I_and_t[(I, t)].append(frame.w)

    # Calculate populations for each I at each t
    for (I, t), weights in frames_by_I_and_t.items():
        t_index = np.searchsorted(ts, t)
        populations[I][t_index] += sum(weights)

    return dict(populations)
