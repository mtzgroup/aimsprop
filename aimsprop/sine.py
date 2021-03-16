import numpy as np


def compute_sine_transform(
    I: np.ndarray,
    s: np.ndarray,
    r: np.ndarray,
    a: float = 0.0,
) -> np.ndarray:
    """Compute the sine transform I(t,r) = \int_0^smax ds I(t,s) * sin(r*s) * exp(-a * s^2)

    Uses the trapezoid rule over provided s (s need not be uniformly spaced).
    The factor of exp(-a * s^2) is a common experimental cutoff factor applied
    to prevent ringing from noise at large s. This tends to smooth the r-domain
    signal, but can limit the resolution.

    Params:
        I: signal to sine transform, e.g.,
            sM(t,s). Time should be on dim 0, s should be on dim 1,
            shape (nt, ns).
        s: s values (can be irregular), shape (ns,)
        r: r values (can be irregular), shape (nr,)
        a: Gaussian exponent in large-s cutoff Gaussian.
    Returns:
        I2: sine transformed signal, shape (nt, nr)
    """

    # Irregular trapezoid weights
    ds = np.diff(s)
    w = np.zeros_like(s)
    w[:-1] += 0.5 * ds
    w[1:] += 0.5 * ds

    V = I * np.outer(np.ones((I.shape[0],)), w * np.exp(-a * s ** 2))
    rr, ss = np.meshgrid(r, s, indexing="ij")
    K = np.sin(rr * ss)
    I2 = np.dot(V, K.T)
    return I2
