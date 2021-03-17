import numpy as np

from .traj import Trajectory


def blur_property(
    traj: Trajectory,
    key: str,
    key2: str,
    R: np.ndarray,
    alpha: float,
) -> Trajectory:

    """Blur a property via Gaussian convolution (blurring applied in the "spatial" coordinate).

    Pararms:
        traj: the Trajectory object to compute the property for (modified in
            place)
        key: the name of the original property
        key2: the name of the blurred property
        R2: the grid of property values to blur to (usually R or
            theta or Q).
        alpha: the Gaussian blurring exponent
    Result/Return:
        traj (Trajectory): reference to the input Trajectory object. The
            property key2 is set to the np.ndarray value of the blurred
            property.

    """

    for frame in traj.frames:
        V = frame.properties[key]
        V2 = np.zeros_like(R)
        if np.array(V).ndim == 0:
            # Yak shave, I hates it!!
            V2 += np.sqrt(alpha / np.pi) * np.exp(-alpha * (R - V) ** 2)
        else:
            for RAB in V:
                V2 += np.sqrt(alpha / np.pi) * np.exp(-alpha * (R - RAB) ** 2)
        frame.properties[key2] = V2
    return traj


def compute_time_blur(
    I: np.ndarray,
    t1: np.ndarray,
    t2: np.ndarray,
    fwhm: float,
) -> np.ndarray:

    """Compute Gaussian blurring in time for an arbitrary property.

    Uses the trapezoid rule over provided t1 (t1/t2 need not be uniformly spaced).

    Params:
        I: signal to Gaussian blur in the 0-th dim, shape (nt1, ...).
        t1: t1 values (can be irregular), shape (nt1,).
        t2: t2 values (can be irregular), shape (nt2,).
        fwhm: the full width at half maximum (FWHM) of the blurring
            kernel.
    Returns:
        I2 - Gaussian blurred signal, shape (nt2, ...).
    """

    # Irregular trapezoid weights in t1
    dt = np.diff(t1)
    w = np.zeros_like(t1)
    w[:-1] += 0.5 * dt
    w[1:] += 0.5 * dt

    # Gaussian exponent corresponding to fwhm
    a = 4.0 * np.log(2.0) / (fwhm ** 2)

    # Blurring kernel
    tt1, tt2 = np.meshgrid(t1, t2, indexing="ij")
    K = (a / np.pi) ** (0.5) * np.exp(-a * (tt1 - tt2) ** 2)

    if I.ndim == 2:
        V = I * np.outer(w, np.ones((I.shape[1],)))
        I2 = np.dot(K, V)
    else:
        raise ValueError("ndim case %d not coded" % I.ndim)

    return I2
