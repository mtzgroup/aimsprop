import numpy as np

from .bundle import Bundle


def compute_pes(
    bundle: Bundle,
    carrier_frequency: float,
    alpha: float,
    eKT: np.ndarray,
) -> Bundle:

    """Compute the simple photoelectron spectroscopy, with Guassian blurring

    User is responsible for calculating and assigning properties to the bundle frames:
            Dyson Orbitals
            Ionization Potential (IP)

    Params:
        bundle: the Bundle object to compute the property for (modified in
            place)
        carrier_frequency: experimental probe pulse carrier frequency (hbar*omega)
        alpha: the Guassian blurring exponent
        eKT: electron energies

    Return:
        bundle: reference to the input Bundle object. The property
            key "pes" is set to computed PES property.
    """

    for frame in bundle.frames:
        IPs = frame.properties["IP"]
        dyson_norms = frame.properties["dyson_norms"]
        pes = np.zeros_like(eKT)
        for ind, (state, IP) in enumerate(IPs):
            dyson_norm = dyson_norms[np.where(dyson_norms[:, 0] == state), 1][0]
            pes += (
                dyson_norm
                * np.sqrt(alpha / np.pi)
                * np.exp(-alpha * (carrier_frequency - IP - eKT) ** 2)
            )
        frame.properties["pes"] = pes

    return bundle
