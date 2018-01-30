import numpy as np

def compute_pes(
    traj,
    carrier_frequency, 
    alpha,
    eKT,
    ):

    
    """ Compute the simple photoelectron spectroscopy, with Guassian blurring

    User is responsible for calculating and assigning properties to the trajectory frames:
            Dyson Orbitals
            Ionization Potential (IP)

    Params:
        traj (Trajectory) - the Trajectory object to compute the property for (modified in
            place)
        key (str) - the name of the property
        carrier_frequency - (float64), experimental probe pulse carrier frequency (hbar*omega)
        alpha (float) - the Guassian blurring exponent
        eKT (np.ndarray) - electron energies

    Result/Return:
        traj - reference to the input Trajectory object. The property
            key is set to computed PES property.
    """

    for frame in traj.frames: 
        IPs = frame.properties['IP']
        dyson_norms = frame.properties['dyson_norms']
        pes = np.zeros_like(eKT)
        for ind, (state, IP) in enumerate(IPs):
            dyson_norm = dyson_norms[np.where(dyson_norms[:,0] == state), 1][0]
            pes += dyson_norm * np.sqrt(alpha / np.pi) * np.exp(-alpha * (carrier_frequency - IP - eKT)**2)
        frame.properties['pes'] = pes

    return traj
