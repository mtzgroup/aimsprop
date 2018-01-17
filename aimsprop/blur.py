import numpy as np

def blur_property(  
    traj,
    key,
    key2,
    R,
    alpha,
    ):

    """ Blur a property via Gaussian convolution.

    Pararms:
        traj (Trajectory) - the Trajectory object to compute the property for (modified in
            place)
        key (str) - the name of the original property
        key2 (str) - the name of the blurred property
        R2 (np.ndarray) - the grid of property values to blur to (usually R or
            theta or Q).
        alpha (float) - the Gaussian blurring exponent
    Result/Return:
        traj (Trajectory) - reference to the input Trajectory object. The
            property key2 is set to the np.ndarray value of the blurred
            property.

    """

    for frame in traj.frames:
        V = frame.properties[key]
        V2 = np.zeros_like(R)
        if np.array(V).ndim == 0:
            # Yak shave, I hates it!!
            V2 += np.sqrt(alpha / np.pi) * np.exp(-alpha * (R - V)**2)
        else: 
            for RAB in V:
                V2 += np.sqrt(alpha / np.pi) * np.exp(-alpha * (R - RAB)**2)
        frame.properties[key2] = V2
    return traj
