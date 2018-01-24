import numpy as np

def rotate_frames(
    traj,
    R,
    ):
    
    """ Rotate xyz coordinates for all frames in trajectory by rotation matrix R. 

    Performs xyz = xyz * R

    Params:
        traj (Trajectory) - the Trajectory object to rotate coordinates in
            place (modified in place).
        R (np.ndarray, shape (3,3)) - rotation or transformation matrix.
    Result/Return:
        traj (Trajectory) - reference to the input Trajectory object. The
            xyz coordinates are overwritten with the transformed xyz
            coordinates.
    """

    for frame in traj.frames:
        frame.xyz = np.dot(frame.xyz, R)
    return traj
