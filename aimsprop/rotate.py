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

def rotate_to_z(
    traj,
    vec,
    ):

    """ Rotate xyz coordinates for all frames in trajectory so that vec is
        rotated onto z. Useful to align molecules to z according to transition
        dipole moment vector.

    Params:
        traj (Trajectory) - the Trajectory object to rotate coordinates in
            place (modified in place).
        vec (np.ndarray, shape (3,)) - vector in current coordinates to rotate
            onto +z
    Result/Return:
        traj (Trajectory) - reference to the input Trajectory object. The
            xyz coordinates are overwritten with the transformed xyz
            coordinates.
    """
    
    n = np.array(vec)
    n /= np.sqrt(np.sum(n**2))
    v = 0.5 * (n + np.array((0.0, 0.0, 1.0)))
    v /= np.sqrt(np.sum(v**2))
    R = np.eye(3)
    R -= 2.0 * np.outer(v, v)
    R = np.dot(R, -np.eye(3)) # preserve chirality
    return rotate_frames(traj, R)
