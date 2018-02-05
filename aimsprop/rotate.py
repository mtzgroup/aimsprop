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

def rotate_frames_angle(
    traj,
    x_angle=0.0,
    y_angle=0.0,
    z_angle=0.0,
    ):
    
    """ Rotate xyz coordinates for all frames in trajectory by rotation matrix R. 

    Performs xyz = xyz * R

    Params:
        traj (Trajectory) - the Trajectory object to rotate coordinates in
            place (modified in place).
        x_angle (angle) - angle (in degrees) of rotation for transformation matrix.
        y_angle (angle) - angle (in degrees) of rotation for transformation matrix.
        z_angle (angle) - angle (in degrees) of rotation for transformation matrix.
    Result/Return:
        traj (Trajectory) - reference to the input Trajectory object. The
            xyz coordinates are overwritten with the transformed xyz
            coordinates.
    """
    
    x_angle = np.radians(x_angle)
    y_angle = np.radians(y_angle)
    z_angle = np.radians(z_angle)
    
    # x rotation matrix (identity if x_angle == 0)
    Rx = np.array([
        [1.0, 0.0, 0.0],
        [0.0, np.cos(x_angle), -np.sin(x_angle)],
        [0.0, np.sin(x_angle), np.cos(x_angle)],
    ])
    
    # y rotation matrix (identity if y_angle == 0)
    Ry = np.array([
        [np.cos(y_angle), 0.0, np.sin(y_angle)],
        [0.0, 1.0, 0.0],
        [-np.sin(y_angle), 0.0, np.cos(y_angle)],
    ])
    
    # z rotation matrix (identity if z_angle == 0)
    Rz = np.array([
        [np.cos(z_angle), -np.sin(z_angle), 0.0],
        [np.sin(z_angle), np.cos(z_angle), 0.0],
        [0.0, 0.0, 1.0],
    ])

    for frame in traj.frames:
        frame.xyz = np.dot(frame.xyz, Rx)
        frame.xyz = np.dot(frame.xyz, Ry)
        frame.xyz = np.dot(frame.xyz, Rz)

    return traj
