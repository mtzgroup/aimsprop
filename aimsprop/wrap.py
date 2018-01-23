def unwrap_property(
    traj,
    key,
    period,
    N=1,
    ):

    """ Unwrap a periodic property to be maximally continuous through time for
        each label. 

    Params:
        traj (Trajectory) - the Trajectory object to compute the property for (modified in
            place)
        key (str) - the name of the original property
        period (float) - the periodicity of the property
        N (int) - the number of images to try in +/- wings
    Result/Return:
        traj (Trajectory) - reference to the input Trajectory object. The
            property key is overwritten with the unwrapped property.
    """

    for label in traj.labels:
        frame_old = None
        for frame in traj.subset_by_label(label).frames:
            if frame_old is None: # Property has just started
                frame_old = frame
                continue
            prop_old = frame_old.properties[key]
            prop = frame.properties[key]
            vals = [prop + k * period for k in range(-N,+N+1)]
            frame.properties[key] = min(vals, key=lambda x : abs(x - prop_old))
            frame_old = frame
    return traj 
                
        

    
