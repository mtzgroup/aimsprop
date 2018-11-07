def unwrap_property(
    traj,
    key,
    period,
    N=1,
    match_all=False,
    ):

    """ Unwrap a periodic property to be maximally continuous through time for
        each label. 

    Params:
        traj (Trajectory) - the Trajectory object to compute the property for (modified in
            place)
        key (str) - the name of the original property
        period (float) - the periodicity of the property
        N (int) - the number of images to try in +/- wings
        match_all (bool) - whether to wrap all labels to same sign
    Result/Return:
        traj (Trajectory) - reference to the input Trajectory object. The
            property key is overwritten with the unwrapped property.
    """

    if match_all:
        ref_frame = traj.subset_by_label(traj.labels[0]).frames[0]
    else:
        ref_frame = None

    for label in traj.labels:
        frame_old = ref_frame
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
                
        

    
