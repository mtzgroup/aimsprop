import copy

import numpy as np

from . import traj


def bootstrap(
    input_traj: traj.Trajectory,
    nsamples: int,
    ICs: list,
    label_ind: int = 0,
):
    """Resample Trajectory Subset According to Bootstrap Algorithm

    Params:
        input_traj: the Trajectory object to resample with replacement
        nsamples: number of resampling sets to conduct
        label_ind: sub label index from which to sample (IC = 0)
    Returns:
        resampled_trajs (list of Trajectories [nsamples]):  re-weighted list of trajectories

    Note:
        Removed parameters: sample_labels - list of trajectory labels from which to resample
    """

    # accumulate samples
    labels = input_traj.labels
    nICs = len(ICs)
    resampled_input_trajs = []
    for ind in range(nsamples):

        # copy of traj to modify
        traj1 = copy.copy(input_traj)

        # generate random numbers for resampling
        samples_inds = np.random.randint(low=0, high=nICs, size=nICs)

        # get new sample set
        samples = []
        new_labels = []
        for ind1, sind in enumerate(samples_inds):
            IC = ICs[sind]
            [samples.append(label) for label in labels if label[label_ind] == IC]
            [new_labels.append(ind1) for label in labels if label[label_ind] == IC]

        # separate new sample set from input_trajectory
        ws = 0.0
        input_trajs = [traj1.subset_by_label(sample) for sample in samples]
        t0 = traj1.ts[0]
        # re-weight trajectories
        for input_traj_t in input_trajs:
            input_trajs_t0 = input_traj_t.subset_by_t(t0)
            for frame in input_trajs_t0.frames:
                ws += frame.w
        ws = 1.0 / ws

        # merge input_trajectories
        traj1 = traj.Trajectory.merge(
            input_trajs, [ws] * len(traj1.frames), labels=new_labels
        )
        resampled_input_trajs.append(traj1)
        print((ind + 1, "/", nsamples, "complete"))

    return resampled_input_trajs


def extract_stats(
    trajs: list,
    key: str,
    diff=False,
):
    """Extract Standard Deviation and Average of Property from Set of Resampled Trajectories

    Params:
        trajs: re-weighted list of trajectories
        key: the property key
        diff: difference property (zeroth-time value subtracted)?
    Returns:
        avg (np.ndarray of shape (ntime, sizeof(prop))): the average property expectation value.
        std (np.ndarray of shape (ntime, sizeof(prop))): the standard deviation of the property expectation value.

    """

    props = []
    for ind, traj in enumerate(trajs):
        PROP = traj.extract_property(key)
        if diff == True:
            PROP -= np.outer(np.ones((PROP.shape[0],)), PROP[0, :])
        props.append(PROP)

    std = np.std(props, axis=0)
    avg = np.average(props, axis=0)

    return avg, std
