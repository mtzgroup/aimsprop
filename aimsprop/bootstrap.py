import copy

import numpy as np

from . import bundle


def bootstrap(
    input_bundle: bundle.Bundle,
    nsamples: int,
    ICs: list,
    label_ind: int = 0,
):
    """Resample Bundle Subset According to Bootstrap Algorithm

    Params:
        input_bundle: the Bundle object to resample with replacement
        nsamples: number of resampling sets to conduct
        label_ind: sub label index from which to sample (IC = 0)
    Returns:
        resampled_bundles (list of Bundles [nsamples]):  re-weighted list of bundles

    Note:
        Removed parameters: sample_labels - list of bundle labels from which to resample
    """

    # accumulate samples
    labels = input_bundle.labels
    nICs = len(ICs)
    resampled_input_bundles = []
    for ind in range(nsamples):

        # copy of bundle to modify
        bundle1 = copy.copy(input_bundle)

        # generate random numbers for resampling
        samples_inds = np.random.randint(low=0, high=nICs, size=nICs)

        # get new sample set
        samples = []
        new_labels = []
        for ind1, sind in enumerate(samples_inds):
            IC = ICs[sind]
            [samples.append(label) for label in labels if label[label_ind] == IC]
            [new_labels.append(ind1) for label in labels if label[label_ind] == IC]

        # separate new sample set from input_bundle
        ws = 0.0
        input_bundles = [bundle1.subset_by_label(sample) for sample in samples]
        t0 = bundle1.ts[0]
        # re-weight bundles
        for input_bundle_t in input_bundles:
            input_bundles_t0 = input_bundle_t.subset_by_t(t0)
            for frame in input_bundles_t0.frames:
                ws += frame.w
        ws = 1.0 / ws

        # merge input_bundles
        bundle1 = bundle.Bundle.merge(
            input_bundles, [ws] * len(bundle1.frames), labels=new_labels
        )
        resampled_input_bundles.append(bundle1)
        print((ind + 1, "/", nsamples, "complete"))

    return resampled_input_bundles


def extract_stats(
    bundles: list,
    key: str,
    diff=False,
):
    """Extract Standard Deviation and Average of Property from Set of Resampled Bundles

    Params:
        bundles: re-weighted list of bundles
        key: the property key
        diff: difference property (zeroth-time value subtracted)?
    Returns:
        avg (np.ndarray of shape (ntime, sizeof(prop))): the average property expectation value.
        std (np.ndarray of shape (ntime, sizeof(prop))): the standard deviation of the property expectation value.

    """

    props = []
    for ind, bundle in enumerate(bundles):
        PROP = bundle.extract_property(key)
        if diff == True:
            PROP -= np.outer(np.ones((PROP.shape[0],)), PROP[0, :])
        props.append(PROP)

    std = np.std(props, axis=0)
    avg = np.average(props, axis=0)

    return avg, std
