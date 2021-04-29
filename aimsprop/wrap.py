from .bundle import Bundle


def unwrap_property(
    bundle: Bundle,
    key: str,
    period: float,
    N: int = 1,
    match_all: bool = False,
):
    """Unwrap a periodic property to be maximally continuous through time for
        each label.

    Params:
        bundle: the Bundle object to compute the property for (modified in
            place)
        key: the name of the original property
        period: the periodicity of the property
        N: the number of images to try in +/- wings
        match_all: whether to wrap all labels to same sign
    Return:
        bundle: reference to the input Bundle object. The
            property key is overwritten with the unwrapped property.
    """

    if match_all:
        ref_frame = bundle.subset_by_label(bundle.labels[0]).frames[0]
    else:
        ref_frame = None

    for label in bundle.labels:
        frame_old = ref_frame
        for frame in bundle.subset_by_label(label).frames:
            if frame_old is None:  # Property has just started
                frame_old = frame
                continue
            prop_old = frame_old.properties[key]
            prop = frame.properties[key]
            vals = [prop + k * period for k in range(-N, +N + 1)]
            frame.properties[key] = min(vals, key=lambda x: abs(x - prop_old))
            frame_old = frame
    return bundle
