import numpy as np


class Frame(object):
    """Class Frame represents one "frame," a molecular geometry with associated time and label."""

    def __init__(
        self,
        label,  # TODO type?
        t: float,
        w: float,
        I: int,
        N: list,
        xyz: np.ndarray,
        properties={},
    ):
        """Verbatim constructor.

        Arguments:
            label: label identifying basis function [e.g., TBF
                index (int), TBF pair index (int, int), etc]
            t: simulation time in au.
            w: weight of frame (e.g., Mulliken population, trajectory
                weight, etc)
            I: electronic state label
            N: list of atomic numbers of molecule
            xyz: coordinates of nuclei in
                Angstrom, shape (natom,3)
            properties: dictionary mapping key (str) to user-defined
                numerical properties (float, complex, or np.ndarray)
        """

        self.label = label
        self.t = t
        self.w = w
        self.I = I
        self.N = N
        self.xyz = xyz
        self.properties = properties.copy()

    def copy(self):
        """Make a copy of self that is sufficiently deep to prevent
        modification of self by modification of the copy's attributes or
        properties keys / references (property values are not deep copied)."""
        return Frame(
            label=self.label,
            t=self.t,
            w=self.w,
            I=self.I,
            N=self.N,
            xyz=self.xyz,
            properties=self.properties,  # __init__ makes a copy of this
        )

    def __lt__(
        self,
        other,
    ):
        """Replacement for python2 __cmp__, to help with sorting : Comparator basis is (t, label, I).

        Old __cmp__: return cmp((self.t, self.label, self.I), (other.t, other.label, other.I))

        Details: https://stackoverflow.com/questions/8276983/why-cant-i-use-the-method-cmp-in-python-3-as-for-python-2
        """
        return (self.t, self.label, self.I) < (other.t, other.label, other.I)

    def __str__(self):
        """Return frame properties"""
        return f"""label: {self.label}
        time: {self.t}
        weight: {self.w}
        state: {self.I}
        atomic number: {self.N}
        xyz: {self.xyz}
        properties: {self.properties}
        """

    def _repr_html_(self):
        """Print nicely"""
        return "<p>" + str(self).replace("\n", "<br>") + "</p>"


class Trajectory(object):
    """Class Trajectory represents a list of Frames, with utility methods to
    extract and merge Trajectories.

    Note that many methods below (e.g., subset_by_*) return shallow copies
    or "views" of the current Trajectory object's frames. If a deeper copy
    is needed, one can easily call Trajectory.copy(), which relies on
    Frame.copy().
    """

    def __init__(
        self,
        frames: list,  # TODO list of frames
    ):
        """Verbatim constructor.

        Params:
            frames: list of Frames in this Trajectory.
        """

        self.frames = frames

    @property
    def labels(self):
        """ Return all unique labels in this trajectory, in sorted order. """
        return list(sorted(set([x.label for x in self.frames])))

    @property
    def ts(self):
        """ Return all unique times in this trajectory, in sorted order. """
        return np.array(sorted(set([x.t for x in self.frames])))

    @property
    def Is(self):
        """ return all unique Is in this trajectory, in sorted order. """
        return list(sorted(set([x.I for x in self.frames])))

    def copy(self):
        """ Return a new Trajectory with frames copied according to Frame.copy()."""
        return Trajectory([frame.copy() for frame in self.frames])

    def subset_by_label(
        self,
        label,
    ):
        """ Return a subset of this trajectory containing all frames with a given label (Frame-sorted) (view) """
        return Trajectory(list(sorted([x for x in self.frames if x.label == label])))

    def subset_by_sublabel(
        self,
        label,
        index,
    ):
        """ Return a subset of this trajectory containing all frames with a given partial label (Frame-sorted) (view) """
        return Trajectory(
            list(sorted([x for x in self.frames if x.label[index] == label]))
        )

    def subset_by_t(
        self,
        t,
        delta=1.0e-11,
    ):
        """Return a subset of this trajectory containing all frames with a given time (Frame-sorted) (view).

        Note that due to possible weirdness with ULP errors in float t
        values, we grab all times within delta absolute error of t
        """
        return Trajectory(
            list(sorted([x for x in self.frames if abs(x.t - t) < delta]))
        )

    def subset_by_I(
        self,
        I,
    ):
        """ Return a subset of this trajectory containing all frames with a given I (Frame-sorted) (view) """
        return Trajectory(list(sorted([x for x in self.frames if x.I == I])))

    def __add__(
        self,
        other,
    ):
        """ Concatenation operator to merge two trajectories (view-based) """
        return Trajectory(self.frames + other.frames)

    def __mul__(
        self,
        w,
    ):
        """Return a new Trajectory with all Frame objects multiplied by
        weight (new copy). This is useful to provide a weight on the initial
        condition due to e.g., oscillator strength and/or conformational population.
        """
        frames = []
        for frame in self.frames:
            frame2 = frame.copy()
            frame2.w *= w
            frames.append(frame2)
        return Trajectory(frames)

    __rmul__ = __mul__

    @staticmethod
    def merge(
        trajs,
        ws,
        labels=None,
    ):
        """Merge a list of trajectories together, with weighting factors (new copy).

        Params:
            trajs (list of Trajectory): trajectories to merge

            ws (list of float): weight factors of trajectories (e.g., from
                oscillator strength or conformational well).
            labels (list of label): if provided, labels are updated to (label,
                frame.label). This is used, e.g., to label the frame by IC as
                well as by any trajectory-specific labeling.
        Returns:
            Trajectory: a single merged Trajectory with updated weights (and labels)
        """

        frames = []
        for I, traj in enumerate(trajs):
            for frame in traj.frames:
                frame2 = frame.copy()
                frame2.w *= ws[I]
                if labels:
                    frame2.label = (labels[I], frame2.label)
                frames.append(frame2)
        return Trajectory(frames)

    def interpolate_nearest(
        self,
        ts,
        delta=1.0e-11,
    ):
        """Return a new trajectory with frame objects interpolated by nearest
            neighbor interpolation if t is inside the range of times of self's
            frames for each label (e.g., no extrapolation is performed). (new
            copy).

            Note that due to possible weirdness with ULP errors in float t
            values, we grab all times within delta absolute error of t

        Params:
            ts (list of float): times to interpolated new Trajectory to
        Returns:
            Trajectory: new trajectory with interpolated frames.
        """

        frames = []
        for label in self.labels:
            traj2 = self.subset_by_label(label)
            t2s = traj2.ts
            for t in ts:
                if t < min(t2s - delta) or t > max(t2s + delta):
                    continue  # Out of range (no extrapolation)
                t2 = min(t2s, key=lambda x: abs(x - t))  # Closest value in t2s
                frames2 = traj2.subset_by_t(t2).copy().frames
                for frame2 in frames2:
                    frame2.t = t
                frames += frames2
        return Trajectory(list(sorted(frames)))

    def interpolate_linear(
        self,
        ts,
    ):
        """Return a new trajectory with frame objects interpolated by linear
            interpolation if t is inside the range of times of self's
            frames for each label (e.g., no extrapolation is performed). (new
            copy).

            Using:
                f(t) = f(t0) + [ (f(t1) - f(t0)) / (t1 - t0) ] * (t - t0)
            Note:
                Not suitable for properties with large second derivative

        Arguments:
            ts (list of float): times to interpolated new Trajectory to
        Returns:
            Trajectory: new trajectory with interpolated frames.
        """

        frames = []
        for label in self.labels:
            traj2 = self.subset_by_label(label)
            t2s = traj2.ts
            for t in ts:
                if t < min(t2s) or t >= max(t2s):
                    continue  # Out of range (no extrapolation)
                # Get points in time to linearly interpolate between
                t3 = t2s[np.where(t2s <= t)[0][-1]]
                t4 = t2s[np.where(t2s > t)[0][0]]
                # Retrieve frames at points of linear interpolation for all labels
                framei = traj2.subset_by_t(t3).copy().frames
                framej = traj2.subset_by_t(t4).copy().frames

                # TODO: Instead of this eliminate duplicate frames at beginning of method?
                if len(framei) >= 1:
                    # print "Warning: duplicate frames present!"
                    framei = framei[0]
                if len(framej) >= 1:
                    # print "Warning: duplicate frames present!"
                    framej = framej[0]

                # linearly interpolate weight
                w = framei.w + ((framej.w - framei.w) / (t4 - t3)) * (t - t3)
                # linearly interpolate position
                xyz = framei.xyz + ((framej.xyz - framei.xyz) / (t4 - t3)) * (t - t3)

                # TODO: create separate option for interpolation of properties (in case different interpolation scheme preferred?
                # also if property is of incompatible type (such as string)

                # linearly interpolate all other properties
                properties = {}
                for key, vali in list(framei.properties.items()):
                    valj = framej.properties[key]
                    valn = vali + ((valj - vali) / (t4 - t3)) * (t - t3)
                    properties[key] = valn

                # create updated frame at new time instace
                nframe = Frame(
                    label=label,
                    t=t,
                    w=w,
                    I=framei.I,
                    N=framei.N,
                    xyz=xyz,
                    properties=properties,
                )
                frames.append(nframe)

        return Trajectory(list(sorted(frames)))

    def extract_property(
        self,
        key: str,
        normalize: bool = True,
        diff: bool = False,
    ):
        """Return a numpy array containing the time-history of a property,
            averaged over all Frames at each time (with frame weight applied).

        Params:
            key: the property key
            normalize: normalize the property by the sum of weights in
                each time?
            diff: difference property (zeroth-time value subtracted)?
        Returns:
            np.ndarray of shape (ntime, sizeof(prop)): the property
            expectation value. Time is always on the rows. If the property is
            scalar, this array will have ndim = 1. If the property is vector,
            this array with have ndim = 2. And so forth.
        """

        Vs = []
        for t in self.ts:
            traj = self.subset_by_t(t)
            V = np.zeros_like(traj.frames[0].properties[key])
            W = 0.0
            for frame in traj.frames:
                V += frame.w * frame.properties[key]
                W += frame.w
            if normalize:
                V /= W
            Vs.append(V)
        if diff:
            R = np.copy(Vs[0])
            for V in Vs:
                V -= R
        return np.array(Vs)

    def remove_duplicates(
        self,
    ):
        """Return a new trajectory with duplicate frame objects removed
            based on frame label and time index

        Returns:
            Trajectory: new trajectory with interpolated frames.
        """

        frames = []
        for ind1, frame1 in enumerate(self.frames):
            unique = True
            for ind2, frame2 in enumerate(self.frames[ind1 + 1 :]):
                if frame1.label == frame2.label and frame1.t == frame2.t:
                    unique = False
            if unique:
                frames.append(frame1)

        return Trajectory(frames)
