import numpy as np
import copy

class Frame(object):

    """ Class Frame represents one "frame," a molecular geometry with associated time and label. 
    """

    def __init__(
        self,
        label,
        t,
        w,
        I,
        N,
        xyz,
        properties={},
        ):

        """ Verbatim constructor.

        Params/Attributes:
            label (hashable) - label identifying basis function [e.g., TBF
                index (int), TBF pair index (int, int), etc]
            t (float) - simulation time in au.
            w (float) - weight of frame (e.g., Mulliken population, trajectory
                weight, etc)
            I (int) - electronic state label
            N (list of int of len natom) - list of atomic numbers of molecule
            xyz [np.ndarray of shape (natom, 3)] - coordinates of nuclei in
                Angstrom.
            properties - dictionary mapping key (str) to user-defined
                numerical properties (float, complex, or np.ndarray)
        """

        self.label = label 
        self.t = t
        self.w = w 
        self.I = I
        self.N = N
        self.xyz = xyz
        self.properties = properties.copy()

    def __cmp__(
        self,
        other,
        ):
        """ Comparator, to help with sorting. """
        return cmp((self.t, self.label, self.I), (other.t, other.label, other.I))

class Trajectory(object):

    """ Class Trajectory represents a list of Frames, with utility methods to
        extract and merge Trajectories. 
    """

    def __init__(
        self,
        frames,
        ): 

        """ Verbatim constructor.

        Params/Attributes:
            frames (list of Frame) - list of Frames in this Trajectory.
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
    
    def subset_by_label(
        self,
        label,
        ):
        """ Return a subset of this trajectory containing all frames with a given label (Frame-sorted) """
        return Trajectory(list(sorted([x for x in self.frames if x.label == label])))

    def subset_by_t(
        self,
        t,
        eps=1.0E-14,
        ):

        """ Return a subset of this trajectory containing all frames with a given time (Frame-sorted). 

            Note that due to possible weirdness with ULP errors in float t
            values, we grab all times within eps relative error of t. 
        """
        return Trajectory(list(sorted([x for x in self.frames if x.t == t or abs(x.t - t) < abs(eps * t)])))

    def subset_by_I(
        self,
        I,
        ):

        """ Return a subset of this trajectory containing all frames with a given I (Frame-sorted) """
        return Trajectory(list(sorted([x for x in self.frames if x.I == I])))
    
    def __add__(
        self,
        other,
        ):
        """ Concatenation operator to merge two trajectories """
        return Trajectory(self.frames + other.frames)

    def __mul__(
        self,
        weight,
        ):

        """ Return a new Trajectory with all Frame objects
        multiplied by weight. This is useful to provide a weight on
        the initial condition due to e.g., oscillator strength
        and/or conformational population. """
        frames = []
        for frame in self.frames:
            frame2 = copy.copy(frame)
            frame2.w *= weight
            frames.append(frame2)
        return 
    
    __rmul__ = __mul__

    def interpolate_nearest(
        self,
        ts,
        ):

        """ Return a new trajectory with frame objects interpolated by nearest
            neighbor interpolation if t is inside the range of times of self's
            frames for each label (e.g., no extrapolation is performed).

        Params:
            ts (list of float) - times to interpolated new Trajectory to
        Returns:
            traj (Trajectory) - new trajectory with interpolated frames.
        """

        frames = []
        for label in self.labels:
            traj2 = self.subset_by_label(label)
            t2s = traj2.ts
            for t in ts:
                if t < min(t2s) or t > max(t2s): continue # Out of range (TODO: eps bounds)
                t2 = min(t2s, key=lambda x:abs(x - t)) # Closest value in t2s
                frames += traj2.subset_by_t(t2).frames
        return Trajectory(list(sorted(frames)))

    # TODO: linear interpolation
                
                 
        
    
         
