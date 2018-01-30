import lebedev
import legendre
import rotation
import numpy as np

def compute_diffraction_moments(
    traj,
    key,
    q,
    factors,
    nlebedev,
    nlebedev2,
    nomega2,
    nlegendre=2,
    print_level=False,
    ):  

    """ Compute the IAM X-Ray Diffraction or UED moments property. See
        aimsprop/notes/xray for details on these moments.

    Notes:
        * All frames for each initial condition (IC) in traj should be aligned so
        that the transition dipole moment from S0 -> Sex at t=0 is on z. This
        is required for proper computation of I2 (I0 is invariant to this). 
        * All frames should be weighted by geometric considerations at the IC
        (e.g., conformational wells, Wigner weights, etc), by the cross
        section for the optical transition at the IC (e.g., oscillator
        strength and excitation energy window), and by the frame weight due
        to non-adiabatic dynamics.

    Params:
        traj (Trajectory) - the Trajectory object to compute the property for (modified in
            place)
        key (str) - the name of the property. 
        q (np.ndarray) - the 1d array of |q| values to collocate the
            diffraction signal to.
        factors (list of AtomicFormFactor) - the list of AtomicFormFactor
            objects. The choice of xray/ued is made by the "mode" field of each
            AtomicFormFactor.
        nlebedev (int) - a Lebedev number for the number of grid points to use
            for angular integration of the Legendre coefficients.
        nlebedev2 (int) - a Lebedev number for the number of grid points to use
            for angular sampling in the Omega angle (surface angle).
        nomega2 (int) - the number of grid points to use for angular sampling
            in the omega angle (surface orientation).
        nlegendre (even int) - the maximum Legendre polynomial order to
            evaluate (usually 0 and 2 are the only Legendre polynomial
            coefficients with any signal).
        print_level (bool) - print progress if true (useful to track long
            property computations)
    Result/Return:
        traj - reference to the input Trajectory object. The properties
            key-l are set where l is [0, 2, ..., nlegendre].
    """
    if nlegendre % 2: raise ValueError('Can only ask for even Legendre functions')

    # Lebedev grid for S(\vec q)
    leb = lebedev.Lebedev.build(nlebedev)
    # Direct product grid for \vec q
    qx = np.outer(q, leb.x)
    qy = np.outer(q, leb.y)
    qz = np.outer(q, leb.z)

    # Rotation quadrature
    Rs, ws = rotation.rotation_quadrature(nomega=nomega2, nlebedev=nlebedev2)

    for find, frame in enumerate(traj.frames):
        if print_level:
            print 'Frame %5d of %5d' % (find, len(traj.frames))
        # Compute N(\vec q) = \sum_{A} f_A (\vec q) * \exp(-1.j * \vec q * \vec r)
        N = np.zeros((len(q), len(leb.x)), dtype=complex)
        for A, factor in enumerate(factors): 
            x = frame.xyz[A,0]
            y = frame.xyz[A,1]
            z = frame.xyz[A,2]
            N += factor.evaluate_N(qx=qx,qy=qy,qz=qz,x=x,y=y,z=z)
        # Compute I(\vec q) = N(\vec q)**2 for this frame
        I = (np.abs(N)**2).real

        # Do angle integration
        Ils = { l: np.zeros((len(q),)) for l in range(0, nlegendre+1, 2) }
        qxyz = leb.xyz
        for R2, w2 in zip(Rs, ws):
            # Account for cos(z)^2 weight
            cos2 = np.sum(R2[:,2])**2
            # Rotate the lebedev grid
            qxyz2 = np.dot(qxyz, R2)
            qx2 = qxyz2[:,0]
            qy2 = qxyz2[:,1]
            qz2 = qxyz2[:,2]
            # Weight by zonal harmonics (evaluated in rotated grid)
            Y = legendre.zonal2(qz2, nlegendre)
            for l in range(0, nlegendre+1, 2):
                Ils[l] += w2 * cos2 * np.einsum('qw,w->q', I, leb.w * Y[l/2, :]) 

        # Assign the properties to frame
        for l in range(0, nlegendre+1, 2):
            frame.properties['%s-%d' % (key, l)] = Ils[l]

    return traj

def compute_diffraction_moment0(
    traj,
    key,
    q,
    factors,
    nlebedev,
    print_level=False,
    ):  

    """ Compute the IAM X-Ray Diffraction or UED moment property for only l=0
        (faster due to lack of rotation quadratures). See aimsprop/notes/xray for
        details on this moment

    Notes:
        * This moment is invariant to the orientation of the frames, so alignment to the IC
        * All frames should be weighted by geometric considerations at the IC
        (e.g., conformational wells, Wigner weights, etc), by the cross
        section for the optical transition at the IC (e.g., oscillator
        strength and excitation energy window), and by the frame weight due
        to non-adiabatic dynamics.

    Params:
        traj (Trajectory) - the Trajectory object to compute the property for (modified in
            place)
        key (str) - the name of the property. 
        q (np.ndarray) - the 1d array of |q| values to collocate the
            diffraction signal to.
        factors (list of AtomicFormFactor) - the list of AtomicFormFactor
            objects. The choice of xray/ued is made by the "mode" field of each
            AtomicFormFactor.
        nlebedev (int) - a Lebedev number for the number of grid points to use
            for angular integration of the Legendre coefficients.
        print_level (bool) - print progress if true (useful to track long
            property computations)
    Result/Return:
        traj - reference to the input Trajectory object. The property
            key-0 is set
    """

    # Lebedev grid for S(\vec q)
    leb = lebedev.Lebedev.build(nlebedev)
    # Direct product grid for \vec q
    qx = np.outer(q, leb.x)
    qy = np.outer(q, leb.y)
    qz = np.outer(q, leb.z)

    for find, frame in enumerate(traj.frames):
        if print_level:
            print 'Frame %5d of %5d' % (find, len(traj.frames))
        # Compute N(\vec q) = \sum_{A} f_A (\vec q) * \exp(-1.j * \vec q * \vec r)
        N = np.zeros((len(q), len(leb.x)), dtype=complex)
        for A, factor in enumerate(factors): 
            x = frame.xyz[A,0]
            y = frame.xyz[A,1]
            z = frame.xyz[A,2]
            N += factor.evaluate_N(qx=qx,qy=qy,qz=qz,x=x,y=y,z=z)
        # Compute I(\vec q) = N(\vec q)**2 for this frame
        I = (np.abs(N)**2).real
        # Compute Y00 (for common normalization)
        Y00 = legendre.zonal2(leb.z, 0)[0]
        # Integrate over S(2)
        I0 = np.einsum('qw,w->q', I, leb.w * Y00)
        # Assign the property to frame
        frame.properties['%s-%d' % (key, 0)] = I0

    return traj

