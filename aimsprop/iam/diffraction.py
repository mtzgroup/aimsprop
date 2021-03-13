import numpy as np

from . import formfactor, rotation


def compute_diffraction(
    traj,
    key,
    s,
    eta,
    L,
    nlebedev=74,
    nomega=12,
    mode="xray",
    form="raw",
    anisotropy="cos2",
    print_level=False,
):

    """Compute the I(s, eta) elastic scattering signal for a Trajectory.
         See aimsprop/notes/ued for details on this property.

    Notes:
        * All frames for each initial condition (IC) in traj should be aligned so
        that the transition dipole moment from S0 -> Sex at t=0 is on z. This
        is required for proper computation of anisotropy.
        * All frames should be weighted by geometric considerations at the IC
        (e.g., conformational wells, Wigner weights, etc), by the cross
        section for the optical transition at the IC (e.g., oscillator
        strength and excitation energy window), and by the frame weight due
        to non-adiabatic dynamics.

    Params:
        traj (Trajectory) - the Trajectory object to compute the property for (modified in
            place)
        key (str) - the name of the property.
        s (np.ndarray) - list of scattering vector norms in Angstrom^-1. The
            relationship between s and theta (scattering angle) is given as,
                s = 4 pi / L * sin(theta / 2).
        eta (np.ndarray) - list of azimuthal scattering angles in radians.
        L (float) - effective wavelength of scattering particle (x-ray
            wavelength or UED deBroglie wavelength) in Angstrom. Used to
            convert through scattering angle theta.
        nlebedev (int) - Lebedev number to use for solid angle orientation
            quadrature.
        nomega (int) - number of uniform quadrature points to use for plane
            orientation quadrature.
        mode (str) - 'xray' or 'ued' for selection of form factors
        form (str) - 'raw' or 'mod' for modified/raw diffraction intensities
            I(s) or M(s).
        anisotropy (str) - 'none' or 'cos2' for isotropic of cos^2 (z)
            anisotropty.
        print_level (bool) - print progress if true (useful to track long
            property computations)
    Result/Return:
        traj - reference to the input Trajectory object. The properties
    """

    # Validity checks
    if mode not in ["xray", "ued"]:
        raise ValueError("Unknown mode: %s" % mode)
    if form not in ["raw", "mod"]:
        raise ValueError("Unknown form: %s" % form)
    if anisotropy not in ["none", "cos2"]:
        raise ValueError("Unknown anisotropy: %s" % anisotropy)

    # Compute scattering angles via Bragg equation
    theta = 2.0 * np.arcsin(s * L / (4.0 * np.pi))
    tt, ee = np.meshgrid(theta, eta, indexing="ij")
    ss, ee = np.meshgrid(s, eta, indexing="ij")
    # Compute scattering vectors
    sx = ss * np.cos(tt / 2.0) * np.sin(ee)
    sy = ss * np.sin(tt / 2.0)
    sz = ss * np.cos(tt / 2.0) * np.cos(ee)

    # Get a rotation quadrature for the orientations of the frames
    if nlebedev == 1 and nomega == 1:
        # Fixed orientation
        Rs = [np.eye(3)]
        ws = [1.0]
    else:
        # Rotation quadrature
        Rs, ws = rotation.rotation_quadrature(nlebedev=nlebedev, nomega=nomega)

    # Get atomic form factors for appropriate x-ray/ued mode
    factors = formfactor.AtomicFormFactor.build_factors(traj.frames[0], mode=mode)

    # Compute atomic scattering Iat
    D = np.zeros_like(sx)
    for A, factor in enumerate(factors):
        F = factor.evaluate_N(qx=sx, qy=sy, qz=sz, x=0.0, y=0.0, z=0.0)
        D += (np.abs(F) ** 2).real

    # Compute IAM scattering, integrating over all orientation angles
    for find, frame in enumerate(traj.frames):
        if print_level:
            print(("Frame %5d of %5d" % (find, len(traj.frames))))
        I = np.zeros_like(sx)
        for R, w in zip(Rs, ws):
            # cos(z)^2 pump anisotropy
            cos2 = R[2, 2] ** 2 if anisotropy == "cos2" else 1.0
            # Rotated molecule
            xyz = np.dot(frame.xyz, R)
            # Compute diffraction
            N = np.zeros_like(I, dtype=complex)
            for A, factor in enumerate(factors):
                x = xyz[A, 0]
                y = xyz[A, 1]
                z = xyz[A, 2]
                N += factor.evaluate_N(qx=sx, qy=sy, qz=sz, x=x, y=y, z=z)
            F = (np.abs(N) ** 2).real
            if form == "mod":
                F = (F - D) / D
            I += w * cos2 * F
        frame.properties[key] = I

    return traj


def compute_diffraction_fast(
    traj,
    key,
    s,
    eta,
    L,
    nlebedev=74,
    nomega=12,
    mode="xray",
    form="raw",
    anisotropy="cos2",
    print_level=False,
):

    """Compute the I(s, eta) elastic scattering signal for a Trajectory.
         See aimsprop/notes/ued for details on this property.

    Notes:
        * All frames for each initial condition (IC) in traj should be aligned so
        that the transition dipole moment from S0 -> Sex at t=0 is on z. This
        is required for proper computation of anisotropy.
        * All frames should be weighted by geometric considerations at the IC
        (e.g., conformational wells, Wigner weights, etc), by the cross
        section for the optical transition at the IC (e.g., oscillator
        strength and excitation energy window), and by the frame weight due
        to non-adiabatic dynamics.

    Params:
        traj (Trajectory) - the Trajectory object to compute the property for (modified in
            place)
        key (str) - the name of the property.
        s (np.ndarray) - list of scattering vector norms in Angstrom^-1. The
            relationship between s and theta (scattering angle) is given as,
                s = 4 pi / L * sin(theta / 2).
        eta (np.ndarray) - list of azimuthal scattering angles in radians.
        L (float) - effective wavelength of scattering particle (x-ray
            wavelength or UED deBroglie wavelength) in Angstrom. Used to
            convert through scattering angle theta.
        nlebedev (int) - Lebedev number to use for solid angle orientation
            quadrature.
        nomega (int) - number of uniform quadrature points to use for plane
            orientation quadrature.
        mode (str) - 'xray' or 'ued' for selection of form factors
        form (str) - 'raw' or 'mod' for modified/raw diffraction intensities
            I(s) or M(s).
        anisotropy (str) - 'none' or 'cos2' for isotropic of cos^2 (z)
            anisotropty.
        print_level (bool) - print progress if true (useful to track long
            property computations)
    Result/Return:
        traj - reference to the input Trajectory object. The properties
    """

    # Validity checks
    if mode not in ["xray", "ued"]:
        raise ValueError("Unknown mode: %s" % mode)
    if form not in ["raw", "mod"]:
        raise ValueError("Unknown form: %s" % form)
    if anisotropy not in ["none", "cos2"]:
        raise ValueError("Unknown anisotropy: %s" % anisotropy)

    # Get a rotation quadrature for the orientations of the frames
    if nlebedev == 1 and nomega == 1:
        # Fixed orientation
        Rs = [np.eye(3)]
        ws = [1.0]
    else:
        # Rotation quadrature
        Rs, ws = rotation.rotation_quadrature(nlebedev=nlebedev, nomega=nomega)

    # Get atomic form factors for appropriate x-ray/ued mode
    factors = formfactor.AtomicFormFactor.build_factors(traj.frames[0], mode=mode)

    import lightspeed as ls

    from . import ext

    s2s = ls.Tensor.array(s)
    eta2s = ls.Tensor.array(eta)

    R2s = ls.Tensor.zeros((len(Rs), 3, 3))
    for (
        Rind,
        R,
    ) in enumerate(Rs):
        R2s[Rind, :, :] = R
    w2s = ls.Tensor.array(ws)

    fA = ls.Tensor.zeros((len(factors), s.size))
    for A, factor in enumerate(factors):
        fA[A, :] = factor.evaluate(qx=0.0, qy=0.0, qz=s)

    # Compute IAM scattering, integrating over all orientation angles
    for find, frame in enumerate(traj.frames):
        if print_level:
            print(("Frame %5d of %5d" % (find, len(traj.frames))))
        xyz = ls.Tensor.array(frame.xyz)
        I = ext.compute_diffraction(
            L,
            s2s,
            eta2s,
            xyz,
            fA,
            R2s,
            w2s,
            True if anisotropy == "cos2" else False,
            True if form == "mod" else False,
        )
        frame.properties[key] = np.array(I)

    return traj


def compute_diffraction_moments_fast(
    traj,
    key,
    s,
    L,
    nlebedev=74,
    nomega=12,
    mode="xray",
    form="raw",
    anisotropy="cos2",
    print_level=False,
):

    """Compute the I(s, eta) elastic scattering moments for a Trajectory.
         See aimsprop/notes/ued for details on this property.

    Notes:
        * All frames for each initial condition (IC) in traj should be aligned so
        that the transition dipole moment from S0 -> Sex at t=0 is on z. This
        is required for proper computation of anisotropy.
        * All frames should be weighted by geometric considerations at the IC
        (e.g., conformational wells, Wigner weights, etc), by the cross
        section for the optical transition at the IC (e.g., oscillator
        strength and excitation energy window), and by the frame weight due
        to non-adiabatic dynamics.

    Params:
        traj (Trajectory) - the Trajectory object to compute the property for (modified in
            place)
        key (str) - the name of the property.
        s (np.ndarray) - list of scattering vector norms in Angstrom^-1. The
            relationship between s and theta (scattering angle) is given as,
                s = 4 pi / L * sin(theta / 2).
        L (float) - effective wavelength of scattering particle (x-ray
            wavelength or UED deBroglie wavelength) in Angstrom. Used to
            convert through scattering angle theta.
        nlebedev (int) - Lebedev number to use for solid angle orientation
            quadrature.
        nomega (int) - number of uniform quadrature points to use for plane
            orientation quadrature.
        mode (str) - 'xray' or 'ued' for selection of form factors
        form (str) - 'raw' or 'mod' for modified/raw diffraction intensities
            I(s) or M(s).
        anisotropy (str) - 'none' or 'cos2' for isotropic of cos^2 (z)
            anisotropty.
        print_level (bool) - print progress if true (useful to track long
            property computations)
    Result/Return:
        traj - reference to the input Trajectory object. The properties "key-0"
            and "key-2" are added to each frame of the Trajectory.
    """

    # Validity checks
    if mode not in ["xray", "ued"]:
        raise ValueError("Unknown mode: %s" % mode)
    if form not in ["raw", "mod"]:
        raise ValueError("Unknown form: %s" % form)
    if anisotropy not in ["none", "cos2"]:
        raise ValueError("Unknown anisotropy: %s" % anisotropy)

    # Special angles to collocate
    eta = np.array([0.0, np.pi / 2.0])

    # Get a rotation quadrature for the orientations of the frames
    if nlebedev == 1 and nomega == 1:
        # Fixed orientation
        Rs = [np.eye(3)]
        ws = [1.0]
    else:
        # Rotation quadrature
        Rs, ws = rotation.rotation_quadrature(nlebedev=nlebedev, nomega=nomega)

    # Get atomic form factors for appropriate x-ray/ued mode
    factors = formfactor.AtomicFormFactor.build_factors(traj.frames[0], mode=mode)

    import lightspeed as ls

    from . import ext

    s2s = ls.Tensor.array(s)
    eta2s = ls.Tensor.array(eta)

    R2s = ls.Tensor.zeros((len(Rs), 3, 3))
    for (
        Rind,
        R,
    ) in enumerate(Rs):
        R2s[Rind, :, :] = R
    w2s = ls.Tensor.array(ws)

    fA = ls.Tensor.zeros((len(factors), s.size))
    for A, factor in enumerate(factors):
        fA[A, :] = factor.evaluate(qx=0.0, qy=0.0, qz=s)

    # Compute IAM scattering, integrating over all orientation angles
    for find, frame in enumerate(traj.frames):
        if print_level:
            print(("Frame %5d of %5d" % (find, len(traj.frames))))
        xyz = ls.Tensor.array(frame.xyz)
        I = ext.compute_diffraction(
            L,
            s2s,
            eta2s,
            xyz,
            fA,
            R2s,
            w2s,
            True if anisotropy == "cos2" else False,
            True if form == "mod" else False,
        )
        # Moment computation
        I0 = 0.5 * (I[:, 0] + I[:, 1])
        I1 = 0.5 * (I[:, 0] - I[:, 1])
        frame.properties["%s-0" % key] = I0
        frame.properties["%s-2" % key] = I1

    return traj


def compute_diffraction_from_moments(
    traj,
    key,
    eta,
):

    for find, frame in enumerate(traj.frames):
        I = np.outer(frame.properties["%s-0" % key], np.cos(0 * eta)) + np.outer(
            frame.properties["%s-2" % key], np.cos(2 * eta)
        )
        frame.properties[key] = I

    return traj


def compute_diffraction_moments_analytical(
    traj,
    key,
    s,
    L,
    mode="xray",
    form="raw",
    anisotropy="perpendicular",
    print_level=False,
):

    """Compute the I(s, eta) elastic scattering moments for a Trajectory.
         See aimsprop/notes/ued for details on this property.

    Notes:
        * All frames for each initial condition (IC) in traj should be aligned so
        that the transition dipole moment from S0 -> Sex at t=0 is on z. This
        is required for proper computation of anisotropy.
        * All frames should be weighted by geometric considerations at the IC
        (e.g., conformational wells, Wigner weights, etc), by the cross
        section for the optical transition at the IC (e.g., oscillator
        strength and excitation energy window), and by the frame weight due
        to non-adiabatic dynamics.

    Params:
        traj (Trajectory) - the Trajectory object to compute the property for (modified in
            place)
        key (str) - the name of the property.
        s (np.ndarray) - list of scattering vector norms in Angstrom^-1. The
            relationship between s and theta (scattering angle) is given as,
                s = 4 pi / L * sin(theta / 2).
        L (float) - effective wavelength of scattering particle (x-ray
            wavelength or UED deBroglie wavelength) in Angstrom. Used to
            convert through scattering angle theta.
        nlebedev (int) - Lebedev number to use for solid angle orientation
            quadrature.
        nomega (int) - number of uniform quadrature points to use for plane
            orientation quadrature.
        mode (str) - 'xray' or 'ued' for selection of form factors
        form (str) - 'raw' or 'mod' for modified/raw diffraction intensities
            I(s) or M(s).
        anisotropy (str) - 'none' or 'perpendicular' or 'parallel'
        print_level (bool) - print progress if true (useful to track long
            property computations)
    Result/Return:
        traj - reference to the input Trajectory object. The properties "key-0"
            and "key-2" are added to each frame of the Trajectory.
    """

    # Validity checks
    if mode not in ["xray", "ued"]:
        raise ValueError("Unknown mode: %s" % mode)
    if form not in ["raw", "mod"]:
        raise ValueError("Unknown form: %s" % form)
    if anisotropy not in ["none", "perpendicular", "parallel"]:
        raise ValueError("Unknown anisotropy: %s" % anisotropy)

    # Compute scattering angles via Bragg equation
    theta = 2.0 * np.arcsin(s * L / (4.0 * np.pi))

    # Get atomic form factors for appropriate x-ray/ued mode
    factors = formfactor.AtomicFormFactor.build_factors(traj.frames[0], mode=mode)

    # Collocate the atomic form factors
    f = np.zeros((len(factors), s.size))
    for A, factor in enumerate(factors):
        f[A, :] = factor.evaluate(qx=0.0, qy=0.0, qz=s)

    # Compute atomic scattering Iat
    D = np.zeros_like(s)
    for A, factor in enumerate(factors):
        D += f[A, :] ** 2

    # Selection fraction
    F = 1.0 if anisotropy == "isotropic" else 1.0 / 3.0

    # Pairs to include
    ABpairs = []
    for A in range(traj.frames[0].xyz.shape[0]):
        for B in range(traj.frames[0].xyz.shape[0]):
            if A >= B:
                continue
            ABpairs.append((A, B))

    # Diffraction moment computation
    for find, frame in enumerate(traj.frames):
        if print_level:
            print(("Frame %5d of %5d" % (find, len(traj.frames))))
        # Geometry
        xyz = frame.xyz
        # Target
        I0 = np.zeros_like(s)
        I2 = np.zeros_like(s)
        for A, B in ABpairs:
            # Geometry
            d = xyz[A, :] - xyz[B, :]
            r2 = np.sum(d ** 2)
            r = np.sqrt(r2)
            sg2 = np.sum(d[:2] ** 2) / r2
            sr = s * r
            # Bessel functions
            J0 = np.sin(sr) / sr
            J0[sr == 0.0] = 1.0
            J1sr = np.sin(sr) / sr ** 3 - np.cos(sr) / sr ** 2
            J1sr[sr == 0.0] = 1.0 / 3.0
            J2 = (3.0 / sr ** 2 - 1.0) * np.sin(sr) / sr - 3.0 * np.cos(sr) / sr ** 2
            J2[sr == 0.0] = 0.0
            # Kernels
            if anisotropy == "isotropic":
                I0 += 2.0 * f[A, :] * f[B, :] * J0
            elif anisotropy == "perpendicular":
                Iz = (
                    2.0
                    * f[A, :]
                    * f[B, :]
                    * (
                        J1sr
                        - (sg2 + (2.0 - 3.0 * sg2) * np.cos(0.5 * theta) ** 2)
                        * J2
                        / 2.0
                    )
                )
                Ix = 2.0 * f[A, :] * f[B, :] * (J1sr - (sg2) * J2 / 2.0)
                I0 += 0.5 * (Iz + Ix)
                I2 += 0.5 * (Iz - Ix)
            elif anisotropy == "parallel":
                I0 += (
                    2.0
                    * f[A, :]
                    * f[B, :]
                    * (
                        J1sr
                        - (sg2 + (2.0 - 3.0 * sg2) * np.sin(0.5 * theta) ** 2)
                        * J2
                        / 2.0
                    )
                )
        # Modified detector pattern
        if form == "raw":
            I0 += F * D
        if form == "mod":
            I0 /= D
            I2 /= D
        # Placement
        frame.properties["%s-0" % (key)] = I0
        frame.properties["%s-2" % (key)] = I2

    return traj


# TODO: These are deprecated, as they are not fully correct for elastic scattering
# def compute_diffraction_moments(
#     traj,
#     key,
#     q,
#     factors,
#     nlebedev,
#     nlebedev2,
#     nomega2,
#     nlegendre=2,
#     print_level=False,
#     ):
#
#     """ Compute the IAM X-Ray Diffraction or UED moments property. See
#         aimsprop/notes/xray for details on these moments.
#
#     Notes:
#         * All frames for each initial condition (IC) in traj should be aligned so
#         that the transition dipole moment from S0 -> Sex at t=0 is on z. This
#         is required for proper computation of I2 (I0 is invariant to this).
#         * All frames should be weighted by geometric considerations at the IC
#         (e.g., conformational wells, Wigner weights, etc), by the cross
#         section for the optical transition at the IC (e.g., oscillator
#         strength and excitation energy window), and by the frame weight due
#         to non-adiabatic dynamics.
#
#     Params:
#         traj (Trajectory) - the Trajectory object to compute the property for (modified in
#             place)
#         key (str) - the name of the property.
#         q (np.ndarray) - the 1d array of |q| values to collocate the
#             diffraction signal to.
#         factors (list of AtomicFormFactor) - the list of AtomicFormFactor
#             objects. The choice of xray/ued is made by the "mode" field of each
#             AtomicFormFactor.
#         nlebedev (int) - a Lebedev number for the number of grid points to use
#             for angular integration of the Legendre coefficients.
#         nlebedev2 (int) - a Lebedev number for the number of grid points to use
#             for angular sampling in the Omega angle (surface angle).
#         nomega2 (int) - the number of grid points to use for angular sampling
#             in the omega angle (surface orientation).
#         nlegendre (even int) - the maximum Legendre polynomial order to
#             evaluate (usually 0 and 2 are the only Legendre polynomial
#             coefficients with any signal).
#         print_level (bool) - print progress if true (useful to track long
#             property computations)
#     Result/Return:
#         traj - reference to the input Trajectory object. The properties
#             key-l are set where l is [0, 2, ..., nlegendre].
#     """
#     if nlegendre % 2: raise ValueError('Can only ask for even Legendre functions')
#
#     # Lebedev grid for S(\vec q)
#     leb = lebedev.Lebedev.build(nlebedev)
#     # Direct product grid for \vec q
#     qx = np.outer(q, leb.x)
#     qy = np.outer(q, leb.y)
#     qz = np.outer(q, leb.z)
#
#     # Rotation quadrature
#     Rs, ws = rotation.rotation_quadrature(nomega=nomega2, nlebedev=nlebedev2)
#
#     for find, frame in enumerate(traj.frames):
#         if print_level:
#             print 'Frame %5d of %5d' % (find, len(traj.frames))
#         # Compute N(\vec q) = \sum_{A} f_A (\vec q) * \exp(-1.j * \vec q * \vec r)
#         N = np.zeros((len(q), len(leb.x)), dtype=complex)
#         for A, factor in enumerate(factors):
#             x = frame.xyz[A,0]
#             y = frame.xyz[A,1]
#             z = frame.xyz[A,2]
#             N += factor.evaluate_N(qx=qx,qy=qy,qz=qz,x=x,y=y,z=z)
#         # Compute I(\vec q) = N(\vec q)**2 for this frame
#         I = (np.abs(N)**2).real
#
#         # Do angle integration
#         Ils = { l: np.zeros((len(q),)) for l in range(0, nlegendre+1, 2) }
#         qxyz = leb.xyz
#         for R2, w2 in zip(Rs, ws):
#             # Account for cos(z)^2 weight
#             cos2 = np.sum(R2[:,2])**2
#             # Rotate the lebedev grid
#             qxyz2 = np.dot(qxyz, R2)
#             qx2 = qxyz2[:,0]
#             qy2 = qxyz2[:,1]
#             qz2 = qxyz2[:,2]
#             # Weight by zonal harmonics (evaluated in rotated grid)
#             Y = legendre.zonal2(qz2, nlegendre)
#             for l in range(0, nlegendre+1, 2):
#                 Ils[l] += w2 * cos2 * np.einsum('qw,w->q', I, leb.w * Y[l/2, :])
#
#         # Assign the properties to frame
#         for l in range(0, nlegendre+1, 2):
#             frame.properties['%s-%d' % (key, l)] = Ils[l]
#
#     return traj
#
# def compute_diffraction_moment0(
#     traj,
#     key,
#     q,
#     factors,
#     nlebedev,
#     print_level=False,
#     ):
#
#     """ Compute the IAM X-Ray Diffraction or UED moment property for only l=0
#         (faster due to lack of rotation quadratures). See aimsprop/notes/xray for
#         details on this moment
#
#     Notes:
#         * This moment is invariant to the orientation of the frames, so alignment to the IC
#         * All frames should be weighted by geometric considerations at the IC
#         (e.g., conformational wells, Wigner weights, etc), by the cross
#         section for the optical transition at the IC (e.g., oscillator
#         strength and excitation energy window), and by the frame weight due
#         to non-adiabatic dynamics.
#
#     Params:
#         traj (Trajectory) - the Trajectory object to compute the property for (modified in
#             place)
#         key (str) - the name of the property.
#         q (np.ndarray) - the 1d array of |q| values to collocate the
#             diffraction signal to.
#         factors (list of AtomicFormFactor) - the list of AtomicFormFactor
#             objects. The choice of xray/ued is made by the "mode" field of each
#             AtomicFormFactor.
#         nlebedev (int) - a Lebedev number for the number of grid points to use
#             for angular integration of the Legendre coefficients.
#         print_level (bool) - print progress if true (useful to track long
#             property computations)
#     Result/Return:
#         traj - reference to the input Trajectory object. The property
#             key-0 is set
#     """
#
#     # Lebedev grid for S(\vec q)
#     leb = lebedev.Lebedev.build(nlebedev)
#     # Direct product grid for \vec q
#     qx = np.outer(q, leb.x)
#     qy = np.outer(q, leb.y)
#     qz = np.outer(q, leb.z)
#
#     for find, frame in enumerate(traj.frames):
#         if print_level:
#             print 'Frame %5d of %5d' % (find, len(traj.frames))
#         # Compute N(\vec q) = \sum_{A} f_A (\vec q) * \exp(-1.j * \vec q * \vec r)
#         N = np.zeros((len(q), len(leb.x)), dtype=complex)
#         for A, factor in enumerate(factors):
#             x = frame.xyz[A,0]
#             y = frame.xyz[A,1]
#             z = frame.xyz[A,2]
#             N += factor.evaluate_N(qx=qx,qy=qy,qz=qz,x=x,y=y,z=z)
#         # Compute I(\vec q) = N(\vec q)**2 for this frame
#         I = (np.abs(N)**2).real
#         # Compute Y00 (for common normalization)
#         Y00 = legendre.zonal2(leb.z, 0)[0]
#         # Integrate over S(2)
#         I0 = np.einsum('qw,w->q', I, leb.w * Y00)
#         # Assign the property to frame
#         frame.properties['%s-%d' % (key, 0)] = I0
#
#     return traj
