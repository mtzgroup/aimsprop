import numpy as np

from . import lebedev


def rotation_quadrature(
    nlebedev,
    nomega,
):

    """Generate a rotation quadrature for SO(3) via a direct product of
        Lebedev grids on S(2) and linear grids on L2(2*pi).

    Params:
        nlebedev (int) - a Lebedev grid number for the number of nodes in the
            solid angle.
        nomega (int) - the number of nodes in the omega angle.
    Returns:
        Rs (list of np.ndarray of shape (3,3)) - Rotation matrices representing
            the quadrature nodes.
        ws (list of float) - the quadrature weights.
    """

    leb = lebedev.Lebedev.build(nlebedev)

    omega = np.linspace(0.0, 2.0 * np.pi, nomega, endpoint=False)
    womega = np.ones_like(omega) / nomega

    Rs = []
    ws = []

    for omega2, womega2 in zip(omega, womega):
        # Rotate about z
        R1 = np.array(
            [
                [np.cos(omega2), +np.sin(omega2), 0.0],
                [-np.sin(omega2), np.cos(omega2), 0.0],
                [0.0, 0.0, 1.0],
            ]
        )
        for U in range(len(leb.x)):
            # Rotate down to Omega from z
            n = np.array((0.0 - leb.x[U], 0.0 - leb.y[U], 1.0 - leb.z[U]))
            if leb.z[U] != 1.0:
                n /= np.sqrt(np.sum(n * n))
            R2 = np.eye(3) - 2.0 * np.outer(n, n)  # Householder reflector
            if leb.z[U] != 1.0:
                R2 *= -1.0  # Parity fix
            Rs.append(np.dot(R1, R2))
            ws.append(
                leb.w[U] * womega2 / (4.0 * np.pi)
            )  # Lebedev grids are normalized to 4*pi

    return Rs, ws
