import numpy as np


def legendre(x, n):

    y = np.zeros((n + 1, len(x)))
    y[0, :] = 1.0
    if n > 0:
        y[1, :] = x
    for k in range(1, n):
        y[k + 1, :] = (
            1.0 / (k + 1.0) * ((2.0 * k + 1.0) * x * y[k, :] - k * y[k - 1, :])
        )

    return y


def legendre2(x, n):

    if n % 2:
        raise ValueError("Can only ask for even Legendre polynomials")

    y = legendre(x, n)

    y2 = np.zeros((n / 2 + 1, len(x)))
    for k in range(0, n + 1, 2):
        y2[k / 2, :] = y[k, :]

    return y2


def zonal(z, n):

    y = legendre(z, n)
    for k in range(0, n + 1):
        y[k, :] *= np.sqrt((2.0 * k + 1.0) / (4.0 * np.pi))
    return y


def zonal2(z, n):

    y = legendre2(z, n)
    for k in range(0, n + 1, 2):
        y[k / 2, :] *= np.sqrt((2.0 * k + 1.0) / (4.0 * np.pi))
    return y
