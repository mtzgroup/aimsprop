import matplotlib
import numpy as np

import aimsprop as ai

matplotlib.use("Agg")
import os

import matplotlib.pyplot as plt


def build_I2_traj(
    tmax=600.0,
    nframe=26,
    dz=0.5,
):

    """Provides a Trajectory representing I2 aligned in z oscillating over one full cycle.

    Params:
        tmax (float) - maximum time
        nframe (int) - total number of frames
        dz (float) - amplitude of oscillation in Angstrom
    Returns:
        Trajectory with I2 oscillating in z
    """

    # I2 in z
    N = [53] * 2
    xyz = np.array(
        [
            [0.000, 0.000, +1.33],
            [0.000, 0.000, -1.33],
        ]
    )

    ts = np.linspace(0.0, tmax, nframe)

    # Oscillate over tmax
    frames = []
    for t in ts:
        xyz2 = xyz + dz * np.sin(2.0 * np.pi * t / tmax) * np.array(
            [
                [0.0, 0.0, +1.0],
                [0.0, 0.0, -1.0],
            ]
        )
        frame = ai.Frame(
            label=0,
            t=t,
            w=1.0,
            I=1,
            N=N,
            xyz=xyz2,
        )
        frames.append(frame)

    traj = ai.Trajectory(frames)

    return traj


def test_ued():

    """Demonstration of UED molecular movie computation.

    Steps:
    1 - get a Trajectory
    2 - setup a grid in s and eta, for an e- beam with deBroglie wavelength L
    3 - compute diffraction moments I0(s, t) and I2(s, t)
    4 - compute diffraction pattern I(s, eta, t)
    5 - make detector movie of I(s, eta, t)
    6 - make plot of I0(s, t) [isotropic signal]
    7 - sine transform I0(s, t) -> R0(r, t)
    8 - make plot of R0(r, t) [molecular movie]
    9 - Gaussian blur R0(r, t) by 100 fs FWHM
    10 - make plot of R0blur(r, t) [molecular movie]

    As is standard practice, difference signals are used throughout. This
    removes background and signal due to unexcited background. Additionally,
    the modified diffraction signal s M(s) is used in place of I(s) throughout.
    Standard perpendicular pump-probe anisotropy (cos^2) is applied.

    The 3.7 MeV e- beam and s range used here is standard for the UED-I
    experiment at SLAC.
    """

    # => 1 - Trajectory <= #

    # Get I2 oscillation trajectory
    traj = build_I2_traj(
        tmax=600.0,
        nframe=151,
        dz=0.5,
    )
    # Unique times in traj
    t = traj.ts

    # => 2 - Scattering Geometry <= #

    # 3.7 MeV electrons
    L = 0.003
    # Minimal s on detector (limited by e- beam hole)
    smin = 3.5
    # Maximal s on detector (limited by detector area and SNR)
    smax = 12.0
    # Sample s in [smin, smax]
    s = np.linspace(smin, smax, 120)
    # Detector azimuthal angles eta
    eta = np.linspace(0.0, 2.0 * np.pi, 61)
    # 2D detector grid
    ss, ee = np.meshgrid(s, eta, indexing="ij")
    # r to sample in sine transform
    r = np.linspace(0.0, 6.0, 120)

    # => 3 - Diffraction Moments <= #

    # Compute diffraction moments I-0 and I-2 (much faster than directly computing I(s, eta)
    traj = ai.iam.compute_diffraction_moments_fast(
        traj,
        "I",
        s=s,
        L=L,
        mode="ued",  # Set up fA(s) for UED
        nlebedev=4334,  # Very large solid angle grid
        nomega=1,  # Special case, on z, omega grid is not needed
        form="mod",  # Return M(s) instead of I(s)
        # print_level=True,
    )

    # => 4 - Diffraction Pattern <= #

    # Compute diffraction pattern I(s, eta) from I-0 and I-2
    ai.iam.compute_diffraction_from_moments(
        traj,
        "I",
        eta,
    )

    # => 5 - Detector Movie <= #

    # Lots of images for this -> place in movie dir
    dirname = "movie"
    if not os.path.exists(dirname):
        os.makedirs(dirname)

    # Get the difference detector patterns M[t,s,eta]
    I = traj.extract_property("I", diff=True)
    # sM version of difference detector patterns
    I *= np.einsum("i,jk->ijk", np.ones((I.shape[0],)), ss)

    # Make some nice levels for colormaps (uniform across t)
    Imax = np.max(np.abs(I))
    levels = np.linspace(-Imax, +Imax, 127)

    # Make a png image for each frame
    plt.clf()
    fig, ax = plt.subplots(subplot_kw=dict(projection="polar"))
    for tind, tval in enumerate(traj.ts):
        plt.cla()
        h = ax.contourf(
            ee,
            ss,
            I[tind, :, :],
            levels=levels,
            cmap=plt.get_cmap("bwr"),
            extend="both",
        )
        th = np.linspace(0.0, 2.0 * np.pi, 1000)
        ax.plot(th, smin * np.ones_like(th), "-k", linewidth=1.0)
        ax.plot(th, smax * np.ones_like(th), "-k", linewidth=1.0)
        ax.set_theta_offset(np.pi / 2.0)
        ax.set_rlim([0.0, smax])
        ax.set_rticks([])
        ax.set_thetagrids([])
        plt.axis("off")
        plt.savefig(
            "%s/I%03d.png" % (dirname, tind), bbox_inches="tight", transparent=True
        )

    # You can go into movie/ and stitch the png files together into a movie, e.g., with ffmpeg (check out movie.sh)

    # => 6 - Isotropic I0(s,t) Plot <= #

    # Get the isotropic part of the difference signal M0[t,s]
    I0 = traj.extract_property("I-0", diff=True)
    # Work with sM(s)
    I0 *= np.outer(np.ones_like(t), s)

    # Plot I0
    I0max = np.max(np.abs(I0))
    I0levels = np.linspace(-I0max, +I0max, 127)
    plt.clf()
    plt.contourf(t, s, I0.T, levels=I0levels, cmap=plt.get_cmap("bwr"))
    plt.xlabel("t [fs]")
    plt.ylabel("s [$\AA{}^{-1}$]")
    plt.colorbar()
    plt.savefig("I0.pdf", bbox_inches="tight")

    # => 7 - Sine Transform s -> r <= #

    # Compute the sine transform I0[t,s] -> R0[t,r]
    R0 = ai.compute_sine_transform(I0, s, r, a=7.0 ** (-2))

    # => 8 - Isotropic R0(r,t) Plot <= #

    R0max = np.max(np.abs(R0))
    R0levels = np.linspace(-R0max, +R0max, 127)
    plt.clf()
    plt.contourf(t, r, R0.T, levels=R0levels, cmap=plt.get_cmap("bwr"))
    plt.xlabel("t [fs]")
    plt.ylabel("r [$\AA{}$]")
    plt.colorbar()
    plt.savefig("R0.pdf", bbox_inches="tight")

    # => 9 - Gaussian Blurring in Time <= #

    # Blur in time to 100 fs FWHM
    R0blur = ai.compute_time_blur(
        R0,
        t,
        t,
        fwhm=100.0,
    )

    # => 10 - Isotropic R0blur(r,t) Plot <= #

    # Original levels
    R0max = np.max(np.abs(R0))
    R0levels = np.linspace(-R0max, +R0max, 127)
    plt.clf()
    plt.contourf(t, r, R0blur.T, levels=R0levels, cmap=plt.get_cmap("bwr"))
    plt.xlabel("t [fs]")
    plt.ylabel("r [$\AA{}$]")
    plt.colorbar()
    plt.savefig("R0blur.pdf", bbox_inches="tight")


if __name__ == "__main__":

    test_ued()
