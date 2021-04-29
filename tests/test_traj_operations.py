# import os

import numpy as np

import aimsprop as ai


def test_major_bundle_operations(tmp_path, bundle):
    # 1. Plot Bond Distances (Spaghetti + Blur)
    # Compute the a bond distance property for all Frames in traj
    ai.compute_bond(bundle, "R01", 0, 1)
    ai.plot_scalar(
        tmp_path / "R.pdf",
        bundle,
        "R01",
        ylabel=r"$R_{CC} [\AA{}]$",
        time_units="fs",
        state_colors=["r", "b"],
        plot_average=True,
    )

    # Blur the bond distance (convolve)
    R = np.linspace(0.5, 3.0, 50)
    ai.blur_property(bundle, "R01", "Rblur", R, alpha=8.0)
    #    Plot the heat map of blurred bond distance
    ai.plot_vector(
        tmp_path / "Rblur.pdf",
        bundle,
        "Rblur",
        y=R,
        ylabel=r"$R [\AA{}]$",
        time_units="fs",
        nlevel=64,
    )

    # 2. Plot of Torsion Angle (Spaghetti + Blur)

    # Compute the a torsion angle property for all Frames in traj
    ai.compute_torsion(bundle, "T0123", 0, 1, 2, 3)
    ai.unwrap_property(bundle, "T0123", 360.0)
    ai.plot_scalar(
        tmp_path / "T.pdf",
        bundle,
        "T0123",
        ylabel=r"$\Theta [^{\circ{}}]$",
        time_units="fs",
        state_colors=["r", "b"],
        plot_average=True,
    )

    # Blur the torsion
    T = np.linspace(-180.0, +180.0, 100)
    ai.blur_property(bundle, "T0123", "Tblur", T, alpha=0.02)
    #    Plot the heat map of blurred torison
    ai.plot_vector(
        tmp_path / "Tblur.pdf",
        bundle,
        "Tblur",
        y=T,
        ylabel=r"$Theta [^{\circ{}}]$",
        time_units="fs",
        nlevel=64,
    )

    # 3. UED Cross Section

    # Compute the "simple" form of the UED cross section in R
    R = np.linspace(1.0, 6.0, 50)
    ai.compute_ued_simple(bundle, "UED", R=R, alpha=8.0)

    # Plot the heat map of the UED cross section detailed above
    ai.plot_vector(
        tmp_path / "UED.pdf",
        bundle,
        "UED",
        y=R,
        ylabel=r"$R [\AA{}]$",
        time_units="fs",
        diff=True,
    )
