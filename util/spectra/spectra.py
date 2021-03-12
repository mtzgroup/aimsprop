#!/usr/bin/env python
import glob
import re

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def read_tddft_outfile(filename):

    lines = open(filename).readlines()

    if len([x for x, y in enumerate(lines) if re.match(r"^\s+Job finished", y)]) == 0:
        return None

    lines = lines[
        (
            [
                x
                for x, y in enumerate(lines)
                if re.match(r"^\s+Final Excited State Results:", y)
            ][0]
            + 4
        ) :
    ]
    vals = []
    for line in lines:
        mobj = re.match(r"^\s*(\d+)\s+\S+\s+(\S+)\s+(\S+)", line)
        if mobj:
            vals.append([float(mobj.group(2)), float(mobj.group(3))])
            assert int(mobj.group(1)) == len(vals)
    return np.array(vals)


def read_fomo_outfile(filename):

    lines = open(filename).readlines()

    if len([x for x, y in enumerate(lines) if re.match(r"^\s+Job finished", y)]) == 0:
        return None

    Es = []
    re_ex = re.compile(r"^\s*(\d+)\s+singlet\s+\S+\s+\S+\s+(\S+)\s+\S+\s*$")
    for line in lines:
        mobj = re.match(re_ex, line)
        if mobj:
            Es.append(float(mobj.group(2)))
            assert int(mobj.group(1)) == len(Es) + 1

    Os = []
    re_osc = re.compile(r"^\s*1\s+->\s+(\d+)\s+\S+\s+\S+\s+\S+\s+\S+\s+(\S+)\s*$")
    for line in lines:
        mobj = re.match(re_osc, line)
        if mobj:
            Os.append(float(mobj.group(2)))
            assert int(mobj.group(1)) == len(Os) + 1
    assert len(Os) == len(Es)

    return np.array([[x, y] for x, y in zip(Es, Os)])


def lorentzian(
    x,
    x0,
    delta,
):

    return 0.5 * delta / np.pi * 1.0 / ((x - x0) ** 2 + (0.5 * delta) ** 2)


def plot_spectra(
    filename,
    spectra,
    E=np.linspace(3.5, 7.0, 1000),
    delta=0.05,
):

    Nspectra = len(spectra)
    Nstate = spectra[0].shape[0]

    Ostates = []
    for state in range(Nstate):
        Ostate = np.zeros_like(E)
        for spectrum in spectra:
            Ostate += spectrum[state, 1] * lorentzian(E, spectrum[state, 0], delta)
        Ostate /= Nspectra
        Ostates.append(Ostate)

    Ototal = np.zeros_like(E)
    for Ostate in Ostates:
        Ototal += Ostate

    cmap = matplotlib.cm.get_cmap("jet")
    colors = [cmap(float(x) / (Nstate - 1)) for x in reversed(list(range(Nstate)))]

    plt.clf()
    for Oind, Ostate in enumerate(Ostates):
        plt.plot(E, Ostate, color=colors[Oind])

    plt.plot(E, Ototal, "-k")
    plt.xlabel(r"$\Delta E$ [eV]")
    plt.ylabel("Cross Section [-]")
    plt.savefig(filename)


filenames = glob.glob("0*/output.dat")
spectra = [read_fomo_outfile(x) for x in filenames]
spectra = [x for x in spectra if x is not None]
plot_spectra("spectra.pdf", spectra, E=np.linspace(5.0, 9.0, 1000), delta=0.10)
