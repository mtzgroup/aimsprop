"""
Written by Rob Parrish

This script reads TeraChem TDDFT and/or FOMO output files and automatically plots absorption spectra.

Each individual excited state is given a color, with the overall spectrum plotted in black on top.

The script includes options to plot with Lorentzian or Gaussian broadening. To change it, modify the ```plot_spectra``` function to change:

```python
Ostate += spectrum[state, 1] * lorentzian(E, spectrum[state, 0], delta)
```

to

```python
Ostate += spectrum[state, 1] * gaussian(E, spectrum[state, 0], delta)
```

It also includes options to plot in eV or nm, with `plot_spectra` or `plot_spectra_nm` respectively.

Modify the following code block with appropriate paths for your system, choice of fomo or tddft for `read_*_output`, and choice of eV of nm for `plot_*`.

```
filenames = glob.glob("/path/to/*/filename*.out")
spectra = [read_fomo_outfile(x) for x in filenames]
spectra = [x for x in spectra if x is not None]
plot_spectra("spectra.pdf", spectra, E=np.linspace(5.0, 9.0, 1000), delta=0.10)
```


"""
import glob
import re

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

matplotlib.use("Agg")


# Python script to automatically extract TeraChem output (TD-DFT or CASCI) and plot
# Lorentzian or Gaussian broadened absorption spectrum.

# Can plot spectra in eV or nm, depending on how you modify the exection.


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


def gaussian(x, x0, delta):
    return np.exp(-np.power(x - x0, 2.0) / (2 * np.power(delta, 2.0)))


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


def plot_spectra_nm(
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

    nm = (1239.84193) / E

    plt.clf()
    for Oind, Ostate in enumerate(Ostates):
        plt.plot(nm, Ostate, color=colors[Oind])

    plt.plot(nm, Ototal, "-k")
    plt.xlabel(r"$\lambda$ [nm]")
    plt.ylabel("Cross Section [-]")
    plt.savefig(filename)


# --- Execution section --- #

if __name__ == "__main__":
    filenames = glob.glob("0*/output.dat")
    spectra = [read_fomo_outfile(x) for x in filenames]
    spectra = [x for x in spectra if x is not None]
    plot_spectra("spectra.pdf", spectra, E=np.linspace(5.0, 9.0, 1000), delta=0.10)
