# aimsprop

To use aimsprop, check out the documentation [here](https://mtzgroup.github.io/aimsprop/)!

To develop / contribute to the aimsprop code,
see [CONTRIBUTING.md](https://mtzgroup.github.io/aimsprop/CONTRIBUTING/) to get your environment and branch set up.

For a tutorial on the basic functionality of aimsprop, go [here](https://mtzgroup.github.io/aimsprop/tutorial/).

## Description

A repository for the representation and manipulation of AIMS-type trajectories, particularly for use in computing
time-dependent properties (bond lengths, angles, torsions, X-ray scattering, UED, photoelectron, UV-Vis, etc) in a
semi-uniform manner.

## Directories

- aimsprop - the Python AIMS trajectory/property code
- docs - mkdocs documentation that is published [here](https://mtzgroup.github.io/aimsprop/)
- tests - unit tests for aimsprop code using pytest
- notes - notes used in implementing
- literature - key papers/book chapters used in implementing certain properties.
- scripts - utility scripts for common ops like plotting absorption spectra, generating difference densities from Molden
  files, etc (no warrenty on these - the usual workflow is to copy these into your own project directories and modify
  until they work)
- examples - examples of use cases of property computations (these are not the most helpful, better to look at
  the [tutorial](https://github.com/mtzgroup/aimsprop/blob/develop/docs/tutorial.md)).

## Todos

- H-Bond coordinates
- Better UED cross sections
- OOP angles
- Update examples
- Adiabatic dynamics are weirdly done in 0.5 fs, non-adiabatic dynamics are in 20 au. These are not the same. We should
  fix/accomodate this.

## Authors

Rob Parrish, Monika Williams, Hayley Weir, Colton Hicks, Alice Walker, Alessio Valentini
