# aimsprop

Read the [documentation](https://mtzgroup.github.io/aimsprop/)!

## Description

A repository for the representation and manipulation of AIMS-type trajectories,
particularly for use in computing time-dependent properties (bond lengths,
angles, torsions, X-ray scattering, UED, photoelectron, UV-Vis, etc) in a
semi-uniform manner.

## Directories

aimsprop - the Python AIMS trajectory/property code
notes - notes used in implementing
literature - key papers/book chapters used in implementing certain properties.
examples - examples of use cases of property computations
util - utility scripts for common ops like plotting absorption spectra,
generating difference densities from Molden files, etc (no warrenty on
these - the usual workflow is to copy these into your own project
directories and modify until they work)

## Todos

- Population/Excited state lifetime plots (incl. spaghetti)
- H-Bond coordinates
- Better UED cross sections
- OOP angles
- Update examples
- Adiabatic dynamics are weirdly done in 0.5 fs, non-adiabatic dynamics are in 20
  au. These are not the same. We should fix/accomodate this.

## Authors

- Rob Parrish
- Monika Williams
- Hayley Weir
- Colton Hicks
