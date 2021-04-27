# aimsprop

## Description

A repository for the representation and manipulation of AIMS-type simulations,
particularly for use in computing time-dependent properties (bond lengths,
angles, torsions, X-ray scattering, UED, photoelectron, UV-Vis, etc) in a
semi-uniform manner.

## Installation

```sh
pip install aimsprop
```

## Basic Information

### Bundle Skeleton

bundle.py - Frame/Bundle classes to represent AIMS simulations

### Parsers

fms90.py - Parser for FMS90 outputs
xyz.py - Parser for XYZ files from adiabatic dynamics

### Properties

geom.py - Easy geometric properties like bond length, bond angle, torsion angle, and out-of-plane angles
pop.py - State populations over time
ued.py - Simple pairwise distance-based I(R) properties for ultrafast electron diffraction (UED)
iam/ - X-ray scattering and ultrafast electron diffraction in the independent atom model (IAM)

### Utility

blur.py - Gaussian blurring of properties (makes a nice projection of the density matrix)
wrap.py - Wrap/unwrap periodic properties
rotate.py - Rotate xyz coordinates by a given rotation matrix

### Plotting

plot.py - Generic interfaces for plotting scalar and vector data associated with Frame objects.

## Authors

- Rob Parrish
- Monika Williams
- Hayley Weir
- Colton Hicks
- Alessio Valentini
- Alice R. Walker
