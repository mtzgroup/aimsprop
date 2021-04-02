# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/), and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [unreleased]

### Added

### Changed

### Removed

## [0.2.0]

### Added

- Extended documentation to include a Code Reference section and much more comprehensive documentation of the main objects.
- GitHub Actions for code quality:
  - `black`
  - `isort`
  - `flake8` (on selected files; adding this functionality file by file due to the size of formatting problems in current repo)
  - `mdl` (markdownlint)
  - Run all tests on each code push
- Extensive documentation of the code
- Tutorial for new users of the package
- CONTRIBUTING.md to guide new users on how to make contributions to the `aimsprop` package
- `mkdocs` to build website from documentation
- Geometry operation tests (`compute_bond`, `compute_angle`, `compute_torsion`, `compute_oop`, `compute_transfer_cood`)

### Changed

- Formatted all code with `black` and `isort`

### Removed

- `pyaims` parser that Monika claimed was no longer needed

## [0.1.0] - 2021-03-05

### Added

- Added basic code tests to help assess that python2->3 transition worked without impacting functionality
- Packaged application using flit

### Changed

- Updated code from python2 -> python3

[unreleased]: https://github.com/mtzgroup/aimsprop/compare/0.2.0...HEAD
[0.2.0]: https://github.com/mtzgroup/aimsprop/releases/tag/0.2.0
[0.1.0]: https://github.com/mtzgroup/aimsprop/releases/tag/0.1.0
