# KDEllipsPy

KDEllipsPy is a Python-based tool for kinematic earthquake rupture inversion using elliptical slip distributions and the Neighbourhood Algorithm (NA).

## Installation

The recommended way to install KDEllipsPy is using a Conda environment. This ensures all scientific dependencies and compilers are correctly managed.

### 1. Create Conda Environment

```bash
conda env create -f environment.yml
conda activate kdellipspy
```

### 2. Install the Package

This will install the Python dependencies and automatically compile the Axitra Fortran engine.

```bash
pip install -e .
```

*Note: If you are on a system where you cannot use Conda, ensure you have `gfortran` and `make` installed on your system before running the pip command.*

## Usage

You can use the CLI to run inversions or use the package as a library in your own scripts/notebooks.

```bash
kdellipspy --help
```

## Features

- Kinematic rupture modeling with elliptical patches.
- Support for multiple stations and 3-component waveforms.
- Integration with Axitra for Green's function computation.
- Neighbourhood Algorithm (NA) and MCMC inversion drivers.
- Comprehensive visualization suite.
