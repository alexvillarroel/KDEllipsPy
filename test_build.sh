#!/bin/bash
set -e
cd kdellipspy/axitra/src
PYTHON=/home/alex/.conda/envs/kdellipspy/bin/python3.11
FC=gfortran
FFLAGS="-Ofast -cpp -fopenmp -fPIC"
CC=gcc
CFLAGS="-O3 -std=gnu89"

echo "Generating f2py wrappers..."
$PYTHON -m numpy.f2py convmPy.F90 parameter.f90 fft2cd.f90 fsource.f90 dll2km.f90 -m convmPy only: moment_conv :

echo "Compiling Fortran wrappers..."
$FC $FFLAGS -c convmPy-f2pywrappers2.f90

echo "Compiling C wrappers..."
F2PY_INC=$( $PYTHON -c "import numpy.f2py; print(numpy.f2py.get_include())" )
$CC $CFLAGS -fPIC -c convmPymodule.c \
    -I$( $PYTHON -c "import numpy; print(numpy.get_include())" ) \
    -I$( $PYTHON -c "import sysconfig; print(sysconfig.get_path('include'))" ) \
    -I$F2PY_INC

$CC $CFLAGS -fPIC -c $F2PY_INC/fortranobject.c \
    -I$( $PYTHON -c "import numpy; print(numpy.get_include())" ) \
    -I$( $PYTHON -c "import sysconfig; print(sysconfig.get_path('include'))" ) \
    -I$F2PY_INC

echo "Linking shared library..."
EXT=$( $PYTHON -c "import sysconfig; print(sysconfig.get_config_var('EXT_SUFFIX'))" )
$FC -shared -o convmPy$EXT \
    convmPymodule.o fortranobject.o convmPy-f2pywrappers2.o \
    parameter.o dimension1.o fft2cd.o fsource.o dll2km.o \
    -lgomp

echo "Done! Generated files:"
ls -la convmPy*
