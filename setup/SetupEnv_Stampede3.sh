#!/bin/bash

module load cmake
module load intel impi

CWD=`pwd`
export COMPILER_ID=Intel
export FC=mpiifort
export CC=mpiicc
export CXX=mpiicpc
export FFTW_PATH=${CWD}/dependencies/fftw-3.3.10
export DECOMP_PATH=${CWD}/dependencies/2decomp_fft
export VTK_IO_PATH=${CWD}/dependencies/Lib_VTK_IO/build
export HDF5_PATH=${CWD}/dependencies/hdf5-1.14.3/build
export FFTPACK_PATH=${CWD}/dependencies/fftpack
export ARCH_OPT_FLAG="-xCORE-AVX512"
