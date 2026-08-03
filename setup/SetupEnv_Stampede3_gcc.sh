#!/bin/bash
module purge
module load gcc
module load cmake/3.31.9
module load openmpi/5.0.8
module load mkl

CWD=`pwd`
export COMPILER_ID=GNU
export CC=mpicc
export CXX=mpicxx
export FC=mpif90
export CFLAGS="-O3 -fno-lto"
export CXXFLAGS="-O3 -fno-lto"
export FFLAGS="-O3 -fno-lto"
export FCFLAGS="-O3 -fno-lto"
export LDFLAGS="-fno-lto"

export FFTW_PATH=${CWD}/dependencies/gcc/fftw-3.3.10
export DECOMP_PATH=${CWD}/dependencies/gcc/2decomp_fft
export HDF5_PATH=${CWD}/dependencies/gcc/hdf5-1.14.3/build
export FFTPACK_PATH=dummy
export ARCH_OPT_FLAG="-march=skylake-avx512"
