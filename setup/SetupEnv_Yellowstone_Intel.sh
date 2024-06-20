#!/bin/bash

module purge
module load cmake
module load intel
module load impi
module load mkl

CWD=`pwd`
export COMPILER_ID=Intel
export FC=mpiifort
export CC=mpiicc
export CXX=mpiicpc
export FFTW_PATH=${CWD}/dependencies/fftw-3.3.10
export DECOMP_PATH=${CWD}/dependencies/2decomp_fft
export VTK_IO_PATH=${CWD}/dependencies/Lib_VTK_IO/build
export HDF5_PATH=${CWD}/dependencies/hdf5-1.8.18
export FFTPACK_PATH=${CWD}/dependencies/fftpack
export ARCH_OPT_FLAG="-axCORE-AVX2,MIC-AVX2 -qopt-zmm-usage=high"

# export ARCH_OPT_FLAG="-xCORE-AVX2 -O3 -axCORE-AVX2"
# export ARCH_OPT_FLAG="-march=core-avx2 -mtune=core-avx2"
