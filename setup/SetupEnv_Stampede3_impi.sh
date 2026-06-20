#!/bin/bash
module purge
module load cmake/3.31.9
module load intel impi
module load fftw3/3.3.10
module load phdf5/2.0.0

CWD=`pwd`
export COMPILER_ID=Intel
export FC=mpiifort
export CC=mpiicc
export CXX=mpiicpc
export FFTW_PATH=${TACC_FFTW3_DIR}
export DECOMP_PATH=${CWD}/dependencies/impi/2decomp_fft
export HDF5_PATH=${TACC_HDF5_DIR}
export FFTPACK_PATH=dummy
export ARCH_OPT_FLAG="-xCORE-AVX512"
