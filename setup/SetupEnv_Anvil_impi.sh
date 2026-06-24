#!/bin/bash
module purge 
module load intel
module load cmake
module load fftw
module load impi
module load hdf5
module list 

export COMPILER_ID=Intel
export FC=mpiifort 
export CC=mpiicc
export CXX=mpiicpc
export FFTW_PATH=${FFTW_HOME}
export HDF5_PATH=${HDF5_HOME}
export DECOMP_PATH=/anvil/projects/x-atm170028/karim/PadeOps/dependencies/impi/2decomp_fft
export FFTPACK_PATH=dummy
export ARCH_OPT_FLAG="-march=core-avx2"
