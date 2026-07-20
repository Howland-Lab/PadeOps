#!/bin/bash
module purge 
module load gcc
module load cmake
module load fftw
module load mvapich2/2.3.6
module load hdf5
module load openblas/0.3.17
module list 

export COMPILER_ID=GNU
export FC=mpif90 
export CC=mpicc
export CXX=mpicpc
export FFTW_PATH=${FFTW_HOME}
export HDF5_PATH=${HDF5_HOME}
export DECOMP_PATH=/anvil/projects/x-atm170028/karim/PadeOps/dependencies/gcc/2decomp_fft
export FFTPACK_PATH=dummy
export ARCH_OPT_FLAG="-march=core-avx2"
