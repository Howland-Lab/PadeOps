#!/bin/bash
module purge 
module load gcc/14.2.0
module load openmpi/5.0.5
module load openblas/0.3.27
module load fftw/3.3.8 
module load hdf5/1.10.7
module load cmake/3.20.0
module list 

export COMPILER_ID=GNU
export FC=mpif90 
export CC=mpicc
export CXX=mpicpc
export FFTW_PATH=${FFTW_HOME}
export HDF5_PATH=${HDF5_HOME}
export DECOMP_PATH=/anvil/projects/x-atm170028/karim/PadeOps/dependencies/gcc/2decomp_fft
export FFTPACK_PATH=dummy
export ARCH_OPT_FLAG="-march=znver3"
