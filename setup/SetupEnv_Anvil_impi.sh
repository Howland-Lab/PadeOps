#!/bin/bash
module purge 
module load intel/19.1.3.304  
module load impi/2019.9.304
module load hdf5/1.10.7
module load fftw/3.3.8
module load cmake/3.20.0
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
