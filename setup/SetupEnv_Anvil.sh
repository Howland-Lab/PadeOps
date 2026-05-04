#!/bin/bash
module purge 
module load intel
module load cmake
module load mvapich2
module load intel-mkl
module list 

export COMPILER_ID=Intel
export FC=mpiifort 
export CC=mpiicc
export CXX=mpiicpc
export FFTW_PATH=/anvil/projects/x-atm170028/karim/PadeOps/dependencies/fftw-3.3.10
export DECOMP_PATH=/anvil/projects/x-atm170028/karim/PadeOps/dependencies/2decomp_fft
export VTK_IO_PATH=/anvil/projects/x-atm170028/padeops_setup/dependencies/Lib_VTK_IO/build
export HDF5_PATH=/anvil/projects/x-atm170028/padeops_setup/dependencies/hdf5-1.8.18
export FFTPACK_PATH=/anvil/projects/x-atm170028/padeops_setup/dependencies/fftpack
export ARCH_OPT_FLAG="-march=core-avx2"
