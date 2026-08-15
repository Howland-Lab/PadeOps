#!/bin/bash
module purge 
module load intel/2024.1 
module load impi/2021.12
module load openblas/0.3.17
module load libszip/2.1.1
module load fftw/3.3.8
module load cmake/3.20.0
module list 

export COMPILER_ID=Intel
export FC=mpiifx
export CC=mpiicx
export CXX=mpiicpx
export I_MPI_F90=ifx
export I_MPI_CC=icx
export I_MPI_CXX=icpx
export FFTW_PATH=${FFTW_HOME}
export HDF5_PATH=/anvil/projects/x-atm170028/karim/PadeOps/dependencies/impi/hdf5-1.14.3/build
export DECOMP_PATH=/anvil/projects/x-atm170028/karim/PadeOps/dependencies/impi/2decomp_fft
export FFTPACK_PATH=dummy
export ARCH_OPT_FLAG="-march=core-avx2"
