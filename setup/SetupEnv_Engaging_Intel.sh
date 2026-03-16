#!/bin/bash
module purge
module load StdEnv
module load community-modules
module load cmake
module load intel-hpc/2025.2.1.44

# Force Intel MPI wrappers to use LLVM-based Intel compilers
export I_MPI_CC=icx
export I_MPI_CXX=icpx
export I_MPI_FC=ifx
export I_MPI_F90=ifx

CWD=$(pwd)

export COMPILER_ID=Intel
export FC=mpiifort
export CC=mpiicc
export CXX=mpiicpc

export FFTW_PATH=${CWD}/dependencies/fftw-3.3.10
export DECOMP_PATH=${CWD}/dependencies/2decomp_fft
export VTK_IO_PATH=${CWD}/dependencies/Lib_VTK_IO/build
export HDF5_PATH=${CWD}/dependencies/hdf5-1.8.18
export FFTPACK_PATH=${CWD}/dependencies/fftpack

export ARCH_OPT_FLAG="-O3 -xHost"