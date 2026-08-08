#!/usr/bin/env bash

# --- Modules ---
module purge 
module load PrgEnv-gnu
# module load PrgEnv-aocc/8.4.0
module load craype-x86-rome
module load cray-libsci
module load cray-fftw
module load load-epcc-module
module load cmake/3.29.4
module load cray-hdf5-parallel
module list

export COMPILER_ID=GNU
export CC=cc
export CXX=CC
export FC=ftn

# --- Project root ---
CWD='/mnt/lustre/a2fs-work3/work/e773/e773/pounds/PadeOps'

# export FFTW_PATH="${CWD}/dependencies/fftw-3.3.10"
export FFTW_PATH=${FFTW_ROOT}

# export HDF5_PATH="${CWD}/dependencies/hdf5-1.14.3/build"
export HDF5_PATH=${HDF5_DIR}

export FFTPACK_PATH="${CWD}/dependencies/fftpack"
export DECOMP_PATH="${CWD}/dependencies/2decomp_fft"
# export VTK_IO_PATH="${CWD}/dependencies/Lib_VTK_IO/build"

export CMAKE_PREFIX_PATH="${HDF5_PATH}:${FFTW_PATH}${CMAKE_PREFIX_PATH:+:${CMAKE_PREFIX_PATH}}"

# --- Architecture flags ---
export ARCH_OPT_FLAG="-march=native"

export OMP_NUM_THREADS=1
