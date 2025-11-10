#!/usr/bin/env bash
# Archer2 GNU + CrayPE environment for building PadeOps

# --- Modules ---
module purge
module load PrgEnv-gnu
module load craype-x86-rome        # target AMD Rome (Zen2) – replaces manual -march
# module load cmake
module load cray-libsci
module load cray-fftw
module load cray-hdf5-parallel
module list

# --- Compilers (use Cray wrappers) ---
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

export CMAKE_PREFIX_PATH="${HDF5_PATH}:${FFTW_PATH}:${VTK_IO_PATH}:${CMAKE_PREFIX_PATH}"

# --- Architecture flags ---
# With craype-x86-rome + wrappers, you usually do NOT need to set -march/-mtune.
# Leave this empty, or only append safe optimisations that won't fight wrappers.
export ARCH_OPT_FLAG=""

# Example of safe extras if you insist:
# export ARCH_OPT_FLAG="-O3 -fopenmp"   # (only if your code uses OpenMP)

# --- Runtime sanity for MPI-only builds ---
export OMP_NUM_THREADS=1
