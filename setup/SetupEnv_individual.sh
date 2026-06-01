#module purge
#module load cmake/3.16.1 intel/19.1.1 impi/19.0.9
#module load cmake/3.8.2 intel/18.0.0
module load cmake/3.17.3 intel/2021.3.0
export COMPILER_ID=Intel
export FC=mpiifort
export CC=mpiicc
export CXX=mpiicpc
export FFTW_PATH=/home/kheck/PadeOps/PadeOps/dependencies/fftw-3.3.5
export HDF5_PATH=/home/kheck/PadeOps/PadeOps/dependencies/hdf5-1.8.18
export FFTPACK_PATH=/home/kheck/PadeOps/PadeOps/dependencies/fftpack
export VTK_IO_PATH=/home/kheck/PadeOps/PadeOps/dependencies/Lib_VTK_IO/build
export DECOMP_PATH=/home/kheck/PadeOps/PadeOps/dependencies/2decomp_fft
export ARCH_OPT_FLAG="-xCOMMON-AVX512 -axCORE-AVX512,MIC-AVX512 -qopt-zmm-usage=high"
