#module purge
#module load cmake/3.16.1 intel/19.1.1 impi/19.0.9
#module load cmake/3.8.2 intel/18.0.0
module load cmake/3.16.1 intel/18.0.0
export COMPILER_ID=Intel
export FC=mpiifort
export CC=mpiicc
export CXX=mpiicpc
export FFTW_PATH=/work/04076/tg833754/stampede2/padeops_dependencies/fftw-3.3.10
export HDF5_PATH=/work/04076/tg833754/stampede2/padeops_dependencies/hdf5-1.8.18
export FFTPACK_PATH=/work/04076/tg833754/stampede2/padeops_dependencies/fftpack
export VTK_IO_PATH=/work/04076/tg833754/stampede2/padeops_dependencies/Lib_VTK_IO/build
export DECOMP_PATH=/work/04076/tg833754/stampede2/padeops_dependencies/2decomp_fft
export ARCH_OPT_FLAG="-xCOMMON-AVX512 -axCORE-AVX512,MIC-AVX512 -qopt-zmm-usage=high"
