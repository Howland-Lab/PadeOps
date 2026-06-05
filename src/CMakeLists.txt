file(GLOB utilities_source_files utilities/*.F90)
file(GLOB derivatives_source_files derivatives/*.F90)
file(GLOB filters_source_files filters/*.F90)
file(GLOB io_source_files io/*.F90)
file(GLOB grid_source_files grid/*.F90)
file(GLOB compressible_source_files compressible/*.F90)
file(GLOB solid_source_files solid/*.F90)
file(GLOB poisson_source_files incompressible/poisson/*.F90)
file(GLOB incompressible_source_files incompressible/*.F90)

# Where to look for FFTW
set(FFTW_LIBRARY_PATH "${FFTW_PATH}/lib")
set(FFTW_INCLUDE_PATH "${FFTW_PATH}/include")

# Include directories
include_directories(${MPI_INCLUDE_PATH} ${FFTW_INCLUDE_PATH} ${DECOMP_INCLUDE_PATH} ${HDF5_INCLUDE_PATH})

# Link directories
link_directories(${FFTW_LIBRARY_PATH} ${DECOMP_LIBRARY_PATH} ${HDF5_LIBRARY_PATH})

# Create a library called PadeOps
add_library(PadeOps STATIC ${utilities_source_files} ${derivatives_source_files} ${filters_source_files} ${io_source_files} )
if ( CMAKE_Fortran_COMPILER_ID MATCHES "Intel" )
   target_link_libraries(PadeOps fftw3 2decomp_fft  ${HDF5_LIBRARY_PATH}/libhdf5hl_fortran.a ${HDF5_LIBRARY_PATH}/libhdf5_hl.a ${HDF5_LIBRARY_PATH}/libhdf5_fortran.a ${HDF5_LIBRARY_PATH}/libhdf5.a ${MPI_LIBRARIES} -lz -ldl -lm -Wl,-rpath -Wl,${HDF5_LIBRARY_PATH})
elseif( CMAKE_Fortran_COMPILER_ID MATCHES "Cray" )
   target_link_libraries(PadeOps fftw3 2decomp_fft ${HDF5_LIBRARY_PATH}/libhdf5hl_fortran.a ${HDF5_LIBRARY_PATH}/libhdf5_hl.a ${HDF5_LIBRARY_PATH}/libhdf5_fortran.a ${HDF5_LIBRARY_PATH}/libhdf5.a ${MPI_LIBRARIES} m z dl)
elseif( CMAKE_Fortran_COMPILER_ID MATCHES "AOCC|AMD" )
   target_link_libraries(PadeOps fftw3 2decomp_fft ${HDF5_LIBRARY_PATH}/libhdf5hl_fortran.a ${HDF5_LIBRARY_PATH}/libhdf5_hl.a ${HDF5_LIBRARY_PATH}/libhdf5_fortran.a ${HDF5_LIBRARY_PATH}/libhdf5.a ${MPI_LIBRARIES} m z dl)
elseif( CMAKE_Fortran_COMPILER_ID MATCHES "GNU" )
   target_link_libraries(PadeOps
      fftw3
      2decomp_fft
      $ENV{CRAY_LIBSCI_PREFIX_DIR}/lib/libsci_gnu_mpi.so
      ${HDF5_LIBRARY_PATH}/libhdf5hl_fortran_parallel_gnu_91.so
      ${HDF5_LIBRARY_PATH}/libhdf5_hl_parallel_gnu_91.so
      ${HDF5_LIBRARY_PATH}/libhdf5_fortran_parallel_gnu_91.so
      ${HDF5_LIBRARY_PATH}/libhdf5_parallel_gnu_91.so
      -lz -ldl -lm -Wl,-rpath -Wl,${HDF5_LIBRARY_PATH}
   )
else()
   # Replace $ENV{CRAY_LIBSCI_PREFIX_DIR}/lib/libsci_gnu.a with local library if needed
   target_link_libraries(PadeOps fftw3 2decomp_fft $ENV{CRAY_LIBSCI_PREFIX_DIR}/lib/libsci_gnu.a  ${HDF5_LIBRARY_PATH}/libhdf5hl_fortran.a ${HDF5_LIBRARY_PATH}/libhdf5_hl.a ${HDF5_LIBRARY_PATH}/libhdf5_fortran.a ${HDF5_LIBRARY_PATH}/libhdf5.a ${MPI_LIBRARIES} -lz -ldl -lm -Wl,-rpath -Wl,${HDF5_LIBRARY_PATH})
endif()
   
if (MPI_Fortran_COMPILER_FLAGS)
    set_target_properties(PadeOps PROPERTIES COMPILE_FLAGS "${CMAKE_Fortran_FLAGS}")
endif()

# Create a library called CompressibleOps
add_library(CompressibleOps STATIC ${grid_source_files} ${compressible_source_files})
target_link_libraries(CompressibleOps PadeOps ${MPI_LIBRARIES})

if (MPI_Fortran_COMPILER_FLAGS)
    set_target_properties(CompressibleOps PROPERTIES COMPILE_FLAGS "${CMAKE_Fortran_FLAGS}")
endif()

# Create a library called SolidOps
add_library(SolidOps STATIC ${grid_source_files} ${solid_source_files})
target_link_libraries(SolidOps CompressibleOps ${MPI_LIBRARIES})

if (MPI_Fortran_COMPILER_FLAGS)
    set_target_properties(SolidOps PROPERTIES COMPILE_FLAGS "${CMAKE_Fortran_FLAGS}")
endif()

# Create a library called IncompressibleOps
add_library(IncompressibleOps STATIC ${grid_source_files} ${poisson_source_files} ${incompressible_source_files})
target_link_libraries(IncompressibleOps PadeOps ${MPI_LIBRARIES})

if (MPI_Fortran_COMPILER_FLAGS)
    set_target_properties(IncompressibleOps PROPERTIES COMPILE_FLAGS "${CMAKE_Fortran_FLAGS}")
endif()
