file(GLOB utilities_source_files utilities/*.F90)
file(GLOB derivatives_source_files derivatives/*.F90)
file(GLOB filters_source_files filters/*.F90)
file(GLOB io_source_files io/*.F90)
file(GLOB grid_source_files grid/*.F90)
file(GLOB compressible_source_files compressible/*.F90)
file(GLOB solid_source_files solid/*.F90)
file(GLOB poisson_source_files incompressible/poisson/*.F90)
file(GLOB incompressible_source_files incompressible/*.F90)

set(PADEOPS_MODULE_DIR "${PadeOps_BINARY_DIR}/src")

# Create a library called PadeOps
add_library(PadeOps STATIC ${utilities_source_files} ${derivatives_source_files} ${filters_source_files} ${io_source_files} )
set_target_properties(PadeOps PROPERTIES Fortran_MODULE_DIRECTORY "${PADEOPS_MODULE_DIR}")
target_include_directories(PadeOps PUBLIC "${PADEOPS_MODULE_DIR}")
target_link_libraries(PadeOps
   PUBLIC
      PadeOps::FortranOptions
      PadeOps::MPI
      PadeOps::FFTW
      PadeOps::2DECOMP
      PadeOps::HDF5
   PRIVATE
      PadeOps::LAPACK
)

# Create a library called CompressibleOps
add_library(CompressibleOps STATIC ${grid_source_files} ${compressible_source_files})
set_target_properties(CompressibleOps PROPERTIES Fortran_MODULE_DIRECTORY "${PADEOPS_MODULE_DIR}")
target_include_directories(CompressibleOps PUBLIC "${PADEOPS_MODULE_DIR}")
target_link_libraries(CompressibleOps PUBLIC PadeOps)

# Create a library called SolidOps
add_library(SolidOps STATIC ${grid_source_files} ${solid_source_files})
set_target_properties(SolidOps PROPERTIES Fortran_MODULE_DIRECTORY "${PADEOPS_MODULE_DIR}")
target_include_directories(SolidOps PUBLIC "${PADEOPS_MODULE_DIR}")
target_link_libraries(SolidOps PUBLIC CompressibleOps)

# Create a library called IncompressibleOps
add_library(IncompressibleOps STATIC ${grid_source_files} ${poisson_source_files} ${incompressible_source_files})
set_target_properties(IncompressibleOps PROPERTIES Fortran_MODULE_DIRECTORY "${PADEOPS_MODULE_DIR}")
target_include_directories(IncompressibleOps PUBLIC "${PADEOPS_MODULE_DIR}")
target_link_libraries(IncompressibleOps PUBLIC PadeOps)
