file(GLOB iproblems_source_files *.F90)

# Include directories
include_directories( ${FFTW_INCLUDE_PATH} ${DECOMP_INCLUDE_PATH} ${PadeOps_BINARY_DIR}/src  )

# Link directories
link_directories( ${FFTW_LIBRARY_PATH} ${DECOMP_LIBRARY_PATH} ${PadeOps_BINARY_DIR}/src )

# Create the incompressible problem executables
foreach ( iproblemfile ${iproblems_source_files} )
    string( REPLACE "${PadeOps_SOURCE_DIR}/problems/incompressible/" "" iproblemfile ${iproblemfile} )
    string( REPLACE ".F90" "" iproblemexec ${iproblemfile} )
  
    add_executable( ${iproblemexec} ${iproblemfile} )
    
    target_link_libraries( ${iproblemexec} IncompressibleOps )

    if (MPI_Fortran_COMPILER_FLAGS)
        set_target_properties(${iproblemexec} PROPERTIES COMPILE_FLAGS "${CMAKE_Fortran_FLAGS}")
    endif()
    
endforeach()