file(GLOB postprocessing_source_files *.F90)

# Include directories
include_directories( ${FFTW_INCLUDE_PATH} ${DECOMP_INCLUDE_PATH} ${PadeOps_BINARY_DIR}/src  )

# Link directories
link_directories( ${FFTW_LIBRARY_PATH} ${DECOMP_LIBRARY_PATH} ${PadeOps_BINARY_DIR}/src )

# Create the compressible problem executables
foreach ( postprocessingfile ${postprocessing_source_files} )
    string( REPLACE "${PadeOps_SOURCE_DIR}/problems/postprocessing/" "" postprocessingfile ${postprocessingfile} )
    string( REPLACE ".F90" "" postprocessingexec ${postprocessingfile} )
  
    add_executable( ${postprocessingexec} ${postprocessingfile} )
    
    target_link_libraries( ${postprocessingexec} PadeOps )

    if (MPI_Fortran_COMPILER_FLAGS)
        set_target_properties(${postprocessingexec} PROPERTIES COMPILE_FLAGS "${CMAKE_Fortran_FLAGS}")
    endif()

    file(MAKE_DIRECTORY ${PadeOps_BINARY_DIR}/problems/compressible/${postprocessingexec}_files)
    file(COPY ${postprocessingexec}_files/input.dat DESTINATION ${PadeOps_BINARY_DIR}/problems/postprocessing/${postprocessingexec}_files) 
endforeach()