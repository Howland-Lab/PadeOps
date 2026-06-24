file(GLOB postprocessing_source_files *.F90)

# Create the compressible problem executables
foreach ( postprocessingfile ${postprocessing_source_files} )
    string( REPLACE "${PadeOps_SOURCE_DIR}/problems/postprocessing/" "" postprocessingfile ${postprocessingfile} )
    string( REPLACE ".F90" "" postprocessingexec ${postprocessingfile} )
  
    add_executable( ${postprocessingexec} ${postprocessingfile} )
    
    target_link_libraries( ${postprocessingexec} PRIVATE PadeOps)

    file(MAKE_DIRECTORY ${PadeOps_BINARY_DIR}/problems/postprocessing/${postprocessingexec}_files)
    file(COPY ${postprocessingexec}_files/input.dat DESTINATION ${PadeOps_BINARY_DIR}/problems/postprocessing/${postprocessingexec}_files) 
endforeach()
