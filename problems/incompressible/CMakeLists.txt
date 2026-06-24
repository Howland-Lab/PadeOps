file(GLOB iproblems_source_files *.F90)

# Create the incompressible problem executables
foreach ( iproblemfile ${iproblems_source_files} )
    string( REPLACE "${PadeOps_SOURCE_DIR}/problems/incompressible/" "" iproblemfile ${iproblemfile} )
    string( REPLACE ".F90" "" iproblemexec ${iproblemfile} )
  
    add_executable( ${iproblemexec} ${iproblemfile} )
    
    target_link_libraries( ${iproblemexec} PRIVATE IncompressibleOps )
    
    file(MAKE_DIRECTORY ${PadeOps_BINARY_DIR}/problems/incompressible/${iproblemexec}_files)
    #    file(COPY ${iproblemexec}_files/*.dat DESTINATION ${PadeOps_BINARY_DIR}/problems/incompressible/${iproblemexec}_files) 
    file(GLOB iproblems_input_files ${iproblemexec}_files/*.dat)
    file(COPY ${iproblems_input_files} DESTINATION ${PadeOps_BINARY_DIR}/problems/incompressible/${iproblemexec}_files) 
endforeach()
