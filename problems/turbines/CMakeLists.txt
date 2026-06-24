file(GLOB turbproblems_source_files *.F90)

# Create the turbines problem executables
foreach ( turbproblemfile ${turbproblems_source_files} )
    string( REPLACE "${PadeOps_SOURCE_DIR}/problems/turbines/" "" turbproblemfile ${turbproblemfile} )
    string( REPLACE ".F90" "" turbproblemexec ${turbproblemfile} )
  
    add_executable( ${turbproblemexec} ${turbproblemfile} )
    
    target_link_libraries( ${turbproblemexec} PRIVATE IncompressibleOps )
    
    file(MAKE_DIRECTORY ${PadeOps_BINARY_DIR}/problems/turbines/${turbproblemexec}_files)
    #    file(COPY ${turbproblemexec}_files/*.dat DESTINATION ${PadeOps_BINARY_DIR}/problems/turbines/${turbproblemexec}_files) 
    file(GLOB turbproblems_input_files ${turbproblemexec}_files/*.dat)
    file(COPY ${turbproblems_input_files} DESTINATION ${PadeOps_BINARY_DIR}/problems/turbines/${turbproblemexec}_files) 
endforeach()
