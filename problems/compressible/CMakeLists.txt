file(GLOB cproblems_source_files *.F90)

# Create the compressible problem executables
foreach ( cproblemfile ${cproblems_source_files} )
    string( REPLACE "${PadeOps_SOURCE_DIR}/problems/compressible/" "" cproblemfile ${cproblemfile} )
    string( REPLACE ".F90" "" cproblemexec ${cproblemfile} )
  
    add_executable( ${cproblemexec} ${cproblemfile} )
    
    target_link_libraries( ${cproblemexec} PRIVATE CompressibleOps )

    file(MAKE_DIRECTORY ${PadeOps_BINARY_DIR}/problems/compressible/${cproblemexec}_files)
    file(COPY ${cproblemexec}_files/input.dat DESTINATION ${PadeOps_BINARY_DIR}/problems/compressible/${cproblemexec}_files) 
endforeach()
