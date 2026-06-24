file(GLOB sproblems_source_files *.F90)

# Create the solids problem executables
foreach ( sproblemfile ${sproblems_source_files} )
    string( REPLACE "${PadeOps_SOURCE_DIR}/problems/solids/" "" sproblemfile ${sproblemfile} )
    string( REPLACE ".F90" "" sproblemexec ${sproblemfile} )
  
    add_executable( ${sproblemexec} ${sproblemfile} )
    
    target_link_libraries( ${sproblemexec} PRIVATE SolidOps )

    file(MAKE_DIRECTORY ${PadeOps_BINARY_DIR}/problems/solids/${sproblemexec}_files)
    file(COPY ${sproblemexec}_files/input.dat DESTINATION ${PadeOps_BINARY_DIR}/problems/solids/${sproblemexec}_files)
    
endforeach()
