file(GLOB tests_source_files *.F90)

# set(tests_source_files "test_cd10.F90;test_cd06.F90")

# Create the test executables
foreach ( testfile ${tests_source_files} )
    string( REPLACE "${PadeOps_SOURCE_DIR}/tests/" "" testfile ${testfile} )
    string( REPLACE ".F90" "" testexec ${testfile} )
  
    add_executable( ${testexec} ${testfile} )
    
    target_link_libraries( ${testexec} PRIVATE IncompressibleOps CompressibleOps )

endforeach()
