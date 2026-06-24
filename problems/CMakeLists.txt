file(GLOB problems_source_files *.F90)
file(GLOB cproblems_source_files compressible/*.F90)
file(GLOB iproblems_source_files incompressible/*.F90)
file(GLOB turbproblems_source_files turbines/*.F90)
file(GLOB pp_source_files postprocessing/*.F90)
file(GLOB pp_igrid_source_files postprocessing_igrid/*.F90)

# Create the problem executables
foreach ( problemfile ${problems_source_files} )
    string( REPLACE "${PadeOps_SOURCE_DIR}/problems/" "" problemfile ${problemfile} )
    string( REPLACE ".F90" "" problemexec ${problemfile} )
  
    add_executable( ${problemexec} ${problemfile} )
    
    target_link_libraries( ${problemexec} PRIVATE PadeOps )
    
endforeach()

add_subdirectory(compressible)
add_subdirectory(incompressible)
add_subdirectory(turbines)
add_subdirectory(postprocessing)
add_subdirectory(postprocessing_igrid)
add_subdirectory(solids)
