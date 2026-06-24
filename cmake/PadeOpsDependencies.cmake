include(CheckFortranSourceCompiles)
include(PadeOpsModuleEnv)

function(_padeops_append_existing_library out_var)
    set(_result "${${out_var}}")
    foreach(_lib IN LISTS ARGN)
        if(_lib AND NOT _lib MATCHES "-NOTFOUND$")
            list(APPEND _result "${_lib}")
        endif()
    endforeach()
    set(${out_var} "${_result}" PARENT_SCOPE)
endfunction()

function(_padeops_filter_found_libraries out_var)
    set(_result)
    foreach(_lib IN LISTS ARGN)
        if(_lib AND NOT _lib MATCHES "-NOTFOUND$")
            list(APPEND _result "${_lib}")
        endif()
    endforeach()
    set(${out_var} "${_result}" PARENT_SCOPE)
endfunction()

padeops_collect_dependency_prefixes(HDF5 PADEOPS_HDF5_PREFIXES)
padeops_collect_dependency_prefixes(FFTW PADEOPS_FFTW_PREFIXES)
padeops_collect_dependency_prefixes(2DECOMP PADEOPS_2DECOMP_PREFIXES)
padeops_collect_dependency_prefixes(SZIP PADEOPS_SZIP_PREFIXES)

foreach(_var HDF5_ROOT HDF5_PATH HDF5_HOME FFTW_ROOT FFTW_PATH FFTW_HOME DECOMP_ROOT DECOMP_PATH DECOMP_HOME SZIP_ROOT SZIP_PATH SZIP_HOME LIBSZIP_ROOT LIBSZIP_PATH LIBSZIP_HOME)
    if(DEFINED ${_var} AND NOT "${${_var}}" STREQUAL "")
        list(APPEND CMAKE_PREFIX_PATH "${${_var}}")
    endif()
    if(DEFINED ENV{${_var}} AND NOT "$ENV{${_var}}" STREQUAL "")
        list(APPEND CMAKE_PREFIX_PATH "$ENV{${_var}}")
    endif()
endforeach()
list(APPEND CMAKE_PREFIX_PATH
    ${PADEOPS_HDF5_PREFIXES}
    ${PADEOPS_FFTW_PREFIXES}
    ${PADEOPS_2DECOMP_PREFIXES}
    ${PADEOPS_SZIP_PREFIXES}
)
list(REMOVE_DUPLICATES CMAKE_PREFIX_PATH)

find_package(MPI REQUIRED COMPONENTS Fortran)
find_package(FFTW REQUIRED)
find_package(2DECOMP REQUIRED)

set(_save_required_includes "${CMAKE_REQUIRED_INCLUDES}")
set(CMAKE_REQUIRED_INCLUDES "${2DECOMP_INCLUDE_DIR}")
check_fortran_source_compiles("
      program padeops_2decomp_compile_check
      use decomp_2d
      end
" PADEOPS_2DECOMP_COMPILES)
set(CMAKE_REQUIRED_INCLUDES "${_save_required_includes}")
if(NOT PADEOPS_2DECOMP_COMPILES)
    message(FATAL_ERROR
        "Found 2DECOMP, but its Fortran module cannot be compiled with the active compiler. Ensure DECOMP_ROOT/DECOMP_PATH points to the 2DECOMP build for this compiler and MPI stack."
    )
endif()

include(PadeOpsHDF5)
include(PadeOpsLAPACK)

if(TARGET MPI::MPI_Fortran)
    add_library(PadeOps::MPI ALIAS MPI::MPI_Fortran)
elseif(NOT TARGET PadeOps::MPI)
    add_library(PadeOps::MPI INTERFACE IMPORTED)
    set_target_properties(PadeOps::MPI PROPERTIES
        INTERFACE_COMPILE_OPTIONS "${MPI_Fortran_COMPILE_OPTIONS}"
        INTERFACE_INCLUDE_DIRECTORIES "${MPI_Fortran_INCLUDE_DIRS}"
        INTERFACE_LINK_LIBRARIES "${MPI_Fortran_LIBRARIES}"
        INTERFACE_LINK_OPTIONS "${MPI_Fortran_LINK_FLAGS}"
    )
endif()

message(STATUS "FFTW library: ${FFTW_LIBRARY}")
message(STATUS "2DECOMP library: ${2DECOMP_LIBRARY}")
