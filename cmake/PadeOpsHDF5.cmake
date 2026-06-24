function(_padeops_check_hdf5_link result_var link_libs)
    set(_save_required_includes "${CMAKE_REQUIRED_INCLUDES}")
    set(_save_required_libs "${CMAKE_REQUIRED_LIBRARIES}")
    set(CMAKE_REQUIRED_INCLUDES "${PADEOPS_HDF5_INCLUDE_DIRS}")
    set(CMAKE_REQUIRED_LIBRARIES "${link_libs}")
    check_fortran_source_compiles("
      program padeops_hdf5_link_check
      use hdf5
      integer ierr
      call h5open_f(ierr)
      call h5close_f(ierr)
      end
" ${result_var})
    set(CMAKE_REQUIRED_INCLUDES "${_save_required_includes}")
    set(CMAKE_REQUIRED_LIBRARIES "${_save_required_libs}")
    set(${result_var} "${${result_var}}" PARENT_SCOPE)
endfunction()

set(HDF5_PREFER_PARALLEL TRUE)
find_package(HDF5 QUIET COMPONENTS Fortran HL)

if(HDF5_FOUND)
    set(PADEOPS_HDF5_INCLUDE_DIRS ${HDF5_INCLUDE_DIRS})
    _padeops_filter_found_libraries(PADEOPS_HDF5_LINK_LIBRARIES ${HDF5_HL_LIBRARIES} ${HDF5_LIBRARIES})
endif()

if(NOT HDF5_FOUND OR NOT PADEOPS_HDF5_LINK_LIBRARIES)
    set(_HDF5_HINTS ${PADEOPS_HDF5_PREFIXES})
    foreach(_var HDF5_ROOT HDF5_PATH HDF5_HOME)
        if(DEFINED ${_var})
            list(APPEND _HDF5_HINTS "${${_var}}")
        endif()
        if(DEFINED ENV{${_var}})
            list(APPEND _HDF5_HINTS "$ENV{${_var}}")
        endif()
    endforeach()
    list(REMOVE_DUPLICATES _HDF5_HINTS)

    find_path(PADEOPS_HDF5_INCLUDE_DIRS
        NAMES hdf5.mod hdf5.h
        HINTS ${_HDF5_HINTS}
        PATH_SUFFIXES include
    )
    foreach(_name hdf5_hl_fortran hdf5hl_fortran hdf5_hl_f90cstub hdf5_fortran hdf5_f90cstub hdf5_hl hdf5)
        find_library(_PADEOPS_HDF5_${_name}_LIBRARY
            NAMES ${_name}
            HINTS ${_HDF5_HINTS}
            PATH_SUFFIXES lib lib64
        )
        _padeops_append_existing_library(PADEOPS_HDF5_LINK_LIBRARIES "${_PADEOPS_HDF5_${_name}_LIBRARY}")
    endforeach()

    if(NOT PADEOPS_HDF5_INCLUDE_DIRS OR NOT PADEOPS_HDF5_LINK_LIBRARIES)
        message(FATAL_ERROR
            "Could not find HDF5 Fortran/HL. Load an HDF5 module or pass -DHDF5_ROOT=/path/to/hdf5."
        )
    endif()
endif()

list(REMOVE_DUPLICATES PADEOPS_HDF5_LINK_LIBRARIES)

_padeops_check_hdf5_link(PADEOPS_HDF5_LINKS "${PADEOPS_HDF5_LINK_LIBRARIES}")

if(NOT PADEOPS_HDF5_LINKS)
    set(_PADEOPS_HDF5_SYSTEM_EXTRAS)
    find_package(ZLIB QUIET)
    if(ZLIB_FOUND)
        list(APPEND _PADEOPS_HDF5_SYSTEM_EXTRAS ZLIB::ZLIB)
    else()
        find_library(PADEOPS_ZLIB_LIBRARY NAMES z)
        _padeops_append_existing_library(_PADEOPS_HDF5_SYSTEM_EXTRAS "${PADEOPS_ZLIB_LIBRARY}")
    endif()

    find_library(PADEOPS_DL_LIBRARY NAMES dl)
    find_library(PADEOPS_M_LIBRARY NAMES m)
    _padeops_append_existing_library(_PADEOPS_HDF5_SYSTEM_EXTRAS "${PADEOPS_DL_LIBRARY}" "${PADEOPS_M_LIBRARY}")

    set(_PADEOPS_HDF5_WITH_SYSTEM ${PADEOPS_HDF5_LINK_LIBRARIES} ${_PADEOPS_HDF5_SYSTEM_EXTRAS})
    _padeops_check_hdf5_link(PADEOPS_HDF5_LINKS_WITH_SYSTEM "${_PADEOPS_HDF5_WITH_SYSTEM}")

    if(PADEOPS_HDF5_LINKS_WITH_SYSTEM)
        set(PADEOPS_HDF5_LINK_LIBRARIES ${_PADEOPS_HDF5_WITH_SYSTEM})
    else()
        find_package(SZIP QUIET)
        if(SZIP_FOUND)
            set(_PADEOPS_HDF5_WITH_SZIP ${PADEOPS_HDF5_LINK_LIBRARIES} PadeOps::SZIP ${_PADEOPS_HDF5_SYSTEM_EXTRAS})
            _padeops_check_hdf5_link(PADEOPS_HDF5_LINKS_WITH_SZIP "${_PADEOPS_HDF5_WITH_SZIP}")
            if(PADEOPS_HDF5_LINKS_WITH_SZIP)
                set(PADEOPS_HDF5_LINK_LIBRARIES ${_PADEOPS_HDF5_WITH_SZIP})
                set(PADEOPS_HDF5_REQUIRES_SZIP TRUE)
            endif()
        endif()
    endif()
endif()

if(NOT PADEOPS_HDF5_LINKS AND NOT PADEOPS_HDF5_LINKS_WITH_SYSTEM AND NOT PADEOPS_HDF5_LINKS_WITH_SZIP)
    message(FATAL_ERROR
        "Found HDF5, but a Fortran HDF5 compile/link test failed. Ensure HDF5 was built with the active Fortran compiler and MPI stack. If this HDF5 was built with SZIP, load/provide SZIP with -DSZIP_ROOT=/path/to/szip or LIBSZIP_HOME."
    )
endif()

if(NOT TARGET PadeOps::HDF5)
    add_library(PadeOps::HDF5 INTERFACE IMPORTED)
    set_target_properties(PadeOps::HDF5 PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${PADEOPS_HDF5_INCLUDE_DIRS}"
        INTERFACE_LINK_LIBRARIES "${PADEOPS_HDF5_LINK_LIBRARIES}"
    )
endif()

message(STATUS "HDF5 include dirs: ${PADEOPS_HDF5_INCLUDE_DIRS}")
message(STATUS "HDF5 link libraries: ${PADEOPS_HDF5_LINK_LIBRARIES}")
if(PADEOPS_HDF5_REQUIRES_SZIP)
    message(STATUS "HDF5 requires SZIP: ${SZIP_LIBRARY}")
endif()
