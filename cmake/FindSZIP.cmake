include(FindPackageHandleStandardArgs)

set(_SZIP_HINTS)
foreach(_var SZIP_ROOT SZIP_PATH SZIP_HOME LIBSZIP_ROOT LIBSZIP_PATH LIBSZIP_HOME)
    if(DEFINED ${_var})
        list(APPEND _SZIP_HINTS "${${_var}}")
    endif()
    if(DEFINED ENV{${_var}})
        list(APPEND _SZIP_HINTS "$ENV{${_var}}")
    endif()
endforeach()
list(REMOVE_DUPLICATES _SZIP_HINTS)

find_library(SZIP_LIBRARY
    NAMES sz szip aec
    HINTS ${_SZIP_HINTS}
    PATH_SUFFIXES lib lib64
)

find_package_handle_standard_args(SZIP
    REQUIRED_VARS SZIP_LIBRARY
)

if(SZIP_FOUND AND NOT TARGET PadeOps::SZIP)
    add_library(PadeOps::SZIP UNKNOWN IMPORTED)
    set_target_properties(PadeOps::SZIP PROPERTIES
        IMPORTED_LOCATION "${SZIP_LIBRARY}"
    )
endif()

mark_as_advanced(SZIP_LIBRARY)
