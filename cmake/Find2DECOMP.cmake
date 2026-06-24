include(FindPackageHandleStandardArgs)

set(_2DECOMP_HINTS)
foreach(_var DECOMP_ROOT DECOMP_PATH DECOMP_HOME)
    if(DEFINED ${_var})
        list(APPEND _2DECOMP_HINTS "${${_var}}")
    endif()
    if(DEFINED ENV{${_var}})
        list(APPEND _2DECOMP_HINTS "$ENV{${_var}}")
    endif()
endforeach()
list(REMOVE_DUPLICATES _2DECOMP_HINTS)

find_path(2DECOMP_INCLUDE_DIR
    NAMES decomp_2d.mod decomp_2d_fft.mod
    HINTS ${_2DECOMP_HINTS}
    PATH_SUFFIXES include modules mod
)

find_library(2DECOMP_LIBRARY
    NAMES 2decomp_fft
    HINTS ${_2DECOMP_HINTS}
    PATH_SUFFIXES lib lib64
)

find_package_handle_standard_args(2DECOMP
    REQUIRED_VARS 2DECOMP_LIBRARY 2DECOMP_INCLUDE_DIR
)

if(2DECOMP_FOUND AND NOT TARGET PadeOps::2DECOMP)
    add_library(PadeOps::2DECOMP INTERFACE IMPORTED)
    set_target_properties(PadeOps::2DECOMP PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${2DECOMP_INCLUDE_DIR}"
        INTERFACE_LINK_LIBRARIES "${2DECOMP_LIBRARY}"
    )
endif()

mark_as_advanced(2DECOMP_INCLUDE_DIR 2DECOMP_LIBRARY)
