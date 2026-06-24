include(FindPackageHandleStandardArgs)

set(_FFTW_HINTS)
foreach(_var FFTW_ROOT FFTW_PATH FFTW_HOME)
    if(DEFINED ${_var})
        list(APPEND _FFTW_HINTS "${${_var}}")
    endif()
    if(DEFINED ENV{${_var}})
        list(APPEND _FFTW_HINTS "$ENV{${_var}}")
    endif()
endforeach()
list(REMOVE_DUPLICATES _FFTW_HINTS)

find_path(FFTW_INCLUDE_DIR
    NAMES fftw3.f fftw3.h
    HINTS ${_FFTW_HINTS}
    PATH_SUFFIXES include
)

find_library(FFTW_LIBRARY
    NAMES fftw3
    HINTS ${_FFTW_HINTS}
    PATH_SUFFIXES lib lib64
)

find_package_handle_standard_args(FFTW
    REQUIRED_VARS FFTW_LIBRARY FFTW_INCLUDE_DIR
)

if(FFTW_FOUND AND NOT TARGET PadeOps::FFTW)
    add_library(PadeOps::FFTW INTERFACE IMPORTED)
    set_target_properties(PadeOps::FFTW PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${FFTW_INCLUDE_DIR}"
        INTERFACE_LINK_LIBRARIES "${FFTW_LIBRARY}"
    )
endif()

mark_as_advanced(FFTW_INCLUDE_DIR FFTW_LIBRARY)
