add_library(PadeOps::FortranOptions INTERFACE IMPORTED)

find_package(OpenMP QUIET COMPONENTS Fortran)

if(DEFINED ENV{ARCH_OPT_FLAG} AND NOT "$ENV{ARCH_OPT_FLAG}" STREQUAL "")
    separate_arguments(PADEOPS_ARCH_OPT_FLAGS NATIVE_COMMAND "$ENV{ARCH_OPT_FLAG}")
else()
    set(PADEOPS_ARCH_OPT_FLAGS)
endif()

set(_PADEOPS_FORTRAN_COMPILE_OPTIONS)
set(_PADEOPS_FORTRAN_LINK_OPTIONS)

if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    if(NOT PADEOPS_ARCH_OPT_FLAGS)
        set(PADEOPS_ARCH_OPT_FLAGS -xhost)
    endif()

    if(CMAKE_BUILD_TYPE MATCHES "Debug")
        list(APPEND _PADEOPS_FORTRAN_COMPILE_OPTIONS
            -g -traceback -heap-arrays 1024 -check all,nouninit
            -check noarg-temp-created -fpe0 -warn -debug extended
            -assume realloc_lhs -fstack-protector -assume protect_parens
            -implicitnone
        )
    else()
        list(APPEND _PADEOPS_FORTRAN_COMPILE_OPTIONS
            -O3 -traceback -heap-arrays 1024 -warn all
            ${PADEOPS_ARCH_OPT_FLAGS}
            -dynamic -qopt-report=2 -qopt-report-phase=vec
        )
    endif()
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    if(NOT PADEOPS_ARCH_OPT_FLAGS)
        set(PADEOPS_ARCH_OPT_FLAGS -march=native)
    endif()

    if(CMAKE_BUILD_TYPE MATCHES "Debug")
        list(APPEND _PADEOPS_FORTRAN_COMPILE_OPTIONS
            -Og -g -fbacktrace -pg -ffree-form -ffree-line-length-none
            -fcheck=all -fbounds-check -ffpe-trap=zero,overflow
            -Wall -Wconversion -Wextra -Waliasing -Wsurprising
        )
    else()
        list(APPEND _PADEOPS_FORTRAN_COMPILE_OPTIONS
            -O3 -Wall -Wconversion -Wextra -Waliasing
            -ffree-form -ffree-line-length-none
            ${PADEOPS_ARCH_OPT_FLAGS}
            -funroll-loops -fallow-argument-mismatch
            -finit-integer=0 -finit-real=zero
        )
    endif()
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Cray")
    if(NOT CMAKE_BUILD_TYPE MATCHES "Debug")
        list(APPEND _PADEOPS_FORTRAN_COMPILE_OPTIONS -hlist=a)
    endif()
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "IBM")
    if(CMAKE_BUILD_TYPE MATCHES "Debug")
        list(APPEND _PADEOPS_FORTRAN_COMPILE_OPTIONS -O0 -qsuffix=cpp=f90 -qxlf2003=polymorphic)
    else()
        list(APPEND _PADEOPS_FORTRAN_COMPILE_OPTIONS -O5 -qsuffix=cpp=f90 -qxlf2003=polymorphic)
    endif()
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "AOCC|AMD")
    if(CMAKE_BUILD_TYPE MATCHES "Debug")
        list(APPEND _PADEOPS_FORTRAN_COMPILE_OPTIONS
            -O0 -g -ffree-form -fcheck=all
            -fbounds-check -ffpe-trap=zero,overflow
        )
    else()
        list(APPEND _PADEOPS_FORTRAN_COMPILE_OPTIONS
            -O3 ${PADEOPS_ARCH_OPT_FLAGS}
            -ffree-form -ffast-math -funroll-loops
        )
    endif()
endif()

if(OpenMP_Fortran_FOUND)
    list(APPEND _PADEOPS_FORTRAN_LINK_OPTIONS OpenMP::OpenMP_Fortran)
    message(STATUS "OpenMP Fortran provider: OpenMP::OpenMP_Fortran")
else()
    message(STATUS "OpenMP Fortran was not found by CMake; building without CMake-managed OpenMP flags")
endif()

if(_PADEOPS_FORTRAN_COMPILE_OPTIONS)
    set_target_properties(PadeOps::FortranOptions PROPERTIES
        INTERFACE_COMPILE_OPTIONS "${_PADEOPS_FORTRAN_COMPILE_OPTIONS}"
    )
endif()

if(_PADEOPS_FORTRAN_LINK_OPTIONS)
    set_target_properties(PadeOps::FortranOptions PROPERTIES
        INTERFACE_LINK_LIBRARIES "${_PADEOPS_FORTRAN_LINK_OPTIONS}"
    )
endif()

message(STATUS "Fortran compiler: ${CMAKE_Fortran_COMPILER_ID} ${CMAKE_Fortran_COMPILER_VERSION}")
