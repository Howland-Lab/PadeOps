option(PADEOPS_REQUIRE_LAPACK
    "Require CMake to find a LAPACK provider at configure time"
    OFF
)

function(_padeops_check_lapack_link result_var link_libs)
    set(_save_required_libs "${CMAKE_REQUIRED_LIBRARIES}")
    set(CMAKE_REQUIRED_LIBRARIES "${link_libs}")
    check_fortran_source_compiles("
      program padeops_lapack_link_check
      double precision a(1,1), b(1,1)
      integer ipiv(1), info
      a(1,1) = 1.0d0
      b(1,1) = 1.0d0
      call dgesv(1, 1, a, 1, ipiv, b, 1, info)
      end
" ${result_var})
    set(CMAKE_REQUIRED_LIBRARIES "${_save_required_libs}")
    set(${result_var} "${${result_var}}" PARENT_SCOPE)
endfunction()

find_package(LAPACK QUIET)

if(LAPACK_FOUND)
    add_library(PadeOps::LAPACK INTERFACE IMPORTED)
    set_target_properties(PadeOps::LAPACK PROPERTIES
        INTERFACE_LINK_LIBRARIES "${LAPACK_LIBRARIES}"
    )
    message(STATUS "LAPACK provider: ${LAPACK_LIBRARIES}")
else()
    _padeops_check_lapack_link(PADEOPS_IMPLICIT_LAPACK_WORKS "")
    if(PADEOPS_IMPLICIT_LAPACK_WORKS)
        add_library(PadeOps::LAPACK INTERFACE IMPORTED)
        message(STATUS "LAPACK provider: active Fortran compiler wrapper/linker defaults")
    elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
        _padeops_check_lapack_link(PADEOPS_INTEL_MKL_FLAG_WORKS "-mkl")
        if(PADEOPS_INTEL_MKL_FLAG_WORKS)
            add_library(PadeOps::LAPACK INTERFACE IMPORTED)
            set_target_properties(PadeOps::LAPACK PROPERTIES
                INTERFACE_LINK_LIBRARIES "-mkl"
            )
            message(STATUS "LAPACK provider: Intel compiler -mkl fallback")
        endif()
    endif()
endif()

if(NOT TARGET PadeOps::LAPACK)
    if(PADEOPS_REQUIRE_LAPACK)
        message(FATAL_ERROR
            "Could not find a working LAPACK provider. Load a LAPACK/MKL module or provide a LAPACK-capable compiler stack."
        )
    endif()

    add_library(PadeOps::LAPACK INTERFACE IMPORTED)
    message(STATUS
        "No LAPACK provider was found. Continuing because PADEOPS_REQUIRE_LAPACK=OFF; LAPACK-using targets will fail to link unless their stack supplies LAPACK implicitly."
    )
endif()
