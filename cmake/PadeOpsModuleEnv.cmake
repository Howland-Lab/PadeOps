set(PADEOPS_MODULES "$ENV{PADEOPS_MODULES}" CACHE STRING
    "Whitespace-separated module names loaded by the setup script for this PadeOps stack"
)

function(_padeops_normalize_prefix out_var candidate)
    set(_result)

    if(candidate AND EXISTS "${candidate}")
        if(IS_DIRECTORY "${candidate}")
            get_filename_component(_name "${candidate}" NAME)
            if(_name MATCHES "^(bin|include|lib|lib64|mod|modules)$")
                get_filename_component(_result "${candidate}" DIRECTORY)
            else()
                set(_result "${candidate}")
            endif()
        else()
            get_filename_component(_dir "${candidate}" DIRECTORY)
            get_filename_component(_name "${_dir}" NAME)
            if(_name MATCHES "^(bin|include|lib|lib64|mod|modules)$")
                get_filename_component(_result "${_dir}" DIRECTORY)
            else()
                set(_result "${_dir}")
            endif()
        endif()
    endif()

    set(${out_var} "${_result}" PARENT_SCOPE)
endfunction()

function(_padeops_append_path_values out_var raw_value)
    set(_result "${${out_var}}")

    if(raw_value)
        string(REPLACE ":" ";" _values "${raw_value}")
        foreach(_value IN LISTS _values)
            if(EXISTS "${_value}")
                list(APPEND _result "${_value}")
            endif()
        endforeach()
    endif()

    set(${out_var} "${_result}" PARENT_SCOPE)
endfunction()

function(_padeops_dependency_patterns dependency name_pattern_var)
    string(TOUPPER "${dependency}" _dep)

    if(_dep STREQUAL "HDF5")
        set(_pattern "HDF5|PHDF5")
    elseif(_dep STREQUAL "FFTW")
        set(_pattern "FFTW")
    elseif(_dep STREQUAL "SZIP")
        set(_pattern "SZIP|LIBSZIP|LIBAEC|AEC")
    elseif(_dep STREQUAL "2DECOMP")
        set(_pattern "2DECOMP|DECOMP")
    else()
        set(_pattern "${_dep}")
    endif()

    set(${name_pattern_var} "${_pattern}" PARENT_SCOPE)
endfunction()

function(padeops_collect_env_candidates dependency out_var)
    _padeops_dependency_patterns("${dependency}" _name_pattern)
    set(_candidates)

    execute_process(
        COMMAND "${CMAKE_COMMAND}" -E environment
        OUTPUT_VARIABLE _environment
        ERROR_QUIET
    )
    string(REPLACE "\n" ";" _environment_lines "${_environment}")

    foreach(_line IN LISTS _environment_lines)
        if(_line MATCHES "^([^=]+)=(.*)$")
            set(_name "${CMAKE_MATCH_1}")
            set(_value "${CMAKE_MATCH_2}")
            string(TOUPPER "${_name}" _name_upper)

            if(_name_upper MATCHES "${_name_pattern}")
                _padeops_append_path_values(_candidates "${_value}")
            endif()
        endif()
    endforeach()

    set(${out_var} "${_candidates}" PARENT_SCOPE)
endfunction()

function(padeops_collect_module_candidates dependency out_var)
    _padeops_dependency_patterns("${dependency}" _name_pattern)
    set(_candidates)

    if(NOT PADEOPS_MODULES)
        set(${out_var} "" PARENT_SCOPE)
        return()
    endif()

    separate_arguments(_modules NATIVE_COMMAND "${PADEOPS_MODULES}")

    foreach(_module IN LISTS _modules)
        string(TOUPPER "${_module}" _module_upper)
        if(NOT _module_upper MATCHES "${_name_pattern}")
            continue()
        endif()

        execute_process(
            COMMAND bash -lc "module show ${_module} 2>&1"
            OUTPUT_VARIABLE _module_show
            ERROR_VARIABLE _module_show_error
            RESULT_VARIABLE _module_show_result
        )

        set(_module_text "${_module_show}\n${_module_show_error}")
        string(REGEX MATCHALL "(/[^ \t\r\n;:\"'()]+)" _paths "${_module_text}")

        foreach(_path IN LISTS _paths)
            if(EXISTS "${_path}")
                list(APPEND _candidates "${_path}")
            endif()
        endforeach()
    endforeach()

    set(${out_var} "${_candidates}" PARENT_SCOPE)
endfunction()

function(padeops_validate_prefixes dependency out_var)
    string(TOUPPER "${dependency}" _dep)
    set(_valid_prefixes)

    foreach(_candidate IN LISTS ARGN)
        _padeops_normalize_prefix(_prefix "${_candidate}")
        if(NOT _prefix)
            continue()
        endif()

        set(_is_valid FALSE)
        if(_dep STREQUAL "HDF5")
            file(GLOB _libs "${_prefix}/lib/libhdf5*" "${_prefix}/lib64/libhdf5*")
            if(_libs AND (EXISTS "${_prefix}/include/hdf5.mod" OR EXISTS "${_prefix}/include/hdf5.h"))
                set(_is_valid TRUE)
            endif()
        elseif(_dep STREQUAL "FFTW")
            file(GLOB _libs "${_prefix}/lib/libfftw3*" "${_prefix}/lib64/libfftw3*")
            if(_libs AND (EXISTS "${_prefix}/include/fftw3.f" OR EXISTS "${_prefix}/include/fftw3.h"))
                set(_is_valid TRUE)
            endif()
        elseif(_dep STREQUAL "SZIP")
            file(GLOB _libs
                "${_prefix}/lib/libsz*"
                "${_prefix}/lib64/libsz*"
                "${_prefix}/lib/libaec*"
                "${_prefix}/lib64/libaec*"
            )
            if(_libs)
                set(_is_valid TRUE)
            endif()
        elseif(_dep STREQUAL "2DECOMP")
            file(GLOB _libs "${_prefix}/lib/lib2decomp_fft*" "${_prefix}/lib64/lib2decomp_fft*")
            if(_libs AND (EXISTS "${_prefix}/include/decomp_2d.mod" OR EXISTS "${_prefix}/modules/decomp_2d.mod" OR EXISTS "${_prefix}/mod/decomp_2d.mod"))
                set(_is_valid TRUE)
            endif()
        endif()

        if(_is_valid)
            list(APPEND _valid_prefixes "${_prefix}")
        endif()
    endforeach()

    if(_valid_prefixes)
        list(REMOVE_DUPLICATES _valid_prefixes)
    endif()

    set(${out_var} "${_valid_prefixes}" PARENT_SCOPE)
endfunction()

function(padeops_collect_dependency_prefixes dependency out_var)
    padeops_collect_env_candidates("${dependency}" _env_candidates)
    padeops_collect_module_candidates("${dependency}" _module_candidates)
    padeops_validate_prefixes("${dependency}" _prefixes ${_env_candidates} ${_module_candidates})

    if(_prefixes)
        message(STATUS "${dependency} module/env prefix candidates: ${_prefixes}")
    endif()

    set(${out_var} "${_prefixes}" PARENT_SCOPE)
endfunction()
