cmake_minimum_required(VERSION 3.16)
project(PadeOps LANGUAGES Fortran C)

list(PREPEND CMAKE_MODULE_PATH "${PadeOps_SOURCE_DIR}/cmake")

if (NOT CMAKE_BUILD_TYPE)
    message(STATUS "No build type selected, default to Release")
    set(CMAKE_BUILD_TYPE "Release")
endif()

include(PadeOpsCompiler)
include(PadeOpsDependencies)

# Add source file directories
add_subdirectory(src)
add_subdirectory(tests)
add_subdirectory(problems)
