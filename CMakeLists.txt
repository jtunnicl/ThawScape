cmake_minimum_required(VERSION 3.1.0)
project(ThawScape LANGUAGES CXX)

# CXX standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# enable testing (disable with -DBUILD_TESTING=OFF)
#include(CTest)

# default to RELEASE build if not specified
if (NOT CMAKE_BUILD_TYPE OR CMAKE_BUILD_TYPE STREQUAL "")
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
endif()

# option to gather profiling data
option(ENABLE_PROFILING "Enable gprof profiling" OFF)
if (ENABLE_PROFILING)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")

    # disable some optimisations as recommended here: https://github.com/jrfonseca/gprof2dot
    if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-omit-frame-pointer -fno-inline-functions -fno-inline-functions-called-once -fno-optimize-sibling-calls")
    endif()
endif()

# Find includes in corresponding build directories
#set(CMAKE_INCLUDE_CURRENT_DIR ON)

# add the executable
add_executable(ThawScape
    main.cpp
    data_structures.h
    streampower.h
    streampower.cpp
    utility.h
    utility.cpp
    priority_flood.hpp
    Array2D.hpp
)

# copy input files to build dir for testing
configure_file(FA.asc ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(topo.asc ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
