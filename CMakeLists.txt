cmake_minimum_required(VERSION 3.1.0)
project(ThawScape LANGUAGES CXX)

# CXX standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# enable testing (disable with -DBUILD_TESTING=OFF)
include(CTest)

# default to RELEASE build if not specified
if (NOT CMAKE_BUILD_TYPE OR CMAKE_BUILD_TYPE STREQUAL "")
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
endif()

# option to use double precision (default is single)
option(DOUBLE_PRECISION "Use double precision" OFF)
if (DOUBLE_PRECISION)
    add_definitions(-DDOUBLE_PRECISION)
endif()

# enable compiler warnings
if (${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU" OR
        ${CMAKE_CXX_COMPILER_ID} STREQUAL "AppleClang" OR
        ${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wpedantic")
elseif (${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
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

# add the library
add_library(ThawScapeLib
    global_defs.h
    timer.hpp
    indexx.hpp
    data_structures.h
    time_fcn.h
    time_fcn.cpp
    streampower.h
    streampower.cpp
    utility.h
    utility.cpp
    priority_flood.hpp
    Array2D.hpp
)
target_include_directories(ThawScapeLib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# add executable
add_executable(ThawScape main.cpp)
target_link_libraries(ThawScape ThawScapeLib)

# copy input files to build dir for testing
configure_file(FA.asc ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(topo.asc ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(ThawScape.ini ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)

# add tests directory
if (BUILD_TESTING)
    add_subdirectory(tests)
endif()
