cmake_minimum_required(VERSION 3.9.0)
project(ThawScape
        DESCRIPTION "ThawScape is an experimental code to simulate the influence of thermokarst dynamics upon ice-cored terrain"
        LANGUAGES CXX)

# CXX standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# enable testing (disable with -DBUILD_TESTING=OFF)
include(CTest)

# default to RELEASE build if not specified
if (NOT CMAKE_BUILD_TYPE OR CMAKE_BUILD_TYPE STREQUAL "")
    message(STATUS "Defaulting to RELWITHDEBINFO build")
    set(CMAKE_BUILD_TYPE "RELWITHDEBINFO" CACHE STRING "" FORCE)
endif()

# option to use double precision (default is single)
option(DOUBLE_PRECISION "Use double precision" OFF)
if (DOUBLE_PRECISION)
    add_definitions(-DDOUBLE_PRECISION)
    message(STATUS "Double precision build")
else()
    message(STATUS "Single precision build")
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

# list of source files
set(ThawScapeSrc
    global_defs.h
    timer.hpp
    raster.h
    raster.cpp
    data_structures.h
    model_time.h
    model_time.cpp
    streampower.h
    streampower.cpp
    utility.h
    utility.cpp
    priority_flood.hpp
    Array2D.hpp
)

# add the library
add_library(ThawScapeLib ${ThawScapeSrc})
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

# add documentation directory
option(BUILD_DOC "Build the documentation" ON)
if (BUILD_DOC)
    find_package(Doxygen REQUIRED dot)
    if (DOXYGEN_FOUND)
        set(DOXYGEN_USE_MDFILE_AS_MAINPAGE README.md)
        set(DOXYGEN_WARN_IF_UNDOCUMENTED NO)
        set(DOXYGEN_BUILTIN_STL_SUPPORT YES)
        set(DOXYGEN_EXTRACT_PRIVATE YES)
        doxygen_add_docs(doc
            ${ThawScapeSrc} README.md
            COMMENT "Build Doxygen documentation"
        )
    else()
        message(WARNING "Doxygen must be installed to build documentation")
    endif()
endif()
