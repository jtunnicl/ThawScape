cmake_minimum_required(VERSION 3.1.0)
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
option(DOUBLE_PRECISION "Use double precision" ON)
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

# OpenMP
option(BUILD_OPENMP "Build with OpenMP support" OFF)
if (BUILD_OPENMP)
    find_package(OpenMP)
    if (NOT OpenMP_CXX_FOUND)
        message(WARNING "Could not find OpenMP")
    endif()
endif()

# list of source files
file(GLOB cpp_src ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)
file(GLOB h_src ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(GLOB hpp_src ${CMAKE_CURRENT_SOURCE_DIR}/*.hpp)
set(ThawScapeSrc ${cpp_src} ${h_src} ${hpp_src})

# add the library
add_library(ThawScapeLib ${ThawScapeSrc})
target_include_directories(ThawScapeLib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
if (OpenMP_CXX_FOUND)
    target_link_libraries(ThawScapeLib OpenMP::OpenMP_CXX)
endif()

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

# add documentation (requires Doxygen and CMake 3.9.0 or later)
if (${CMAKE_VERSION} VERSION_LESS "3.9.0")
    message(STATUS "Building documentation requires CMake >= 3.9.0")
else()
    find_package(Doxygen)
    if (DOXYGEN_FOUND)
        set(DOXYGEN_USE_MDFILE_AS_MAINPAGE README.md)
        set(DOXYGEN_WARN_IF_UNDOCUMENTED NO)
        set(DOXYGEN_BUILTIN_STL_SUPPORT YES)
        set(DOXYGEN_EXTRACT_PRIVATE YES)
        set(DOXYGEN_GENERATE_XML YES)
        set(DOXYGEN_EXTRACT_ALL YES)
        doxygen_add_docs(doc
            ${ThawScapeSrc} README.md
            COMMENT "Build Doxygen documentation"
        )
    else()
        message(STATUS "Building documentation requires Doxygen")
    endif()
endif()
