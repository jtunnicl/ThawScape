# catch
include("${CMAKE_CURRENT_SOURCE_DIR}/catch2/cmake/Catch.cmake")

# test ThawScape
set(TEST_SRC
    tests.cpp
    test_indexx.cpp
    test_streampower.cpp
)
add_executable(ThawScapeTests ${TEST_SRC})
target_link_libraries(ThawScapeTests ThawScapeLib)
catch_discover_tests(ThawScapeTests)

# copy test input files
configure_file(ThawScapeTestInit.ini ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
