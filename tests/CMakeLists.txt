# catch
include("${CMAKE_CURRENT_SOURCE_DIR}/catch2/cmake/Catch.cmake")

# test ThawScape
set(TEST_SRC
    tests.cpp
    test_indexx.cpp
)
add_executable(ThawScapeTests ${TEST_SRC})
target_link_libraries(ThawScapeTests ThawScapeLib)
target_include_directories(ThawScapeTests PRIVATE ${PROJECT_SOURCE_DIR})
catch_discover_tests(ThawScapeTests)
