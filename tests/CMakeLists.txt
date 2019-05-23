# test indexx.hpp
add_executable(test_indexx test_indexx.cpp)
target_link_libraries(test_indexx ThawScapeLib)
add_test(NAME test_indexx COMMAND test_indexx)
