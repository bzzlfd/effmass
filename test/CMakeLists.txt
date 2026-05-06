# Test executable for NCPP UPF reader
add_executable(test_ncpp_upf test_ncpp_upf.cpp)
target_link_libraries(test_ncpp_upf PRIVATE pseudo)

add_test(NAME test_ncpp_upf COMMAND test_ncpp_upf WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
