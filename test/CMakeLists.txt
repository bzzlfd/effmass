# IO reader tests
add_executable(test_io_ncpp test_io_ncpp.cpp)
target_link_libraries(test_io_ncpp PRIVATE pseudo)
add_test(NAME test_io_ncpp COMMAND test_io_ncpp WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})

# Potential operator tests
add_executable(test_ncpp test_ncpp.cpp)
target_link_libraries(test_ncpp PRIVATE pseudo)
add_test(NAME test_ncpp COMMAND test_ncpp WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})

# Lattice tests
add_executable(test_lattice test_lattice.cpp)
target_link_libraries(test_lattice PRIVATE io)
add_test(NAME test_lattice COMMAND test_lattice WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})

# EIGEN tests
add_executable(test_eigen test_eigen.cpp)
target_link_libraries(test_eigen PRIVATE io)
add_test(NAME test_eigen COMMAND test_eigen WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
