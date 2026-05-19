# IO reader tests
add_executable(test_io_ncpp test_io_ncpp.cpp)
target_link_libraries(test_io_ncpp PRIVATE pseudo)
add_test(NAME test_io_ncpp COMMAND test_io_ncpp WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})

add_executable(test_io_eigen test_io_eigen.cpp)
target_link_libraries(test_io_eigen PRIVATE io)
add_test(NAME test_io_eigen COMMAND test_io_eigen WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})

add_executable(test_io_rho test_io_rho.cpp)
target_link_libraries(test_io_rho PRIVATE io)
add_test(NAME test_io_rho COMMAND test_io_rho WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})

add_executable(test_io_gkk test_io_gkk.cpp)
target_link_libraries(test_io_gkk PRIVATE io)
add_test(NAME test_io_gkk COMMAND test_io_gkk WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})

add_executable(test_io_wg test_io_wg.cpp)
target_link_libraries(test_io_wg PRIVATE io)
add_test(NAME test_io_wg COMMAND test_io_wg WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})

add_executable(test_io_lattice test_io_lattice.cpp)
target_link_libraries(test_io_lattice PRIVATE io)
add_test(NAME test_io_lattice COMMAND test_io_lattice WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})

# Potential operator tests
add_executable(test_pseudo_ncpp test_pseudo_ncpp.cpp)
target_link_libraries(test_pseudo_ncpp PRIVATE pseudo)
add_test(NAME test_pseudo_ncpp COMMAND test_pseudo_ncpp WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})

# Math tests (FFT, spherical harmonics, spherical Bessel)
add_executable(test_math test_math.cpp)
target_link_libraries(test_math PRIVATE math)
add_test(NAME test_math COMMAND test_math WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
