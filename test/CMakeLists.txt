# IO reader tests
add_executable(test-io-eigen io/test_eigen_occ.cpp)
target_link_libraries(test-io-eigen PRIVATE io)
add_test(NAME test-io-eigen COMMAND test-io-eigen WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})

add_executable(test-io-atom io/test_atom.cpp)
target_link_libraries(test-io-atom PRIVATE io)
add_test(NAME test-io-atom COMMAND test-io-atom WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})

add_executable(test-io-rho io/test_rho.cpp)
target_link_libraries(test-io-rho PRIVATE io)
add_test(NAME test-io-rho COMMAND test-io-rho WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})

add_executable(test-io-gkk io/test_gkk.cpp)
target_link_libraries(test-io-gkk PRIVATE io)
add_test(NAME test-io-gkk COMMAND test-io-gkk WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})

add_executable(test-io-wg io/test_wg.cpp)
target_link_libraries(test-io-wg PRIVATE io)
add_test(NAME test-io-wg COMMAND test-io-wg WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})

add_executable(test-io-lattice io/test_lattice.cpp)
target_link_libraries(test-io-lattice PRIVATE io)
add_test(NAME test-io-lattice COMMAND test-io-lattice WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})

# Pseudopotential tests
add_executable(test-pseudo-io-upf pseudo/io/test_ncpp_upf.cpp)
target_link_libraries(test-pseudo-io-upf PRIVATE pseudo)
add_test(NAME test-pseudo-io-upf COMMAND test-pseudo-io-upf WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})

add_executable(test-pseudo-ncpp pseudo/test_ncpp.cpp)
target_link_libraries(test-pseudo-ncpp PRIVATE pseudo)
add_test(NAME test-pseudo-ncpp COMMAND test-pseudo-ncpp WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})

add_executable(test-pseudo-io-soc-upf pseudo/io/test_soc_upf.cpp)
target_link_libraries(test-pseudo-io-soc-upf PRIVATE pseudo)
add_test(NAME test-pseudo-io-soc-upf COMMAND test-pseudo-io-soc-upf WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})

add_executable(test-pseudo-ncpp-upf2upfSO pseudo/test_ncpp_upf2upfSO.cpp)
target_link_libraries(test-pseudo-ncpp-upf2upfSO PRIVATE pseudo)
add_test(NAME test-pseudo-ncpp-upf2upfSO COMMAND test-pseudo-ncpp-upf2upfSO WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})

# Math tests
add_executable(test-math-sph-bessel math/test_sph_bessel.cpp)
target_link_libraries(test-math-sph-bessel PRIVATE math)
add_test(NAME test-math-sph-bessel COMMAND test-math-sph-bessel WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})

add_executable(test-math-sph-harmonics math/test_sph_harmonics.cpp)
target_link_libraries(test-math-sph-harmonics PRIVATE math)
add_test(NAME test-math-sph-harmonics COMMAND test-math-sph-harmonics WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})

add_executable(test-math-fft math/test_fft.cpp)
target_link_libraries(test-math-fft PRIVATE math)
add_test(NAME test-math-fft COMMAND test-math-fft WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})

add_executable(test-math-fourier-bessel math/test_fourier_bessel.cpp)
target_link_libraries(test-math-fourier-bessel PRIVATE math pseudo)
add_test(NAME test-math-fourier-bessel COMMAND test-math-fourier-bessel WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})

add_executable(test-math-linalg math/test_linalg.cpp)
target_link_libraries(test-math-linalg PRIVATE math)
add_test(NAME test-math-linalg COMMAND test-math-linalg WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})

# Cross-module tests
add_executable(test-wg-fft test_wg_fft.cpp)
target_link_libraries(test-wg-fft PRIVATE io math)
add_test(NAME test-wg-fft COMMAND test-wg-fft WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})

# Hamiltonian tests
add_executable(test-hamiltonian H_psi/test_hamiltonian.cpp)
target_link_libraries(test-hamiltonian PRIVATE H_psi)
add_test(NAME test-hamiltonian COMMAND test-hamiltonian WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})

# StructureFactor tests
add_executable(test-structure-factor H_psi/test_structure_factor.cpp)
target_link_libraries(test-structure-factor PRIVATE H_psi)
add_test(NAME test-structure-factor COMMAND test-structure-factor WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})

# H|ψ⟩ = E|ψ⟩ eigenstate test
add_executable(test-hpsi-eigen H_psi/test_hpsi_eigen.cpp)
target_link_libraries(test-hpsi-eigen PRIVATE H_psi)
add_test(NAME test-hpsi-eigen COMMAND test-hpsi-eigen WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
