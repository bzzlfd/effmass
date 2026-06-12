cmake_minimum_required(VERSION 3.31)

set(CMAKE_CXX_COMPILER "clang++")
set(CMAKE_CXX_FLAGS "-stdlib=libc++")

# Required for experimental import std in CMake
set(CMAKE_EXPERIMENTAL_CXX_IMPORT_STD "d0edc3af-4c50-42ea-a356-e2862fe7a444")
set(CMAKE_CXX_MODULE_STD ON)
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project(effmass LANGUAGES CXX)

set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# pugixml (used by UPF reader)
add_subdirectory(vendor/pugixml)

# Library target for io module
add_library(io)
target_sources(io
    PUBLIC
        FILE_SET cxx_modules TYPE CXX_MODULES FILES
            src/io.cppm
            src/io/GKK.cppm
            src/io/WG.cppm
            src/io/ATOM.cppm
            src/io/EIGEN.cppm
            src/io/OCC.cppm
            src/io/RHO.cppm
            src/io/VR.cppm
            src/io/lattice.cppm
)
target_link_libraries(io PRIVATE utils)

# Utility types
add_library(utils)
target_sources(utils
    PUBLIC
        FILE_SET cxx_modules TYPE CXX_MODULES FILES
            src/utils/array2d.cppm
            src/utils/vector3d.cppm
            src/utils/logger.cppm
)

# Library target for pseudo-potential module
add_library(pseudo)
target_sources(pseudo
    PUBLIC
        FILE_SET cxx_modules TYPE CXX_MODULES FILES
            src/pseudo.cppm
            src/pseudo/io/upf.cppm
            src/pseudo/ncpp.cppm
)
target_link_libraries(pseudo PRIVATE pugixml::pugixml utils math)

# Library target for math module
add_library(math)
target_sources(math
    PUBLIC
        FILE_SET cxx_modules TYPE CXX_MODULES FILES
            src/math.cppm
            src/math/sph_bessel.cppm
            src/math/sph_harmonics.cppm
            src/math/fft.cppm
            src/math/fourier_bessel.cppm
            src/math/linalg.cppm
    PRIVATE
        vendor/pocketfft/pocketfft.c
)
set_source_files_properties(vendor/pocketfft/pocketfft.c PROPERTIES
    LANGUAGE CXX
    COMPILE_DEFINITIONS "restrict=__restrict__"
)
target_include_directories(math PRIVATE vendor/pocketfft)
target_link_libraries(math PRIVATE utils)

# Library target for H_psi module
add_library(H_psi)
set(HPSI_DEBUG 0 CACHE STRING "H|psi> debug level: 0=off, 1=summary, 2=detailed")
if(HPSI_DEBUG)
    target_compile_definitions(H_psi PRIVATE HPSI_DEBUG=${HPSI_DEBUG})
endif()
target_sources(H_psi
    PUBLIC
        FILE_SET cxx_modules TYPE CXX_MODULES FILES
            src/H_psi.cppm
            src/H_psi/hamiltonian.cppm
            src/H_psi/structure_factor.cppm
            src/support/density.cppm
    PRIVATE
        src/H_psi/hamiltonian.cpp
        src/H_psi/hamiltonian_callable.cpp
        src/H_psi/hamiltonian_gradient.cpp
        src/H_psi/hamiltonian_hessian.cpp
)
target_link_libraries(H_psi PRIVATE io math pseudo utils)

# Tests
enable_testing()
add_subdirectory(test)

