cmake_minimum_required(VERSION 3.31)

set(CMAKE_CXX_COMPILER "clang++")
set(CMAKE_CXX_FLAGS "-stdlib=libc++")

# Required for experimental import std in CMake
set(CMAKE_EXPERIMENTAL_CXX_IMPORT_STD "d0edc3af-4c50-42ea-a356-e2862fe7a444")
set(CMAKE_CXX_MODULE_STD ON)
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project(effmass LANGUAGES CXX)

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
            src/io/EIGEN.cppm
            src/io/RHO.cppm
            src/io/VR.cppm
            src/io/lattice.cppm
)

# Library target for pseudo-potential module
add_library(pseudo)
target_sources(pseudo
    PUBLIC
        FILE_SET cxx_modules TYPE CXX_MODULES FILES
            src/pseudo.cppm
            src/pseudo/io/ncpp.upf.cppm
            src/pseudo/io/uspp.upf.cppm
            src/pseudo/io/paw.upf.cppm
            src/pseudo/ncpp.cppm
)
target_link_libraries(pseudo PRIVATE pugixml::pugixml)

option(EFFMASS_USE_POCKETFFT "Use PocketFFT for FFT (falls back to custom implementation)" ON)

# Library target for transform module
add_library(transform)
target_sources(transform
    PUBLIC
        FILE_SET cxx_modules TYPE CXX_MODULES FILES
            src/transform.cppm
            src/transform/sph_bessel.cppm
            src/transform/sph_harmonics.cppm
            src/transform/fft3d.cppm
)

if(EFFMASS_USE_POCKETFFT)
    target_compile_definitions(transform PRIVATE EFFMASS_USE_POCKETFFT)
    target_sources(transform PRIVATE vendor/pocketfft/pocketfft.cpp)
    set_source_files_properties(vendor/pocketfft/pocketfft.cpp PROPERTIES COMPILE_DEFINITIONS "restrict=__restrict__")
    target_include_directories(transform PRIVATE vendor/pocketfft)
endif()

# Tests
enable_testing()
add_subdirectory(test)

