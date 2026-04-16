---
name: cpp23-modules-cmake
description: Configure and build C++23 code with modules using CMake and Ninja. Use when setting up C++23 projects with `import std`, custom modules (.cppm), or mixed header/module code. Covers clang++ with libc++ configuration for experimental `import std` support.
---

# C++23 Modules with CMake/Ninja

## Quick Setup

### Basic CMakeLists.txt (import std)

```cmake
cmake_minimum_required(VERSION 3.31)

set(CMAKE_CXX_COMPILER "clang++")
set(CMAKE_CXX_FLAGS "-stdlib=libc++")

# Required for experimental import std
set(CMAKE_EXPERIMENTAL_CXX_IMPORT_STD "d0edc3af-4c50-42ea-a356-e2862fe7a444")

project(my_project LANGUAGES CXX)

set(CMAKE_CXX_MODULE_STD ON)
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(my_app "main.cpp")
```

### With Custom Module

```cmake
add_executable(my_app "main.cpp")

target_sources(my_app
    PUBLIC
        FILE_SET cxx_modules TYPE CXX_MODULES FILES "mymodule.cppm"
)
```

### Mixed: Headers + Custom Module (no import std)

```cmake
cmake_minimum_required(VERSION 3.31)

set(CMAKE_CXX_COMPILER "clang++")

project(my_project LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(my_app "main.cpp")
target_sources(my_app
    PUBLIC
        FILE_SET cxx_modules TYPE CXX_MODULES FILES "mymodule.cppm"
)
```

## Build Commands

```bash
# Configure
cmake -B build -G Ninja

# Build
cmake --build build

# Or with Ninja directly
ninja -C build
```

## Key Source Patterns

### import std only
```cpp
import std;
int main() { std::cout << "Hello\n"; }
```

### Custom module (cm.cppm)
```cpp
export module cm;
export int sum(int, int);
int sum(int a, int b) { return a+b; }
```

### Use custom module
```cpp
import std;
import cm;
int main() { std::cout << sum(1, 2); }
```

## Requirements

| Component | Minimum Version |
|-----------|-----------------|
| CMake | 3.31 |
| Ninja | 1.11 |
| Clang | 18+ |
| libc++ | 18+ |

## CXX Version Requirement

**C++23 is REQUIRED.** Do NOT use C++20 or lower.

```cmake
# CORRECT
set(CMAKE_CXX_STANDARD 23)

# WRONG - C++20 modules are incompatible
set(CMAKE_CXX_STANDARD 20)
```

**Why:** C++23 significantly changed modules semantics:
- `import std` only works in C++23
- Module interface syntax differs between C++20 and C++23
- CMake's `CXX_MODULES` support targets C++23

**Build failure without C++23:**
```
error: 'import' requires a C++23 compatible compiler and standard library
error: unknown type name 'module'
```

If C++23 is unavailable, use traditional headers (`#include`) instead of modules.

## Troubleshooting

| Error | Solution |
|-------|----------|
| `import std` not found | Ensure `CMAKE_CXX_MODULE_STD=ON` and libc++ is installed |
| Module interface errors | Check `.cppm` files are in `FILE_SET cxx_modules` |
| Cache errors | Delete `build/` and `CMakeCache.txt`, reconfigure |

## Experimental UUID Reference

The UUID `d0edc3af-4c50-42ea-a356-e2862fe7a444` enables experimental `import std` support in CMake 3.31. This is required until the feature stabilizes.
