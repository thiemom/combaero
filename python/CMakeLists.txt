cmake_minimum_required(VERSION 3.15)

# Wrapper CMake project used when building the Python wheel via scikit-build-core.
# It pulls in the main C++ project from the parent directory.

project(combaero_python LANGUAGES CXX)

# Disable C++ tests when building the Python extension; they are not needed here
set(COMBAERO_ENABLE_TESTS OFF CACHE BOOL "" FORCE)

# Add the main project (one directory up) as a subdirectory.
# The second argument is a binary-dir name for that subproject.
add_subdirectory(.. combaero-core)
