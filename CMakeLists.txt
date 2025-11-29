cmake_minimum_required(VERSION 3.15)

# Use scikit-build-core provided project name/version when available, with
# sensible fallbacks for plain CMake builds.
set(_proj_name "${SKBUILD_PROJECT_NAME}")
if(NOT _proj_name)
    set(_proj_name "combaero")
endif()

set(_proj_version "${SKBUILD_PROJECT_VERSION}")
if(NOT _proj_version)
    set(_proj_version "0.0.1")
endif()

project(${_proj_name} LANGUAGES CXX VERSION ${_proj_version})

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Include directories
include_directories(include)

# Core library sources
set(SOURCES
    src/thermo.cpp
    src/transport.cpp
    src/combustion.cpp
    src/humidair.cpp
    src/equilibrium.cpp
    src/utils.cpp
)

# Core library
add_library(combaero_lib STATIC ${SOURCES})
target_include_directories(combaero_lib PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

# Create main executable
add_executable(combaero_main src/main.cpp)
if (MSVC)
    target_compile_options(combaero_main PRIVATE /W4)
else()
    target_compile_options(combaero_main PRIVATE -Wall -Wextra)
endif()
target_link_libraries(combaero_main PRIVATE combaero_lib)

# Print message when main build is complete
add_custom_command(TARGET combaero_main POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E echo "Main executable built! Run with: ./combaero_main"
)

# Examples
add_executable(combaero_example examples/thermo_example.cpp)
add_executable(humidair_example examples/humidair_example.cpp)

# Add compiler warnings for examples
if (MSVC)
    target_compile_options(combaero_example PRIVATE /W4)
    target_compile_options(humidair_example PRIVATE /W4)
else()
    target_compile_options(combaero_example PRIVATE -Wall -Wextra)
    target_compile_options(humidair_example PRIVATE -Wall -Wextra)
endif()

target_link_libraries(combaero_example PRIVATE combaero_lib)
target_link_libraries(humidair_example PRIVATE combaero_lib)

# Print message when examples are built
add_custom_command(TARGET combaero_example POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E echo "CombAero example built! Run with: ./combaero_example"
)

add_custom_command(TARGET humidair_example POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E echo "Humid air example built! Run with: ./humidair_example"
)

# Optional Python bindings (enabled explicitly or via scikit-build-core)
option(COMBAERO_BUILD_PYTHON "Build Python extension" OFF)

if (COMBAERO_BUILD_PYTHON)
    # Let pybind11 find and configure Python
    set(PYBIND11_FINDPYTHON ON)
    find_package(pybind11 CONFIG REQUIRED)

    # Define the extension module via pybind11_add_module
    pybind11_add_module(_core python/combaero/_core.cpp)

    target_link_libraries(_core PRIVATE combaero_lib)
    target_include_directories(_core PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/include
    )

    # Install the extension into the Python package namespace; scikit-build-core
    # will map this into the wheel under the combaero package.
    install(TARGETS _core DESTINATION ${_proj_name})
endif()

# Testing setup (can be disabled, e.g. for Python wheel builds)
option(COMBAERO_ENABLE_TESTS "Enable building C++ tests" ON)

if (COMBAERO_ENABLE_TESTS)
    enable_testing()

    # Include FetchContent for downloading Google Test
    include(FetchContent)

    # Download and configure Google Test
    FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG release-1.12.1
    )

    # For Windows: Prevent overriding the parent project's compiler/linker settings
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(googletest)

    # Add test directory
    add_subdirectory(tests)
endif()

# Optional Python tests that exercise the installed combaero wheel.
option(COMBAERO_ENABLE_PYTHON_TESTS "Enable Python-level tests (requires pytest and installed wheel)" OFF)

if (COMBAERO_ENABLE_PYTHON_TESTS)
    find_package(Python COMPONENTS Interpreter REQUIRED)

    add_test(
        NAME python_import_tests
        COMMAND ${Python_EXECUTABLE} -m pytest
                ${CMAKE_SOURCE_DIR}/python/tests/test_import.py
    )
endif()
