cmake_minimum_required(VERSION 3.15)

# Use scikit-build-core provided project name/version when available, with
# sensible fallbacks for plain CMake builds.
set(_proj_name "${SKBUILD_PROJECT_NAME}")
if(NOT _proj_name)
    set(_proj_name "combaero")
endif()

set(_proj_version "${SKBUILD_PROJECT_VERSION}")
if(NOT _proj_version)
    set(_proj_version "0.0.2")
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
    src/compressible.cpp
    src/incompressible.cpp
    src/friction.cpp
    src/geometry.cpp
    src/acoustics.cpp
    src/can_annular_solvers.cpp
    src/heat_transfer.cpp
    src/humidair.cpp
    src/equilibrium.cpp
    src/state.cpp
    src/utils.cpp
    src/orifice.cpp
    src/pipe_flow.cpp
    src/materials.cpp
    src/cooling_correlations.cpp
)

# Core library
add_library(combaero_lib STATIC ${SOURCES})
target_include_directories(combaero_lib PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)
# Enable position-independent code for shared library linking
set_target_properties(combaero_lib PROPERTIES POSITION_INDEPENDENT_CODE ON)

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
add_executable(combustion_example examples/combustion_example.cpp)
add_executable(cooling_example examples/cooling_example.cpp)

# Add compiler warnings for examples
if (MSVC)
    target_compile_options(combaero_example PRIVATE /W4)
    target_compile_options(humidair_example PRIVATE /W4)
    target_compile_options(combustion_example PRIVATE /W4)
    target_compile_options(cooling_example PRIVATE /W4)
else()
    target_compile_options(combaero_example PRIVATE -Wall -Wextra)
    target_compile_options(humidair_example PRIVATE -Wall -Wextra)
    target_compile_options(combustion_example PRIVATE -Wall -Wextra)
    target_compile_options(cooling_example PRIVATE -Wall -Wextra)
endif()

target_link_libraries(combaero_example PRIVATE combaero_lib)
target_link_libraries(humidair_example PRIVATE combaero_lib)
target_link_libraries(combustion_example PRIVATE combaero_lib)
target_link_libraries(cooling_example PRIVATE combaero_lib)

# Print message when examples are built
add_custom_command(TARGET combaero_example POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E echo "CombAero example built! Run with: ./combaero_example"
)

add_custom_command(TARGET humidair_example POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E echo "Humid air example built! Run with: ./humidair_example"
)

add_custom_command(TARGET combustion_example POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E echo "Combustion example built! Run with: ./combustion_example"
)

add_custom_command(TARGET cooling_example POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E echo "Cooling example built! Run with: ./cooling_example"
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

# Pre-build source code checks (non-ASCII detection)
option(COMBAERO_ENABLE_SOURCE_CHECKS "Run source code quality checks before build" ON)

if (COMBAERO_ENABLE_SOURCE_CHECKS AND NOT WIN32)
    find_program(BASH_EXECUTABLE bash)
    if (BASH_EXECUTABLE)
        add_custom_target(check_source_style ALL
            COMMAND ${BASH_EXECUTABLE} ${CMAKE_SOURCE_DIR}/scripts/check-source-style.sh
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
            COMMENT "Checking source code style (non-ASCII, block comments)..."
        )
        # Make library depend on the check
        add_dependencies(combaero_lib check_source_style)
    endif()

    # Regenerate UNITS.md from units_data.h
    find_program(PYTHON_EXECUTABLE python3 python)
    if (PYTHON_EXECUTABLE AND EXISTS "${CMAKE_SOURCE_DIR}/scripts/generate_units_md.py")
        add_custom_target(generate_units_md ALL
            COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_SOURCE_DIR}/scripts/generate_units_md.py
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
            COMMENT "Regenerating docs/UNITS.md..."
            DEPENDS ${CMAKE_SOURCE_DIR}/include/units_data.h
        )
        add_dependencies(combaero_lib generate_units_md)
    endif()
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
