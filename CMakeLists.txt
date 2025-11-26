cmake_minimum_required(VERSION 3.10)
project(combaero)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Include directories
include_directories(include)

# Add source files
set(SOURCES
    src/thermo_transport.cpp
    src/humidair.cpp
)

# Create main executable
add_executable(combaero_main ${SOURCES} src/main.cpp)
target_compile_options(combaero_main PRIVATE -Wall -Wextra)

# Print message when main build is complete
add_custom_command(TARGET combaero_main POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E echo "Main executable built! Run with: ./combaero_main"
)

# Examples
add_executable(combaero_example ${SOURCES} examples/thermo_example.cpp)
add_executable(humidair_example ${SOURCES} examples/humidair_example.cpp)

# Add compiler flags for examples
target_compile_options(combaero_example PRIVATE -Wall -Wextra)
target_compile_options(humidair_example PRIVATE -Wall -Wextra)

# Print message when examples are built
add_custom_command(TARGET combaero_example POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E echo "CombAero example built! Run with: ./combaero_example"
)

add_custom_command(TARGET humidair_example POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E echo "Humid air example built! Run with: ./humidair_example"
)

# Testing setup
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
