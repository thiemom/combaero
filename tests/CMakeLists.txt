# Tests CMakeLists.txt

# Create a library target for the test sources
set(TEST_SOURCES
    test_thermo_transport.cpp  # Core thermodynamic property tests
    test_humidair.cpp          # Humid air calculation tests
)

# First create a library from our source files
add_library(combaero_lib STATIC
    ${CMAKE_SOURCE_DIR}/src/thermo_transport.cpp
    ${CMAKE_SOURCE_DIR}/src/humidair.cpp
)

# Link against Google Test and our library
add_executable(combaero_tests ${TEST_SOURCES})
target_link_libraries(combaero_tests PRIVATE 
    gtest 
    gtest_main
    combaero_lib
)

# Add tests to CTest
include(GoogleTest)
gtest_discover_tests(combaero_tests)

# Print message when tests are built
add_custom_command(TARGET combaero_tests POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E echo "Tests built! Run with: ./tests/combaero_tests"
)

# Add standalone test executables (these have their own main())
add_executable(test_fraction_functions test_fraction_functions.cpp)
target_link_libraries(test_fraction_functions PRIVATE combaero_lib)

add_executable(test_ice_equation_accuracy test_ice_equation_accuracy.cpp)
target_link_libraries(test_ice_equation_accuracy PRIVATE combaero_lib)

add_executable(test_water_equation_accuracy test_water_equation_accuracy.cpp)
target_link_libraries(test_water_equation_accuracy PRIVATE combaero_lib)

# Add custom targets to run the accuracy tests
add_custom_target(run_accuracy_tests
    COMMAND ${CMAKE_COMMAND} -E echo "\nRunning Ice Equation Accuracy Test..."
    COMMAND $<TARGET_FILE:test_ice_equation_accuracy>
    COMMAND ${CMAKE_COMMAND} -E echo "\nRunning Water Equation Accuracy Test..."
    COMMAND $<TARGET_FILE:test_water_equation_accuracy>
    COMMENT "Running accuracy tests against Hyland-Wexler reference values"
    DEPENDS test_ice_equation_accuracy test_water_equation_accuracy
)
