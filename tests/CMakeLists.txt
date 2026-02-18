# Tests CMakeLists.txt

# Create a library target for the test sources
set(TEST_SOURCES
    test_thermo_transport.cpp  # Core thermodynamic property tests
    test_humidair.cpp          # Humid air calculation tests
    test_incompressible.cpp    # Incompressible flow tests (zeta, Cd, pressure loss)
    test_orifice.cpp           # Orifice Cd correlation tests
    test_acoustic_adapter.cpp  # Geometry-to-acoustics adapter tests
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
