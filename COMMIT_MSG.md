feat(network): Step 3 - Add FlowNetwork.add_wall() and thermal coupling toggle

- Add walls: dict[str, WallConnection] to FlowNetwork.__init__
- Add thermal_coupling_enabled: bool = True toggle for debug/fast-solve
- Add add_wall() method with validation:
  * Checks for duplicate wall IDs
  * Validates element_a and element_b exist in network
- Update to_dict() to serialize walls and thermal_coupling_enabled
- Update from_dict() to deserialize walls and restore thermal toggle
- Add comprehensive test suite (test_flow_network_walls.py) with 9 test cases:
  * Default values
  * Successful wall addition
  * Error handling (duplicate ID, unknown elements)
  * Thermal coupling toggle
  * Serialization/deserialization with walls
  * Default thermal coupling restoration
- All tests pass (31/31 total: 22 existing + 9 new)
- Error handling works correctly with descriptive messages

This enables thermal wall coupling management in the network,
preparing for solver integration in Step 4.
