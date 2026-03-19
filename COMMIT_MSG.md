feat(network): Step 2 - Add htc_and_T() to element hierarchy

- Add htc_and_T() method to NetworkElement base class (returns None by default)
- Add surface parameter to PipeElement constructor with ConvectiveSurface default
- Override htc_and_T() in PipeElement to compute HTC using surface
- Handle T_wall=None by passing math.nan to C++ functions
- Add comprehensive test suite (test_pipe_element_htc.py) with 4 test cases:
  * Default pipe (area=0 returns None)
  * Pipe with custom surface (computes HTC)
  * Pipe with specified wall temperature
  * NetworkElement base class default behavior
- All tests pass (22/22 total: 18 existing + 4 new)
- Example script continues to work

This enables elements to compute convective heat transfer as specified
in the design document, preparing for full solver integration.
