fix: propagate channel jacobians through pipe and walls

- return ChannelResult from ConvectiveSurface/PipeElement and expose jacobians in tests
- extend solver wall coupling to store channel results and relay full chain rule
- update regression tests plus wall connection coverage to new API
- ruff-format and refresh docs/comments for wall evaluation
