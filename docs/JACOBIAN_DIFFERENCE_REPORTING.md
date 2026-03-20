# Jacobian Difference Reporting Flag

## Overview

A new environment variable `PRINT_JACOBIAN_DIFFERENCES` has been added to the C++ Jacobian tests to collect and analyze numerical accuracy patterns across all derivative calculations.

## Usage

### Basic Usage
```bash
# Run tests normally (no detailed output, just summary at end)
cd build && ./test_solver_jacobians

# Run with detailed difference reporting
cd build && PRINT_JACOBIAN_DIFFERENCES=1 ./test_solver_jacobians

# Run specific tests with detailed reporting
cd build && PRINT_JACOBIAN_DIFFERENCES=1 ./test_solver_jacobians --gtest_filter="*Orifice*"
```

### File Output (NEW)
```bash
# Write differences to CSV file for persistent records
cd build && JACOBIAN_OUTPUT_FILE=jacobian_differences.csv PRINT_JACOBIAN_DIFFERENCES=1 ./test_solver_jacobians

# File output without console display
cd build && JACOBIAN_OUTPUT_FILE=jacobian_differences.csv ./test_solver_jacobians

# Append to existing file (useful for collecting data over time)
cd build && JACOBIAN_OUTPUT_FILE=jacobian_history.csv PRINT_JACOBIAN_DIFFERENCES=1 ./test_solver_jacobians --gtest_filter="*Orifice*"
```

### Environment Variables
- `PRINT_JACOBIAN_DIFFERENCES=1` or `true`: Enable detailed console reporting
- `JACOBIAN_OUTPUT_FILE=filename.csv`: Write differences to CSV file
- Both can be used independently or together

## Output Format

### Detailed Output (when flag is enabled)
```
JACOBIAN_DIFF: TestName | derivative_name | analytical=X.XXXXXX | fd=X.XXXXXX | abs_diff=X.XXXXXX | rel_diff=X.XXXXXX
```

### CSV File Output (when JACOBIAN_OUTPUT_FILE is set)
```csv
test_name,derivative,analytical,finite_diff,abs_diff,rel_diff
OrificeDerivatives,d_mdot_dP_total_up,0.000004745402,0.000004745343,0.000000000059,0.000012499853
OrificeDerivatives,d_mdot_dP_static_down,-0.000004745402,-0.000004745461,0.000000000059,0.000012500164
...
```

### Summary Output (always printed at end)
```
=== JACOBIAN DIFFERENCE SUMMARY ===
                     Test          Derivative       Abs Diff       Rel Diff
---------------------------------------------------------------------
OrificeDerivatives  d_mdot_dT_up          0.000          0.000
PipeDerivatives     d_dP_dT_up            0.001          0.000
...
---------------------------------------------------------------------
Maximum absolute difference: 1.234e-10
Maximum relative difference: 5.678e-06
Total comparisons: 42
=== END SUMMARY ===
```

## Implementation Details

### Adding to New Tests
To add difference reporting to a new Jacobian test:

```cpp
// For each derivative comparison
double analytical = d_mdot_dT_up;      // From function
double finite_diff = fd_dT_up;         // From finite difference

// Report the difference
report_jacobian_difference("TestName", "d_mdot_dT_up", analytical, finite_diff);

// Continue with normal test
EXPECT_NEAR(analytical, finite_diff, tolerance);
```

### Data Collection
- **Global storage**: All differences are stored in `g_jacobian_differences` vector
- **Automatic summary**: Printed at test suite completion via test environment
- **Statistics**: Maximum absolute/relative differences and total comparisons

## Pattern Analysis Use Cases

### 1. Identify Problematic Functions
```bash
PRINT_JACOBIAN_DIFFERENCES=1 ./test_solver_jacobians | grep "JACOBIAN_DIFF" | sort -k6 -n
```
Shows derivatives with largest absolute differences.

### 2. Compare Relative vs Absolute Accuracy
The summary helps distinguish between:
- Large absolute differences (might be acceptable for large-magnitude derivatives)
- Large relative differences (problematic even for small magnitudes)

### 3. Track Regression
Run before and after changes to detect accuracy regressions:
```bash
# Before changes
PRINT_JACOBIAN_DIFFERENCES=1 ./test_solver_jacobians > before.txt

# After changes
PRINT_JACOBIAN_DIFFERENCES=1 ./test_solver_jacobians > after.txt

# Compare
diff before.txt after.txt
```

### 4. Optimize Tolerances
Use the actual differences to set appropriate test tolerances:
- If max relative difference is 1e-8, set tolerance to 1e-6 for safety margin
- Identify outliers that might need special handling

## Current Test Coverage

### Tests with Reporting Enabled
- `OrificeDerivatives` - 4 derivatives
- `AdiabaticWallDerivatives` - 1 derivative

### Tests Needing Updates
The following tests need `report_jacobian_difference()` calls added:
- `PipeDerivatives` - 4 derivatives
- `OrificeCompressibleDerivatives` - 3 derivatives
- `PipeCompressibleMdotDerivatives` - 3 derivatives
- `MomentumChamberDerivatives` - 2 derivatives
- `AdiabaticTCompleteDerivatives` - 1 derivative
- `AdiabaticTEquilibriumDerivatives` - 2 derivatives
- `T0FromStaticDerivatives` - 1 derivative
- `P0FromStaticDerivatives` - 1 derivative
- `MachNumberDerivatives` - 1 derivative

## Examples

### Example 1: Check All Orifice Tests
```bash
PRINT_JACOBIAN_DIFFERENCES=1 ./test_solver_jacobians --gtest_filter="*Orifice*"
```

### Example 2: Find Largest Differences
```bash
PRINT_JACOBIAN_DIFFERENCES=1 ./test_solver_jacobians 2>&1 | grep "Maximum"
```

### Example 3: Export to CSV File Directly
```bash
# Write directly to CSV file (recommended)
JACOBIAN_OUTPUT_FILE=jacobian_differences.csv ./test_solver_jacobians

# With detailed console output too
JACOBIAN_OUTPUT_FILE=jacobian_differences.csv PRINT_JACOBIAN_DIFFERENCES=1 ./test_solver_jacobians
```

### Example 4: Time Series Data Collection
```bash
# Create timestamped files for tracking changes over time
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
JACOBIAN_OUTPUT_FILE=jacobian_history_${TIMESTAMP}.csv ./test_solver_jacobians

# Append to daily log file
JACOBIAN_OUTPUT_FILE=jacobian_$(date +%Y%m%d).csv PRINT_JACOBIAN_DIFFERENCES=1 ./test_solver_jacobians --gtest_filter="*Orifice*"
```

### Example 5: Export to CSV for Analysis
```bash
# Manual CSV extraction (alternative method)
PRINT_JACOBIAN_DIFFERENCES=1 ./test_solver_jacobians 2>&1 | \
  grep "JACOBIAN_DIFF" | \
  sed 's/.*| //' | \
  sed 's/ | /,/g' > manual_jacobian_differences.csv
```

## Benefits

1. **Pattern Recognition**: Identify which types of derivatives are consistently less accurate
2. **Tolerance Optimization**: Set tolerances based on actual numerical behavior
3. **Regression Detection**: Catch accuracy regressions in mathematical functions
4. **Debugging Aid**: Quickly identify problematic derivative calculations
5. **Quality Assurance**: Ensure consistent numerical accuracy across the codebase
6. **Persistent Records**: CSV files create historical data for long-term analysis
7. **Time Series Analysis**: Track accuracy changes over git history and code modifications
8. **Statistical Analysis**: Export to pandas/Excel/R for deeper statistical analysis
9. **Automated Reporting**: Integrate with CI/CD pipelines for automated accuracy tracking

## Future Enhancements

1. **Automatic Tolerance Setting**: Suggest optimal tolerances based on collected data
2. **Trend Analysis**: Track accuracy changes over git history
3. **Integration with CI**: Fail builds if differences exceed thresholds
4. **Statistical Analysis**: Identify outliers and statistical patterns
5. **Visualization**: Generate plots of accuracy distributions
