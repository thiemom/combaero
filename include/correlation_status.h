#ifndef CORRELATION_STATUS_H
#define CORRELATION_STATUS_H

#include <cmath>
#include <functional>
#include <string>

namespace combaero {

// Status returned by correlations when called outside their validated range.
// Valid        : input within the correlation's fit range.
// Extrapolated : input outside fit range; result is a smooth power-law
//                extrapolation and is still usable in Newton/CasADi solvers.
enum class CorrelationStatus : uint8_t {
    Valid        = 0,
    Extrapolated = 1,
};

// Warning handler type.  Default implementation writes to std::cerr.
using WarningHandler = std::function<void(const std::string&)>;

// Replace the global warning handler.  Pass an empty function to suppress all
// warnings.  Not thread-safe — call once at program startup or use
// suppress_warnings() in Python for batch runs.
void set_warning_handler(WarningHandler fn);

// Retrieve the current handler (useful for save/restore patterns).
WarningHandler get_warning_handler();

// Emit a warning through the current handler.
void warn(const std::string& msg);

// Returns true when a correlation result is safe to use in a Jacobian-based
// solver: finite, strictly positive (> lo), and below the overflow guard (< hi).
// lo defaults to 1e-12 (covers Nu, h, f — all must be positive).
// hi defaults to 1e15 (well below double overflow).
inline bool is_well_behaved(double v, double lo = 1e-12, double hi = 1e15) {
    return std::isfinite(v) && v > lo && v < hi;
}

}  // namespace combaero

#endif  // CORRELATION_STATUS_H
