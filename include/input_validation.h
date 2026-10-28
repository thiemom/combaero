#pragma once

// NaN-safe input validation helpers (#769).
//
// A check written "if (x <= 0) throw ..." lets NaN through: every comparison
// with NaN is false, so the function answers with a silent NaN or, after a
// clamp, a plausible number. The predicates below are written so that NaN
// and +/-inf fail them, and the require_* helpers throw when they fail.
//
// The exception type is a template parameter so that each function keeps the
// type it already throws (std::invalid_argument -> Python ValueError,
// std::runtime_error -> Python RuntimeError, std::out_of_range -> Python
// IndexError). The default is std::invalid_argument.
//
// The require_* helpers take a C string so that a passing check builds no
// std::string. A message that formats the offending value should test the
// predicate and build the message only on failure.

#include <cmath>
#include <stdexcept>

namespace combaero::validation {

// x is finite (not NaN, not +/-inf).
inline bool is_finite(double x) { return std::isfinite(x); }

// x is finite and x > 0.
inline bool is_finite_positive(double x) { return std::isfinite(x) && x > 0.0; }

// x is finite and x >= 0.
inline bool is_finite_nonnegative(double x) {
  return std::isfinite(x) && x >= 0.0;
}

// x is finite and lo <= x <= hi (closed interval).
inline bool is_finite_in_range(double x, double lo, double hi) {
  return std::isfinite(x) && x >= lo && x <= hi;
}

// Throw E(message) unless x is finite.
template <class E = std::invalid_argument>
inline void require_finite(double x, const char *message) {
  if (!is_finite(x)) {
    throw E(message);
  }
}

// Throw E(message) unless x is finite and x > 0.
template <class E = std::invalid_argument>
inline void require_finite_positive(double x, const char *message) {
  if (!is_finite_positive(x)) {
    throw E(message);
  }
}

// Throw E(message) unless x is finite and x >= 0.
template <class E = std::invalid_argument>
inline void require_finite_nonnegative(double x, const char *message) {
  if (!is_finite_nonnegative(x)) {
    throw E(message);
  }
}

// Throw E(message) unless x is finite and lo <= x <= hi.
template <class E = std::invalid_argument>
inline void require_finite_in_range(double x, double lo, double hi,
                                    const char *message) {
  if (!is_finite_in_range(x, lo, hi)) {
    throw E(message);
  }
}

} // namespace combaero::validation
