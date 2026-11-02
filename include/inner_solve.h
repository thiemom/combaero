#pragma once

// Exit handling shared by the kernels' iterative inner solves (#868): a loop
// that reaches its iteration cap without meeting its tolerance reports it
// instead of returning the last iterate as if it were a root.

#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <string>

namespace combaero::inner_solve {

// A loop that reaches its cap still returns its last iterate when the last
// step (or its bracket) is within this fraction of the iterate: a tolerance
// set at 1e-13..1e-15 can cycle in the last bits of a root it has found
// (the Jobson expansibility fixed point does, at 3e-14). A wider exit is not
// a root and throws.
constexpr double kCapAcceptRel = 1e-10;

// Iteration budgets of the solves that report a cap exit through this header
// (the audit table, docs/archive/AUDIT_868_INNER_SOLVES.md, lists them all).
// Values unchanged from the literals they replace.
constexpr int kStaticStateMaxIter = 100;          // chamber_ / sonic_static_state
constexpr int kSetVHMaxIter = 20;                 // State::set_VH_mass
constexpr double kSetVHTolK = 1e-4;               // State::set_VH_mass step [K]
constexpr int kMergeFaceMaxIter = 100;            // side-stream merge face fixed point
constexpr int kCrossflowCdMaxIter = 100;          // crossflow Cd, falling branch
constexpr int kAreaChangeMaxIter = 200;           // area_change Mach / flux solves
constexpr int kPipeEquivalentCMaxIter = 100;      // orifice pipe-equivalent C
constexpr int kJobsonExpansibilityMaxIter = 200;  // Jobson expansibility fixed point

// True when a cap exit is a root to rounding: last step or bracket width
// within kCapAcceptRel of the iterate's magnitude (false for a NaN).
inline bool cap_exit_is_root(double last_step, double x) {
  return std::abs(last_step) <= kCapAcceptRel * std::abs(x);
}

// std::runtime_error naming the solve and its last residual (or step).
// Iterations is any integer type (the budgets are int or std::size_t).
template <class N>
[[noreturn]] void throw_no_convergence(const char* fn, N iterations, const char* what,
                                       double value) {
  char buf[64];
  std::snprintf(buf, sizeof buf, "%.6e", value);
  throw std::runtime_error(std::string(fn) + ": no convergence in " +
                           std::to_string(iterations) + " iterations (last " + what + " = " +
                           buf + ")");
}

}  // namespace combaero::inner_solve
