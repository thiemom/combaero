#pragma once
// Numerical integration utilities — header-only.
//
// composite_simpson: Composite Simpson's 1/3 rule.
// Error: O(h^4), where h = (b-a)/n.  For smooth integrands with n=64,
// absolute error is typically < 1e-9.

template <typename F>
double composite_simpson(F &&f, double a, double b, int n = 64) {
  // n must be even
  const double h = (b - a) / n;
  double sum = f(a) + f(b);
  for (int i = 1; i < n; i += 2)
    sum += 4.0 * f(a + i * h);
  for (int i = 2; i < n - 1; i += 2)
    sum += 2.0 * f(a + i * h);
  return (h / 3.0) * sum;
}
