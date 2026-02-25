#ifndef COMBAERO_MATH_CONSTANTS_H
#define COMBAERO_MATH_CONSTANTS_H

// Portable math constants for cross-platform compatibility.
// MSVC does not define M_PI by default (it's a POSIX extension).

#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifndef M_PI_2
#define M_PI_2 1.57079632679489661923
#endif

#ifndef M_PI_4
#define M_PI_4 0.78539816339744830962
#endif

#ifndef M_SQRT2
#define M_SQRT2 1.41421356237309504880
#endif

#ifndef M_LN10
#define M_LN10 2.30258509299404568402
#endif

#endif // COMBAERO_MATH_CONSTANTS_H
