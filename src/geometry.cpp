#include "../include/geometry.h"
#include "../include/math_constants.h"
#include <stdexcept>

double hydraulic_diameter(double A, double P_wetted) {
    if (A <= 0.0 || P_wetted <= 0.0) {
        throw std::invalid_argument("hydraulic_diameter: A and P_wetted must be positive");
    }
    return 4.0 * A / P_wetted;
}

double hydraulic_diameter_rect(double a, double b) {
    if (a <= 0.0 || b <= 0.0) {
        throw std::invalid_argument("hydraulic_diameter_rect: a and b must be positive");
    }
    return 2.0 * a * b / (a + b);
}

double hydraulic_diameter_annulus(double D_outer, double D_inner) {
    if (D_outer <= D_inner || D_inner < 0.0) {
        throw std::invalid_argument("hydraulic_diameter_annulus: D_outer > D_inner >= 0 required");
    }
    return D_outer - D_inner;
}

// -------------------------------------------------------------
// Tube methods
// -------------------------------------------------------------

double Tube::area() const {
    return M_PI * D * D / 4.0;
}

double Tube::volume() const {
    return area() * L;
}

double Tube::perimeter() const {
    return M_PI * D;
}

// -------------------------------------------------------------
// Annulus methods
// -------------------------------------------------------------

double Annulus::D_mean() const {
    return (D_inner + D_outer) / 2.0;
}

double Annulus::gap() const {
    return (D_outer - D_inner) / 2.0;
}

double Annulus::area() const {
    return M_PI * (D_outer * D_outer - D_inner * D_inner) / 4.0;
}

double Annulus::volume() const {
    return area() * L;
}

double Annulus::circumference() const {
    return M_PI * D_mean();
}
