#include "../include/geometry.h"
#include "../include/input_validation.h"
#include "../include/math_constants.h"
#include <stdexcept>

namespace combaero {

using validation::is_finite_nonnegative;
using validation::is_finite_positive;
using validation::require_finite_nonnegative;
using validation::require_finite_positive;

namespace {

// D_outer > D_inner >= 0, both finite. A swapped annulus used to return a
// negative area.
void require_annulus_diameters(double D_outer, double D_inner,
                               const char *message) {
  if (!is_finite_nonnegative(D_inner) || !is_finite_positive(D_outer) ||
      !(D_outer > D_inner)) {
    throw std::invalid_argument(message);
  }
}

} // namespace

double hydraulic_diameter(double A, double P_wetted) {
    if (!is_finite_positive(A) || !is_finite_positive(P_wetted)) {
        throw std::invalid_argument("hydraulic_diameter: A and P_wetted must be positive");
    }
    return 4.0 * A / P_wetted;
}

double hydraulic_diameter_rect(double a, double b) {
    if (!is_finite_positive(a) || !is_finite_positive(b)) {
        throw std::invalid_argument("hydraulic_diameter_rect: a and b must be positive");
    }
    return 2.0 * a * b / (a + b);
}

double hydraulic_diameter_annulus(double D_outer, double D_inner) {
    require_annulus_diameters(
        D_outer, D_inner,
        "hydraulic_diameter_annulus: D_outer > D_inner >= 0 required");
    return D_outer - D_inner;
}

// -------------------------------------------------------------
// Simple Geometry Helpers
// -------------------------------------------------------------

double circular_area(double D) {
    require_finite_positive(D, "circular_area: D must be positive");
    return M_PI * D * D / 4.0;
}

double annular_area(double D_outer, double D_inner) {
    require_annulus_diameters(D_outer, D_inner,
                              "annular_area: D_outer > D_inner >= 0 required");
    return M_PI * (D_outer * D_outer - D_inner * D_inner) / 4.0;
}

double cylinder_volume(double D, double L) {
    if (!is_finite_positive(D) || !is_finite_positive(L)) {
        throw std::invalid_argument("cylinder_volume: D and L must be positive");
    }
    return M_PI * D * D * L / 4.0;
}

// -------------------------------------------------------------
// Tube methods
// -------------------------------------------------------------

double Tube::area() const {
    require_finite_positive(D, "Tube: D must be positive");
    return M_PI * D * D / 4.0;
}

double Tube::volume() const {
    require_finite_positive(L, "Tube: L must be positive");
    return area() * L;
}

double Tube::perimeter() const {
    require_finite_positive(D, "Tube: D must be positive");
    return M_PI * D;
}

// -------------------------------------------------------------
// Annulus methods
// -------------------------------------------------------------

static constexpr const char *kAnnulusDiameters =
    "Annulus: D_outer > D_inner >= 0 required";

double Annulus::D_mean() const {
    require_annulus_diameters(D_outer, D_inner, kAnnulusDiameters);
    return (D_inner + D_outer) / 2.0;
}

double Annulus::gap() const {
    require_annulus_diameters(D_outer, D_inner, kAnnulusDiameters);
    return (D_outer - D_inner) / 2.0;
}

double Annulus::area() const {
    require_annulus_diameters(D_outer, D_inner, kAnnulusDiameters);
    return M_PI * (D_outer * D_outer - D_inner * D_inner) / 4.0;
}

double Annulus::volume() const {
    require_finite_positive(L, "Annulus: L must be positive");
    return area() * L;
}

double Annulus::circumference() const {
    return M_PI * D_mean();
}

// -------------------------------------------------------------
// CanAnnularFlowGeometry methods
// -------------------------------------------------------------

static constexpr const char *kCanAnnularDiameters =
    "CanAnnularFlowGeometry: D_outer > D_inner >= 0 required";

double CanAnnularFlowGeometry::D_mean() const {
    require_annulus_diameters(D_outer, D_inner, kCanAnnularDiameters);
    return (D_inner + D_outer) / 2.0;
}

double CanAnnularFlowGeometry::gap() const {
    require_annulus_diameters(D_outer, D_inner, kCanAnnularDiameters);
    return (D_outer - D_inner) / 2.0;
}

double CanAnnularFlowGeometry::area() const {
    require_annulus_diameters(D_outer, D_inner, kCanAnnularDiameters);
    return M_PI * (D_outer * D_outer - D_inner * D_inner) / 4.0;
}

double CanAnnularFlowGeometry::volume() const {
    require_finite_positive(L, "CanAnnularFlowGeometry: L must be positive");
    return area() * L;
}

double CanAnnularFlowGeometry::circumference() const {
    return M_PI * D_mean();
}

double CanAnnularFlowGeometry::area_primary() const {
    require_finite_positive(D_primary,
                            "CanAnnularFlowGeometry: D_primary must be positive");
    return M_PI * D_primary * D_primary / 4.0;
}

double CanAnnularFlowGeometry::volume_primary() const {
    require_finite_nonnegative(
        L_primary, "CanAnnularFlowGeometry: L_primary must be >= 0");
    return area_primary() * L_primary;
}

double CanAnnularFlowGeometry::volume_transition() const {
    require_finite_nonnegative(
        L_transition, "CanAnnularFlowGeometry: L_transition must be >= 0");
    return 0.5 * (area_primary() + area()) * L_transition;
}

double CanAnnularFlowGeometry::volume_total() const {
    return volume_primary() + volume_transition() + volume();
}

double CanAnnularFlowGeometry::length_total() const {
    require_finite_positive(L, "CanAnnularFlowGeometry: L must be positive");
    require_finite_nonnegative(
        L_primary, "CanAnnularFlowGeometry: L_primary must be >= 0");
    require_finite_nonnegative(
        L_transition, "CanAnnularFlowGeometry: L_transition must be >= 0");
    return L_primary + L_transition + L;
}

// -------------------------------------------------------------
// Residence Time
// -------------------------------------------------------------

double residence_time(double V, double Q) {
    if (!is_finite_positive(V) || !is_finite_positive(Q)) {
        throw std::invalid_argument("residence_time: V and Q must be positive");
    }
    return V / Q;
}

double residence_time_mdot(double V, double mdot, double rho) {
    if (!is_finite_positive(V) || !is_finite_positive(mdot) ||
        !is_finite_positive(rho)) {
        throw std::invalid_argument("residence_time_mdot: V, mdot, and rho must be positive");
    }
    return V * rho / mdot;
}

double residence_time(const Tube& tube, double Q) {
    return residence_time(tube.volume(), Q);
}

double residence_time(const Annulus& annulus, double Q) {
    return residence_time(annulus.volume(), Q);
}

double residence_time(const CanAnnularFlowGeometry& geom, double Q) {
    return residence_time(geom.volume_total(), Q);
}

double residence_time_mdot(const Tube& tube, double mdot, double rho) {
    return residence_time_mdot(tube.volume(), mdot, rho);
}

double residence_time_mdot(const Annulus& annulus, double mdot, double rho) {
    return residence_time_mdot(annulus.volume(), mdot, rho);
}

double residence_time_mdot(const CanAnnularFlowGeometry& geom, double mdot, double rho) {
    return residence_time_mdot(geom.volume_total(), mdot, rho);
}

double space_velocity(double Q, double V) {
    if (!is_finite_positive(Q) || !is_finite_positive(V)) {
        throw std::invalid_argument("space_velocity: Q and V must be positive");
    }
    return Q / V;
}

} // namespace combaero
