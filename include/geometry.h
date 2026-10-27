#ifndef GEOMETRY_H
#define GEOMETRY_H

namespace combaero {

// -------------------------------------------------------------
// Geometric utilities for flow calculations
// -------------------------------------------------------------
// These functions compute geometric properties used across
// incompressible, compressible, and heat transfer calculations.

// Hydraulic diameter for arbitrary cross-section
// Dh = 4 * A / P_wetted
//
// Parameters:
//   A        : cross-sectional area [m^2]
//   P_wetted : wetted perimeter [m]
// Returns: hydraulic diameter [m]
double hydraulic_diameter(double A, double P_wetted);

// Hydraulic diameter for rectangular duct
// Dh = 2*a*b / (a + b)
//
// Parameters:
//   a, b : side lengths [m]
// Returns: hydraulic diameter [m]
double hydraulic_diameter_rect(double a, double b);

// Hydraulic diameter for annulus (concentric channels)
// Dh = D_outer - D_inner
//
// Parameters:
//   D_outer : outer diameter [m]
//   D_inner : inner diameter [m]
// Returns: hydraulic diameter [m]
double hydraulic_diameter_annulus(double D_outer, double D_inner);

// -------------------------------------------------------------
// Simple Geometry Helpers
// -------------------------------------------------------------
// Common geometric calculations for channel flow

// Circular channel cross-sectional area
// A = pi * (D/2)^2
//
// Parameters:
//   D : diameter [m]
// Returns: area [m^2]
double circular_area(double D);


// Annular cross-sectional area
// A = pi * ((D_outer/2)^2 - (D_inner/2)^2)
//
// Parameters:
//   D_outer : outer diameter [m]
//   D_inner : inner diameter [m]
// Returns: area [m^2]
double annular_area(double D_outer, double D_inner);

// Cylindrical channel volume
// V = pi * (D/2)^2 * L
//
// Parameters:
//   D : diameter [m]
//   L : length [m]
// Returns: volume [m^3]
double cylinder_volume(double D, double L);


// -------------------------------------------------------------
// Acoustic Geometry Primitives
// -------------------------------------------------------------
// Standard shapes for acoustic mode calculations.
// All dimensions in SI units (meters).

// Cylindrical tube (channel, duct, combustor liner)
struct Tube {
    double L;   // Length [m]
    double D;   // Diameter [m]

    double area() const;           // Cross-sectional area [m^2]
    double volume() const;         // Volume [m^3]
    double perimeter() const;      // Circumference [m]
};

// Annular duct (gas turbine combustor, annular gaps)
// Thin annulus approximation valid when (D_outer - D_inner) << D_mean
struct Annulus {
    double L;        // Length [m]
    double D_inner;  // Inner diameter [m]
    double D_outer;  // Outer diameter [m]

    // Acoustic solver compatibility
    int n_azimuthal_max = 5;  // Maximum azimuthal mode number to search

    // Diameter-based accessors (flow/geometry convention)
    double D_mean() const;         // Mean diameter [m]: (D_inner + D_outer) / 2
    double gap() const;            // Annular gap [m]: (D_outer - D_inner) / 2
    double area() const;           // Cross-sectional area [m^2]
    double volume() const;         // Volume [m^3]
    double circumference() const;  // Mean circumference [m]: pi * D_mean

    // Radius-based accessors (acoustic solver convention)
    double radius_inner() const { return D_inner / 2.0; }  // Inner radius [m]
    double radius_outer() const { return D_outer / 2.0; }  // Outer radius [m]
};

// Parameterized can-annular combustor geometry.
// Represents a circular primary zone, a circular-to-annular transition,
// and a downstream annular section.
//
// Notes:
// - Annular members (L, D_inner, D_outer) match Annulus naming/semantics.
// - Transition volume uses linear area interpolation:
//     V_transition = 0.5 * (A_primary + A_annulus) * L_transition
struct CanAnnularFlowGeometry {
    // Annular section (same members as Annulus)
    double L;        // Annular section length [m]
    double D_inner;  // Annular inner diameter [m]
    double D_outer;  // Annular outer diameter [m]

    // Dedicated can parameters
    double L_primary;     // Primary circular zone length [m]
    double D_primary;     // Primary circular zone diameter [m]
    double L_transition;  // Circular-to-annular transition length [m]

    // Annular-section helpers (same meaning as Annulus)
    double D_mean() const;         // Mean annular diameter [m]
    double gap() const;            // Annular gap [m]
    double area() const;           // Annular cross-sectional area [m^2]
    double volume() const;         // Annular section volume [m^3]
    double circumference() const;  // Mean annular circumference [m]

    // Can-specific helpers
    double area_primary() const;       // Primary circular area [m^2]
    double volume_primary() const;     // Primary circular volume [m^3]
    double volume_transition() const;  // Transition volume [m^3]
    double volume_total() const;       // Total volume [m^3]
    double length_total() const;       // Total length [m]
};

// -------------------------------------------------------------
// Residence Time
// -------------------------------------------------------------
// Time for fluid to pass through a volume.
// tau = V / Q_dot  [s]
//
// Applications:
// - Reactor design (Damkoehler number: Da = tau_res / tau_reaction)
// - Combustor sizing
// - Mixing time estimates

// Generic: tau = V / Q_dot
// V : volume [m^3]
// Q : volumetric flow rate [m^3/s]
double residence_time(double V, double Q);

// From mass flow: tau = V*rho / m_dot
// V    : volume [m^3]
// mdot : mass flow rate [kg/s]
// rho  : density [kg/m^3]
double residence_time_mdot(double V, double mdot, double rho);

// Convenience overloads for geometry structs
double residence_time(const Tube& tube, double Q);
double residence_time(const Annulus& annulus, double Q);
double residence_time(const CanAnnularFlowGeometry& geom, double Q);
double residence_time_mdot(const Tube& tube, double mdot, double rho);
double residence_time_mdot(const Annulus& annulus, double mdot, double rho);
double residence_time_mdot(const CanAnnularFlowGeometry& geom, double mdot, double rho);

// Space velocity (inverse of residence time)
// SV = Q_dot / V = 1/tau  [1/s]
double space_velocity(double Q, double V);

} // namespace combaero

#endif // GEOMETRY_H
