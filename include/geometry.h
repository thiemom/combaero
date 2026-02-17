#ifndef GEOMETRY_H
#define GEOMETRY_H

// -------------------------------------------------------------
// Geometric utilities for flow calculations
// -------------------------------------------------------------
// These functions compute geometric properties used across
// incompressible, compressible, and heat transfer calculations.

// Hydraulic diameter for arbitrary cross-section
// Dh = 4 * A / P_wetted
//
// Parameters:
//   A        : cross-sectional area [m²]
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

// Hydraulic diameter for annulus (concentric pipes)
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
// Common geometric calculations for pipe flow

// Circular pipe cross-sectional area
// A = π * (D/2)²
//
// Parameters:
//   D : diameter [m]
// Returns: area [m²]
double pipe_area(double D);

// Annular cross-sectional area
// A = π * ((D_outer/2)² - (D_inner/2)²)
//
// Parameters:
//   D_outer : outer diameter [m]
//   D_inner : inner diameter [m]
// Returns: area [m²]
double annular_area(double D_outer, double D_inner);

// Cylindrical pipe volume
// V = π * (D/2)² * L
//
// Parameters:
//   D : diameter [m]
//   L : length [m]
// Returns: volume [m³]
double pipe_volume(double D, double L);

// -------------------------------------------------------------
// Acoustic Geometry Primitives
// -------------------------------------------------------------
// Standard shapes for acoustic mode calculations.
// All dimensions in SI units (meters).

// Cylindrical tube (pipe, duct, combustor liner)
struct Tube {
    double L;   // Length [m]
    double D;   // Diameter [m]

    double area() const;           // Cross-sectional area [m²]
    double volume() const;         // Volume [m³]
    double perimeter() const;      // Circumference [m]
};

// Annular duct (gas turbine combustor, annular gaps)
// Thin annulus approximation valid when (D_outer - D_inner) << D_mean
struct Annulus {
    double L;        // Length [m]
    double D_inner;  // Inner diameter [m]
    double D_outer;  // Outer diameter [m]

    double D_mean() const;         // Mean diameter [m]: (D_inner + D_outer) / 2
    double gap() const;            // Annular gap [m]: (D_outer - D_inner) / 2
    double area() const;           // Cross-sectional area [m²]
    double volume() const;         // Volume [m³]
    double circumference() const;  // Mean circumference [m]: π * D_mean
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
    double area() const;           // Annular cross-sectional area [m²]
    double volume() const;         // Annular section volume [m³]
    double circumference() const;  // Mean annular circumference [m]

    // Can-specific helpers
    double area_primary() const;       // Primary circular area [m²]
    double volume_primary() const;     // Primary circular volume [m³]
    double volume_transition() const;  // Transition volume [m³]
    double volume_total() const;       // Total volume [m³]
    double length_total() const;       // Total length [m]
};

// -------------------------------------------------------------
// Residence Time
// -------------------------------------------------------------
// Time for fluid to pass through a volume.
// τ = V / Q̇  [s]
//
// Applications:
// - Reactor design (Damköhler number: Da = τ_res / τ_reaction)
// - Combustor sizing
// - Mixing time estimates

// Generic: τ = V / Q̇
// V : volume [m³]
// Q : volumetric flow rate [m³/s]
double residence_time(double V, double Q);

// From mass flow: τ = V·ρ / ṁ
// V    : volume [m³]
// mdot : mass flow rate [kg/s]
// rho  : density [kg/m³]
double residence_time_mdot(double V, double mdot, double rho);

// Convenience overloads for geometry structs
double residence_time(const Tube& tube, double Q);
double residence_time(const Annulus& annulus, double Q);
double residence_time(const CanAnnularFlowGeometry& geom, double Q);
double residence_time_mdot(const Tube& tube, double mdot, double rho);
double residence_time_mdot(const Annulus& annulus, double mdot, double rho);
double residence_time_mdot(const CanAnnularFlowGeometry& geom, double mdot, double rho);

// Space velocity (inverse of residence time)
// SV = Q̇ / V = 1/τ  [1/s]
double space_velocity(double Q, double V);

#endif // GEOMETRY_H
