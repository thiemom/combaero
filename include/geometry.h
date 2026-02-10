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
double residence_time_mdot(const Tube& tube, double mdot, double rho);
double residence_time_mdot(const Annulus& annulus, double mdot, double rho);

// Space velocity (inverse of residence time)
// SV = Q̇ / V = 1/τ  [1/s]
double space_velocity(double Q, double V);

#endif // GEOMETRY_H
