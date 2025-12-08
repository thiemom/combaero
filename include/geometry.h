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

#endif // GEOMETRY_H
