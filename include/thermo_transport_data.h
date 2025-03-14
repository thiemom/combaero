#ifndef THERMO_TRANSPORT_DATA_H
#define THERMO_TRANSPORT_DATA_H

#include <vector>
#include <string>
#include <unordered_map>

struct NASA_Coeffs {
    double T_low;
    double T_mid;
    double T_high;
    std::vector<double> low_coeffs;
    std::vector<double> high_coeffs;
};

struct Transport_Props {
    std::string geometry;
    double well_depth;
    double diameter;
    double polarizability;
};

struct Molecular_Structure {
    int C;  // Number of carbon atoms
    int H;  // Number of hydrogen atoms
    int O;  // Number of oxygen atoms
    int N;  // Number of nitrogen atoms
};

const std::vector<std::string> species_names = {"H2", "O2", "CH4", "N2", "Ar", "CO2", "H2O"};

const std::unordered_map<std::string, int> species_index = {
    {"H2", 0},
    {"O2", 1},
    {"CH4", 2},
    {"N2", 3},
    {"Ar", 4},
    {"CO2", 5},
    {"H2O", 6}
};

const std::vector<double> molar_masses = {2.016, 32.0, 16.043, 28.014, 39.948, 44.01, 18.015};

const std::vector<NASA_Coeffs> nasa_coeffs = {
// H2
{200.0, 1000.0, 6000.0, {2.34433112, 0.00798052075, -1.9478151e-05, 2.01572094e-08, -7.37611761e-12, -917.935173, 0.683010238}, {2.93286579, 0.000826607967, -1.46402335e-07, 1.54100359e-11, -6.88804432e-16, -813.065597, -1.02432887}},
// O2
{200.0, 1000.0, 6000.0, {3.78245636, -0.00299673415, 9.847302e-06, -9.68129508e-09, 3.24372836e-12, -1063.94356, 3.65767573}, {3.66096083, 0.000656365523, -1.41149485e-07, 2.05797658e-11, -1.29913248e-15, -1215.97725, 3.41536184}},
// CH4
{200.0, 1000.0, 6000.0, {5.14987613, -0.0136709788, 4.91800599e-05, -4.84743026e-08, 1.66693956e-11, -10246.6476, -4.64130376}, {1.63552643, 0.0100842795, -3.36916254e-06, 5.34958667e-10, -3.15518833e-14, -10005.6455, 9.99313326}},
// N2
{300.0, 1000.0, 5000.0, {3.298677, 0.0014082404, -3.963222e-06, 5.641515e-09, -2.444854e-12, -1020.8999, 3.950372}, {2.92664, 0.0014879768, -5.68476e-07, 1.0097038e-10, -6.753351e-15, -922.7977, 5.980528}},
// Ar
{300.0, 1000.0, 5000.0, {2.5, 0.0, 0.0, 0.0, 0.0, -745.375, 4.366}, {2.5, 0.0, 0.0, 0.0, 0.0, -745.375, 4.366}},
// CO2
{200.0, 1000.0, 6000.0, {2.35677352, 0.00898459677, -7.12356269e-06, 2.45919022e-09, -1.43699548e-13, -48371.9697, 9.90105222}, {4.63659493, 0.00274131991, -9.95828531e-07, 1.60373011e-10, -9.16103468e-15, -49024.9341, -1.93534855}},
// H2O
{200.0, 1000.0, 6000.0, {4.19864056, -0.00203643410, 6.52040211e-06, -5.48797062e-09, 1.77197817e-12, -30293.7267, -0.849032208}, {2.67703787, 0.00297318329, -7.73769690e-07, 9.44336689e-11, -4.26900959e-15, -29885.8938, 6.88255571}}
};

const std::vector<Transport_Props> transport_props = {
{"linear", 38.0, 2.92, 0.79},      // H2
{"linear", 107.4, 3.46, 1.6},      // O2
{"nonlinear", 141.4, 3.75, 2.6},   // CH4
{"linear", 97.53, 3.62, 1.76},     // N2
{"atom", 136.5, 3.33, 1.64},       // Ar
{"linear", 244.0, 3.76, 2.65},     // CO2
{"nonlinear", 572.4, 2.605, 1.84}  // H2O
};

// Molecular structure data for each species (H2, O2, CH4, N2, Ar, CO2, H2O)
const std::vector<Molecular_Structure> molecular_structures = {
{0, 2, 0, 0},  // H2:  H=2
{0, 0, 2, 0},  // O2:  O=2
{1, 4, 0, 0},  // CH4: C=1, H=4
{0, 0, 0, 2},  // N2:  N=2
{0, 0, 0, 0},  // Ar:  (noble gas)
{1, 0, 2, 0},  // CO2: C=1, O=2
{0, 2, 1, 0}   // H2O: H=2, O=1
};

#endif // THERMO_TRANSPORT_DATA_H