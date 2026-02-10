import pandas as pd
import math

# This script generates thermo_transport_data.h from species_data.csv to keep
# all thermo and transport data in the C++ thermo project consistent and typo-free.


def generate_cpp_header(csv_file, output_header):
    df = pd.read_csv(csv_file)

    header_content = """#ifndef THERMO_TRANSPORT_DATA_H
#define THERMO_TRANSPORT_DATA_H

#include <vector>
#include <string>
#include <unordered_map>
#include <limits>

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
    std::size_t C;
    std::size_t H;
    std::size_t O;
    std::size_t N;
};

const std::vector<std::string> species_names = {"""

    species_list: list[str] = []
    nasa_coeffs_list: list[str] = []
    transport_props_list: list[str] = []
    molar_masses_list: list[str] = []
    molecular_structures_list: list[str] = []

    def format_double(value):
        """Format a double value, converting NaN to C++ quiet_NaN()"""
        if pd.isna(value) or (isinstance(value, float) and math.isnan(value)):
            return "std::numeric_limits<double>::quiet_NaN()"
        return str(value)

    def get_int_or_zero(val):
        """Convert a possibly-NaN numeric value to int, defaulting to 0."""
        if pd.isna(val):
            return 0
        return int(val)

    def species_sort_key(row):
        """Return a tuple used to sort species into human-friendly groups.

        Groups (by primary key):
        0: air species [N2, O2, AR, CO2, H2O] in that order
        1: inert species
        2: hydrocarbons (C>0, H>0, O==0, N==0) ordered by C, then name
        3: other fuel species (C>0 but not hydrocarbons) alphabetical
        4: remaining species alphabetical
        """

        name = str(row["species_name"]).upper()

        air_order = ["N2", "O2", "AR", "CO2", "H2O"]
        if name in air_order:
            return (0, air_order.index(name), name)

        C = get_int_or_zero(row.get("C"))
        H = get_int_or_zero(row.get("H"))
        O_count = get_int_or_zero(row.get("O"))
        N = get_int_or_zero(row.get("N"))

        # Inert species: no C, H, O; allow N or others, but exclude air species
        if C == 0 and H == 0 and O_count == 0:
            return (1, name)

        # Hydrocarbons: C>0, H>0, no O or N
        if C > 0 and H > 0 and O_count == 0 and N == 0:
            return (2, C, name)

        # Other fuels: contain carbon but are not pure hydrocarbons
        if C > 0:
            return (3, name)

        # Everything else: group 4, alphabetical
        return (4, name)

    # Sort rows according to the custom species ordering rules
    sorted_rows = sorted((row for _, row in df.iterrows()), key=species_sort_key)

    for row in sorted_rows:
        species = str(row["species_name"])
        species_list.append(species)

        nasa_low = str(row["NASA_low"]).strip("[]").split(", ")
        nasa_high = str(row["NASA_high"]).strip("[]").split(", ")

        nasa_coeffs_list.append(
            f"{{{row['T_low']}, {row['T_mid']}, {row['T_high']}, "
            f"{{{', '.join(nasa_low)}}}, "
            f"{{{', '.join(nasa_high)}}}}}"
        )

        polarizability_str = format_double(row["polarizability"])
        transport_props_list.append(
            f"{{\"{row['geometry']}\", {row['well-depth']}, "
            f"{row['diameter']}, {polarizability_str}}}"
        )

        molar_masses_list.append(str(row["molar_mass"]))

        molecular_structures_list.append(
            f"{{{get_int_or_zero(row['C'])}, {get_int_or_zero(row['H'])}, "
            f"{get_int_or_zero(row['O'])}, {get_int_or_zero(row['N'])}}}"
        )

    header_content += ", ".join(f'"{s}"' for s in species_list) + "};\n\n"

    header_content += "const std::unordered_map<std::string, int> species_index = {\n"
    header_content += ",\n".join([f'    {{"{s}", {i}}}' for i, s in enumerate(species_list)])
    header_content += "\n};\n\n"

    header_content += (
        "const std::vector<double> molar_masses = {" + ", ".join(molar_masses_list) + "};\n\n"
    )

    header_content += (
        "const std::vector<NASA_Coeffs> nasa_coeffs = {\n"
        + ",\n".join(nasa_coeffs_list)
        + "\n};\n\n"
    )

    header_content += (
        "const std::vector<Transport_Props> transport_props = {\n"
        + ",\n".join(transport_props_list)
        + "\n};\n\n"
    )

    header_content += (
        "const std::vector<Molecular_Structure> molecular_structures = {\n"
        + ",\n".join(molecular_structures_list)
        + "\n};\n\n"
    )

    header_content += "#endif // THERMO_TRANSPORT_DATA_H" ""

    with open(output_header, "w") as f:
        f.write(header_content)

    print(f"Header file '{output_header}' generated successfully.")


if __name__ == "__main__":
    csv_file = "species_data.csv"
    output_header = "thermo_transport_data.h"
    generate_cpp_header(csv_file, output_header)
