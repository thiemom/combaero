import pandas as pd
import math
import os

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
    int C;
    int H;
    int O;
    int N;
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
    
    for _, row in df.iterrows():
        species = str(row['species_name'])
        species_list.append(species)
        
        nasa_low = str(row['NASA_low']).strip('[]').split(', ')
        nasa_high = str(row['NASA_high']).strip('[]').split(', ')
        
        nasa_coeffs_list.append(
            f"{{{row['T_low']}, {row['T_mid']}, {row['T_high']}, "
            f"{{{', '.join(nasa_low)}}}, "
            f"{{{', '.join(nasa_high)}}}}}"
        )
        
        polarizability_str = format_double(row['polarizability'])
        transport_props_list.append(
            f"{{\"{row['geometry']}\", {row['well-depth']}, {row['diameter']}, {polarizability_str}}}"
        )
        
        molar_masses_list.append(str(row['molar_mass']))

        molecular_structures_list.append(
            f"{{{int(row['C'])}, {int(row['H'])}, {int(row['O'])}, {int(row['N'])}}}"
        )
    
    header_content += ", ".join(f'"{s}"' for s in species_list) + "};\n\n"
    
    header_content += "const std::unordered_map<std::string, int> species_index = {\n"
    header_content += ",\n".join([f'    {{"{s}", {i}}}' for i, s in enumerate(species_list)])
    header_content += "\n};\n\n"
    
    header_content += "const std::vector<double> molar_masses = {" + ", ".join(molar_masses_list) + "};\n\n"
    
    header_content += "const std::vector<NASA_Coeffs> nasa_coeffs = {\n" + ",\n".join(nasa_coeffs_list) + "\n};\n\n"
    
    header_content += "const std::vector<Transport_Props> transport_props = {\n" + ",\n".join(transport_props_list) + "\n};\n\n"

    header_content += "const std::vector<Molecular_Structure> molecular_structures = {\n" + ",\n".join(molecular_structures_list) + "\n};\n\n"
    
    header_content += "#endif // THERMO_TRANSPORT_DATA_H"""
    
    with open(output_header, 'w') as f:
        f.write(header_content)
    
    print(f"Header file '{output_header}' generated successfully.")

if __name__ == "__main__":
    csv_file = "species_data.csv"
    output_header = "thermo_transport_data.h"
    generate_cpp_header(csv_file, output_header)
