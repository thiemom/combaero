import yaml
import csv

def extract_species_data(yaml_file, output_file, selected_species):
    with open(yaml_file, 'r') as f:
        data = yaml.safe_load(f)
    
    species_data = []
    headers = [
        "species_name", "molar_mass", "P_ref", "C", "H", "O", "N", "T_low", "T_mid", "T_high", 
        "NASA_low", "NASA_high", "model", "geometry", "well-depth", "diameter", "polarizability", "rotational-relaxation",
    ]
    extra_headers = set()

    # Normalize selected species names to uppercase so matching is case-insensitive
    selected_species_upper = {name.upper() for name in selected_species}
    
    for species in data.get("species", []):
        species_name = species["name"]
        if species_name.upper() in selected_species_upper:
            comp = species.get("composition", {})
            thermo = species.get("thermo", {})
            transport = species.get("transport", {})
            
            entry = {
                "species_name": species_name.upper(),
                "molar_mass": sum(comp.get(el, 0) * get_atomic_mass(el) for el in comp),
                "P_ref": data.get("phases", [{}])[0].get("state", {}).get("P", ""),
                "C": comp.get("C", 0),
                "H": comp.get("H", 0),
                "O": comp.get("O", 0),
                "N": comp.get("N", 0),
                "T_low": thermo.get("temperature-ranges", [None, None, None])[0],
                "T_mid": thermo.get("temperature-ranges", [None, None, None])[1],
                "T_high": thermo.get("temperature-ranges", [None, None, None])[2],
                "NASA_low": thermo.get("data", [[], []])[0],
                "NASA_high": thermo.get("data", [[], []])[1],
                "model": thermo.get("model", ""),
                "geometry": transport.get("geometry", ""),
                "well-depth": transport.get("well-depth", ""),
                "diameter": transport.get("diameter", ""),
                "polarizability": transport.get("polarizability", "")
            }
            
            for key, value in species.items():
                if key not in entry:
                    entry[key] = value
                    extra_headers.add(key)
            
            species_data.append(entry)
    
    headers.extend(sorted(extra_headers))
    
    with open(output_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=headers)
        writer.writeheader()
        writer.writerows(species_data)

def get_atomic_mass(element):
    atomic_masses = {
        "C": 12.011,
        "H": 1.008,
        "O": 16.00,
        "N": 14.007,
        "AR": 39.948,
        "HE": 4.002602,
    }

    key = element.upper()
    return atomic_masses.get(key, 0)

if __name__ == "__main__":
    mechanism_file = "JetSurf2.yaml"  # Input mechanism file in yaml format
    output_file = "species_data.csv" 
    selected_species = ["O2", "N2", "AR", "CO2", "H2O", "CH4", "C2H6", "C3H8", "iC4H10", "NC4H10", "iC5H12", "NC5H12", "NC6H14", "NC7H16", "H2", "CO"]  # User-defined species
    extract_species_data(mechanism_file, output_file, selected_species)
    print(f"Data extracted to {output_file}")
