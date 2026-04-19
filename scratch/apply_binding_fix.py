import sys

path = 'python/combaero/_core.cpp'
with open(path, 'r') as f:
    lines = f.readlines()

# Targeting line 3931 (0-indexed: 3930)
target_idx = 3930
if 'm.def("Cd_orifice", &Cd' in lines[target_idx]:
    new_block = [
        '  m.def(\n',
        '      "Cd_orifice",\n',
        '      [](const OrificeGeometry &geom, const OrificeState &state,\n',
        '         CdCorrelation correlation) {\n',
        '        auto corr = make_correlation(correlation);\n',
        '        return corr ? corr->Cd(geom, state) : Cd(geom, state);\n',
        '      },\n',
        '      py::arg("geom"), py::arg("state"),\n',
        '      py::arg("correlation") = CdCorrelation::ReaderHarrisGallagher,\n'
    ]
    # Replace lines 3931 only (keeping the docstring which starts at 3932)
    lines[target_idx] = "".join(new_block)

    with open(path, 'w') as f:
        f.writelines(lines)
    print("Success")
else:
    print(f"Failed: line {target_idx+1} does not match expected content: {lines[target_idx].strip()}")
