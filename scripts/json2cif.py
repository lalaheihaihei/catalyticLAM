import json
import numpy as np
from ase import Atoms
from ase.io import write

def decode_ndarray(obj):
    if isinstance(obj, dict) and '__ndarray__' in obj:
        shape, dtype, data = obj['__ndarray__']
        return np.array(data, dtype=dtype).reshape(shape)
    return obj

def json_line_to_atoms(json_line):
    entry = json.loads(json_line, object_hook=decode_ndarray)

    if "structure" in entry:
        structure = entry["structure"]
        
        cell = structure["lattice"]["matrix"]
        
        positions = []
        symbols = []
        for site in structure["sites"]:
            positions.append(site["xyz"])
            symbols.append(site["label"])

        atoms = Atoms(symbols=symbols, positions=positions, cell=cell, pbc=True)
        return atoms
    
    return None

def convert_json_lines_to_cifs(json_file, output_dir):
    with open(json_file, 'r') as file:
        lines = file.readlines()
    
    for i, line in enumerate(lines):
        try:
            atoms = json_line_to_atoms(line)
            if atoms is not None:
                cif_file = f"{output_dir}/structure-{i + 1}.cif"
                write(cif_file, atoms)
            else:
                print(f"No valid structure found in line {i + 1}.")
        except json.JSONDecodeError as e:
            print(f"JSON decoding error in line {i + 1}: {e}")

json_file = 'db.json'
output_dir = 'cif'
import os
if not os.path.exists(output_dir):
    os.makedirs(output_dir)
convert_json_lines_to_cifs(json_file, output_dir)
