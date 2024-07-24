import json
import numpy as np
from ase import Atoms
from ase.io import write
import os
import re
from collections import defaultdict
import argparse

def decode_ndarray(obj):
    """Custom decoder for handling numpy arrays stored in JSON."""
    if isinstance(obj, dict) and '__ndarray__' in obj:
        shape, dtype, data = obj['__ndarray__']
        return np.array(data, dtype=dtype).reshape(shape)
    return obj

def json_line_to_atoms(json_line):
    """Convert a single JSON line to an ASE Atoms object and return the formula_pretty."""
    entry = json.loads(json_line, object_hook=decode_ndarray)

    formula_pretty = entry.get("formula_pretty", "unknown")
    
    if "structure" in entry:
        structure = entry["structure"]
        
        cell = structure["lattice"]["matrix"]
        
        positions = [site["xyz"] for site in structure["sites"]]
        symbols = [site["label"] for site in structure["sites"]]

        atoms = Atoms(symbols=symbols, positions=positions, cell=cell, pbc=True)
        return atoms, formula_pretty
    
    return None, formula_pretty

def normalize_formula(formula):
    def parse_formula(formula):
        pattern = re.compile(r'([A-Za-z][a-z]?)(\d*)')
        stack = []
        current = []
        i = 0
        while i < len(formula):
            if formula[i] == '(':
                stack.append(current)
                current = []
                i += 1
            elif formula[i] == ')':
                i += 1
                count = ''
                while i < len(formula) and formula[i].isdigit():
                    count += formula[i]
                    i += 1
                count = int(count or 1)
                previous = stack.pop()
                for elem, num in current:
                    previous.append((elem, num * count))
                current = previous
            else:
                match = pattern.match(formula[i:])
                if match:
                    element = match.group(1)
                    count = int(match.group(2) or 1)
                    current.append((element, count))
                    i += len(match.group(0))
                else:
                    i += 1
        return current

    parsed = parse_formula(formula)
    element_count = defaultdict(int)
    for element, count in parsed:
        element_count[element] += count
    
    order = []
    seen = set()
    for element, _ in parsed:
        if element not in seen:
            seen.add(element)
            order.append(element)
    
    return ''.join(f"{el}{element_count[el] if element_count[el] > 1 else ''}" for el in order)

def get_unique_filename(output_dir, base_name, count):
    if count == 0:
        return os.path.join(output_dir, f"{base_name}.cif")
    else:
        return os.path.join(output_dir, f"{base_name}-{count + 1}.cif")

def convert_json_lines_to_cifs(json_file, output_dir):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    formula_count = defaultdict(int)

    with open(json_file, 'r') as file:
        lines = file.readlines()
    
    for i, line in enumerate(lines):
        try:
            atoms, formula_pretty = json_line_to_atoms(line)
            if atoms is not None:
                safe_formula_pretty = normalize_formula(formula_pretty)
                
                count = formula_count[safe_formula_pretty]
                
                cif_file = get_unique_filename(output_dir, safe_formula_pretty, count)
                
                write(cif_file, atoms)
                
                formula_count[safe_formula_pretty] += 1
            else:
                print(f"No valid structure found in line {i + 1}.")
        except json.JSONDecodeError as e:
            print(f"JSON decoding error in line {i + 1}: {e}")

def main():
    parser = argparse.ArgumentParser(description='Convert JSON lines to CIF files.')
    parser.add_argument('--file', type=str, required=True, help='Path to the JSON file.')
    parser.add_argument('--outdir', type=str, default='CIF', help='Output directory for CIF files.')
    
    args = parser.parse_args()
    
    convert_json_lines_to_cifs(args.file, args.outdir)

if __name__ == '__main__':
    main()
