import os
import numpy as np

def read_xyz(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
        num_atoms = int(lines[0].strip())
        atomic_symbols = []
        coordinates = []
        for line in lines[2:2 + num_atoms]:
            parts = line.split()
            atomic_symbols.append(parts[0])
            coordinates.append([float(x) for x in parts[1:4]])
        return atomic_symbols, np.array(coordinates)

def construct_unit_cell(coordinates, padding=5.0):
    min_coords = np.min(coordinates, axis=0) - padding
    max_coords = np.max(coordinates, axis=0) + padding
    cell_vectors = np.diag(max_coords - min_coords)
    return min_coords, cell_vectors

def write_poscar(file_path, atomic_symbols, coordinates, min_coords, cell_vectors):
    unique_symbols = sorted(set(atomic_symbols))
    num_atoms_per_type = [atomic_symbols.count(symbol) for symbol in unique_symbols]
    
    with open(file_path, 'w') as file:
        file.write("Generated POSCAR\n")
        file.write("1.0\n")
        for vector in cell_vectors:
            file.write(f" {vector[0]} {vector[1]} {vector[2]}\n")
        
        file.write(" ".join(unique_symbols) + "\n")
        file.write(" ".join(map(str, num_atoms_per_type)) + "\n")
        file.write("Cartesian\n")
        
        for symbol in unique_symbols:
            for i, atom_symbol in enumerate(atomic_symbols):
                if atom_symbol == symbol:
                    shifted_coord = coordinates[i] - min_coords
                    file.write(f" {shifted_coord[0]} {shifted_coord[1]} {shifted_coord[2]}\n")

def convert_xyz_to_poscar(input_dir, output_dir, padding=5.0):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    for file_name in os.listdir(input_dir):
        if file_name.endswith('.xyz'):
            input_file_path = os.path.join(input_dir, file_name)
            
            base_name = os.path.splitext(file_name)[0]
            
            output_file_name = f"POSCAR-{base_name}"
            output_file_path = os.path.join(output_dir, output_file_name)
            
            atomic_symbols, coordinates = read_xyz(input_file_path)
            min_coords, cell_vectors = construct_unit_cell(coordinates, padding)
            write_poscar(output_file_path, atomic_symbols, coordinates, min_coords, cell_vectors)

input_directory = './'
output_directory = './'
convert_xyz_to_poscar(input_directory, output_directory)
