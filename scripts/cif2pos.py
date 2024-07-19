import os
from pymatgen.core import Structure as Str

input_folder = "./"

for i in os.listdir(input_folder):
    if i.endswith(".cif"):
        cif_path = os.path.join(input_folder, i)
        
        structure = Str.from_file(cif_path)
        
        poscar_filename = f"POSCAR-{i.replace('.cif', '')}"
        poscar_path = os.path.join(input_folder, poscar_filename)
        
        structure.to(fmt="poscar", filename=poscar_path)
        
        print(f"Converted {i} to {poscar_path}")

print("Conversion complete.")
