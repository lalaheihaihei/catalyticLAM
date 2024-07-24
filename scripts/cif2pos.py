import os
import argparse
from pymatgen.core import Structure as Str

def convert_cif_to_poscar(input_dir, output_dir):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    for filename in os.listdir(input_dir):
        if filename.endswith(".cif"):
            cif_path = os.path.join(input_dir, filename)
            
            try:
                structure = Str.from_file(cif_path)
                
                poscar_filename = f"POSCAR-{filename.replace('.cif', '')}"
                poscar_path = os.path.join(output_dir, poscar_filename)
                
                structure.to(fmt="poscar", filename=poscar_path)
                
                print(f"Converted {filename} to {poscar_path}")
            
            except Exception as e:
                print(f"Failed to convert {filename}: {e}")

def main():
    parser = argparse.ArgumentParser(description='Convert CIF files to POSCAR files.')
    parser.add_argument('--input_dir', type=str, required=True, help='Path to the input directory containing CIF files.')
    parser.add_argument('--output_dir', type=str, default='POSCAR', help='Path to the output directory for POSCAR files.')
    
    args = parser.parse_args()
    
    convert_cif_to_poscar(args.input_dir, args.output_dir)

if __name__ == '__main__':
    main()
