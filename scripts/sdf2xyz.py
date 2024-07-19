from rdkit import Chem
from rdkit.Chem import AllChem
import os

def sdf_to_xyz(sdf_file, output_dir):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    suppl = Chem.SDMolSupplier(sdf_file, removeHs=False)
    valid_molecules = 0

    for idx, mol in enumerate(suppl):
        if mol is None:
            print(f"Warning: Molecule at index {idx} could not be parsed.")
            continue

        try:
            mol_name = f'mol-{idx + 1}'
            xyz_file = os.path.join(output_dir, f'{mol_name}.xyz')

            with open(xyz_file, 'w') as f:
                atoms = mol.GetAtoms()
                num_atoms = len(atoms)
                f.write(f"{num_atoms}\n")
                f.write(f"{mol_name}\n")
                for atom in atoms:
                    pos = mol.GetConformer().GetAtomPosition(atom.GetIdx())
                    f.write(f"{atom.GetSymbol()} {pos.x:.6f} {pos.y:.6f} {pos.z:.6f}\n")
            
            valid_molecules += 1
        except Exception as e:
            print(f"Error processing molecule at index {idx}: {e}")

    print(f"Converted {valid_molecules} valid molecules from {sdf_file} to XYZ format and saved in '{output_dir}' directory.")

sdf_file = 'gdb9.sdf'
output_dir = 'molxyz'
sdf_to_xyz(sdf_file, output_dir)
