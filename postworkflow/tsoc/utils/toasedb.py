import argparse
import os
from ase.db import connect
import dpdata as dp

def convert_outcar_to_ase_db(outcar_path, db_path):
    """
    Convert OUTCAR file to ASE database format.
    
    Parameters:
    outcar_path (str): Path to the OUTCAR file.
    db_path (str): Path where the ASE database will be saved.
    """
    # Ensure the database directory exists
    db_dir = os.path.dirname(db_path)
    if not os.path.exists(db_dir):
        os.makedirs(db_dir)

    # Load the data from OUTCAR
    try:
        data = dp.LabeledSystem(outcar_path, fmt='vasp/outcar')
    except Exception as e:
        print(f"Error loading OUTCAR file: {e}")
        return

    # Open the ASE database for writing
    db = connect(db_path)
    
    # Iterate over the data and add each frame to the database
    for i, structure in enumerate(data):
        # Convert dpdata structure to ASE Atoms object(s)
        ase_structures = structure.to_ase_structure()
        
        # Handle case where to_ase_structure() returns a list of Atoms objects
        if isinstance(ase_structures, list):
            for j, ase_atoms in enumerate(ase_structures):
                # Add each Atoms object to the database
                db.write(ase_atoms, key_value_pairs={"step": i, "frame": j})
        else:
            # If it's a single Atoms object, write it to the database
            db.write(ase_structures, key_value_pairs={"step": i})

    print(f"Data successfully written to {db_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert OUTCAR file to ASE database format.")
    parser.add_argument("outcar_path", type=str, help="Path to the OUTCAR file.")
    parser.add_argument("db_path", type=str, help="Path where the ASE database will be saved.")

    args = parser.parse_args()

    convert_outcar_to_ase_db(args.outcar_path, args.db_path)
