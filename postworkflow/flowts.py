import os
from ase import Atoms, io
from deepmd.calculator import DP
from ase.optimize import FIRE
from ase.neb import NEB
from ase.io import read, write
from ase.constraints import FixAtoms
from ase.neb import IDPP
import argparse

def neb_calculation(initial_structure, final_structure, model_path, n_images=5, fmax=0.02, interpolation='linear', spring_constant=1.0):
    # Read initial and final structures
    initial = read(initial_structure)
    final = read(final_structure)

    # Identify the atoms to be fixed (bottom 48 atoms)
    z_positions = initial.positions[:, 2]  # Get z-coordinates of atoms
    sorted_indices = z_positions.argsort()  # Sort indices by z-coordinate
    fixed_indices = sorted_indices[:48]  # Indices of the bottom 48 atoms

    # Generate intermediate images
    images = [initial]
    images += [initial.copy() for i in range(n_images)]
    images += [final]

    # Interpolate intermediate images
    neb = NEB(images, k=spring_constant, climb=True)
    
    if interpolation == 'idpp':
        neb.idpp_interpolate()
    else:
        neb.interpolate()

    # Output the initial interpolated structures as POSCAR
    for i, image in enumerate(images):
        write(f'POSCAR{str(i).zfill(2)}', image)

    # Set the calculator for each image
    for image in images:
        image.set_calculator(DP(model=model_path))

    # Constrain the bottom 48 atoms for each image
    for image in images:
        image.set_constraint(FixAtoms(indices=fixed_indices))

    # Run NEB calculation
    optimizer = FIRE(neb)
    optimizer.run(fmax=fmax)

    # Output the final NEB structures as CONTCAR
    for i, image in enumerate(images):
        write(f'CONTCAR{str(i).zfill(2)}', image)
    
    return images

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run NEB calculation using ASE and DeepMD.")
    parser.add_argument('initial_structure', type=str, help='Path to the initial structure (POSCAR format).')
    parser.add_argument('final_structure', type=str, help='Path to the final structure (POSCAR format).')
    parser.add_argument('model_path', type=str, help='Path to the DeepMD model.')
    parser.add_argument('--n_images', type=int, default=5, help='Number of intermediate images.')
    parser.add_argument('--fmax', type=float, default=1.1, help='Maximum force criteria for optimization.')
    parser.add_argument('--interpolation', choices=['linear', 'idpp'], default='linear', help='Interpolation method for generating intermediate images.')
    parser.add_argument('--spring_constant', type=float, default=8.0, help='Spring constant for NEB calculations.')

    args = parser.parse_args()
    
    # Run NEB calculation
    images = neb_calculation(
        args.initial_structure,
        args.final_structure,
        args.model_path,
        args.n_images,
        args.fmax,
        args.interpolation,
        args.spring_constant
    )
