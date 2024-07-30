import os
import argparse
import matplotlib.pyplot as plt
from ase import io
from ase.mep import NEB, NEBTools
from ase.optimize import FIRE
from ase.constraints import FixAtoms
from fairchem.core import OCPCalculator

def tag_atoms(slab):
    # Set tags: Ru atoms to 1, others to 0
    tags = [1 if atom.symbol == 'Ru' else 0 for atom in slab]
    slab.set_tags(tags)
    print(slab.get_tags())
    return slab

def neb_calculation(initial_structure, final_structure, checkpoint_path, n_images=5, fmax=0.02, interpolation='linear', spring_constant=1.0, apply_constraint=False):
    # Read initial and final structures
    initial = io.read(initial_structure)
    final = io.read(final_structure)

    # Identify the atoms to be fixed (bottom 48 atoms)
    z_positions = initial.positions[:, 2]  # Get z-coordinates of atoms
    sorted_indices = z_positions.argsort()  # Sort indices by z-coordinate
    fixed_indices = sorted_indices[:48]  # Indices of the bottom 48 atoms

    # Tag atoms in the initial and final structures
    initial = tag_atoms(initial)
    final = tag_atoms(final)

    # Generate intermediate images
    images = [initial]
    images += [initial.copy() for _ in range(n_images)]
    images += [final]

    # Interpolate intermediate images
    neb = NEB(images, k=spring_constant, climb=True)
    
    if interpolation == 'idpp':
        print("IDPP interpolation is not available in the current ASE version.")
    else:
        neb.interpolate()

    # Set the calculator for each image
    for image in images:
        image.calc = OCPCalculator(checkpoint_path=checkpoint_path)

    # Constrain the bottom 48 atoms for each image
    for image in images:
        image.set_constraint(FixAtoms(indices=fixed_indices))

    # Run NEB calculation
    optimizer = FIRE(neb, trajectory='neb.traj')
    optimizer.run(fmax=fmax)

    return images

def read_vasp_energy(folder_path, n_images):
    energy = []
    specific_energy = []  # To store energy from the first frame
    for i in range(n_images + 2):  # +2 for the initial and final images
        folder = os.path.join(folder_path, f"{str(i).zfill(2)}")
        outcar_file = os.path.join(folder, "OUTCAR")
        if os.path.exists(outcar_file):
            with open(outcar_file, 'r') as file:
                lines = file.readlines()
                # Extract the last free energy entry with 'free  energy   TOTEN'
                for line in reversed(lines):
                    if 'free  energy   TOTEN' in line:
                        energy_value = float(line.split()[4])  # Extract the 5th column
                        energy.append(energy_value)
                        break
                for line in lines:
                    if 'free  energy   TOTEN' in line:
                        energy_value = float(line.split()[4])  # Extract the 5th column
                        specific_energy.append(energy_value)
                        break

        else:
            print(f"Warning: OUTCAR file not found at {outcar_file}")
    # Set the first frame energy to 0
    if energy:
        energy = [e - energy[0] for e in energy]
        specific_energy = [e - specific_energy[0] for e in specific_energy]  # Normalize specific energy
    return energy, specific_energy

def read_vasp_forces(folder_path, n_images):
    forces = []
    specific_forces = []  # To store forces from the first frame
    for i in range(n_images + 2):  # +2 for the initial and final images
        folder = os.path.join(folder_path, f"{str(i).zfill(2)}")
        outcar_file = os.path.join(folder, "OUTCAR")
        if os.path.exists(outcar_file):
            with open(outcar_file, 'r') as file:
                lines = file.readlines()
                # Extract the last RMS force entry
                for line in reversed(lines):
                    if 'RMS' in line:
                        force_value = float(line.split()[5])  # Modify index based on file format
                        forces.append(force_value)
                        # Collect force from the first frame
                        if i == 0:
                            specific_forces.append(force_value)
                        break
        else:
            print(f"Warning: OUTCAR file not found at {outcar_file}")
    return forces, specific_forces

def plot_comparison(ase_images, vasp_energy, vasp_forces, vasp_specific_energy, vasp_specific_forces):
    # Extract ASE energies
    ase_energies = [image.get_potential_energy() for image in ase_images]
    
    # Extract ASE forces
    ase_forces = [image.get_forces().max() for image in ase_images]

    # Plot NEB path
    fig1, ax1 = plt.subplots(figsize=(10, 6))
    nebtools = NEBTools(ase_images)
    nebtools.plot_band(ax1)
    ax1.set_xlabel('Image Index')
    ax1.set_ylabel('Energy (eV)')
    ax1.set_title('NEB Path')
    
    # Plot VASP energies for comparison
    ax1.plot(range(len(vasp_energy)), vasp_energy, label='VASP Energies final', marker='x')
    ax1.plot(range(len(vasp_specific_energy)), vasp_specific_energy, label='VASP Energies initial', marker='x')
    ax1.legend()
    fig1.savefig('neb_path.png')  # Save the NEB path plot as a PNG file
    plt.close(fig1)

    # Plot forces
    fig2, ax2 = plt.subplots(figsize=(10, 6))
    ax2.plot(range(len(ase_forces)), ase_forces, label='ASE Forces', marker='o')
    ax2.plot(range(len(vasp_forces)), vasp_forces, label='VASP Forces', marker='x')
    ax2.set_xlabel('Image Index')
    ax2.set_ylabel('Force (eV/Å)')
    ax2.set_title('NEB Force Comparison')
    ax2.legend()
    fig2.savefig('neb_forces.png')  # Save the forces comparison plot as a PNG file
    plt.close(fig2)

def main():
    parser = argparse.ArgumentParser(description="Run NEB calculation, plot results, and compare with VASP results.")
    parser.add_argument('initial_structure', type=str, help='Path to the initial structure (POSCAR format).')
    parser.add_argument('final_structure', type=str, help='Path to the final structure (POSCAR format).')
    parser.add_argument('checkpoint_path', type=str, help='Path to the OCP model checkpoint.')
    parser.add_argument('--n_images', type=int, default=4, help='Number of intermediate images.')
    parser.add_argument('--fmax', type=float, default=1.0, help='Maximum force criteria for optimization.')
    parser.add_argument('--interpolation', choices=['linear'], default='linear', help='Interpolation method for generating intermediate images.')
    parser.add_argument('--spring_constant', type=float, default=1.0, help='Spring constant for NEB calculations.')
    parser.add_argument('--apply_constraint', type=bool, default=False, help='Whether to apply constraints during interpolation.')
    parser.add_argument('--vasp_results_path', type=str, required=True, help='Path to the VASP results folder.')

    args = parser.parse_args()

    # Calculate NEB
    images = neb_calculation(
        args.initial_structure,
        args.final_structure,
        args.checkpoint_path,
        n_images=args.n_images,
        fmax=args.fmax,
        interpolation=args.interpolation,
        spring_constant=args.spring_constant,
        apply_constraint=args.apply_constraint
    )

    # Read VASP results
    vasp_energy, vasp_specific_energy = read_vasp_energy(args.vasp_results_path, args.n_images)
    vasp_forces, vasp_specific_forces = read_vasp_forces(args.vasp_results_path, args.n_images)

    # Plot comparison
    plot_comparison(images, vasp_energy, vasp_forces, vasp_specific_energy, vasp_specific_forces)

if __name__ == "__main__":
    main()
