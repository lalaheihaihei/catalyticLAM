import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
from ase import io
from fairchem.core import OCPCalculator

def tag_atoms(slab):
    # Set tags: Ru atoms to 1, others to 0
    tags = [1 if atom.symbol == 'Ru' else 0 for atom in slab]
    slab.set_tags(tags)
    return slab

def calculate_energies_forces(structure_paths, checkpoint_path):
    energies = []
    forces = []

    for path in structure_paths:
        atoms = io.read(path)
        atoms = tag_atoms(atoms)
        atoms.calc = OCPCalculator(checkpoint_path=checkpoint_path)

        energy = atoms.get_potential_energy()
        force = atoms.get_forces().max()

        energies.append(energy)
        forces.append(force)

    return energies, forces

def read_vasp_energy_forces(folder_path, n_images):
    energy = []
    forces = []
    structure_paths = []

    for i in range(n_images + 2):  # +2 for the initial and final images
        if i == 0:
            structure_paths.append(os.path.join(folder_path, "00", "POSCAR"))
        elif i == n_images + 1:
            structure_paths.append(os.path.join(folder_path, f"{str(i).zfill(2)}", "POSCAR"))
        else:
            structure_paths.append(os.path.join(folder_path, f"{str(i).zfill(2)}", "CONTCAR"))

        outcar_file = os.path.join(folder_path, f"{str(i).zfill(2)}", "OUTCAR")
        if os.path.exists(outcar_file):
            with open(outcar_file, 'r') as file:
                lines = file.readlines()
                # Extract the last free energy entry with 'free  energy   TOTEN'
                for line in reversed(lines):
                    if 'free  energy   TOTEN' in line:
                        energy_value = float(line.split()[4])  # Extract the 5th column
                        energy.append(energy_value)
                        break
                # Extract the last RMS force entry
                for line in reversed(lines):
                    if 'RMS' in line:
                        force_value = float(line.split()[5])  # Modify index based on file format
                        forces.append(force_value)
                        break
        else:
            print(f"Warning: OUTCAR file not found at {outcar_file}")

    return structure_paths, energy, forces

def plot_comparison(checkpoint_energies, checkpoint_forces, vasp_energy, vasp_forces):
    # Define colors for the checkpoints
    colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red']
    
    # Plot energies
    fig, ax1 = plt.subplots(figsize=(10, 6))

    for idx, energies in enumerate(checkpoint_energies):
        ax1.plot(range(len(energies)), energies, label=f'Checkpoint {idx + 1} Energies', marker='o', color=colors[idx])

    ax1.plot(range(len(vasp_energy)), vasp_energy, label='VASP Energies', marker='x', color='black', linewidth=2)
    ax1.set_xlabel('Image Index')
    ax1.set_ylabel('Energy (eV)')
    ax1.set_title('Energy Comparison')
    ax1.legend()
    fig.savefig('energy_comparison.png')
    plt.close(fig)

    # Plot forces
    fig, ax2 = plt.subplots(figsize=(10, 6))

    for idx, forces in enumerate(checkpoint_forces):
        ax2.plot(range(len(forces)), forces, label=f'Checkpoint {idx + 1} Forces', marker='o', color=colors[idx])

    ax2.plot(range(len(vasp_forces)), vasp_forces, label='VASP Forces', marker='x', color='black', linewidth=2)
    ax2.set_xlabel('Image Index')
    ax2.set_ylabel('Force (eV/Å)')
    ax2.set_title('Force Comparison')
    ax2.legend()
    fig.savefig('force_comparison.png')
    plt.close(fig)

def calculate_pointwise_errors(checkpoint_energies, vasp_raw_energy):
    mse_list = []
    rmse_list = []

    for i in range(len(vasp_raw_energy)):
        errors = [(energies[i] - vasp_raw_energy[i]) ** 2 for energies in checkpoint_energies]
        mse = np.mean(errors)
        rmse = np.sqrt(mse)
        mse_list.append(mse)
        rmse_list.append(rmse)

    return mse_list, rmse_list

def main():
    parser = argparse.ArgumentParser(description="Calculate energies and forces with multiple checkpoints, plot results, and compare with VASP results.")
    parser.add_argument('checkpoint_paths', type=str, nargs=4, help='Paths to the four OCP model checkpoints.')
    parser.add_argument('--vasp_results_path', type=str, required=True, help='Path to the VASP results folder.')
    parser.add_argument('--n_images', type=int, default=4, help='Number of intermediate images.')

    args = parser.parse_args()

    # Read VASP results
    structure_paths, vasp_raw_energy, vasp_forces = read_vasp_energy_forces(args.vasp_results_path, args.n_images)

    # Normalize VASP energies to the first frame
    vasp_energy = [e - vasp_raw_energy[0] for e in vasp_raw_energy]

    # Initialize lists to store results
    checkpoint_energies = []
    checkpoint_forces = []
    checkpoint_energies_raw = []

    # Calculate first frame energies for normalization
    first_frame_energies = []

    for checkpoint_path in args.checkpoint_paths:
        energies, _ = calculate_energies_forces([structure_paths[0]], checkpoint_path)
        first_frame_energies.append(energies[0])

    avg_first_frame_energy = sum(first_frame_energies) / len(first_frame_energies)

    # Calculate energies and forces for each checkpoint
    for checkpoint_path in args.checkpoint_paths:
        energies, forces = calculate_energies_forces(structure_paths, checkpoint_path)
        checkpoint_energies_raw.append(energies.copy())  # Store raw energies for error calculation
        # Normalize energies to the average first frame energy
        energies = [e - avg_first_frame_energy for e in energies]
        checkpoint_energies.append(energies)
        checkpoint_forces.append(forces)

    # Plot comparison
    plot_comparison(checkpoint_energies, checkpoint_forces, vasp_energy, vasp_forces)

    # Calculate pointwise errors
    mse_list, rmse_list = calculate_pointwise_errors(checkpoint_energies_raw, vasp_raw_energy)
    print("Pointwise MSE:", mse_list)
    print("Pointwise RMSE:", rmse_list)

if __name__ == "__main__":
    main()

