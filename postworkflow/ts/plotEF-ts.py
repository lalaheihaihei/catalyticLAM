import os
import matplotlib.pyplot as plt

def extract_energies(oszicar_path):
    energies = []
    if os.path.exists(oszicar_path):
        with open(oszicar_path, 'r') as file:
            for line in file:
                if 'E0=' in line:
                    energy = float(line.split()[4])
                    energies.append(energy)
    return energies

def extract_forces(outcar_path):
    forces = []
    if os.path.exists(outcar_path):
        with open(outcar_path, 'r') as file:
            for line in file:
                if 'RMS' in line:
                    force = float(line.split()[4])
                    forces.append(force)
    return forces

def main():
    num_iterations = 3  # Adjust this to the number of ts iterations you have
    n_images = 4  # Adjust this to the number of images you have

    all_energies_list = []
    all_forces_list = []
    all_comparison_energies_list = []
    all_comparison_forces_list = []

    # Process each image from 01 to n_images
    for j in range(1, n_images + 1):
        all_energies = []
        all_forces = []

        for i in range(1, num_iterations + 1):
            folder = f'ts{i}/{str(j).zfill(2)}'
            oszicar_path = os.path.join(folder, 'OSZICAR')
            outcar_path = os.path.join(folder, 'OUTCAR')

            energies = extract_energies(oszicar_path)
            forces = extract_forces(outcar_path)

            all_energies.extend(energies)
            all_forces.extend(forces)

        # Include the final ts folder
        folder = f'tsfinal/{str(j).zfill(2)}'
        oszicar_path = os.path.join(folder, 'OSZICAR')
        outcar_path = os.path.join(folder, 'OUTCAR')

        final_energies = extract_energies(oszicar_path)
        final_forces = extract_forces(outcar_path)

        all_energies.extend(final_energies)
        all_forces.extend(final_forces)

        all_energies_list.append(all_energies)
        all_forces_list.append(all_forces)

        # Extract data from the tsall folder for comparison
        folder = f'tsall/{str(j).zfill(2)}'
        oszicar_path = os.path.join(folder, 'OSZICAR')
        outcar_path = os.path.join(folder, 'OUTCAR')

        all_comparison_energies = extract_energies(oszicar_path)
        all_comparison_forces = extract_forces(outcar_path)

        all_comparison_energies_list.append(all_comparison_energies)
        all_comparison_forces_list.append(all_comparison_forces)

    # Plot all the energy data together in a 2-column layout
    plt.figure(figsize=(12, 8), dpi=300)  # Increase dpi for higher resolution
    for j in range(1, n_images + 1):
        plt.subplot((n_images + 1) // 2, 2, j)
        plt.plot(range(len(all_energies_list[j - 1])), all_energies_list[j - 1], marker='o', markersize=3, linestyle='-', color='b', label=f'ts DP+DFT iterations')
        plt.plot(range(len(all_comparison_energies_list[j - 1])), all_comparison_energies_list[j - 1], marker='x', markersize=3, linestyle='--', color='g', label=f'ts DFT')
        plt.xlabel('Step')
        plt.ylabel('Energy (eV)')
        plt.title(f'Energy Change for {str(j).zfill(2)}')
        plt.legend()
        plt.grid(True)

    plt.tight_layout()
    plt.savefig('ts_energy_change_comparison.png')
    plt.show()

    # Plot all the force data together in a 2-column layout
    plt.figure(figsize=(12, 8), dpi=300)  # Increase dpi for higher resolution
    for j in range(1, n_images + 1):
        plt.subplot((n_images + 1) // 2, 2, j)
        plt.plot(range(len(all_forces_list[j - 1])), all_forces_list[j - 1], marker='o', markersize=3, linestyle='-', color='r', label=f'ts DP+DFT iterations')
        plt.plot(range(len(all_comparison_forces_list[j - 1])), all_comparison_forces_list[j - 1], marker='x', markersize=3, linestyle='--', color='m', label=f'ts DFT')
        plt.xlabel('Step')
        plt.ylabel('Max Force (eV/Å)')
        plt.title(f'Max Force Change for {str(j).zfill(2)}')
        plt.legend()
        plt.grid(True)

    plt.tight_layout()
    plt.savefig('ts_force_change_comparison.png')
    plt.show()

if __name__ == '__main__':
    main()


