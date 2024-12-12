import os
import matplotlib.pyplot as plt
import numpy as np
from ase.io import read
from ase.visualize.plot import plot_atoms
from matplotlib.ticker import FormatStrFormatter

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

def extract_time(outcar_path):
    time = None
    if os.path.exists(outcar_path):
        with open(outcar_path, 'r') as file:
            for line in file:
                if 'Elapsed time (sec):' in line:
                    time = float(line.split()[-1])
    return time

def extract_finetune_time(out_path):
    time = None
    if os.path.exists(out_path):
        with open(out_path, 'r') as file:
            for line in file:
                if 'Total time taken:' in line:
                    time = float(line.split()[-1])
    return time

def plot_structure(file_path, ax, title, rotation='0x,0y,0z'):
    if os.path.exists(file_path):
        atoms = read(file_path)
        # Translate all atoms by 1/2 of the lattice vector in the x-direction
        cell = atoms.get_cell()
        translation = 0.1 * cell[0]  # 1/2 of the a-vector (x-direction)
        atoms.translate(translation)
        
        # Wrap atoms back into the unit cell
        atoms.wrap()

        plot_atoms(atoms, ax, rotation=rotation, show_unit_cell=2)
    ax.set_title(title)
    ax.axis('off')

def main():
    folders = [
        '11', '12', '13', '14', '21', '71', '23', '24',
    '31', '32', '33', '34', '41', '42', '43', '44',
    '51', '52', '53', '54', '61', '62', '63', '64'
    ]
    times_dict = {
        'opt1': [],
        'opt2': [],
        'opt3': [],
        'optfinal': [],
        'finetune1': [],
        'finetune2': [],
        'finetune3': []
    }
    
    num_folders = len(folders)
    iterations = 3
    energy_rows = (num_folders + 3) // 4
    force_rows = (num_folders + 3) // 4
    structure_rows = (num_folders * 2 + 7) // 8  # 2 plots per folder, 8 columns layout
    
    fig_energy, axs_energy = plt.subplots(energy_rows, 4, figsize=(20, 5 * energy_rows))
    fig_force, axs_force = plt.subplots(force_rows, 4, figsize=(20, 5 * force_rows))
    fig_structure, axs_structure = plt.subplots(structure_rows, 8, figsize=(24, 3 * structure_rows))

    deltas_0 = []
    deltas_5 = []
    deltas_10 = []
    deltas_15 = []
    
    for idx, folder in enumerate(folders):
        print(folder)
        all_energies = []
        all_forces = []
        steps = []
        opt_lengths = []  # Record the length of energies for each opt step

        # Extract data from the iterative optimization folders (opt1, opt2, ...)
        for i in range(1, iterations + 1):
            oszicar_path = os.path.join(folder, f'opt{i}', 'OSZICAR')
            outcar_path = os.path.join(folder, f'opt{i}', 'OUTCAR')
            energies = extract_energies(oszicar_path)
            print(energies)
            forces = extract_forces(outcar_path)
            time = extract_time(outcar_path)
            
            times_dict[f'opt{i}'].append(time)

            all_energies.extend(energies)
            all_forces.extend(forces)
            steps.extend([i] * len(energies))
            opt_lengths.append(len(energies))  # Record the length of this step's energies

        # Include the final optimization data
        final_oszicar_path = os.path.join(folder, 'optfinal', 'OSZICAR')
        final_outcar_path = os.path.join(folder, 'optfinal', 'OUTCAR')

        final_energies = extract_energies(final_oszicar_path)
        final_forces = extract_forces(final_outcar_path)
        final_time = extract_time(final_outcar_path)

        times_dict['optfinal'].append(final_time)

        all_energies.extend(final_energies)
        all_forces.extend(final_forces)
        steps.extend(['final'] * len(final_energies))


        energy_row, energy_col = divmod(idx, 4)
        structure_row, structure_col = divmod(idx * 2, 8)

        # Plotting the energy change
        axs_energy[energy_row, energy_col].plot(range(len(all_energies)), all_energies, marker='o', linestyle='-', color='b', label='opt CLAM iterations')
        axs_energy[energy_row, energy_col].yaxis.set_major_formatter(FormatStrFormatter('%.2f'))  

        if len(all_energies) > 1:
            axs_energy[energy_row, energy_col].scatter(0, all_energies[0], color='black')
            delta_0 = all_energies[0] - all_energies[-1]
            deltas_0.append(delta_0)
            axs_energy[energy_row, energy_col].text(0, all_energies[0], f'ΔE={delta_0:.2f}eV', color='black')
        # Use accumulated length for step positions
        if len(opt_lengths) > 0:
            axs_energy[energy_row, energy_col].scatter(opt_lengths[0], all_energies[opt_lengths[0]], color='red')
            delta_5 = all_energies[opt_lengths[0]] - all_energies[-1]
            deltas_5.append(delta_5)
            axs_energy[energy_row, energy_col].text(opt_lengths[0], all_energies[opt_lengths[0]], f'ΔE={delta_5:.2f}eV', color='red')
        if len(opt_lengths) > 1:
            step_10 = opt_lengths[0] + opt_lengths[1]
            axs_energy[energy_row, energy_col].scatter(step_10, all_energies[step_10], color='orange')
            delta_10 = all_energies[step_10] - all_energies[-1]
            deltas_10.append(delta_10)
            axs_energy[energy_row, energy_col].text(step_10, all_energies[step_10], f'ΔE={delta_10:.2f}eV', color='orange')
        if len(opt_lengths) > 2:
            step_15 = opt_lengths[0] + opt_lengths[1] + opt_lengths[2]
            axs_energy[energy_row, energy_col].scatter(step_15, all_energies[step_15], color='purple')
            delta_15 = all_energies[step_15] - all_energies[-1]
            deltas_15.append(delta_15)
            axs_energy[energy_row, energy_col].text(step_15, all_energies[step_15], f'ΔE={delta_15:.2f}eV', color='purple')


        # Extract finetune times
        finetune_times = []
        for i in range(1, 4):
            finetune_out_path = os.path.join(folder, f'finetune{i}', 'out')
            finetune_time = extract_finetune_time(finetune_out_path)
            times_dict[f'finetune{i}'].append(finetune_time)
            finetune_times.append(finetune_time)

        time_text = (f"opt1: {times_dict['opt1'][-1]:.2f}s, opt2: {times_dict['opt2'][-1]:.2f}s, opt3: {times_dict['opt3'][-1]:.2f}s,\n"
                     f"optfinal: {times_dict['optfinal'][-1]:.2f}s\n"
                     f"ft1: {times_dict['finetune1'][-1]:.2f}s, ft2: {times_dict['finetune2'][-1]:.2f}s, ft3: {times_dict['finetune3'][-1]:.2f}s")
        axs_energy[energy_row, energy_col].text(0.05, 0.95, time_text, transform=axs_energy[energy_row, energy_col].transAxes, fontsize=10, verticalalignment='top')
        
        axs_energy[energy_row, energy_col].set_xlabel('Step')
        axs_energy[energy_row, energy_col].set_ylabel('Energy (eV)')
        axs_energy[energy_row, energy_col].set_title(f'Energy Change - {folder}')
        axs_energy[energy_row, energy_col].legend()
        axs_energy[energy_row, energy_col].grid(True)

        # Plotting the max force change
        axs_force[energy_row, energy_col].plot(range(len(all_forces)), all_forces, marker='o', linestyle='-', color='r', label='opt CLAM iterations')
        #axs_force[energy_row, energy_col].plot(all_comparison_steps, all_comparison_forces, marker='x', linestyle='--', color='m', label='opt DFT')
        axs_force[energy_row, energy_col].set_xlabel('Step')
        axs_force[energy_row, energy_col].set_ylabel('Max Force (eV/Å)')
        axs_force[energy_row, energy_col].set_title(f'Max Force Change - {folder}')
        axs_force[energy_row, energy_col].legend()
        axs_force[energy_row, energy_col].grid(True)

        # Plotting the initial and final structures
        poscar_path = os.path.join(folder, 'POSCAR')
        contcar_path = os.path.join(folder, 'CONTCAR-ase-final')
        
        plot_structure(poscar_path, axs_structure[structure_row, structure_col], title=f'{folder} - Initial Top View', rotation='0x,0y,0z')
        plot_structure(contcar_path, axs_structure[structure_row, structure_col + 1], title=f'{folder} - Final Top View', rotation='0x,0y,0z')

    plt.tight_layout()
    fig_energy.savefig('energy_change_comparison.png')
    fig_force.savefig('force_change_comparison.png')
    fig_structure.savefig('structure_comparison.png')
    plt.show()

    # Define bin edges with a fixed interval of 0.05
    bin_edges = np.arange(min(deltas_0 + deltas_5 + deltas_10 + deltas_15), max(deltas_0 + deltas_5 + deltas_10 + deltas_15) + 0.4, 0.05)
    
    mae_0 = np.mean(np.abs(deltas_0))
    mae_5 = np.mean(np.abs(deltas_5))
    mae_10 = np.mean(np.abs(deltas_10))
    mae_15 = np.mean(np.abs(deltas_15))

    print("mae: ",mae_0, mae_5, mae_10, mae_15)
    # Plot the distribution of energy differences with consistent bin edges
    fig, axs = plt.subplots(4, 1, figsize=(10, 20))
    
    font_size = 25 

    axs[0].hist(deltas_0, bins=bin_edges, alpha=0.5, label='ΔE (0th step)')
    axs[0].set_xlabel('Energy Difference (eV)', fontsize=font_size)
    axs[0].set_ylabel('Frequency', fontsize=font_size)
    axs[0].set_title('Distribution of Energy Differences at 0th Step', fontsize=font_size)
    axs[0].legend(fontsize=font_size)
    axs[0].grid(True)
    axs[0].tick_params(axis='both', which='major', labelsize=font_size)

    axs[1].hist(deltas_5, bins=bin_edges, alpha=0.5, label='ΔE (5th step)')
    axs[1].set_xlabel('Energy Difference (eV)', fontsize=font_size)
    axs[1].set_ylabel('Frequency', fontsize=font_size)
    axs[1].set_title('Distribution of Energy Differences at 5th Step', fontsize=font_size)
    axs[1].legend(fontsize=font_size)
    axs[1].grid(True)
    axs[1].tick_params(axis='both', which='major', labelsize=font_size)

    axs[2].hist(deltas_10, bins=bin_edges, alpha=0.5, label='ΔE (10th step)')
    axs[2].set_xlabel('Energy Difference (eV)', fontsize=font_size)
    axs[2].set_ylabel('Frequency', fontsize=font_size)
    axs[2].set_title('Distribution of Energy Differences at 10th Step', fontsize=font_size)
    axs[2].legend(fontsize=font_size)
    axs[2].grid(True)
    axs[2].tick_params(axis='both', which='major', labelsize=font_size)

    axs[3].hist(deltas_15, bins=bin_edges, alpha=0.5, label='ΔE (15th step)')
    axs[3].set_xlabel('Energy Difference (eV)', fontsize=font_size)
    axs[3].set_ylabel('Frequency', fontsize=font_size)
    axs[3].set_title('Distribution of Energy Differences at 15th Step', fontsize=font_size)
    axs[3].legend(fontsize=font_size)
    axs[3].grid(True)
    axs[3].tick_params(axis='both', which='major', labelsize=font_size)
    
    plt.tight_layout()
    plt.savefig('energy_differences_distribution.png')
    plt.show()
    
    # Print average times
    for key in times_dict:
        avg_time = np.mean([t for t in times_dict[key] if t is not None])
        print(f"Average time for {key}: {avg_time:.2f} seconds")

if __name__ == '__main__':
    main()

