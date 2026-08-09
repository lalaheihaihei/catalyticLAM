import numpy as np
from ase import Atoms
import matplotlib.pyplot as plt
import torch
from fairchem.core import OCPCalculator
from fairchem.core.datasets.oc22_lmdb_dataset import OC22LmdbDataset
from tqdm import tqdm
import argparse

parser = argparse.ArgumentParser(description='make data test')
parser.add_argument('--checkpoint', type=str, help='the path of the checkpoint file')
parser.add_argument('--dataset', type=str, help='the path of the dataset')
args = parser.parse_args()

model = OCPCalculator(checkpoint_path=args.checkpoint, cpu=not torch.cuda.is_available())
data = OC22LmdbDataset({"src": f'{args.dataset}'})

def torch_geometric_to_ase(data):
    positions = data.pos.numpy()
    sys = data.atomic_numbers.numpy()
    cell = data.cell.numpy().reshape(3, 3)
    numbers = data.natoms
    atoms = Atoms(positions=positions, cell=cell, symbols=sys, pbc=True)
    atoms.set_tags(data.tags)
    return atoms

def calculate_error(prediction, target, metric_func):
    error = metric_func(target - prediction)
    return error.mean().item(), error.sum().item()

def calculate_metrics(energy_dft, energy_pred, forces_dft, forces_pred, numel):
    rmse_energy = np.sqrt(calculate_error(torch.tensor(energy_pred), torch.tensor(energy_dft), torch.abs)[0] / numel)
    mae_energy = calculate_error(torch.tensor(energy_pred), torch.tensor(energy_dft), torch.abs)[0] / numel
    rmse_atom_energy = np.sqrt(calculate_error(torch.tensor(energy_pred), torch.tensor(energy_dft), torch.abs)[0] / (numel * torch.tensor(energy_dft).size(0)))
    mae_atom_energy = calculate_error(torch.tensor(energy_pred), torch.tensor(energy_dft), torch.abs)[0] / (numel * torch.tensor(energy_dft).size(0))
    return rmse_energy, mae_energy, rmse_atom_energy, mae_atom_energy

def calculate_forces_metrics(forces_dft, forces_pred, numel):
    rmse_forces = np.sqrt(calculate_error(forces_pred, forces_dft, torch.abs)[0] / numel)
    mae_forces = calculate_error(forces_pred, forces_dft, torch.abs)[0] / numel
    return rmse_forces, mae_forces

def main():
    energy_dft = []
    energy_pred = []
    forces_dft = []
    forces_pred = []
    numel = 0

    for i in tqdm(data):
        ase_data = torch_geometric_to_ase(i)
        ase_data.set_calculator(calc=model)
        energy_dft.append(i.energy)
        energy_pred.append(ase_data.get_potential_energy())
        forces_dft.append(i.forces)
        forces_pred.append(ase_data.get_forces())
        numel += 1

    energy_dft = torch.tensor(energy_dft)
    energy_pred = torch.tensor(energy_pred)
    forces_dft = torch.tensor(forces_dft)
    forces_pred = torch.tensor(forces_pred)

    rmse_energy, mae_energy, rmse_atom_energy, mae_atom_energy = calculate_metrics(energy_dft, energy_pred, forces_dft, forces_pred, numel)
    rmse_forces, mae_forces = calculate_forces_metrics(forces_dft, forces_pred, numel)

    print(f'rmse_energy: {rmse_energy}')
    print(f'rmse_atom_energy: {rmse_atom_energy}')
    print(f'rmse_forces: {rmse_forces}')
    print(f'mae_energy: {mae_energy}')
    print(f'mae_atom_energy: {mae_atom_energy}')
    print(f'mae_forces: {mae_forces}')

    plt.figure(figsize=(12, 8))
    plt.subplot(2, 2, 1)
    plt.scatter(energy_dft, energy_pred)
    plt.text(np.mean(energy_dft), np.mean(energy_pred), f'RMSE = {rmse_energy:.3f} eV')
    plt.xlabel('DFT energy (eV)')
    plt.ylabel('Predicted energy (eV)')

    for i, label in enumerate(['x', 'y', 'z']):
        plt.subplot(2, 2, i + 2)
        plt.scatter(forces_dft[:, i], forces_pred[:, i])
        plt.text(torch.mean(forces_dft[:, i]), torch.mean(forces_pred[:, i]), f'RMSE = {rmse_forces[i]:.3f} eV/Å')
        plt.xlabel(f'DFT force_{label} (eV/Å)')
        plt.ylabel(f'Predicted force_{label} (eV/Å)')

    plt.tight_layout()
    plt.savefig('test.png')

if __name__ == '__main__':
    main()