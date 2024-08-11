import numpy as np
from ase import Atoms
import matplotlib.pyplot as plt
import torch
from fairchem.core import OCPCalculator
from fairchem.core.datasets.oc22_lmdb_dataset import OC22LmdbDataset
from tqdm import tqdm
import argparse

parser = argparse.ArgumentParser(description='make data test')

parser.add_argument('--checkpoint', type=str,
                    help='the path of the checkpoint file')

parser.add_argument('--dataset', type=str,
                    help='the path of the dataset')

args = parser.parse_args()

# checkpoint = '/home/ljcgroup/ljc/oc-finetune/checkpoints/2024-08-03-23-10-56/checkpoint.pt'
model = OCPCalculator(checkpoint_path=args.checkpoint, cpu= not torch.cuda.is_available())
# prefix = '/home/ljcgroup/wzh/data/lmdb/val/bulk'
data = OC22LmdbDataset({"src":f'{args.dataset}'})

def torch_geometric_to_ase(data):
    positions = data.pos.numpy()  # get atom positions
    sys = data.atomic_numbers.numpy()     # get atomic numbers
    cell = data.cell.numpy().reshape(3,3)
    numbers = data.natoms
    atoms = Atoms( positions=positions, cell=cell,symbols=sys, pbc=True)
    atoms.set_tags(data.tags)
    return atoms

def mae(
    prediction: dict[str, torch.Tensor],
    target: dict[str, torch.Tensor],
) -> dict[str, float | int]:
    error = torch.abs(target - prediction)
    return {
        "metric": torch.mean(error).item(),
        "total": torch.sum(error).item(),
    }

def rmse(
    prediction: dict[str, torch.tensor],
    target: dict[str, torch.tensor],
) -> dict[str, float | int]:
    error = (target - prediction) ** 2
    return {
        "metric": torch.mean(error).item(),
        "total": torch.sum(error).item()
    }

def forcesx_rmse(
    prediction: dict[str, torch.tensor],
    target: dict[str, torch.tensor]
):
    return rmse(prediction[:, 0], target[:, 0])
def forcesx_mae(
    prediction: dict[str, torch.Tensor],
    target: dict[str, torch.Tensor],
):
    return mae(prediction[:, 0], target[:, 0])

def forcesy_rmse(
    prediction: dict[str, torch.tensor],
    target: dict[str, torch.tensor]
):
    return rmse(prediction[:, 1], target[:, 1])
def forcesy_mae(
    prediction: dict[str, torch.Tensor],
    target: dict[str, torch.Tensor],
):
    return mae(prediction[:, 1], target[:, 1])

def forcesz_rmse(
    prediction: dict[str, torch.tensor],
    target: dict[str, torch.tensor]
):
    return rmse(prediction[:, 2], target[:, 2])

def forcesz_mae(
    prediction: dict[str, torch.Tensor],
    target: dict[str, torch.Tensor],
):
    return mae(prediction[:, 2], target[:, 2])

metrics = {
                "rmse_energy_total": 0,
                "rmse_atom_energy_total": 0,
                "rmse_force_x_total": 0,
                "rmse_force_y_total": 0,
                "rmse_force_z_total": 0,
                "mae_energy_total": 0,
                "mae_atom_energy_total": 0,
                "mae_force_x_total": 0,
                "mae_force_y_total": 0,
                "mae_force_z_total": 0,
                "numel": 0,
            }
energy_dft = []
energy_pred = []
forces_x_dft = torch.tensor([])
forces_x_pred = torch.tensor([])
forces_y_dft = torch.tensor([])
forces_y_pred = torch.tensor([])
forces_z_dft = torch.tensor([])
forces_z_pred = torch.tensor([])
for i in tqdm(data):
    ase_data = torch_geometric_to_ase(i)
    ase_data.set_calculator(calc=model)
    metrics['rmse_energy_total'] += rmse(torch.tensor(ase_data.get_potential_energy()), torch.tensor(i.energy))['metric']
    metrics['mae_energy_total'] += mae(torch.tensor(ase_data.get_potential_energy()), torch.tensor(i.energy))['metric']
    metrics['rmse_atom_energy_total'] += rmse(torch.tensor(ase_data.get_potential_energy()), torch.tensor(i.energy))['metric']/torch.tensor(i.natoms) # i.energy.to_tensor())['metric']
    metrics['mae_atom_energy_total'] += mae(torch.tensor(ase_data.get_potential_energy()), torch.tensor(i.energy))['metric']/torch.tensor(i.natoms)
    energy_dft.append(i.energy)
    energy_pred.append(ase_data.get_potential_energy())
    metrics['rmse_force_x_total'] += forcesx_rmse(ase_data.get_forces(), i.forces)['metric']
    metrics['mae_force_x_total'] += forcesx_mae(ase_data.get_forces(), i.forces)['metric']
    forces_x_dft = torch.cat((forces_x_dft, i.forces[:,0]))
    forces_x_pred = torch.cat((forces_x_pred, torch.tensor(ase_data.get_forces()[:,0])))
    metrics['rmse_force_y_total'] += forcesy_rmse(ase_data.get_forces(), i.forces)['metric']
    metrics['mae_force_y_total'] += forcesy_mae(ase_data.get_forces(), i.forces)['metric']
    forces_y_dft = torch.cat((forces_y_dft, i.forces[:,1]))
    forces_y_pred = torch.cat((forces_y_pred, torch.tensor(ase_data.get_forces()[:,1])))
    metrics['rmse_force_z_total'] += forcesz_rmse(ase_data.get_forces(), i.forces)['metric']
    metrics['mae_force_z_total'] += forcesz_mae(ase_data.get_forces(), i.forces)['metric']
    forces_z_dft = torch.cat((forces_z_dft, i.forces[:,2]))
    forces_z_pred = torch.cat((forces_z_pred, torch.tensor(ase_data.get_forces()[:,2])))
    metrics['numel'] += 1

rmse_energy = np.sqrt(metrics['rmse_energy_total']/metrics['numel'])
mae_energy = metrics['mae_energy_total']/metrics['numel']
rmse_atom_energy = np.sqrt(metrics['rmse_atom_energy_total']/metrics['numel'])
mae_atom_energy = metrics['mae_atom_energy_total']/metrics['numel']
rmse_forces_x = np.sqrt(metrics['rmse_force_x_total']/metrics['numel'])
mae_forces_x = metrics['mae_force_x_total']/metrics['numel']
rmse_forces_y = np.sqrt(metrics['rmse_force_y_total']/metrics['numel'])
mae_forces_y = metrics['mae_force_y_total']/metrics['numel']
rmse_forces_z = np.sqrt(metrics['rmse_force_z_total']/metrics['numel'])
mae_forces_z = metrics['mae_force_z_total']/metrics['numel']
print('rmse_energy:', rmse_energy)
print('rmse_atom_energy:', rmse_atom_energy)
print('rmse_forces_x:', rmse_forces_x)
print('rmse_forces_y:', rmse_forces_y)
print('rmse_forces_z:', rmse_forces_z)
print('rmse_force_mean:', np.mean([rmse_forces_x, rmse_forces_y, rmse_forces_z]))

print('mae_energy:', mae_energy)
print('mae_atom_energy:', mae_atom_energy)
print('mae_forces_x:', mae_forces_x)
print('mae_forces_y:', mae_forces_y)
print('mae_forces_z:', mae_forces_z)
print('mae_force_mean:', np.mean([mae_forces_x, mae_forces_y, mae_forces_z]))

#draw figure with rmse
plt.figure(figsize=(12, 8))
plt.subplot(2, 2, 1)
plt.scatter(energy_dft, energy_pred)
plt.text(np.mean(energy_dft), np.mean(energy_pred), f'RMSE = {rmse_energy:.3f} eV')
plt.xlabel('DFT energy (eV)')
plt.ylabel('Predicted energy (eV)')
plt.subplot(2, 2, 2)
plt.scatter(forces_x_dft, forces_x_pred)
plt.text(torch.mean(forces_x_dft), torch.mean(forces_x_pred), f'RMSE = {rmse_forces_x:.3f} eV/Å')
plt.xlabel('DFT force_x (eV/Å)')
plt.ylabel('Predicted force_x (eV/Å)')
plt.subplot(2, 2, 3)
plt.scatter(forces_y_dft, forces_y_pred)
plt.text(torch.mean(forces_y_dft), torch.mean(forces_y_pred), f'RMSE = {rmse_forces_y:.3f} eV/Å')
plt.xlabel('DFT force_y (eV/Å)')
plt.ylabel('Predicted force_y (eV/Å)')
plt.subplot(2, 2, 4)
plt.scatter(forces_z_dft, forces_z_pred)
plt.xlabel('DFT force_z (eV/Å)')
plt.ylabel('Predicted force_z (eV/Å)')
plt.text(torch.mean(forces_z_dft), torch.mean(forces_z_pred), f'RMSE = {rmse_forces_z:.3f} eV/Å')
plt.tight_layout()
plt.savefig('test.png')

