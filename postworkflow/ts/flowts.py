import os
import shutil
import json
import time
import subprocess
import argparse
from ase import Atoms, io
from deepmd.calculator import DP
from ase.optimize import FIRE
from ase.neb import NEB
from ase.io import read, write
from ase.constraints import FixAtoms
from ase.neb import IDPP
import dpdata

def write_record(record_file, step):
    with open(record_file, 'a') as f:
        f.write(step + '\n')

def read_last_step(record_file):
    if not os.path.exists(record_file):
        return None
    with open(record_file, 'r') as f:
        lines = f.readlines()
    return lines[-1].strip() if lines else None

def neb_calculation(initial_structure, final_structure, model_path, n_images=5, fmax=0.02, interpolation='linear', spring_constant=1.0, apply_constraint=False):
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
        neb.interpolate(apply_constraint=apply_constraint)

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

def prepare_vasp_calculation(n_images, ts_dir, initial_outcar, final_outcar, skip_first_neb=False):
    os.makedirs(ts_dir, exist_ok=True)
    for i in range(n_images + 2):
        subdir = os.path.join(ts_dir, str(i).zfill(2))
        os.makedirs(subdir, exist_ok=True)
        if skip_first_neb:
            shutil.copy(f'POSCAR{str(i).zfill(2)}', os.path.join(subdir, 'POSCAR'))
        else:
            shutil.copy(f'CONTCAR{str(i).zfill(2)}', os.path.join(subdir, 'POSCAR'))

    shutil.copy('./utils/INCAR', ts_dir)
    shutil.copy('./utils/KPOINTS', ts_dir)
    shutil.copy('./utils/POTCAR', ts_dir)
    shutil.copy('./utils/sub.vasp', ts_dir)

    # Change to the new directory and generate the POTCAR file using vaspkit
    os.chdir(ts_dir)
    # os.system('echo -e "103\n" | vaspkit')

    # Submit the VASP job and retrieve the job ID
    result = subprocess.run(['sbatch', 'sub.vasp'], stdout=subprocess.PIPE)
    job_id = result.stdout.decode().strip().split()[-1]

    # Return to the parent directory
    os.chdir('..')

    # Copy initial_outcar and final_outcar to the appropriate subdirectories
    shutil.copy(initial_outcar, os.path.join(ts_dir, '00', 'OUTCAR'))
    shutil.copy(final_outcar, os.path.join(ts_dir, str(n_images + 1).zfill(2), 'OUTCAR'))

    return job_id

def prepare_final_vasp_calculation(n_images, ts_dir, initial_outcar, final_outcar):
    os.makedirs(ts_dir, exist_ok=True)
    for i in range(n_images + 2):
        subdir = os.path.join(ts_dir, str(i).zfill(2))
        os.makedirs(subdir, exist_ok=True)
        shutil.copy(f'CONTCAR{str(i).zfill(2)}', os.path.join(subdir, 'POSCAR'))

    shutil.copy('./utils/INCAR', ts_dir)
    shutil.copy('./utils/KPOINTS', ts_dir)
    shutil.copy('./utils/POTCAR', ts_dir)
    shutil.copy('./utils/sub.vasp', ts_dir)

    # Modify INCAR to set NSW to 300
    incar_path = os.path.join(ts_dir, 'INCAR')
    with open(incar_path, 'r') as file:
        lines = file.readlines()
    with open(incar_path, 'w') as file:
        for line in lines:
            if 'NSW' in line:
                file.write('NSW = 300\n')
            else:
                file.write(line)

    # Change to the new directory and generate the POTCAR file using vaspkit
    os.chdir(ts_dir)
    # os.system('echo -e "103\n" | vaspkit')

    # Submit the VASP job and retrieve the job ID
    result = subprocess.run(['sbatch', 'sub.vasp'], stdout=subprocess.PIPE)
    job_id = result.stdout.decode().strip().split()[-1]

    # Return to the parent directory
    os.chdir('..')

    # Copy initial_outcar and final_outcar to the appropriate subdirectories
    shutil.copy(initial_outcar, os.path.join(ts_dir, '00', 'OUTCAR'))
    shutil.copy(final_outcar, os.path.join(ts_dir, str(n_images + 1).zfill(2), 'OUTCAR'))

    return job_id

def check_job_status(job_id):
    result = subprocess.run(['squeue', '--job', job_id], stdout=subprocess.PIPE)
    return job_id in result.stdout.decode()

def wait_for_job_completion(job_id, record_file, step):
    while check_job_status(job_id):
        print(f"Waiting for job {job_id} to complete...")
        time.sleep(10)
    print(f"Job {job_id} has completed.")
    write_record(record_file, step)

def generate_new_dataset(n_images, ts_dir):
    all_dsys = dpdata.LabeledSystem()  # Create an empty LabeledSystem object to merge all subfolder OUTCAR data

    for i in range(1, n_images + 1):  # Ignore 00 and the last folder
        outcar_path = os.path.join(ts_dir, str(i).zfill(2), 'OUTCAR')
        if os.path.exists(outcar_path):
            dsys = dpdata.LabeledSystem(outcar_path)
            all_dsys.append(dsys)  # Append the current subfolder LabeledSystem data to all_dsys
        else:
            print(f"Warning: OUTCAR file not found at {outcar_path}")
    
    # Save the merged data to the deepmd_data directory
    all_dsys.to("deepmd/npy", os.path.join(ts_dir, 'deepmd_data'), set_size=all_dsys.get_nframes())

def update_finetune_json(step, previous_finetune_file, steps_per_iteration):
    finetune_dir = f'finetune{step}'
    os.makedirs(finetune_dir, exist_ok=True)

    shutil.copy('./utils/finetune1.json', f'{finetune_dir}/finetune1.json')
    shutil.copy('./utils/sub.dp', f'{finetune_dir}/sub.dp')

    finetune_file = f'{finetune_dir}/finetune1.json'

    with open(finetune_file, 'r') as f:
        data = json.load(f)
    
    if previous_finetune_file:
        with open(previous_finetune_file, 'r') as f_prev:
            prev_data = json.load(f_prev)
        data["training"]["training_data"]["systems"] = prev_data["training"]["training_data"]["systems"]

    new_data_path = f"../ts{step}/deepmd_data/"
    data["training"]["training_data"]["systems"].append(new_data_path)
    data["training"]["validation_data"]["systems"].append(new_data_path)

    data["training"]["numb_steps"] = steps_per_iteration * step

    with open(finetune_file, 'w') as f:
        json.dump(data, f, indent=4)

def finetune_model(step, previous_finetune_file, steps_per_iteration):
    update_finetune_json(step, previous_finetune_file, steps_per_iteration)
    finetune_dir = f'finetune{step}'

    os.chdir(finetune_dir)
    result = subprocess.run(['sbatch', 'sub.dp'], stdout=subprocess.PIPE)
    job_id = result.stdout.decode().strip().split()[-1]
    os.chdir('..')

    return job_id

def freeze(iteration, record_file):
    step = f'freeze{iteration}'
    if read_last_step(record_file) == step:
        return

    finetune_dir = f'finetune{iteration}'
    os.chdir(finetune_dir)
    os.system('dp --pt freeze')
    os.chdir('..')
    write_record(record_file, step)

def main():
    parser = argparse.ArgumentParser(description="Run iterative NEB and finetuning.")
    parser.add_argument('initial_structure', type=str, help='Path to the initial structure (POSCAR format).')
    parser.add_argument('final_structure', type=str, help='Path to the final structure (POSCAR format).')
    parser.add_argument('model_path', type=str, help='Path to the DeepMD model.')
    parser.add_argument('initial_outcar', type=str, help='Path to the initial OUTCAR file.')
    parser.add_argument('final_outcar', type=str, help='Path to the final OUTCAR file.')
    parser.add_argument('--num_steps', type=int, default=1, help='Number of steps to run.')
    parser.add_argument('--n_images', type=int, default=4, help='Number of intermediate images.')
    parser.add_argument('--fmax', type=float, default=0.5, help='Maximum force criteria for optimization.')
    parser.add_argument('--interpolation', choices=['linear', 'idpp'], default='linear', help='Interpolation method for generating intermediate images.')
    parser.add_argument('--spring_constant', type=float, default=1.0, help='Spring constant for NEB calculations.')
    parser.add_argument('--steps_per_iteration', type=int, default=1000, help='Number of steps per iteration.')
    parser.add_argument('--apply_constraint', type=bool, default=False, help='Whether to apply constraints during interpolation.')
    parser.add_argument('--skip_first_neb', action='store_true', help='Skip the first ASE NEB optimization.')

    args = parser.parse_args()

    record_file = 'record.txt'
    last_step = read_last_step(record_file)
    print(f"Last step in record file: {last_step}")
    start_step = 1
    if last_step:
        try:
            start_step = int(''.join(filter(str.isdigit, last_step))) + 1
        except (ValueError, IndexError):
            print(f"Error parsing the last step: {last_step}")
            return

    for step in range(start_step, args.num_steps + 1):
        ts_dir = f'ts{step}'
        images = neb_calculation(
            args.initial_structure,
            args.final_structure,
            args.model_path,
            args.n_images,
            args.fmax,
            args.interpolation,
            args.spring_constant,
            args.apply_constraint
        )

        job_id = prepare_vasp_calculation(args.n_images, ts_dir, args.initial_outcar, args.final_outcar, skip_first_neb=args.skip_first_neb and step == 1)
        wait_for_job_completion(job_id, record_file, f'vaspopt{step}')
        
        generate_new_dataset(args.n_images, ts_dir)
        previous_finetune_file = f'finetune{step-1}/finetune1.json' if step > 1 else None
        job_id = finetune_model(step, previous_finetune_file, args.steps_per_iteration)
        wait_for_job_completion(job_id, record_file, f'finetune{step}')

        # Update model path to use the newly generated model
        args.model_path = f'finetune{step}/frozen_model.pth'

        # Run the freeze step
        freeze(step, record_file)

    # Run NEB calculation for the final step
    args.model_path = f'finetune{args.num_steps}/frozen_model.pth'
    final_ts_dir = 'tsfinal'
    images = neb_calculation(
        args.initial_structure,
        args.final_structure,
        args.model_path,
        args.n_images,
        args.fmax,
        args.interpolation,
        args.spring_constant,
        args.apply_constraint
    )

    # Prepare final VASP calculation in tsfinal directory with NSW set to 300
    job_id = prepare_final_vasp_calculation(args.n_images, final_ts_dir, args.initial_outcar, args.final_outcar)
    wait_for_job_completion(job_id, record_file, 'vaspopt-final')

if __name__ == "__main__":
    main()


