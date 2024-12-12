import os
import shutil
import json
import time
import subprocess
import argparse
from ase import Atoms, io
from ase.optimize import FIRE, BFGS
from ase.neb import NEB
from ase.io import read, write
from ase.constraints import FixAtoms
from fairchem.core import OCPCalculator
from ase.db import connect
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

def tag_atoms(slab):
    tags = [1 if atom.symbol not in ['H', 'C', 'N', 'O', 'S', 'Cl', 'Br', 'I', 'F'] else 2 for atom in slab]
    slab.set_tags(tags)
    return slab

def neb_calculation(initial_structure, final_structure, checkpoint_path, fixed_atoms, n_images=5, fmax=0.02, interpolation='linear', spring_constant=1.0):
    initial = read(initial_structure)
    final = read(final_structure)

    z_positions = initial.positions[:, 2]
    sorted_indices = z_positions.argsort()
    fixed_indices = sorted_indices[:fixed_atoms]

    initial = tag_atoms(initial)
    final = tag_atoms(final)

    images = [initial] + [initial.copy() for _ in range(n_images)] + [final]

    neb = NEB(images, k=spring_constant, climb=True)
    if interpolation == 'idpp':
        neb.idpp_interpolate()
    else:
        neb.interpolate()

    for i, image in enumerate(images):
        write(f'POSCAR{str(i).zfill(2)}', image)

    for image in images:
        image.set_calculator(OCPCalculator(checkpoint_path=checkpoint_path))
        image.set_constraint(FixAtoms(indices=fixed_indices))

    optimizer = FIRE(neb)
    optimizer.run(fmax=fmax, steps=50)

    for i, image in enumerate(images):
        write(f'CONTCAR{str(i).zfill(2)}', image)
    
    return images

def convert_outcar_to_ase_db(outcar_path, db_path):
    db_dir = os.path.dirname(db_path)
    if not os.path.exists(db_dir):
        os.makedirs(db_dir)

    try:
        data = dpdata.LabeledSystem(outcar_path)
    except Exception as e:
        print(f"Error loading OUTCAR file: {e}")
        return

    db = connect(db_path)
    for i, structure in enumerate(data):
        ase_structures = structure.to_ase_structure()
        if isinstance(ase_structures, list):
            for j, ase_atoms in enumerate(ase_structures):
                db.write(ase_atoms, key_value_pairs={"step": i, "frame": j})
        else:
            db.write(ase_structures, key_value_pairs={"step": i})

    print(f"Data successfully written to {db_path}")

def prepare_vasp_calculation(step, n_images, ts_dir, initial_outcar, final_outcar, skip_first_neb=False):
    os.makedirs(ts_dir, exist_ok=True)
    for i in range(n_images + 2):
        subdir = os.path.join(ts_dir, str(i).zfill(2))
        os.makedirs(subdir, exist_ok=True)
        if skip_first_neb and step == 1:
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

def check_for_checkpoints_and_handle(job_id, iteration, record_file, step):
    checkpoints_dir = f'finetune{iteration}/checkpoints'
    while check_job_status(job_id):
        print(f"Waiting for job {job_id} to complete...")
        time.sleep(10)
    print(f"Job {job_id} has completed.")

    found_checkpoint = False
    for root, dirs, files in os.walk(checkpoints_dir):
        if 'best_checkpoint.pt' in files:
            found_checkpoint = True
            checkpoint_path = os.path.join(root, 'best_checkpoint.pt')
            break

    if found_checkpoint:
        shutil.copy(checkpoint_path, f'finetune{iteration}/checkpoint.pt')
        print(f"Checkpoint found and copied from {checkpoint_path}.")
    else:
        print("No checkpoint.pt file found in the subdirectories.")

    write_record(record_file, step)

def wait_for_oc_job_checkpoint(job_id, record_file, step, iteration):
    check_for_checkpoints_and_handle(job_id, iteration, record_file, step)

def wait_for_job_completion(job_id, record_file, step):
    while check_job_status(job_id):
        print(f"Waiting for job {job_id} to complete...")
        time.sleep(10)
    print(f"Job {job_id} has completed.")
    write_record(record_file, step)

def generate_new_dataset(iteration, n_images, record_file):
    step = f'dpdata{iteration}'
    if read_last_step(record_file) == step:
        return

    db_path = f'./output_database.db'
    for i in range(1, n_images + 1):
        outcar_path = f'ts{iteration}/{str(i).zfill(2)}/OUTCAR'
        if os.path.exists(outcar_path):
            convert_outcar_to_ase_db(outcar_path, db_path)
        else:
            print(f"Warning: OUTCAR file not found at {outcar_path}")

    write_record(record_file, step)

def update_finetune_json(iteration, epoch_per_iteration):
    finetune_dir = f'finetune{iteration}'
    os.makedirs(finetune_dir, exist_ok=True)

    for file in ['finetune1.yml', 'base.yml', 'main.py', 'sub.oc']:
        shutil.copy(f'./utils/{file}', f'{finetune_dir}/{file}')

def finetune_model(iteration, epoch_per_iteration):
    update_finetune_json(iteration, epoch_per_iteration)
    finetune_dir = f'finetune{iteration}'

    os.chdir(finetune_dir)
    result = subprocess.run(['sbatch', 'sub.oc'], stdout=subprocess.PIPE)
    job_id = result.stdout.decode().strip().split()[-1]
    os.chdir('..')

    return job_id

def freeze(iteration, record_file):
    step = f'freeze{iteration}'
    write_record(record_file, step)

def clean():
    # Remove directories and files as specified
    os.system('rm -rf ts1 ts2 ts3 ts4 ts5 ts6 ts7 ts8 tsfinal finetune* CONTCAR* POSCAR0* *.png output_database.db nohup.out is fs record.txt')

def main():
    parser = argparse.ArgumentParser(description="Run iterative optimization and finetuning.")
    parser.add_argument('--initial_structure', type=str, default='POSCARis', help='Path to the initial structure (POSCAR format).')
    parser.add_argument('--final_structure', type=str, default='POSCARfs', help='Path to the final structure (POSCAR format).')
    parser.add_argument('--model_path', type=str, default='./utils/best_checkpoint.pt', help='Path to the OCP model checkpoint.')
    parser.add_argument('--initial_outcar', type=str, default='OUTCARis', help='Path to the initial OUTCAR file.')
    parser.add_argument('--final_outcar', type=str, default='OUTCARfs', help='Path to the final OUTCAR file.')
    parser.add_argument('--fixed_atoms', type=int, default=48, help='Fixed Atoms.')
    parser.add_argument('--num_steps', type=int, default=1, help='Number of steps to run.')
    parser.add_argument('--n_images', type=int, default=4, help='Number of intermediate images.')
    parser.add_argument('--fmax', type=float, default=0.5, help='Maximum force criteria for optimization.')
    parser.add_argument('--interpolation', choices=['linear', 'idpp'], default='linear', help='Interpolation method for generating intermediate images.')
    parser.add_argument('--spring_constant', type=float, default=1.0, help='Spring constant for NEB calculations.')
    parser.add_argument('--epoch_per_iteration', type=int, default=30, help='Number of steps per iteration.')
    parser.add_argument('--apply_constraint', type=bool, default=False, help='Whether to apply constraints during interpolation.')
    parser.add_argument('--skip_first_neb', action='store_true', help='Whether to skip the first ase cineb step.')
    parser.add_argument('--clean', action='store_true', help='Clean up specified directories and files.')
    args = parser.parse_args()
    # If the clean flag is set, call the clean function and exit
    if args.clean:
        clean()
        return  # Exit after cleaning
        
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
    
    if os.path.basename(args.initial_outcar) != 'OUTCAR':
        os.makedirs('is', exist_ok=True)
        shutil.copy(args.initial_outcar, './is/OUTCAR')
        args.initial_outcar = './is/OUTCAR'

    if os.path.basename(args.final_outcar) != 'OUTCAR':
        os.makedirs('fs', exist_ok=True)
        shutil.copy(args.final_outcar, './fs/OUTCAR')
        args.final_outcar = './fs/OUTCAR'
    
    if start_step == 1:
        convert_outcar_to_ase_db(args.initial_outcar, './output_database.db')
        convert_outcar_to_ase_db(args.final_outcar, './output_database.db')
        if not args.skip_first_neb: # if skip the first ase-neb calc, there is no need to do finetune0
            update_finetune_json(0, args.epoch_per_iteration)
            job_id = finetune_model(0, args.epoch_per_iteration)
            wait_for_oc_job_checkpoint(job_id, record_file, 'finetune0', 0)
            args.model_path = f'finetune0/checkpoint.pt'
    
    if args.num_steps == 0:

        if args.skip_first_neb: # if skip the first ase-neb calc, only to 1 step images to get POSCAR0*
            images = neb_calculation(
                args.initial_structure,
                args.final_structure,
                args.model_path,
                args.fixed_atoms,
                args.n_images,
                10.0,
                args.interpolation,
                args.spring_constant)
        else:
            images = neb_calculation(
                args.initial_structure,
                args.final_structure,
                args.model_path,
                args.fixed_atoms,
                args.n_images,
                args.fmax,
                args.interpolation,
                args.spring_constant)

        return
    
    for step in range(start_step, args.num_steps + 1):
        ts_dir = f'ts{step}'

        images = neb_calculation(
            args.initial_structure,
            args.final_structure,
            args.model_path,
            args.fixed_atoms,
            args.n_images,
            args.fmax,
            args.interpolation,
            args.spring_constant
        )
    
        job_id = prepare_vasp_calculation(step, args.n_images, ts_dir, args.initial_outcar, args.final_outcar, args.skip_first_neb)
        wait_for_job_completion(job_id, record_file, f'vaspopt{step}')
        
        generate_new_dataset(step, args.n_images, record_file)
        job_id = finetune_model(step, args.epoch_per_iteration)
        wait_for_oc_job_checkpoint(job_id, record_file, f'finetune{step}', step)
        freeze(step, record_file)
        args.model_path = f'finetune{step}/checkpoint.pt'
    
    args.model_path = f'finetune{args.num_steps}/checkpoint.pt'
    final_ts_dir = 'tsfinal'
    images = neb_calculation(
        args.initial_structure,
        args.final_structure,
        args.model_path,
        args.fixed_atoms,
        args.n_images,
        args.fmax,
        args.interpolation,
        args.spring_constant
    )
    
    job_id = prepare_final_vasp_calculation(args.n_images, final_ts_dir, args.initial_outcar, args.final_outcar)
    wait_for_job_completion(job_id, record_file, 'vaspopt-final')

if __name__ == "__main__":
    main()
