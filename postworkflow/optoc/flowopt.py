import os
import shutil
import json
import time
import subprocess
import argparse
from ase import Atoms
from ase import io
from ase.constraints import FixAtoms  
from fairchem.core import OCPCalculator
from ase.optimize import BFGS
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

def convert_outcar_to_ase_db(outcar_path, db_path):
    # Ensure the database directory exists
    db_dir = os.path.dirname(db_path)
    if not os.path.exists(db_dir):
        os.makedirs(db_dir)

    # Load the data from OUTCAR
    try:
        data = dpdata.LabeledSystem(outcar_path)
    except Exception as e:
        print(f"Error loading OUTCAR file: {e}")
        return

    # Open the ASE database for writing
    db = connect(db_path)
    
    # Iterate over the data and add each frame to the database
    for i, structure in enumerate(data):
        # Convert dpdata structure to ASE Atoms object(s)
        ase_structures = structure.to_ase_structure()
        
        # Handle case where `to_ase_structure()` returns a list of Atoms objects
        if isinstance(ase_structures, list):
            for j, ase_atoms in enumerate(ase_structures):
                # Add each Atoms object to the database
                db.write(ase_atoms, key_value_pairs={"step": i, "frame": j})
        else:
            # If it's a single Atoms object, write it to the database
            db.write(ase_structures, key_value_pairs={"step": i})

    print(f"Data successfully written to {db_path}")

def optimize_structure(iteration, model_path, record_file, fixed_atoms, fmax, final=False):
    step = f'aseopt{iteration}' if not final else 'aseopt-final'
    if read_last_step(record_file) == step:
        return
    
    # Determine the input file for this iteration
    if final:
        input_file = f'opt{iteration}/CONTCAR'
    else:
        input_file = 'POSCAR' if iteration == 1 else f'opt{iteration-1}/CONTCAR'
    
    # Read atomic structure from the input file
    slab = io.read(input_file)
    tags = []
    for atom in slab:
        if atom.symbol == 'Ru':
            tags.append(1)
        else:
            tags.append(0)
    slab.set_tags(tags)

    # Fix the positions of the lowest fixed_atoms atoms
    sorted_indices = sorted(range(len(slab)), key=lambda i: slab.positions[i][2])
    fixed_indices = sorted_indices[:fixed_atoms]
    constraints = FixAtoms(indices=fixed_indices)
    slab.set_constraint(constraints)

    # Set the calculator for the Atoms object using OCPCalculator
    calculator = OCPCalculator(checkpoint_path=model_path)
    slab.set_calculator(calculator)

    # Optimize the structure
    dyn = BFGS(slab)
    dyn.run(fmax=fmax,steps=30)
    
    # Save the optimized structure
    output_file = f'CONTCAR-ase-{iteration}' if not final else 'CONTCAR-ase-final'
    io.write(output_file, slab)
    
    write_record(record_file, step)

def prepare_vasp_calculation(iteration, final=False):
    if final:
        opt_dir = 'optfinal'
        os.makedirs(opt_dir, exist_ok=True)

        # Copy the final optimized structure to the new directory as POSCAR
        shutil.copy('CONTCAR-ase-final', f'{opt_dir}/POSCAR')
    else:
        opt_dir = f'opt{iteration}'
        os.makedirs(opt_dir, exist_ok=True)

        # Copy the optimized structure to the new directory as POSCAR
        shutil.copy(f'CONTCAR-ase-{iteration}', f'{opt_dir}/POSCAR')

    # Copy necessary VASP input files from the utils directory to the new directory
    shutil.copy('./utils/INCAR', f'{opt_dir}/INCAR')
    shutil.copy('./utils/KPOINTS', f'{opt_dir}/KPOINTS')
    shutil.copy('./utils/sub.vasp', f'{opt_dir}/sub.vasp')

    # If this is not the first iteration, copy WAVECAR and CHGCAR from the previous iteration
    if iteration > 1 and not final:
        prev_opt_dir = f'opt{iteration-1}'
        shutil.copy(f'{prev_opt_dir}/WAVECAR', f'{opt_dir}/WAVECAR')
        shutil.copy(f'{prev_opt_dir}/CHGCAR', f'{opt_dir}/CHGCAR')
    
    # Modify INCAR file for the final calculation
    if final:
        incar_path = f'{opt_dir}/INCAR'
        prev_opt_dir = f'opt{iteration-1}'
        shutil.copy(f'{prev_opt_dir}/WAVECAR', f'{opt_dir}/WAVECAR')
        shutil.copy(f'{prev_opt_dir}/CHGCAR', f'{opt_dir}/CHGCAR')

        with open(incar_path, 'r') as file:
            lines = file.readlines()
        with open(incar_path, 'w') as file:
            for line in lines:
                if 'NSW' not in line:
                    file.write(line)
            file.write('NSW = 300\n')
    
    # Change to the new directory and generate the POTCAR file using vaspkit
    os.chdir(opt_dir)
    os.system('echo -e "103\n" | vaspkit')

    # Submit the VASP job and retrieve the job ID
    result = subprocess.run(['sbatch', 'sub.vasp'], stdout=subprocess.PIPE)
    job_id = result.stdout.decode().strip().split()[-1]

    # Return to the parent directory
    os.chdir('..')

    return job_id

def check_job_status(job_id):
    # Check the status of the job using the job ID
    result = subprocess.run(['squeue', '--job', job_id], stdout=subprocess.PIPE)
    return job_id in result.stdout.decode()

def check_for_checkpoints_and_handle(job_id, iteration, record_file, step):
    checkpoints_dir = f'finetune{iteration}/checkpoints'

    # Wait for the job to complete
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
        # Copy the found checkpoint to the main finetune directory
        shutil.copy(checkpoint_path, f'finetune{iteration}/checkpoint.pt')
        print(f"Checkpoint found and copied from {checkpoint_path}.")
    else:
        print("No checkpoint.pt file found in the subdirectories.")

    # Record the step
    write_record(record_file, step)

# Update wait_for_job_completion to use the new function for checking and handling checkpoints
def wait_for_oc_job_checkpoint(job_id, record_file, step, iteration):
    check_for_checkpoints_and_handle(job_id, iteration, record_file, step)

def wait_for_job_completion(job_id, record_file, step):
    # Wait for the job to complete by polling its status every 10 seconds
    while check_job_status(job_id):
        print(f"Waiting for job {job_id} to complete...")
        time.sleep(10)
    print(f"Job {job_id} has completed.")
    write_record(record_file, step)

def generate_new_dataset(iteration, record_file):
    step = f'dpdata{iteration}'
    if read_last_step(record_file) == step:
        return
    
    # Generate a new dataset from the VASP output file OUTCAR
    outcar_path = f'opt{iteration}/OUTCAR'
    db_path = f'./output_database.db'
    convert_outcar_to_ase_db(outcar_path, db_path)
    
    write_record(record_file, step)

def update_finetune_json(iteration, previous_finetune_file, steps_per_iteration):
    finetune_dir = f'finetune{iteration}'
    os.makedirs(finetune_dir, exist_ok=True)

    # Copy the finetune1.json and sub.dp files to the new finetune directory
    shutil.copy('./utils/finetune1.yml', f'{finetune_dir}/finetune1.yml')
    shutil.copy('./utils/base.yml', f'{finetune_dir}/base.yml')
    shutil.copy('./utils/main.py', f'{finetune_dir}/main.py')
    shutil.copy('./utils/sub.oc', f'{finetune_dir}/sub.oc')

    finetune_file = f'{finetune_dir}/finetune1.yml'

def finetune_model(iteration, previous_finetune_file, steps_per_iteration):
    update_finetune_json(iteration, previous_finetune_file, steps_per_iteration)
    finetune_dir = f'finetune{iteration}'

    # Submit the finetuning job and retrieve the job ID
    os.chdir(finetune_dir)
    result = subprocess.run(['sbatch', 'sub.oc'], stdout=subprocess.PIPE)
    job_id = result.stdout.decode().strip().split()[-1]
    os.chdir('..')

    return job_id

def freeze(iteration, record_file):
    step = f'freeze{iteration}'
    #if read_last_step(record_file) == step:
    #    return
    
    #finetune_dir = f'finetune{iteration}'
    #os.chdir(finetune_dir)
    #os.system('dp --pt freeze')
    #os.chdir('..')
    
    write_record(record_file, step)

def main():
    parser = argparse.ArgumentParser(description="Run iterative optimization and finetuning.")
    parser.add_argument('--num_iterations', type=int, default=4, help='Number of iterations to run.')
    parser.add_argument('--steps_per_iteration', type=int, default=100, help='Number of steps per iteration.')
    parser.add_argument('--fixed_atoms', type=int, default=0, help='Number of atoms to fix in position.')
    parser.add_argument('--iffinal', type=bool, default=True, help='Whether to perform the final optimization step.')
    parser.add_argument('--skip_first_opt', type=bool, default=False, help='Whether to skip the first optimization step.')
    parser.add_argument('--fmax', type=float, default=0.2, help='Maximum force criteria for optimization.')
    args = parser.parse_args()

    record_file = 'record.txt'
    num_iterations = args.num_iterations
    steps_per_iteration = args.steps_per_iteration
    fixed_atoms = args.fixed_atoms
    fmax = args.fmax

    last_step = read_last_step(record_file)
    print(f"Last step in record file: {last_step}")
    start_iteration = 1
    if last_step:
        # Determine which iteration to start from based on the last completed step
        try:
            start_iteration = int(''.join(filter(str.isdigit, last_step))) + 1
            model_path = f'finetune{int(start_iteration-1)}/checkpoint.pt'
        except (ValueError, IndexError):
            print(f"Error parsing the last step: {last_step}")
            return
    else:
        model_path = "./utils/gnoc_oc22_oc20_all_s2ef.pt"

    for i in range(start_iteration, num_iterations + 1):
        if not (i == 1 and args.skip_first_opt):
            # model_path = f'finetune{i-1}/frozen_model.pth'
            optimize_structure(i, model_path, record_file, fixed_atoms, fmax)
        else:
            # Directly copy the initial POSCAR to opt1 for VASP calculation
            os.makedirs(f'opt{i}', exist_ok=True)
            shutil.copy('POSCAR', f'opt{i}/POSCAR')
            shutil.copy('POSCAR', f'./CONTCAR-ase-1')
        job_id = prepare_vasp_calculation(i)
        wait_for_job_completion(job_id, record_file, f'vaspopt{i}')
        
        generate_new_dataset(i, record_file)
        previous_finetune_file = f'finetune{i-1}/finetune1.yml' if i > 1 else None
        job_id = finetune_model(i, previous_finetune_file, steps_per_iteration)
        wait_for_oc_job_checkpoint(job_id, record_file, f'finetune{i}', i)
        #wait_for_job_completion(job_id, record_file, f'finetune{i}')
        freeze(i, record_file)
        # Update model path to use the newly generated model
        model_path = f'finetune{i}/checkpoint.pt'

    # Run the final optimization step if specified
    if args.iffinal:
        optimize_structure(num_iterations, model_path, record_file, fixed_atoms, fmax, final=True)
        job_id = prepare_vasp_calculation(num_iterations, final=True)
        wait_for_job_completion(job_id, record_file, 'vaspopt-final')

if __name__ == "__main__":
    main()

