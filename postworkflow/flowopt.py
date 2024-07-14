import os
import shutil
import json
import time
import subprocess
import argparse
from ase import Atoms
from ase import io
from deepmd.calculator import DP
from ase.optimize import BFGS
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

def optimize_structure(iteration, model_path, record_file, final=False):
    step = f'aseopt{iteration}' if not final else 'aseopt-final'
    if read_last_step(record_file) == step:
        return
    
    # Determine the input file for this iteration
    if final:
        input_file = f'opt{iteration}/CONTCAR'
    else:
        input_file = 'POSCAR' if iteration == 1 else f'opt{iteration-1}/CONTCAR'
    
    # Read atomic structure from the input file
    water = io.read(input_file)

    # Set the calculator for the Atoms object using DeepMD model
    water.set_calculator(DP(model=model_path))

    # Optimize the structure
    dyn = BFGS(water)
    dyn.run(fmax=3e-2)
    
    # Save the optimized structure
    output_file = f'CONTCAR-ase-{iteration}' if not final else 'CONTCAR-ase-final'
    io.write(output_file, water)
    
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
    dsys = dpdata.LabeledSystem(f'opt{iteration}/OUTCAR')
    dsys.to("deepmd/npy", f'opt{iteration}/deepmd_data', set_size=dsys.get_nframes())
    
    write_record(record_file, step)

def update_finetune_json(iteration, previous_finetune_file, steps_per_iteration):
    finetune_dir = f'finetune{iteration}'
    os.makedirs(finetune_dir, exist_ok=True)

    # Copy the finetune1.json and sub.dp files to the new finetune directory
    shutil.copy('./utils/finetune1.json', f'{finetune_dir}/finetune1.json')
    shutil.copy('./utils/sub.dp', f'{finetune_dir}/sub.dp')

    finetune_file = f'{finetune_dir}/finetune1.json'

    # Read the existing finetune1.json file and the previous one to append data paths
    with open(finetune_file, 'r') as f:
        data = json.load(f)
    
    if previous_finetune_file:
        with open(previous_finetune_file, 'r') as f_prev:
            prev_data = json.load(f_prev)
        data["training"]["training_data"]["systems"] = prev_data["training"]["training_data"]["systems"]

    # Update the training_data and validation paths
    new_data_path = f"../opt{iteration}/deepmd_data/"
    data["training"]["training_data"]["systems"].append(new_data_path)
    data["training"]["validation_data"]["systems"].append(new_data_path)

    # Update numb_steps based on the current iteration
    data["training"]["numb_steps"] = steps_per_iteration * iteration

    # Write the updated content back to the finetune1.json file
    with open(finetune_file, 'w') as f:
        json.dump(data, f, indent=4)

def finetune_model(iteration, previous_finetune_file, steps_per_iteration):
    update_finetune_json(iteration, previous_finetune_file, steps_per_iteration)
    finetune_dir = f'finetune{iteration}'

    # Submit the finetuning job and retrieve the job ID
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
    parser = argparse.ArgumentParser(description="Run iterative optimization and finetuning.")
    parser.add_argument('--num_iterations', type=int, default=5, help='Number of iterations to run.')
    parser.add_argument('--steps_per_iteration', type=int, default=100, help='Number of steps per iteration.')
    args = parser.parse_args()

    record_file = 'record.txt'
    num_iterations = args.num_iterations
    steps_per_iteration = args.steps_per_iteration
    model_path = "./frozen_model.pth"

    last_step = read_last_step(record_file)
    print(f"Last step in record file: {last_step}")
    start_iteration = 1
    if last_step:
        # Determine which iteration to start from based on the last completed step
        try:
            start_iteration = int(''.join(filter(str.isdigit, last_step))) + 1
        except (ValueError, IndexError):
            print(f"Error parsing the last step: {last_step}")
            return
    
    for i in range(start_iteration, num_iterations + 1):
        optimize_structure(i, model_path, record_file)
        job_id = prepare_vasp_calculation(i)
        wait_for_job_completion(job_id, record_file, f'vaspopt{i}')
        
        generate_new_dataset(i, record_file)
        previous_finetune_file = f'finetune{i-1}/finetune1.json' if i > 1 else None
        job_id = finetune_model(i, previous_finetune_file, steps_per_iteration)
        wait_for_job_completion(job_id, record_file, f'finetune{i}')
        freeze(i, record_file)
        # Update model path to use the newly generated model
        model_path = f'finetune{i}/frozen_model.pth'

    # Run the final optimization step
    optimize_structure(num_iterations, model_path, record_file, final=True)
    job_id = prepare_vasp_calculation(num_iterations, final=True)
    wait_for_job_completion(job_id, record_file, 'vaspopt-final')

if __name__ == "__main__":
    main()
