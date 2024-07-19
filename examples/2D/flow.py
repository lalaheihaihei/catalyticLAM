import argparse
import os
import shutil
import time
from glob import glob
import subprocess
#from dpdata import LabeledSystem, MultiSystems

class VASPManager:
    def __init__(self, work_path, vasp_sub_path, structre_db_path, max_jobs, sleep_time, prefix):
        self.work_path = work_path
        self.vasp_sub_path = vasp_sub_path
        self.structre_db_path = structre_db_path
        self.max_jobs = max_jobs
        self.sleep_time = sleep_time
        self.files_to_copy = glob(os.path.join(self.structre_db_path, f"*{prefix}*"))

    def copy_and_enter_folder(self, prefix):
        files_to_copy = glob(os.path.join(self.structre_db_path, f"*{prefix}*"))
        if not files_to_copy:
            print(f"No files found with prefix '{prefix}' in '{self.structre_db_path}'.")
            return False
        target_folder = os.path.join(os.getcwd(), prefix)
        os.makedirs(target_folder, exist_ok=True)
        for file_path in files_to_copy:
            shutil.copy(file_path, target_folder)
        os.chdir(target_folder)
        print(f"Entered folder '{target_folder}' and perform operations")
        return True

class VASPJobManager:
    def __init__(self, manager, atomic_structure_manager, opt_incar_path, md_incar_path):
        self.manager = manager
        self.atomic_structure_manager = atomic_structure_manager
        self.opt_incar_path = opt_incar_path  # Store the path for use in operations
        self.md_incar_path = md_incar_path

    def change_magmom_incar(self, directory):
        poscar_path = os.path.join(directory, "POSCAR")
        incar_path = os.path.join(directory, "INCAR")
        
        try:
            with open(poscar_path, "r") as file:
                lines = file.readlines()
        except FileNotFoundError:
            print(f"{poscar_path} not found")
            return
    
        if len(lines) < 7:
            print("POSCAR file format is incorrect")
            return
    
        element_names = lines[5].strip().split()
        atom_counts = lines[6].strip().split()
        mag_list = []
    
        for element in element_names:
            if element in ["Fe", "Co", "Ni", "Mn"]:
                mag_list.append("3.0")  # Example value, adjust based on your simulation needs
            else:
                mag_list.append("0.0")  # Default magnetic moment for non-magnetic elements
    
        if '_mag' in poscar_path:  # Check if the POSCAR filename indicates magnetic calculations
            magmom_line = "MAGMOM = " + " ".join(f"{count}*{mag}" for count, mag in zip(atom_counts, mag_list))
            with open(incar_path, "a") as incar_file:
                incar_file.write("\nAMIX = 0.2\n")
                incar_file.write("BMIX = 0.0001\n")
                incar_file.write("AMIX_MAG = 0.8\n")
                incar_file.write("BMIX_MAG = 0.0001\n")
                incar_file.write("LMAXMIX = 4\n")  # Assuming LMAXMIX value is needed, replace 4 with appropriate value
                incar_file.write("ISPIN = 2\n")
                incar_file.write(magmom_line + "\n")
            print(f"MAGMOM and other magnetic parameters updated in {incar_path}")
    
    def opt_operation(self):
        current_directory = os.getcwd()
        record_file_path = os.path.join(self.manager.work_path, "record.txt")  # Path for the log file
    
        files = [f for f in os.listdir(current_directory) if f.startswith("POSCAR") and os.path.isfile(os.path.join(current_directory, f))] # not include dirctory f, because of restart calculation 
        for file in files:
            poscar_path = os.path.join(current_directory, file)
            new_folder = os.path.join(current_directory, file + "d")
            opt_folder = os.path.join(new_folder, "OPT")
    
            if not self.atomic_structure_manager.check_atom_number(poscar_path):
                with open(record_file_path, "a") as record_file:
                    record_file.write(f"Skipped {file} due to excessive atoms (>300) in {opt_folder}\n")
                continue
    
            if os.path.exists(opt_folder):
                if not ConvergenceChecker().check_scf_convergence(opt_folder):
                    with open(record_file_path, "a") as record_file:
                        record_file.write(f"Re-running OPT due to non-convergence for {file} in {opt_folder}\n")
                    self.prepare_opt_folder(new_folder, opt_folder, poscar_path)
                else:
                    print(f"OPT already converged in {opt_folder}. Skipping re-calculation.")
            else:
                os.makedirs(opt_folder, exist_ok=True)
                self.prepare_opt_folder(new_folder, opt_folder, poscar_path)
    
    def prepare_opt_folder(self, new_folder, opt_folder, poscar_path):
        shutil.move(poscar_path, os.path.join(new_folder, "POSCAR"))
        shutil.copy(os.path.join(new_folder, "POSCAR"), opt_folder)
        os.system(f'echo -e "102\n2\n0" | vaspkit')
        shutil.copy(self.opt_incar_path, os.path.join(opt_folder, "INCAR"))
        self.change_magmom_incar(opt_folder)  # Update MAGMOM before submitting the job
        self.wait_and_submit(opt_folder, os.path.join(self.manager.work_path, "record.txt"))
    
    def md_operation(self):
        current_directory = os.getcwd()
        record_file_path = os.path.join(self.manager.work_path, "record.txt")  # For logging
        folders = [os.path.join(current_directory, f) for f in os.listdir(current_directory) if os.path.isdir(f) and f.endswith('d')]
    
        for folder in folders:
            poscar_path = os.path.join(folder, 'POSCAR')
            if not os.path.exists(poscar_path):
                print(f"POSCAR file not found in {folder}, skipping...")
                continue
    
            with open(poscar_path, 'r') as poscar_file:
                poscar_lines = poscar_file.readlines()
                if len(poscar_lines) >= 7:
                    element_names = poscar_lines[5].strip().split()
                    if 'Ni' in element_names or 'Co' in element_names:
                        print(f"Skipping MD operation for {folder} due to presence of Ni or Co in POSCAR.")
                        with open(record_file_path, "a") as record_file:
                            record_file.write(f"Skipped MD operation for {folder} due to presence of Ni or Co.\n")
                        continue
    
            opt_folder = os.path.join(folder, 'OPT')
            md_folder_path = os.path.join(folder, 'MD')
            if os.path.exists(opt_folder):
                if os.path.exists(md_folder_path):
                    if self.is_md_complete(md_folder_path):
                        print(f"MD calculations already complete at {md_folder_path}.")
                    else:
                        print(f"MD folder already exists at {md_folder_path}, but job is incomplete. Re-running MD job.")
                        self.prepare_and_submit_md_job(opt_folder, md_folder_path, record_file_path)
                else:
                    if ConvergenceChecker().check_scf_convergence(opt_folder):
                        self.prepare_and_submit_md_job(opt_folder, md_folder_path, record_file_path)
                    else:
                        with open(record_file_path, "a") as record_file:
                            record_file.write(f"OPT results in {opt_folder} are not converged.\n")
                        print(f"OPT results in {opt_folder} are not converged.")
                    
    def is_md_complete(self, md_folder_path):
        outcar_path = os.path.join(md_folder_path, 'OUTCAR')
        try:
            with open(outcar_path, 'r') as file:
                for line in file:
                    if 'Total CPU time used' in line:
                        return True
        except FileNotFoundError:
            return False
        return False

    def prepare_and_submit_md_job(self, opt_folder, md_folder_path, record_file_path):
        os.makedirs(md_folder_path, exist_ok=True)
        shutil.copy(os.path.join(opt_folder, 'CONTCAR'), os.path.join(md_folder_path, 'POSCAR'))
        self.atomic_structure_manager.remove_zero_velocities(os.path.join(md_folder_path, 'POSCAR'))
        shutil.copy(os.path.join(opt_folder, 'KPOINTS'), md_folder_path)
        shutil.copy(os.path.join(opt_folder, 'POTCAR'), md_folder_path)
        shutil.copy(os.path.join(opt_folder, 'sub.vasp'), md_folder_path)
        
        if os.path.exists(self.md_incar_path):
            shutil.copy(self.md_incar_path, os.path.join(md_folder_path, 'INCAR'))
            self.change_magmom_incar(md_folder_path)  # Update MAGMOM before submitting the job
        else:
            incar_path = os.path.join(opt_folder, 'INCAR')
            new_incar_path = os.path.join(md_folder_path, 'INCAR')
            with open(incar_path, 'r') as file:
                lines = file.readlines()

            with open(new_incar_path, 'w') as file:
                for line in lines:
                    if 'IBRION' not in line and 'POTIM' not in line and 'NSW' not in line:
                        file.write(line)

            # Append MD-specific settings
            with open(new_incar_path, 'a') as file:
                file.write('NSW = 300\n')
                file.write('IBRION = 0\n')
                file.write('POTIM = 1\n')
                file.write('SMASS = 0\n')
                file.write('MDALGO = 2\n')
                file.write('TEBEG = 500\n')
                file.write('TEEND = 500\n')
                file.write('NBLOCK = 1\n')

        self.wait_and_submit(md_folder_path, record_file_path)

    def wait_and_submit(self, job_directory, record_file_path):
        """ Submits a job only when job slots are available under the user's limit. """
        os.chdir(job_directory)
        if not os.path.exists('sub.vasp'):
            shutil.copy(self.manager.vasp_sub_path, ".")
        while True:
            try:
                result = subprocess.run(['squeue', "-u", "gengzi", '--format=%A'], capture_output=True, text=True)
                num_jobs = len(result.stdout.strip().split('\n')) - 1
                if num_jobs < int(self.manager.max_jobs):
                    subprocess.run(['sbatch', 'sub.vasp'])
                    print(f"Job submitted in directory: {os.getcwd()}")
                    with open(record_file_path, "a") as record_file:
                        record_file.write(f"Submitted job in {os.getcwd()}\n")
                    break
                else:
                    print(f"Job queue full, waiting {self.manager.sleep_time} seconds...")
                    time.sleep(int(self.manager.sleep_time))
            except subprocess.CalledProcessError as e:
                print(f"Error during job submission in {os.getcwd()}: {e}")
                with open(record_file_path, "a") as record_file:
                    record_file.write(f"Failed to submit job in {os.getcwd()} due to error: {e}\n")
                break


class AtomicStructureManager:
    def __init__(self):
        pass

    def check_atom_number(self, filepath):
        try:
            with open(filepath, 'r') as f:
                lines = f.readlines()
                if len(lines) >= 7:
                    seventh_line = lines[6]
                    numbers = [int(num) for num in seventh_line.split()]
                    if sum(numbers) > 300:
                        print(f"Skipped file with more than 300 atoms: {filepath}")
                        return False
        except FileNotFoundError:
            print(f"File {filepath} not found for atomic number check.")
        return True

    def remove_zero_velocities(self, poscar_path):
        """ Removes lines containing '0 0 0' from the POSCAR file, assuming these are velocity entries. """
        try:
            with open(poscar_path, 'r') as file:
                lines = file.readlines()

            new_lines = [line for line in lines if not line.strip().endswith('0.00000000E+00  0.00000000E+00  0.00000000E+00')]

            with open(poscar_path, 'w') as file:
                file.writelines(new_lines)

            print(f"Removed zero velocities from {poscar_path}")
        except FileNotFoundError:
            print(f"File not found: {poscar_path}")
        except Exception as e:
            print(f"Error processing {poscar_path}: {e}")

class ConvergenceChecker:
    def check_convergence(self, directory):
        outcar_path = os.path.join(directory, 'OUTCAR')
        try:
            with open(outcar_path, 'r') as file:
                if 'reached required accuracy' in file.read():
                    return True
        except FileNotFoundError:
            print(f"OUTCAR file not found in {directory}")
        return False

    def check_scf_convergence(self, directory):
        outcar_path = os.path.join(directory, 'OUTCAR')
        oszicar_path = os.path.join(directory, 'OSZICAR')
        scf_converged = False
        try:
            with open(outcar_path, 'r') as file:
                for line in file:
                    if 'EDIFF is reached' in line:
                        scf_converged = True
                        print(f"SCF converged in {directory}")
                        return scf_converged
        except FileNotFoundError:
            print(f"OUTCAR file not found in {directory}")
            return False  # Return immediately if file is not found
    
        if not scf_converged:
            # Check OSZICAR for additional clues
            try:
                with open(oszicar_path, 'r') as file:
                    lines = file.readlines()
                    if lines:
                        last_line = lines[-1]
                        if 'E0=' in last_line:
                            print(f"SCF did not converge in {directory}")
                        else:
                            print(f"Potential error in SCF process, check OSZICAR in {directory}")
                    else:
                        print(f"OSZICAR file is empty, indicating an early termination or error in {directory}")
            except FileNotFoundError:
                print(f"OSZICAR file not found in {directory}")
            return False
    
        return scf_converged

    def md_check(self, md_folder):
        outcar_path = os.path.join(md_folder, 'OUTCAR')
        record_file_path = os.path.join(md_folder, "record.txt")  # Assuming each folder may have its own record file
        try:
            with open(outcar_path, 'r') as file:
                lines = file.readlines()
                total_cpu_time_used = any('Total CPU time used' in line for line in lines)
                scf_cycles = sum('EDIFF is reached' in line for line in lines)

                if total_cpu_time_used:
                    print(f"MD calculations completed at {md_folder}. SCF cycles converged: {scf_cycles}")
                else:
                    print(f"MD calculations incomplete at {md_folder}. SCF cycles converged: {scf_cycles}")

                with open(record_file_path, "a") as record_file:
                    record_file.write(f"MD status checked at {md_folder}, SCF cycles: {scf_cycles}\n")
        except FileNotFoundError:
            print(f"{outcar_path} not found in {md_folder}")
            with open(record_file_path, "a") as record_file:
                record_file.write(f"{outcar_path} file not found in {md_folder}\n")

class DPDataProcessor:
    def process(self, prefix):
        outcar_files = glob("**/OUTCAR", recursive=True)
        ms = MultiSystems()
        for outcar_file in outcar_files:
            try:
                ls = LabeledSystem(outcar_file)
                ms.append(ls)
            except Exception as e:
                print(f"Failed to load {outcar_file}: {e}")
        ms.to_deepmd_raw("dpdata-" + prefix)
        ms.to_deepmd_npy("dpdata-" + prefix)
        print(f"You have converted {len(outcar_files)} OUTCARs into deepmd-readable format")

def parse_arguments():
    parser = argparse.ArgumentParser(description="Manage VASP calculations and processing with various commands.")
    parser.add_argument("prefix", help="Prefix for files to process.")
    parser.add_argument("operation", choices=["opt", "md", "convcheck", "mdcheck", "dpdata"], help="Specify the operation to perform.")
    parser.add_argument("-r", "--restart", choices=["restart"], default=None, help="Restart tasks if needed.")
    return parser.parse_args()

def read_parameters_from_file(file_path="input"):
    parameters = {}
    with open(file_path, 'r') as file:
        for line in file:
            if line.strip():
                key, value = line.strip().split(":")
                parameters[key.strip()] = value.strip()
    return (
        parameters.get('work_path'),
        parameters.get('opt_INCAR_path'),
        parameters.get('md_INCAR_path'),
        parameters.get('vasp_sub_path'),
        parameters.get('structre_db_path'),
        parameters.get('max_jobs'),
        parameters.get('sleep_time'),
    )

def main():
    args = parse_arguments()
    work_path, opt_incar_path, md_incar_path, vasp_sub_path, structre_db_path, max_jobs, sleep_time = read_parameters_from_file()
    manager = VASPManager(work_path, vasp_sub_path, structre_db_path, int(max_jobs), int(sleep_time), args.prefix)
    atomic_manager = AtomicStructureManager()  
    job_manager = VASPJobManager(manager, atomic_manager, opt_incar_path, md_incar_path)
    convergence_checker = ConvergenceChecker()

    if args.operation == "opt" and manager.copy_and_enter_folder(args.prefix):
        job_manager.opt_operation()
    elif args.operation == "md" and manager.copy_and_enter_folder(args.prefix):
        job_manager.md_operation()
    elif args.operation == "convcheck" and manager.copy_and_enter_folder(args.prefix):
        current_directory = os.getcwd()
        folders = [os.path.join(current_directory, f) for f in os.listdir(current_directory) if os.path.isdir(f) and f.endswith('d')]
        for folder in folders:
            opt_folder = os.path.join(folder, 'OPT')
            if os.path.exists(opt_folder):
                convergence_checker.check_scf_convergence(opt_folder)
    elif args.operation == "mdcheck" and manager.copy_and_enter_folder(args.prefix):
        current_directory = os.getcwd()
        folders = [os.path.join(current_directory, f) for f in os.listdir(current_directory) if os.path.isdir(f) and f.endswith('d')]
        for folder in folders:
            md_folder = os.path.join(folder, 'MD')
            if os.path.exists(md_folder):
                convergence_checker.md_check(md_folder)
    
if __name__ == "__main__":
    main()

