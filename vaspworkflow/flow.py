import time
import argparse
import os
import shutil
import time
from glob import glob
import subprocess
import dpdata as dp
from ase.io import read
from tqdm import tqdm

class VASPManager:
    """
    Class to manage VASP calculations

    Args:
    work_path: path to the working directory
    vasp_sub_path: path to the VASP submission script
    structre_db_path: path to the structure database
    max_jobs: maximum number of jobs to run at once
    sleep_time: time to sleep between checks
    user_name: user name for the VASP submission script
    prefix: prefix for the structure files
    """
    def __init__(self, work_path, vasp_sub_path, structre_db_path, max_jobs, sleep_time, user_name, prefix):
        self.work_path = work_path
        self.vasp_sub_path = vasp_sub_path
        self.structre_db_path = structre_db_path
        self.max_jobs = max_jobs
        self.sleep_time = sleep_time
        self.user_name = user_name
        self.files_to_copy = glob(os.path.join(self.structre_db_path, f"*{prefix}*"))

    def copy_and_enter_folder(self, prefix):
        """ Copy files from the structure database to the current working directory and enter the folder """
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
    """
    Class to manage VASP jobs

    Args:
    manager: VASPManager object
    atomic_structure_manager: AtomicStructureManager object
    opt_incar_path: path to the OPT INCAR file
    md_incar_path: path to the MD INCAR file
    user_name: user name for the VASP submission script 
    """
    def __init__(self, manager, atomic_structure_manager, opt_incar_path, md_incar_path, user_name):
        self.manager = manager
        self.atomic_structure_manager = atomic_structure_manager
        self.opt_incar_path = opt_incar_path  # Store the path for use in operations
        self.md_incar_path = md_incar_path
        self.user_name = user_name

    def change_magmom_incar(self, directory):
        """ Change the MAGMOM and other magnetic parameters in the INCAR file """
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
        """ Perform OPT operation """
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
        """ Prepare the OPT folder """
        shutil.move(poscar_path, os.path.join(new_folder, "POSCAR"))
        shutil.copy(os.path.join(new_folder, "POSCAR"), opt_folder)
        time.sleep(0.2)
        os.chdir(opt_folder)
        os.system(f'pwd')
        os.system(f'echo -e "102\n2\n0" | vaspkit')
        os.chdir("../..")
        shutil.copy(self.opt_incar_path, os.path.join(opt_folder, "INCAR"))
        self.change_magmom_incar(opt_folder)  # Update MAGMOM before submitting the job
        self.wait_and_submit(opt_folder, os.path.join(self.manager.work_path, "record.txt"))
    
    def md_operation(self):
        """ Perform MD operation """
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
        """ Check if MD calculations are complete"""
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
        """ Prepare and submit MD job """
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
                result = subprocess.run(['squeue', "-u", self.user_name, '--format=%A'], capture_output=True, text=True)
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
    """
    Class for managing atomic structures.
    """
    def __init__(self):
        pass

    def check_atom_number(self, filepath):
        """ Check if the number of atoms in the POSCAR file is greater than 300. """
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
    """ Class for checking convergence of calculations."""
    def check_convergence(self, directory):
        """ Check if the calculation has converged. """
        outcar_path = os.path.join(directory, 'OUTCAR')
        try:
            with open(outcar_path, 'r') as file:
                if 'reached required accuracy' in file.read():
                    return True
        except FileNotFoundError:
            print(f"OUTCAR file not found in {directory}")
        return False

    def check_scf_convergence(self, directory):
        """check if scf converged """
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
                            print(f"Unkown error in SCF process, check OSZICAR in {directory}")
                    else:
                        print(f"OSZICAR file is empty, indicating an early termination or error in {directory}")
            except FileNotFoundError:
                print(f"OSZICAR file not found in {directory}")
            return False
    
        return scf_converged

    def md_check(self, md_folder):
        """ Check if the MD calculations have completed. """
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
    """
    Class for processing DP data.

    Attributes:
        step (int): Step size for data processing.
        test_size (float): Percentage of data to be used for testing.
        prefix (str): Prefix for data files.
        fmt (str): Format of data files.
    """
    def __init__(self,step=1,test_size=0.1,prefix='./data',fmt='npy_single') -> None:
        self.step=step
        self.test_size=test_size
        self.prefix=prefix
        self.fmt=fmt
    def process(self):
        """ Process DP data."""
        def find_folders_with_outcar(search_path='./POSCAR',label='OUTCAR'):
            """ Find folders containing 'OUTCAR' files. """
            seen_folders = []  
           # Traverse all 'OUTCAR' files under the specified search path.
            for outcar_path in glob(os.path.join(search_path, '**/MD/', label)):  
             # Get the folder path where the OUTCAR file resides.
                folder_path = os.path.dirname(outcar_path)  
            # If the folder path is not one that we have already processed (to avoid duplicate prints). 
                if folder_path not in seen_folders:  
                  # Add the folder path to the collection of folders that have been processed.
                    seen_folders.append(folder_path)
            return seen_folders
        
        import re
        def extract_nonzero_elements(compound_string):
            """ Extract non-zero elements from a compound string."""
            elements = re.findall(r'([A-Z][a-z]*)(\d+)', compound_string)
            nonzero_elements = [(element, int(quantity)) for element, quantity in elements if int(quantity) != 0]
            formula = ''.join(f"{element}{quantity}" if quantity != 1 else element for element, quantity in nonzero_elements)         
            return formula
        
        d=find_folders_with_outcar()
        ms=dp.MultiSystems()
        for i in tqdm(d):
            try:
                data=dp.LabeledSystem(f'{i}/OUTCAR',fmt='vasp/outcar',step=int(self.step))
            except Exception as e:
                print(f"system {i} failed to load,please check the OUTCAR file")
                continue
            ms.append(data)
        # split data into train and validation
        dataset=ms.train_test_split(test_size=float(self.test_size))
        if not os.path.exists(self.prefix):
            os.makedirs(self.prefix)
        else:
            pass
        if not os.path.exists(os.path.join(self.prefix,'train')):
            os.makedirs(os.path.join(self.prefix,'train'))
        if not os.path.exists(os.path.join(self.prefix,'val')):
            os.makedirs(os.path.join(self.prefix,'val'))

        if self.fmt=='npy_single':
            for i in dataset[0]:
                name=extract_nonzero_elements(i.formula)
                i.to_deepmd_npy(self.prefix + "/train"+f"/{name}")
                
            for i in dataset[1]:
                name=extract_nonzero_elements(i.formula)
                i.to_deepmd_npy(self.prefix + "/val"+f"/{name}")

        if self.fmt=='npy_mix':
            dataset[0].to_deepmd_npy_mixed(self.prefix + "/train")
            dataset[1].to_deepmd_npy_mixed(self.prefix + "/val")

        elif self.fmt=='lmdb':

            from fairchem.core.preprocessing import AtomsToGraphs
            import lmdb
            import pickle
            import torch
            import numpy as np
            def tolmdb(file,data):
                """ Convert data to lmdb format. note: this format is justed for s2ef task with fairchem package."""
                a2g = AtomsToGraphs(
                max_neigh=50,
                radius=6,
                r_energy=True,    # False for test data
                r_forces=True,    # False for test data
                r_distances=False,
                r_fixed=True,
            )
                # tags = raw_data[0].get_tags()

                for i in data:
                    name=extract_nonzero_elements(i.formula)
                    ase_data=i.to_ase_structure()

                # simple devide atom according to z coordinate, 1 for subsurface, 2 for surface
                    b=ase_data[0].positions[:,2]
                    # b=a[1]['pos'][:,2]
                    c1= b<np.max(b)-2 
                    c2= np.min(b)+2< b
                    c=c1 & c2
                    tag=np.where(c,1,2)


                    # ase_data[0].set_tags(tag)
                    data_objects = a2g.convert_all(ase_data, disable_tqdm=True)
                    # tags = ase_data[0].get_tags()
                    db = lmdb.open(
                    f'{self.prefix}/{file}/{name}.lmdb',
                    map_size=1099511627776 * 2,
                    subdir=False,
                    meminit=False,
                    map_async=True,
                    )
                    for fid, data1 in tqdm(enumerate(data_objects), total=len(data_objects)):
                    #assign fid
                        data1.fid = torch.LongTensor([fid])
                        data1.sid = torch.LongTensor([0])
                        data1.tags = torch.LongTensor(tag)

                        txn = db.begin(write=True)
                        txn.put(f"{fid}".encode("ascii"), pickle.dumps(data1, protocol=-1))
                        txn.commit()

                        txn = db.begin(write=True)
                        txn.put(f"length".encode("ascii"), pickle.dumps(len(data_objects), protocol=-1))
                        txn.commit()

                    db.sync()
                    db.close()
            tolmdb(file='train',data=dataset[0])
            tolmdb(file='val',data=dataset[1])


def parse_arguments():
    """ Parse command-line arguments. """
    parser = argparse.ArgumentParser(description="Manage VASP calculations and processing with various commands.")
    parser.add_argument("prefix", help="Prefix for files to process.")
    parser.add_argument("operation", choices=["opt", "md", "optcheck", "mdcheck", "dpdata",'plot'], help="Specify the operation to perform.")
    parser.add_argument("-r", "--restart", choices=["restart"], default=None, help="Restart tasks if needed.")
    return parser.parse_args()

def read_parameters_from_file(file_path="input"):
    """ Read parameters from a file. """
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
        parameters.get('user_name'),
        parameters.get('step_data',None),
        parameters.get('test_size',None),
        parameters.get('dataset_prefix',None),
        parameters.get('dataset_fmt',None),
        parameters.get('plot_name',None)

    )

def main():
    """ Main function. the entry for the program."""
    args = parse_arguments()
    work_path, opt_incar_path, md_incar_path, vasp_sub_path, structre_db_path, max_jobs, sleep_time, user_name, step,test_size, dataset_prefix, dataset_fmt, plot_name = read_parameters_from_file()
    manager = VASPManager(work_path, vasp_sub_path, structre_db_path, int(max_jobs), int(sleep_time), user_name, args.prefix)
    atomic_manager = AtomicStructureManager()  
    job_manager = VASPJobManager(manager, atomic_manager, opt_incar_path, md_incar_path, user_name)
    convergence_checker = ConvergenceChecker()

    if args.operation == "opt" and manager.copy_and_enter_folder(args.prefix):
        job_manager.opt_operation()
    elif args.operation == "md" and manager.copy_and_enter_folder(args.prefix):
        job_manager.md_operation()
    elif args.operation == "optcheck" and manager.copy_and_enter_folder(args.prefix):
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
    elif args.operation == "dpdata":
        dpdata=DPDataProcessor(step=step,test_size=test_size,prefix=dataset_prefix,fmt=dataset_fmt)
        dpdata.process()
    
    elif args.operation == "plot" and os.path.isdir(os.path.join(dataset_prefix,'train')):
        from utils.plot_weight import PeriodicTableWeightsVisualizer
        visualizer = PeriodicTableWeightsVisualizer(poscar_dir = os.path.join(dataset_prefix,'train'),fmt=dataset_fmt)
        visualizer.plot_2d(name = plot_name)
    
if __name__ == "__main__":
    main()

