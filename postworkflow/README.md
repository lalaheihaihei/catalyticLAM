# postworkflow

The purpose of `flowopt.py` and `flowts.py` is to reduce the number of ionic steps required for structure optimization and transition state (TS) calculations by using a "local fine-tuning" algorithm, thereby minimizing the computational cost of self-consistent field (SCF) steps. Our approach involves selectively incorporating a minimal amount of new DFT SCF data to fine-tune the model, enhancing its accuracy on the local potential energy surface (PES) specifically relevant to optimization or TS tasks. 

![Schematic representation of local fine-tuning on the PES. The bold black line represents the real PES, while the thin blue line indicates the uncertainty in the model-predicted PES.](../docs/Schematic.png)

## Table of Contents

- [0. Prerequisites](#0-prerequisites)
- [1. flowopt.py based on GemNet-OC Pretrained Model](#1-flowoptpy-based-on-gemnet-oc-pretrained-model)
  - [1.1 flowopt.py File Structure](#11-flowoptpy-file-structure)
  - [1.2 flowopt.py Usage](#12-flowoptpy-usage)
  - [1.3 Workflow Description for flowopt.py](#13-workflow-description-for-flowoptpy)
- [2. flowts.py based on GemNet-OC Pretrained Model](#2-flowtspy-based-on-gemnet-oc-pretrained-model)
- [3. flowopt.py based on DeePMD-kit Pretrained Model](#3-flowoptpy)
  - [3.1 flowopt.py File Structure](#31-flowoptpy-file-structure)
  - [3.2 flowopt.py Usage](#32-flowoptpy-usage)
    - [3.2.1 Command-Line Arguments](#321-command-line-arguments)
    - [3.2.2 Example Command](#322-example-command)
  - [3.3 Workflow Description for flowopt.py](#33-workflow-description-for-flowoptpy)
  - [3.4 Resuming from Checkpoints](#34-resuming-from-checkpoints)
  - [3.5 Notes for flowopt.ts](#35-notes-for-flowoptts)
- [4. flowts.py based on DeePMD-kit Pretrained Model](#4-flowtspy)
  - [4.1 flowts.py File Structure](#41-flowtspy-file-structure)
  - [4.2 flowts.py Usage](#42-flowtspy-usage)
    - [4.2.1 Command-Line Arguments](#421-command-line-arguments)
    - [4.2.2 Example Command](#422-example-command)
  - [4.3 Workflow Description for flowts.py](#43-workflow-description-for-flowtspy)
  - [4.4 Resuming from Checkpoints](#44-resuming-from-checkpoints)
- [5. License](#5-license)
- [6. Acknowledgements](#6-acknowledgements)

## 0. Prerequisites

- Python 3 (version > 3.10)
- ASE (Atomic Simulation Environment, version > 3.22)
- VASP (Vienna Ab initio Simulation Package)
- fairchem (version 1.0.0) or DeePMD-kit (version > 3.0.0a1)
- dpdata (optional, version > 0.2.18)
- SLURM (for job scheduling)

## 1. flowopt.py based on GemNet-OC Pretrained Model

This section describes how to use the `./optoc/flowopt.py` script based on the [Gemnet-OC](https://github.com/FAIR-Chem/fairchem) pretrained model.

When you use `./optoc/flowopt.py`, please download the FairChem model checkpoints, such as: [gnoc_oc22_oc20_all_s2ef.pt](https://dl.fbaipublicfiles.com/opencatalystproject/models/2022_09/oc22/s2ef/gnoc_oc22_oc20_all_s2ef.pt) or our pretrained model [best_checkpoint.pt] check points from release(Coming soon) and put it into `utils` directory.

`flowopt.py` script performs an iterative optimization and fine-tuning workflow using ASE (Atomic Simulation Environment), VASP (Vienna Ab initio Simulation Package), and [Gemnet-OC](https://github.com/FAIR-Chem/fairchem). The workflow involves optimizing atomic structures, running VASP calculations, generating new datasets for finetune, and fine-tuning machine learning models iteratively. After completing the specified number of iterations, a final VASP optimization step (optional) is performed to get accurate configuration and energy.

**<div style="display: flex; justify-content: space-between;">**
<img src="../docs/opt_energy_change_comparison.png" alt="Energy per ionic steps calculated by VASP and CLAM+VASP for OPT calculation" width="47%">
<img src="../docs/opt_force_change_comparison.png" alt="Force per ionic steps calculated by VASP and CLAM+VASP for CI-NEB calculation" width="45%">

</div>

## 1.1 flowopt.py File Structure

```
.
├── flowopt.py # Main script for the workflow
├── POSCAR # Initial structure for opt
├── utils/
│ ├── INCAR # VASP input file for OPT
│ ├── POTCAR # VASP pseudopotentials file
│ ├── KPOINTS # VASP k-points file
│ ├── sub.vasp # VASP submission script for slurm
│ ├── best_checkpoint.pt # pretrained checkpoint of Gemnet-OC
│ ├── finetune1.yml # Fine-tuning input file of Gemnet-OC
│ ├── base.yml # Fine-tuning input file of Gemnet-OC
│ ├── main.py # Fine-tuning input file of Gemnet-OC
│ └── sub.oc # Gemnet-OC submission script for slurm
└── record.txt # optional. Record the progress of the workflow for restart
```

## 1.2 flowopt.py Usage

1. Ensure all prerequisites are installed and configured correctly.
2. Prepare the necessary input files (`INCAR`, `KPOINTS`, `POTCAR `,  `finetune1.yml`, `base.yml`, `main.py`, `sub.vasp`, `sub.oc`, `best_checkpoint.pt`) in the `./utils` directory.
3. Place the initial atomic structure file (`POSCAR`) in the working directory.
4. Run the script with the desired parameters.

### flowopt.py Runing Command

To run the script, use the following command:

```sh
python flowopt.py [--num_iterations NUM] [--fixed_atoms NUM] [--iffinal BOOL] [--do_first_aseopt BOOL] [--fmax NUM] [--nsw NUM]
```

Usually, this script needs to be run in the background using `nohup` comand, for example:

```bash
nohup python flowopt.py --num_iterations 2 --fixed_atoms 72 --iffinal true --do_first_aseopt true --fmax 0.05 &
```

This will return a job number when the script runs in the background. You can find this number using the `ps -ef` command. If you want to stop the job, use `kill [job number]`. Note that, currently, this script only supports the SLURM queue system.

### Command-Line Arguments

This script supports the following command-line arguments:

* `--num_iterations` (type: int, default: 4): The number of iterations to run. Each iteration consists of structure optimization using ASE, DFT calculations using VASP, and fine-tuning with Gemnet-OC. For simple systems, 1 iteration is often sufficient.
* `--fixed_atoms` (type: int, default: 0): The number of atoms to keep fixed during the structure optimization step. The atoms with the lowest z-coordinates will be fixed in place.
  `--iffinal` (type: bool, default: true): Determines whether to perform a final DFT optimization step. This step is carried out after the last iteration and includes an additional VASP-based optimization for the final structure. It is generally recommended to keep `--iffinal` set to `true` to obtain accurate final energy and configuration values. Setting it to `false` allows a faster, rougher structure.
  `--do_first_aseopt` (type: bool, default: true): Specifies whether to perform an initial ASE optimization step using the pretrained model (`utils/best_checkpoint.pt`). If set to `false`, the script will skip the initial optimization and directly use the initial POSCAR file for the first VASP calculation. When the pretrained model includes structures similar to POSCAR, setting this to `true` can help reach a near-minimum configuration more quickly. If POSCAR represents a completely new system, it is recommended to set this to `false` to start with VASP optimization and collect data for fine-tuning.
* `--fmax` (type: float, default: 0.2): The maximum force criterion for structure optimization using ASE with the pretrained model. The optimization stops when the forces on all atoms are below this threshold or if the ASE optimization reaches 30 steps, preventing unreasonable results.
* `--nsw` (type: int, default: 5): The `NSW` parameter in the INCAR file, which sets the number of DFT ionic optimization steps to perform in each iteration. Typically, a value of 5 is a good choice.

If you want to remove all output files, use the following command:

```bash
python ./flowopt.py --clean
```

## 1.3 Workflow Description for flowopt.py

```mermaid
graph TD
A[start] -->|do_first_aseopt?| X{do first ASE?}
X -.yes.-> B[ase optimization]
X -.no.-> C[VASP optimization, NSW=3]
B --get CONTCAR--- C[VASP optimization by NSW step]
C  -->D[collect data]
D --> E[finetune new data by Gemnet-OC]
E --> F{num_iterations?}
F -.yes.-> H[final ASE optimization]
F -.no.-> B
H --get CONTCAR--> I(if_final?)
I -.yes.-> K[final VASP optimization]
I -.no.-> J(END)
K --> J
```

The script follows these main steps:

1. **Optimization of Atomic Structures by ASE** :
   * Optimize the atomic structure using pretrained or finetuned model by ASE.
   * Save the optimized structure as `CONTCAR-ase-{iteration}`.
2. **Preparation for VASP Calculation** :
   * Create a new directory for each iteration (`opt{iteration}`).
   * Copy the optimized structure and necessary VASP input files to the new directory.
   * If not the first iteration, copy `WAVECAR` and `CHGCAR` from the previous iteration.
3. **Run VASP Calculation** :
   * Submit the VASP job and retrieve the job ID.
   * Wait for the VASP job to complete.
4. **Generate New Dataset** :
   * Generate new dataset from the VASP output file `OUTCAR`, and add it to `./output_database.db`.
5. **Fine-Tuning of Models** :
   * Update the fine-tuning input file (`finetune1.yml`), use `./output_database.db` dataset.
   * Submit the fine-tuning job and retrieve the job ID.
   * Wait for the fine-tuning job to complete.
6. **Final Optimization Step** :
   * After completing all iterations, run a final ASE optimization step.
   * Create a new directory (`optfinal`) and copy the final optimized structure in it.
   * Modify the `INCAR` file to set `NSW = 300`.
   * if `--iffinal true`Submit the final VASP job and wait for its completion.

## 2. flowts.py based on GemNet-OC Pretrained Model

## 3. flowopt.py

`flowopt.py` script performs an iterative optimization and fine-tuning workflow using ASE (Atomic Simulation Environment) and VASP (Vienna Ab initio Simulation Package) with DeePMD-kit. The workflow involves optimizing atomic structures, running VASP calculations, generating new datasets for finetune, and fine-tuning machine learning models iteratively. After completing the specified number of iterations, a final VASP optimization step is performed to get accurate configuration and energy.

## 3.1 flowopt.py File Structure

```
.
├── flowopt.py # Main script for the workflow
├── POSCAR # Initial structure for opt
├── frozen_model.pth # Pretrained model needed when --skip_first_opt false
├── utils
│ ├── INCAR # VASP input file for OPT `NSW = 3`
│ ├── POTCAR # VASP pseudopotentials file
│ ├── KPOINTS # VASP k-points file
│ ├── sub.vasp # VASP submission script for slurm
│ ├── finetune1.json # Fine-tuning configuration template
│ └── sub.dp # DeePMD-kit submission script for slurm
└── record.txt # File to record the progress of the workflow for restart
```

## 3.2 flowopt.py Usage

1. Ensure all prerequisites are installed and configured correctly.
2. Prepare the necessary input files (`INCAR`, `KPOINTS`, `sub.vasp`, `finetune1.json`, `sub.dp`, `POTCAR`) in the `utils` directory.
3. Place the initial atomic structure file (`POSCAR`) in the working directory.
4. Run the script with the desired parameters.

### 3.2.1 Command-Line Arguments

This script supports the following command-line arguments:

- `--num_iterations` (type: `int`, default: `4`): The number of optimization and finetuning iterations to run. Each iteration consists of optimizing the structure using ASE, running a VASP calculation, generating a new dataset, and finetuning the model.
- `--steps_per_iteration` (type: `int`, default: `100`): The number of training, i.e. "numb_steps" parameter in finetune.json, steps to perform during each finetuning iteration. This value will be multiplied by the iteration number to determine the total number of training steps for each iteration. For example: if --steps_per_iteration is 100, numb_steps will be 100, 200, 300, 400 in finetune1, finetune2, finetune3, and finetune4, respectively.
- `--fixed_atoms` (type: `int`, default: `0`): The number of atoms to fix in position during the structure optimization step. The lowest atoms based on their z-coordinate will be fixed.
- `--iffinal` (type: `bool`, default: `True`): Whether to perform the final optimization step. This step is carried out after the last iteration and includes an additional structure optimization and VASP calculation. Usually, it is necessary to keep `--iffinal` `True` to get an accurate final energy and configuration.
- `--skip_first_opt` (type: `bool`, default: `False`): Whether to skip the first optimization step. If set to `True`, the script will directly use the initial `POSCAR` file for the first VASP calculation without performing the initial structure optimization by ASE with DP potential. Usually, if the training set of large atomic model includes similar structure of POSCAR, it recommands to use `--skip_first_opt` `False` to get faster result. If the POSCAR is a new system, it recommands to use `--skip_first_opt` `True` to do VASP optimization of several steps to get some dataset for finetune.
- `--fmax` (type: `float`, default: `0.2`): The maximum force criterion for structure optimization by ASE with DP potential. The optimization will stop when the forces on all atoms are below this threshold.

### 3.2.2 Example Command

To run the script, use the following command:

```sh
python flowopt.py [--num_iterations NUM] [--steps_per_iteration NUM] [--fixed_atoms NUM] [--iffinal BOOL] [--skip_first_opt BOOL]
```

```bash
for in-domain (ID) system:
nohup python flowopt.py --num_iterations 3 --steps_per_iteration 200 --fixed_atoms 0 --iffinal true --skip_first_opt false --fmax 0.1 &
for out-of-domain (OOD) system:
nohup python script.py --num_iterations 5 --steps_per_iteration 200 --fixed_atoms 0 --iffinal true --skip_first_opt true --fmax 0.2 &
```

## 3.3 Workflow Description for flowopt.py

```mermaid
graph TD
A[start] -->|skip_first_opt?| X{skip first?}
X -.no.-> B[ase optimization with DP potential]
X -.yes.-> C[VASP optimization, NSW=3]
B --get CONTCAR--- C[VASP optimization by NSW=3 step]
C  -->D[collect data by dpdata]
D --> E[finetune with new data by DeePMD-kit]
E --> F{num_steps?}
F -.yes.-> H[final ASE optimization with DP]
F -.no.-> B
H --get CONTCAR--> I(if final?)
I -.yes.-> K[final VASP optimization with DFT]
I -.no.-> J(END)
K --> J
```

The script follows these main steps:

1. **Optimization of Atomic Structures by ASE** :
   * Optimize the atomic structure using catalyticLAM DeePMD-kit models.
   * Save the optimized structure as `CONTCAR-ase-{iteration}`.
2. **Preparation for VASP Calculation** :
   * Create a new directory for each iteration (`opt{iteration}`).
   * Copy the optimized structure and necessary VASP input files to the new directory.
   * If not the first iteration, copy `WAVECAR` and `CHGCAR` from the previous iteration.
3. **Run VASP Calculation** :
   * Submit the VASP job and retrieve the job ID.
   * Wait for the VASP job to complete.
4. **Generate New Dataset** :
   * Generate a new dataset from the VASP output file `OUTCAR` by dpdata.
   * Save the dataset in the specified directory.
5. **Fine-Tuning of Models** :
   * Update the fine-tuning configuration (`finetune1.json`) with new data paths and iteration-specific parameters.
   * Submit the fine-tuning job and retrieve the job ID.
   * Wait for the fine-tuning job to complete.
6. **Freezing the Model** :
   * Freeze the fine-tuned model using DeePMD-kit.
7. **Final Optimization Step** :
   * After completing all iterations, run a final optimization step.
   * Create a new directory (`optfinal`) and copy the final optimized structure.
   * Modify the `INCAR` file to set `NSW = 300`.
   * Submit the final VASP job and wait for its completion.

## 3.4 Resuming from Checkpoints

The script records the progress of each step in `record.txt`. If the script is interrupted, it can resume from the last completed step by reading the `record.txt` file.
examples of `record.txt`

```
aseopt1
vaspopt1
dpdata1
finetune1
freeze1
aseopt2
vaspopt2
dpdata2
finetune2
freeze2
```

## 3.5 Notes for flowopt.ts

* Ensure the VASP and DeePMD-kit executables are accessible in your environment.
* Customize the `INCAR`, `KPOINTS`, `POTCAR`, and submission scripts (`sub.vasp`, `sub.dp`) as needed for your specific system and requirements.
* Adjust the polling interval in `wait_for_job_completion` if necessary.

## 4. flowts.py

`flowts.py` script performs an iterative NEB calculation and fine-tuning workflow using ASE (Atomic Simulation Environment) and VASP (Vienna Ab initio Simulation Package) with DeePMD-kit. The workflow involves ASE-CINEB calculations, running VASP-CINEB calculations, generating new datasets for finetune, and fine-tuning machine learning models iteratively. After completing the specified number of iterations, a final VASP-CINEB step is performed to get accurate transition states configuration and energy.

![Energy per ionic steps calculated by VASP and CLAM+VASP for CI-NEB calulation](../docs/ts_energy_change_comparison.png)
![Force per ionic steps calculated by VASP and CLAM+VASP for CI-NEB calulation](../docs/ts_force_change_comparison.png)

## 4.1 flowts.py File Structure

```
.
├── flowts.py # Main script for the workflow
├── POSCARis # Initial structure for NEB calculation (POSCAR format)
├── POSCARfs # Final structure for NEB calculation (POSCAR format)
├── OUTCARis # Initial structure's opt OUTCAR
├── OUTCARfs # Final structure's opt OUTCAR
├── frozen_model.pth # Pretrained model
├── utils
│ ├── INCAR # VASP input file for CINEB NSW = 3
│ ├── POTCAR # VASP pseudopotentials file
│ ├── KPOINTS # VASP k-points file
│ ├── sub.vasp # VASP submission script for slurm
│ ├── finetune1.json # Fine-tuning configuration template
│ └── sub.dp # DeePMD-kit submission script for slurm
└── record.txt # File to record the progress of the workflow for restart
```

## 4.2 flowts.py Usage

1. Ensure all prerequisites are installed and configured correctly.
2. Prepare the necessary input files (`INCAR`, `KPOINTS`, `sub.vasp`, `finetune1.json`, `sub.dp`, `POTCAR`) in the `utils` directory.
3. Place the initial and final atomic structure files (`initial_structure`, `final_structure`) in the working directory.
4. Place the initial and final opt OUTCAR (`initial_opt_OUTCAR`, `final_opt_OUTCAR`) in the working directory.
5. Run the script with the desired parameters.

### 4.2.1 Command-Line Arguments

This script supports the following command-line arguments:

- `initial_structure` (type: `str`): Path to the initial structure (POSCAR format).
- `final_structure` (type: `str`): Path to the final structure (POSCAR format).
- `model_path` (type: `str`): Path to the DeePMD model, ex: frozen_model.pth.
- `initial_outcar` (type: `str`): Path to the initial OUTCAR file.
- `final_outcar` (type: `str`): Path to the final OUTCAR file.
- `--num_steps` (type: `int`, default: `1`): Number of loop steps to run.
- `--n_images` (type: `int`, default: `4`): Number of intermediate images.
- `--fmax` (type: `float`, default: `0.5`): Maximum force criteria for optimization by ASE.
- `--interpolation` (choices: `['linear', 'idpp']`, default: `linear`): Interpolation method for generating intermediate images by ASE.
- `--spring_constant` (type: `float`, default: `1.0`): Spring constant for ASE CI-NEB calculations.
- `--steps_per_iteration` (type: `int`, default: `1000`): The number of training, i.e. "numb_steps" parameter in finetune1.json, steps to perform during each finetuning iteration. This value will be multiplied by the iteration number to determine the total number of training steps for each iteration. For example: if `--steps_per_iteration` is 1000, "numb_steps" will be 1000, 2000, 3000, 4000 in finetune1, finetune2, finetune3, and finetune4, respectively.
- `--apply_constraint` (type: `bool`, default: `False`): Whether to apply constraints during interpolation.
- `--skip_first_opt` (type: `bool`, default: `True`): Whether to skip the first ASE CI-NEB step. If set to `True`, the script will still to do first ASE CI-NEB step based on DP potential, but the obtained CONTCARs will not be used to the initial VASP ci-neb calculation. If set to `True`, the `POSCARis` and `POSCARfs` files will be used to generate first VASP calculation. Usually, it recommands to keep `--skip_first_opt` `True` to do VASP calculation of first step to get some dataset for finetune.

### 4.2.2 Example Command

To run the script, use the following command:

```sh
python flowts.py initial_structure final_structure model_path initial_outcar final_outcar [--num_steps NUM] [--n_images NUM] [--fmax NUM] [--interpolation METHOD] [--spring_constant NUM] [--steps_per_iteration NUM] [--apply_constraint BOOL] [--skip_first_neb BOOL]
```

```sh
nohup python ./flowts.py POSCARis POSCARfs ./frozen_model.pth OUTCARis OUTCARfs &
```

## 4.3 Workflow Description for flowts.py

```mermaid
graph TD
A[start] -->|skip_first_neb?| X{skip first?}
X -.no.-> B[ASE CINEB optimization]
X -.yes.-> C[VASP CINEB, NSW=3]
B -->|get CONTCAR0*| C
C --> D[collect data by dpdata]
D --> E[finetune with new data by DeePMD-kit]
E --> F{num_steps?}
F -.yes.-> H[final VASP CINEB optimization]
F -.no.-> B
H --> I[end]
```

The script follows these main steps:

1. **NEB Calculation (ASE)** :
   * Generates intermediate images using the initial and final structures.
   * Optimizes the NEB path using DeePMD-kit models.
   * Saves the optimized structures as `CONTCAR{n}`.
2. **Preparation for VASP Calculation** :
   * Creates a new directory for each step (`ts{step}`).
   * Copies the optimized structures and necessary VASP input files to the new directory.
3. **Run VASP NEB Calculation** :
   * Generates the `POTCAR` file using VASPKIT.
   * Submits the VASP job and retrieves the job ID.
   * Waits for the VASP job to complete.
4. **Generate New Dataset** :
   * Generates a new dataset from the VASP output files `OUTCAR`.
   * Saves the dataset in the specified directory.
5. **Fine-Tuning of Models** :
   * Updates the fine-tuning configuration (`finetune1.json`) with new data paths and iteration-specific parameters.
   * Submits the fine-tuning job and retrieves the job ID.
   * Waits for the fine-tuning job to complete.
6. **Freezing the Model** :
   * Freezes the fine-tuned model using DeePMD-kit.
7. **Final NEB Calculation** :
   * After completing all steps, runs a final NEB calculation.
   * Creates a new directory (`tsfinal`) and copies the final optimized structures.
   * Modifies the `INCAR` file to set `NSW = 300`.
   * Submits the final VASP job and waits for its completion.

## 4.4 Resuming from Checkpoints

The script records the progress of each step in `record.txt`. If the script is interrupted, it can resume from the last completed step by reading the `record.txt` file.
examples of `record.txt`

```
vaspopt1
finetune1
freeze1
vaspopt-final
```

## 5. License

This project is licensed under the LGPL-3.0 License.

## 6. Acknowledgements

* ASE: [https://wiki.fysik.dtu.dk/ase/](https://wiki.fysik.dtu.dk/ase/)
* VASP: [https://www.vasp.at/](https://www.vasp.at/)
* DeePMD-kit: [https://github.com/deepmodeling/deepmd-kit](https://github.com/deepmodeling/deepmd-kit)
* dpdata: [https://github.com/deepmodeling/dpdata](https://github.com/deepmodeling/dpdata)
* VASPKIT: [https://vaspkit.com/](https://vaspkit.com/)
* SLURM: [https://slurm.schedmd.com/](https://slurm.schedmd.com/)
  

