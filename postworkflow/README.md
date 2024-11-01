# postworkflow

Aims of flowopt.py and flowts.py are reducing the ionic steps of optimization and transition state calculations, to avoid computational cost of SCF steps

## Table of Contents

- [0. Prerequisites](#0-prerequisites)
- [1. flowopt.py](#1-flowoptpy)
  - [1.1 flowopt.py File Structure](#11-flowoptpy-file-structure)
  - [1.2 flowopt.py Usage](#12-flowoptpy-usage)
    - [1.2.1 Command-Line Arguments](#121-command-line-arguments)
    - [1.2.2 Example Command](#122-example-command)
  - [1.3 Workflow Description for flowopt.py](#13-workflow-description-for-flowoptpy)
  - [1.4 Resuming from Checkpoints](#14-resuming-from-checkpoints)
  - [1.5 Notes for flowopt.ts](#15-notes-for-flowoptts)
- [2. flowts.py](#2-flowtspy)
  - [2.1 flowts.py File Structure](#21-flowtspy-file-structure)
  - [2.2 flowts.py Usage](#22-flowtspy-usage)
    - [2.2.1 Command-Line Arguments](#221-command-line-arguments)
    - [2.2.2 Example Command](#222-example-command)
  - [2.3 Workflow Description for flowts.py](#23-workflow-description-for-flowtspy)
  - [2.4 Resuming from Checkpoints](#24-resuming-from-checkpoints)
- [3. flowopt.py based on GemNet-OC Pretrained Model](#3-flowoptpy-based-on-gemnet-oc-pretrained-model)
- [4. flowts.py based on GemNet-OC Pretrained Model](#3-flowtspy-based-on-gemnet-oc-pretrained-model)
- [5. License](#4-license)
- [6. Acknowledgements](#5-acknowledgements)


## 0. Prerequisites

- Python 3 (version > 3.10)
- ASE (Atomic Simulation Environment, version > 3.22)
- VASP (Vienna Ab initio Simulation Package)
- DeePMD-kit (version > 3.0.0a1)
- dpdata (version > 0.2.18)
- VASPKIT
- SLURM (for job scheduling)


## 1. flowopt.py

`flowopt.py` script performs an iterative optimization and fine-tuning workflow using ASE (Atomic Simulation Environment) and VASP (Vienna Ab initio Simulation Package) with DeePMD-kit. The workflow involves optimizing atomic structures, running VASP calculations, generating new datasets for finetune, and fine-tuning machine learning models iteratively. After completing the specified number of iterations, a final VASP optimization step is performed to get accurate configuration and energy.

![Energy per ionic steps calculated by VASP and CLAM+VASP for OPT calulation](../docs/opt_energy_change_comparison.png)
![Force per ionic steps calculated by VASP and CLAM+VASP for CI-NEB calulation](../docs/opt_force_change_comparison.png)

## 1.1 flowopt.py File Structure

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

## 1.2 flowopt.py Usage

1. Ensure all prerequisites are installed and configured correctly.
2. Prepare the necessary input files (`INCAR`, `KPOINTS`, `sub.vasp`, `finetune1.json`, `sub.dp`, `POTCAR`) in the `utils` directory.
3. Place the initial atomic structure file (`POSCAR`) in the working directory.
4. Run the script with the desired parameters.

### 1.2.1 Command-Line Arguments

This script supports the following command-line arguments:

- `--num_iterations` (type: `int`, default: `4`): The number of optimization and finetuning iterations to run. Each iteration consists of optimizing the structure using ASE, running a VASP calculation, generating a new dataset, and finetuning the model.
- `--steps_per_iteration` (type: `int`, default: `100`): The number of training, i.e. "numb_steps" parameter in finetune.json, steps to perform during each finetuning iteration. This value will be multiplied by the iteration number to determine the total number of training steps for each iteration. For example: if --steps_per_iteration is 100, numb_steps will be 100, 200, 300, 400 in finetune1, finetune2, finetune3, and finetune4, respectively.
- `--fixed_atoms` (type: `int`, default: `0`): The number of atoms to fix in position during the structure optimization step. The lowest atoms based on their z-coordinate will be fixed.
- `--iffinal` (type: `bool`, default: `True`): Whether to perform the final optimization step. This step is carried out after the last iteration and includes an additional structure optimization and VASP calculation. Usually, it is necessary to keep `--iffinal` `True` to get an accurate final energy and configuration.
- `--skip_first_opt` (type: `bool`, default: `False`): Whether to skip the first optimization step. If set to `True`, the script will directly use the initial `POSCAR` file for the first VASP calculation without performing the initial structure optimization by ASE with DP potential. Usually, if the training set of large atomic model includes similar structure of POSCAR, it recommands to use `--skip_first_opt` `False` to get faster result. If the POSCAR is a new system, it recommands to use `--skip_first_opt` `True` to do VASP optimization of several steps to get some dataset for finetune.
- `--fmax` (type: `float`, default: `0.2`): The maximum force criterion for structure optimization by ASE with DP potential. The optimization will stop when the forces on all atoms are below this threshold.

### 1.2.2 Example Command

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

## 1.3 Workflow Description for flowopt.py
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

## 1.4 Resuming from Checkpoints

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

## 1.5 Notes for flowopt.ts

* Ensure the VASP and DeePMD-kit executables are accessible in your environment.
* Customize the `INCAR`, `KPOINTS`, `POTCAR`, and submission scripts (`sub.vasp`, `sub.dp`) as needed for your specific system and requirements.
* Adjust the polling interval in `wait_for_job_completion` if necessary.


## 2. flowts.py

`flowts.py` script performs an iterative NEB calculation and fine-tuning workflow using ASE (Atomic Simulation Environment) and VASP (Vienna Ab initio Simulation Package) with DeePMD-kit. The workflow involves ASE-CINEB calculations, running VASP-CINEB calculations, generating new datasets for finetune, and fine-tuning machine learning models iteratively. After completing the specified number of iterations, a final VASP-CINEB step is performed to get accurate transition states configuration and energy.

![Energy per ionic steps calculated by VASP and CLAM+VASP for CI-NEB calulation](../docs/ts_energy_change_comparison.png)
![Force per ionic steps calculated by VASP and CLAM+VASP for CI-NEB calulation](../docs/ts_force_change_comparison.png)


## 2.1 flowts.py File Structure
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

## 2.2 flowts.py Usage

1. Ensure all prerequisites are installed and configured correctly.
2. Prepare the necessary input files (`INCAR`, `KPOINTS`, `sub.vasp`, `finetune1.json`, `sub.dp`, `POTCAR`) in the `utils` directory.
3. Place the initial and final atomic structure files (`initial_structure`, `final_structure`) in the working directory.
4. Place the initial and final opt OUTCAR (`initial_opt_OUTCAR`, `final_opt_OUTCAR`) in the working directory.
5. Run the script with the desired parameters.

### 2.2.1 Command-Line Arguments

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

### 2.2.2 Example Command

To run the script, use the following command:

```sh
python flowts.py initial_structure final_structure model_path initial_outcar final_outcar [--num_steps NUM] [--n_images NUM] [--fmax NUM] [--interpolation METHOD] [--spring_constant NUM] [--steps_per_iteration NUM] [--apply_constraint BOOL] [--skip_first_neb BOOL]
```
```sh
nohup python ./flowts.py POSCARis POSCARfs ./frozen_model.pth OUTCARis OUTCARfs &
```


## 2.3 Workflow Description for flowts.py 
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

## 2.4 Resuming from Checkpoints

The script records the progress of each step in `record.txt`. If the script is interrupted, it can resume from the last completed step by reading the `record.txt` file.
examples of `record.txt`
```
vaspopt1
finetune1
freeze1
vaspopt-final

```

## 3. flowopt.py based on GemNet-OC Pretrained Model

This section describes how to use `flowopt.py` based on the GemNet-OC pretrained model.

You can find more details and download the model checkpoints from the [FairChem model checkpoints page](https://fair-chem.github.io/core/model_checkpoints.html).

When you use `optoc/flowopt.py`, please download the FairChem model checkpoints, such as:[gnoc_oc22_oc20_all_s2ef.pt](https://dl.fbaipublicfiles.com/opencatalystproject/models/2022_09/oc22/s2ef/gnoc_oc22_oc20_all_s2ef.pt), and put it into `utils` directory.

## 4. flowts.py based on GemNet-OC Pretrained Model


## 5. License

This project is licensed under the LGPL-3.0 License.

## 6. Acknowledgements

* ASE: [https://wiki.fysik.dtu.dk/ase/](https://wiki.fysik.dtu.dk/ase/)
* VASP: [https://www.vasp.at/](https://www.vasp.at/)
* DeePMD-kit: [https://github.com/deepmodeling/deepmd-kit](https://github.com/deepmodeling/deepmd-kit)
* dpdata: [https://github.com/deepmodeling/dpdata](https://github.com/deepmodeling/dpdata)
* VASPKIT: [https://vaspkit.com/](https://vaspkit.com/)
* SLURM: [https://slurm.schedmd.com/](https://slurm.schedmd.com/)
