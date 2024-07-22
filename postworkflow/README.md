# postworkflow

Aims of flowopt.py and flowts.py are reducing the ionic steps of optimization and transition state calculations, to avoid computational cost of SCF step,

## flowopt.py

flowopt.py script performs an iterative optimization and fine-tuning workflow using ASE (Atomic Simulation Environment) and VASP (Vienna Ab initio Simulation Package) with DeepMD-kit. The workflow involves optimizing atomic structures, running VASP calculations, generating new datasets for finetune, and fine-tuning machine learning models iteratively. After completing the specified number of iterations, a final VASP optimization step is performed to get accurate configuration and energy.

![Energy per ionic steps calculated by VASP and CLAM+VASP for OPT calulation](../docs/opt_energy_change_comparison.png)
![Force per ionic steps calculated by VASP and CLAM+VASP for CI-NEB calulation](../docs/opt_force_change_comparison.png)

## flowts.py

flowts.py script performs an iterative NEB calculation and fine-tuning workflow using ASE (Atomic Simulation Environment) and VASP (Vienna Ab initio Simulation Package) with DeepMD-kit. The workflow involves ase-CINEB calculations, running VASP-CINEB calculations, generating new datasets for finetune, and fine-tuning machine learning models iteratively. After completing the specified number of iterations, a final VASP-CINEB step is performed to get accurate transition states configuration and energy.

![Energy per ionic steps calculated by VASP and CLAM+VASP for CI-NEB calulation](../docs/ts_energy_change_comparison.png)
![Force per ionic steps calculated by VASP and CLAM+VASP for CI-NEB calulation](../docs/ts_force_change_comparison.png)

## Prerequisites

- Python 3 (version > 3.10)
- ASE (Atomic Simulation Environment, version > 3.22)
- VASP (Vienna Ab initio Simulation Package)
- DeepMD-kit (version > 3.0.0a1)
- dpdata (version > 0.2.18)
- VASPkit
- SLURM (for job scheduling)

## flowopt.py File Structure
.
├── flowopt.py # Main script for the workflow
├── POSCAR # initial structure for opt
├── utils
│ ├── INCAR # VASP input file, with NSW = 3
│ ├── POTCAR # VASP pseudopotentials file
│ ├── KPOINTS # VASP k-points file
│ ├── sub.vasp # VASP submission script, slurm
│ ├── finetune1.json # Fine-tuning configuration template
│ └── sub.dp # DeepMD-kit submission script, slurm
└── record.txt # File to record the progress of the workflow for restart

## flowopt.py Usage

1. Ensure all prerequisites are installed and configured correctly.
2. Prepare the necessary input files (`INCAR`, `KPOINTS`, `sub.vasp`, `finetune1.json`, `sub.dp`, `POTCAR`) in the `utils` directory.
3. Place the initial atomic structure file (`POSCAR`) in the working directory.
4. Run the script with the desired parameters.

### Command-Line Arguments

This script supports the following command-line arguments:

- `--num_iterations` (type: `int`, default: `4`): The number of optimization and finetuning iterations to run. Each iteration consists of optimizing the structure using ASE, running a VASP calculation, generating a new dataset, and finetuning the model.
- `--steps_per_iteration` (type: `int`, default: `100`): The number of training, i.e. "numb_steps" parameter in finetune.json, steps to perform during each finetuning iteration. This value will be multiplied by the iteration number to determine the total number of training steps for each iteration. for example: if --steps_per_iteration 100, numb_steps will be 100, 200, 300, 400 in finetune1, finetune2, finetune3, finetune4, respectively.
- `--fixed_atoms` (type: `int`, default: `0`): The number of atoms to fix in position during the structure optimization step. The lowest atoms based on their z-coordinate will be fixed.
- `--iffinal` (type: `bool`, default: `True`): Whether to perform the final optimization step. This step is carried out after the last iteration and includes an additional structure optimization and VASP calculation. Usually, it is necessary to keep --iffinal true to get an accurate final energy and configuration.
- `--skip_first_opt` (type: `bool`, default: `False`): Whether to skip the first optimization step. If set to `True`, the script will directly use the initial `POSCAR` file for the first VASP calculation without performing the initial structure optimization by ASE with DP potential. Usually, if the training set of large atomic model include similar structureof POSCAR, it recommand to use --skip_first_opt false to get faster result. if the POSCAR is a new system, it recommand to use --skip_first_opt true to do VASP optimization of several step to get some dataset for finetune.
- `--fmax` (type: `float`, default: `0.2`): The maximum force criterion for structure optimization by ase with DP potential. The optimization will stop when the forces on all atoms are below this threshold.

### Example Command

To run the script, use the following command:

```sh
python script_name.py [--num_iterations NUM] [--steps_per_iteration NUM] [--fixed_atoms NUM] [--iffinal BOOL] [--skip_first_opt BOOL]
```

```bash
for in-domain system:
nohup python script.py --num_iterations 3 --steps_per_iteration 200 --fixed_atoms 0 --iffinal true --skip_first_opt false --fmax 0.1 &
for out-of-domain system:
nohup python script.py --num_iterations 5 --steps_per_iteration 200 --fixed_atoms 0 --iffinal true --skip_first_opt true --fmax 0.2 &
```

## Workflow Description

```mermaid
graph 
A(start) --> X{if skip_first}
X -.no.-> B[ase optimization with DP potential]
X -.yes.-> C[VASP optimization by NSW=3 step]
B --get CONTCAR--- C[VASP optimization by NSW=3 step]
C  -->D[collect data by dpdata]
D --> E[finetune with new data by DeePMD-kit]
E --> F{if num_iterations}
F -.yes.-> G{if final}
F -.no.-> B
G-.yes.-> H(final VASP optimization with DFT)
G-.no.-> I(end)
H --> I
```


The script follows these main steps:

1. **Optimization of Atomic Structures by ASE** :
   * Optimizes the atomic structure using catalyticLAM DeepMD-kit models.
   * Saves the optimized structure as `CONTCAR-ase-{iteration}`.
2. **Preparation for VASP Calculation** :
   * Creates a new directory for each iteration (`opt{iteration}`).
   * Copies the optimized structure and necessary VASP input files to the new directory.
   * If not the first iteration, copies `WAVECAR` and `CHGCAR` from the previous iteration.
3. **Run VASP Calculation** :
   * Submits the VASP job and retrieves the job ID.
   * Waits for the VASP job to complete.
4. **Generate New Dataset** :
   * Generates a new dataset from the VASP output file `OUTCAR` by dpdata.
   * Saves the dataset in the specified directory.
5. **Fine-Tuning of Models** :
   * Updates the fine-tuning configuration (`finetune1.json`) with new data paths and iteration-specific parameters.
   * Submits the fine-tuning job and retrieves the job ID.
   * Waits for the fine-tuning job to complete.
6. **Freezing the Model** :
   * Freezes the fine-tuned model using DeepMD-kit.
7. **Final Optimization Step** :
   * After completing all iterations, runs a final optimization step.
   * Creates a new directory (`optfinal`) and copies the final optimized structure.
   * Modifies the `INCAR` file to set `NSW = 300`.
   * Submits the final VASP job and waits for its completion.

## Resuming from Checkpoints

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

## Notes

* Ensure the VASP and DeepMD-kit executables are accessible in your environment.
* Customize the `INCAR`, `KPOINTS`, `POTCAR`, and submission scripts (`sub.vasp`, `sub.dp`) as needed for your specific system and requirements.
* Adjust the polling interval in `wait_for_job_completion` if necessary.

## License

This project is licensed under the LGPL-3.0 License.

## Acknowledgements

* ASE: [https://wiki.fysik.dtu.dk/ase/]()
* VASP: [https://www.vasp.at/](https://www.vasp.at/)
* DeepMD-kit: [https://github.com/deepmodeling/deepmd-kit](https://github.com/deepmodeling/deepmd-kit)
* dpdata: [https://github.com/deepmodeling/dpdata](https://github.com/deepmodeling/dpdata)
* VASPkit: [https://vaspkit.com/](https://vaspkit.com/)
* SLURM: [https://slurm.schedmd.com/]()

