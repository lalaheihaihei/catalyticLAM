# catalyticLAM
Machine-Learning-Based Interatomic Potentials for Catalysis: an Universal Catalytic Large Atomic Model

## Directory Structure

```
├── pretrainedCLAM 
├── datasetworkflow
│   ├── flow.py
│   ├── input
│   ├── structure_db
│   └── utils
├── generation
│   ├── alloy_db.py
│   ├── make-bulk.py
│   ├── make-slab.py
│   ├── material_db.py
│   ├── metal_db.py
│   └── molecule_db.py
├── LICENSE
├── postworkflow
│   ├── flowopt.py
│   ├── flowts.py
│   ├── POSCAR
│   └── utils
└── README.md
```

## Overview

### 1. Generation

This section is responsible for generating structures, including bulk and slab structures for VASP calculation. Please read README.md in generation for details.

- `make-bulk.py`: Generates bulk structures.
- `make-slab.py`: Generates slab structures.
- `material_db.py`: Information for generation database.
- `molecule_db.py`: Molecular structures database.

### 2. DataSetWorkflow

This section manages VASP tasks and workflows, as well as collects data for dpdata.

- `flow.py`: Main workflow script for managing VASP calculations and data processing.
- `input`: Input parameter file.
- `POSCAR`: Directory containing various POSCAR files and their corresponding VASP calculation results.
- `structure_db`: Stores the structure database.
- `utils`: Contains configuration files and scripts required for VASP calculations.

### 3. PostWorkflow

This section is responsible for model-accelerated structure optimization, transition state search, and catalytic reaction network construction.

- `flowopt.py`: Workflow script for structure optimization.
- `flowts.py`: Workflow script for transition state search.
- `POSCAR`: Initial structure file.
- `utils`: Contains configuration files and scripts for optimization and transition state search.

### 4. Pretrained Catalytic Large Atomic Model for PostWorkflow

## Installation and Usage

### Prerequisites

Ensure that your system has the following software installed:

- Python 3 (version > 3.10)
- ASE (version > 3.22)
- DeepMD-kit (version > 3.0.0a1)
- dpdata (version > 0.2.18)
- VASP and VASPkit
- SLURM (for job scheduling)
- tqdm (for progress bars)

### Installation

1. Clone the repository:

```bash
git clone https://github.com/lalaheihaihei/catalyticLAM.git
cd catalyticLAM
```

### Usage

#### Generate Structures

Navigate to the `generation` directory and run the appropriate script to generate the desired structures:

```
cd generation
python ./make-bulk.py --task search --ificsd --bulktype alloy --elementNumber 2
python ./make-slab.py --up-down UUUUDDDD --plot --element Ru --molecule-type all
```

#### Run VASP Workflow

Navigate to the `datasetworkflow` directory, edit the `input` file as needed, and run `flow.py`:

```
cd datasetworkflow
nohup python flow.py POSCAR opt &
python flow.py POSCAR optcheck
nohup python flow.py POSCAR md &
python flow.py POSCAR mdcheck
python flow.py POSCAR dpdata
python flow.py POSCAR plot
```

#### Run Structure Optimization and Transition State Search

Navigate to the `postworkflow` directory, edit the configuration files, and run the relevant scripts:

```
cd postworkflow
python flowopt.py
```

## License

This project is licensed under the GNU General Public License v3.0. See the [LICENSE](LICENSE) file for details.

## Acknowledgements

* ASE: [https://wiki.fysik.dtu.dk/ase/](https://markdown.lovejade.cn/)
* VASP: [https://www.vasp.at/](https://www.vasp.at/)
* DeepMD-kit: [https://github.com/deepmodeling/deepmd-kit](https://github.com/deepmodeling/deepmd-kit)
* dpdata: [https://github.com/deepmodeling/dpdata](https://github.com/deepmodeling/dpdata)
* VASPkit: [https://vaspkit.com/](https://vaspkit.com/)
* fairchem: [https://github.com/FAIR-Chem/fairchem](https://github.com/FAIR-Chem/fairchem)
* SLURM: [https://slurm.schedmd.com/](https://markdown.lovejade.cn/)

