# CatalyticLAM

Machine-Learning-Based Interatomic Potentials for Catalysis: an Universal Catalytic Large Atomic Model

## Overview

### 1. Generation

This section is responsible for generating structures, including bulk and slab structures for VASP calculation. Please read [README.md](./generation/README.md) in generation for details.

- `get-bulk.py`: Generates bulk structures.
- `get-slab.py`: Generates slab structures.
- `element_list.json`: Metal and alloy elements for bulk generation.
- `material.json`: Information for database generation.
- `molecule.json`: Molecular structures database.

### 2. VASPWorkflow

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

This section contains the pretrained CLAM model for postworkflow and its training, finetune, checkpoint files.

- `finetune.json`: input file for fine-tuning the CLAM model.
- `input.json`: input file for training the CLAM model.
- `model.ckpt-10000000.pt`: The checkpoint of pretrained CLAM model.
- `lcurve`: output file containing learning curves.
- `dataset`: Directory containing the dataset for finetune operation.

### 5. Scripts

This section contains some useful scripts to generate cluster structures and convert the format of files.

- `cif2pos.py`: Convert CIF file to POSCAR.
- `get-cluster.py`: Generate the structures of metal clusters in xyz format.
- `json2cif.py`: Convert JSON file to CIF file.
- `xyz2pos.py`: Convert XYZ file to POSCAR.
- `sim_model.py`: A script for delete the unnecessary keys in checkpoint files(oc22).

### 6. Structure_db

This section contains the initial structure files and POSCAR files for VASP calculation.

- `2D.tgz`: The total 6351 POSCAR files of 2D materials for VASP calculation.
- `2D-raw.tgz`: The initial json file containing the information of 2D materials and the corresponding cif files.
- `bulk.tgz`: The POSCAR files of metals and alloys for VASP calculation.
- `cluster.tgz`: The POSCAR files of clusters for VASP calculation.
- `cluster-raw.tgz`: The initial xyz files of clusters.
- `molecule.tgz`: The total POSCAR files of molecules for VASP calculation.
- `molecule-raw.tgz`: The initial xyz files of molecules.
- `slab.tgz`: The POSCAR files of slabs for VASP calculation.

## Installation and Usage

### Prerequisites

Ensure that your system has the following software installed:

- Python 3 (version > 3.10)
- ASE (version > 3.22)
- DeepMD-kit (version > 3.0.0a1)
- dpdata (version > 0.2.18)
- VASP (version > 5.4.4)
- VASPKIT
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
python get-bulk.py --task search --ificsd --bulktype alloy --elementNumber 2
python get-slab.py --up-down UUD --plot --element Au --molecule-type CO --type type1
python get-slab.py --up-down UUUUDDDD --plot --element Ru --molecule-type all --type type3
```

#### Run VASP Workflow

Navigate to the `datasetworkflow` directory, edit the `input` file as needed, and run `flow.py`:

```
cd vaspworkflow
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
cd optdp or optoc
python flowopt.py
```

#### Construct Pretrained Catalytic Large Atomic Model for PostWorkflow

Navigate to the `script` directory, edit the input files, and run the relevant scripts:
more information please refer to Deepmd-kit official website(below).

```
dp --pt train input.json > out
dp --pt train --finetune model.ckpt.10000000.pt --model-branch <head> finetune.json > out (At present, the head name is only supported for oc22 and metal)
```

#### Run script in Scripts folder

Navigate to the `train` directory and run the appropriate script to generate the cluster structures or convert file formats.

1. cif2pos.py:

```
python cif2pos.py --input_dir CIF --output_dir POSCAR   # Convert the cif files in the CIF folder to corresponding POSCAR files in the POSCAR folder
python cif2pos.py --input_dir CIF
# --input_dir represents the path to the input directory containing CIF files
# --output_dir represents the path to the output directory for POSCAR files, default is 'POSCAR'
```

2. get-cluster.py:

```
python get-cluster.py   # generate all cluster structures in xyz format
python get-cluster.py --type all --element all   # generate all cluster structures in xyz format
python get-cluster.py --type all --element Cu   # generate all Cu-related clusters
python get-cluster.py --type s_fcc --element all   # generate all supported metal clusters
python get-cluster.py --type fcc --element Au   # generate Cu-related clusters in FCC lattice
python get-cluster.py --type hcp --element Zn   # generate Zn-related clusters in HCP lattice
# --type represents which type of structure to generate, including 's_fcc', 'fcc', 'bcc', 'hcp', and 'all', default is 'all', where 's_fcc' reprents the supported cluster
# --type represents which element or metal to use, default is 'all'
```

3. json2cif.py:

```
python json2cif.py --file input.json --outdir CIF   # Convert the json file to corresponding cif files in the CIF folder
python json2cif.py --file input.json
# --file represents the path to the JSON file
# --outdir represents the path to the output directory for CIF files, default is 'CIF'
```

4. xyz2pos.py:

```
python xyz2pos.py --input_dir XYZ --output_dir POSCAR --padding 5.0
python xyz2pos.py --input_dir XYZ
# --input_dir represents the path to the input directory containing XYZ files
# --output_dir represents the path to the output directory for POSCAR files, default is 'POSCAR'
# --padding reprensts the distance of the atoms from the box boundary in three directions, default is '5.0'
```

5. sim_model.py:
```
python sim_model.py
# You should revise the path of checkpoint files and specify the keys you want to delete by yourself.
```

## License

This project is licensed under the GNU General Public License v3.0. See the [LICENSE](LICENSE) file for details.

## Acknowledgements

* ASE: [https://wiki.fysik.dtu.dk/ase/](https://markdown.lovejade.cn/)
* VASP: [https://www.vasp.at/](https://www.vasp.at/)
* DeepMD-kit: [https://github.com/deepmodeling/deepmd-kit](https://github.com/deepmodeling/deepmd-kit)
* dpdata: [https://github.com/deepmodeling/dpdata](https://github.com/deepmodeling/dpdata)
* VASPKIT: [https://vaspkit.com/](https://vaspkit.com/)
* fairchem: [https://github.com/FAIR-Chem/fairchem](https://github.com/FAIR-Chem/fairchem)
* SLURM: [https://slurm.schedmd.com/](https://markdown.lovejade.cn/)
