# CatalyticLAM

Machine-Learning-Based Interatomic Potentials for Catalysis: a Universal Catalytic Large Atomic Model

* [Paper](https://chemrxiv.org/engage/chemrxiv/article-details/66c56756a4e53c48760ee61c)

0. [Pre-trained Model](#0-pre-trained-model)
1. [Overview](#1-overview)
   - [1.1 Generation](#11-generation)
   - [1.2 VASP Workflow](#12-vasp-workflow)
   - [1.3 Post Workflow](#13-post-workflow)
   - [1.4 Pretrained CLAM for Post Workflow](#14-pretrained-clam-for-post-workflow)
   - [1.5 Scripts](#15-scripts)
   - [1.6 Initial Structures](#16-initial-structures)
2. [Installation](#2-installation)
   - [2.1 Prerequisites](#21-prerequisites)
   - [2.2 Installation](#22-installation)
3. [Quick Usage](#3-quick-usage)
   - [3.1 Generate Structures](#31-generate-structures)
   - [3.2 Run VASP Workflow](#32-run-vasp-workflow)
   - [3.3 Run Structure Optimization and Transition State Search](#33-run-structure-optimization-and-transition-state-search)
   - [3.4 Run Reaction Network generation](#34-run-reaction-network-generation)
   - [3.5 Construct Pretrained CLAM for Post Workflow](#35-construct-pretrained-clam-for-post-workflow)
   - [3.6 Other Scripts](#36-other-scripts)
   - [3.7 Initial Structures](#37-initial-structures)
4. [License](#4-license)
5. [Acknowledgements](#5-acknowledgements)
6. [Citation](#6-citation)

## 0. Pre-trained Model
Our Pre-trained Models can be obtained in the following configurations saved in google driver:

| Model   | training strategy     | Download            | val force MAE(meV/A) on metal system |val force MAE(meV/A) on qm9 | val force MAE(meV/A) on 2d system| val energy MAE(meV/atom) on metal system | val energy MAE(meV/atom) on qm9 |val energy MAE(meV/atom) on 2d system |
|----------|------------|---------------------|-------------------|---------|---------|---------|---------|---------|
| Gemnet-OC (1.0.0)  | finetuned based on [GemNet-OC-S2EFS-OC20+OC22](https://fair-chem.github.io/catalysts/models.html) for 5 epoch without quad_interaction | [checkpoint](https://drive.google.com/file/d/17U-dW6PvADtQFePky64fsN0PgFmdYxrp/view?usp=drive_link)<br>[config](https://drive.google.com/file/d/1R-Wi7YJ8ddZsp_lPZ2QMriWrjXZ3sO1V/view?usp=drive_link)    | 34.5 | 15.7| 42.0  |  3.9| 2.4| 7.0   |
| Gemnet-OC (1.0.0)  | finetuned based on [GemNet-OC-S2EFS-OC20+OC22](https://fair-chem.github.io/catalysts/models.html) for 50 epoch with quad_interaction | [checkpoint](https://drive.google.com/file/d/1-eDLr-ZMnMXLxdfPmFugzdZrFXQc5dwt/view?usp=drive_link)<br>[config](https://drive.google.com/file/d/16aNrEzUTzPbQ9HCykNMUiDLEnJeBzLAD/view?usp=drive_link)    | 28.5 | 21.5| 27.1  |  1.8 | 2.9| 3.2   |
| DPA2 (v3.0.0b4) | finetuned based on [DPA2_medium_28_10M_beta4.pt](https://www.aissquare.com/models/detail?pageType=models&name=DPA-2.3.0-v3.0.0b4) for 10000000 steps| [model.ckpt-10000000.pt](https://drive.google.com/file/d/1dyPteWsJr4wXCFn1dB3cq0XoUhSyZiGA/view?usp=drive_link)<br>[config](https://drive.google.com/file/d/1i04b0IFnjW7PX5GCWhVnaqR58v9r_0qs/view?usp=drive_link) | 72.2  | 78.0  | 84.5 | 5.8 | 7.7 | 5.5 | 

## 1. Overview

### 1.1 Generation

This section is responsible for generating structures, including bulk and slab structures for VASP calculations, and generating initial adsorption structures. We preset a series of commonly used small adsorbates and place them on slab models with 1, 2, or 4 molecules.

### 1.2 VASP Workflow

This section manages VASP tasks and workflows, as well as collects data for `dpdata`. It can perform high-throughput optimization and molecular dynamics (MD) jobs based on pre-generated structures. It can also check the convergence of SCF steps in optimization and MD, and perform high-throughput conversion to LMDB or NPY format for further training.

### 1.3 Post Workflow

This section is responsible for model-accelerated structure optimization, transition state search, and catalytic reaction network construction. The optimization and transition state search are based on a local fine-tune method, which involves a Labeling, Fine-tuning, and Inference loop to accelerate the optimization and MD process. It can also automatically construct reaction networks to generate possible intermediates and transition state structures.

### 1.4 Pretrained CLAM for Post Workflow

This section contains the pretrained CLAM model for the post workflow, including its training, fine-tuning, and checkpoint files.

### 1.5 Scripts

This section contains some useful scripts to generate cluster structures and convert the format of files.

### 1.6 Initial Structures

This section contains the initial structure files and POSCAR files for VASP optimization and MD calculations to generate datasets.

## 2. Installation

### 2.1 Prerequisites

Ensure that your system has the following software installed:

- Python 3 (version > 3.10)
- ASE (version > 3.22)
- Pymatgen (version > 2023.3.23)
- DeePMD-kit (version > 3.0.0a1)
- dpdata (version > 0.2.18)
- fairchem
- VASP (version > 5.4.4)
- VASPKIT
- SLURM (for job scheduling)
- tqdm (for progress bars)

### 2.2 Installation

1. Clone the repository:

```bash
git clone https://github.com/lalaheihaihei/catalyticLAM.git
```

2. Enter the folder:

```
cd catalyticLAM
```

## 3. Quick Usage

### 3.1 Generate Structures

Navigate to the [generation](./generation/) directory and run the appropriate script to generate the desired structures:

- `get-bulk.py`: Generates bulk structures.
- `get-slab.py`: Generates slab structures.
- `element_list.json`: Metal and alloy elements for bulk generation.
- `material.json`: Information for database generation.
- `molecule.json`: Molecular structures database.

```
cd generation
python get-bulk.py --api-key Your-Api-Key --bulktype metal --elementNumber 1 --task search --ificsd
python get-bulk.py --plot --api-key Your-Api-Key --min-lw 10.0 --task generate
python get-slab.py --plot --api-key Your-Api-Key --molecule-type CO --up-down UUD --element Au --type type1
python get-slab.py --plot --api-key Your-Api-Key --molecule-type all --up-down UUUUDDDD --element Pd --type type3
```

Detailed usages are seen in [README.md](./generation/README.md).

### 3.2 Run VASP Workflow

Navigate to the [vaspworkflow](./vaspworkflow/) directory, edit the `input` file as needed, and run `flow.py`:

- `flow.py`: Main workflow script for managing VASP calculations and data processing.
- `input`: Input parameter file.
- `POSCAR`: Directory containing various POSCAR files and their corresponding VASP calculation results.
- `structure_db`: Stores the structure database.
- `utils`: Contains configuration files and scripts required for VASP calculations.

```
cd vaspworkflow
nohup python flow.py POSCAR opt &
python flow.py POSCAR optcheck
nohup python flow.py POSCAR md &
python flow.py POSCAR mdcheck
python flow.py POSCAR dpdata
python flow.py POSCAR plot
```

Detailed usages are seen in [README.md](./vaspworkflow/README.md).

### 3.3 Run Structure Optimization and Transition State Search

Navigate to the [postworkflow](./postworkflow) directory, prepare input files, and run the relevant scripts:

- `flowopt.py`: Workflow script for structure optimization.
- `flowts.py`: Workflow script for transition state search.
- `POSCAR`: Initial structure file.
- `utils`: Contains configuration files and scripts for optimization and transition state search.

```
cd postworkflow
cd optdp or optoc
nohup python ./flowopt.py --num_iterations 3 --steps_per_iteration 200 --fixed_atoms 0 --iffinal true --fmax 0.1 &
cd tsdp or tsoc
nohup python ./flowts.py POSCARis POSCARfs ./frozen_model.pth OUTCARis OUTCARfs &
```

Detailed usages are seen in [README.md](./postworkflow/README.md).

### 3.4 Run Reaction Network generation

Navigate to the [postworkflow/RNET](./postworkflow/RNET) directory, prepare input files, and run the relevant scripts:

- `RNet.py`:  Genarate  reaction network diagram.
- `MakeSlab.py`: Construct all possible structures for intermediats adsorption on metal surfaces.
- `plot_all.py`: Plot the energy changes and energy differences MAE.

```
cd postworkflow/RNET
python RNet.py 1 2 --layout spring
python MakeSlab.py --element Pt --max-index 1
```

Detailed usages are seen in [README.md](./postworkflow/RNET/README.md).

### 3.5 Construct Pretrained CLAM for Post Workflow

Navigate to the [train](./train) directory, edit the input files, and run the training or fine-tuning jobs.
Details of CLAM are in [README.md](./train/README.md)

```
dp --pt train input.json > out
dp --pt train --finetune model.ckpt.10000000.pt --model-branch <head> finetune.json > out (At present, the head name is only supported for oc22, qm and metal)
```

```
python main.py --mode train --config-yml finetune1.yml --print-every 1000 >> out
python main.py --mode train --config-yml finetune1.yml --checkpoint gnoc_oc22_oc20_all_s2ef.pt --print-every 1000 >> out
```

Detailed usages are seen in [README.md](./train/README.md).

More information please refer to [Deepmd-kit official website](https://github.com/deepmodeling/deepmd-kit) and [fairchem official website](https://github.com/FAIR-Chem/fairchem).

### 3.6 Other Scripts

Navigate to the [scripts](./scripts) directory and run the appropriate script to generate the cluster structures or convert file formats.

* `cif2pos.py`: Convert CIF file to POSCAR.
* `get-cluster.py`: Generate the structures of metal clusters in xyz format.
* `json2cif.py`: Convert JSON file to CIF file.
* `xyz2pos.py`: Convert XYZ file to POSCAR.
* `sim_model.py`: For deleting the unnecessary keys in checkpoint files (oc22).
* `cal_nframes.py`: Calculate the number of frames in a dataset with dp (deepmd-kit) format.
* `make_test.py`: Make a dataset test with lmdb format.

More details are seen in [README.md](./scripts/README.md).

### 3.7 Initial structures

Navigate to the [structure_db](./structure_db) directory, you can find compressed files, which containing the initial structures.

- `2D.tgz`: The total 6351 POSCAR files of 2D materials for VASP calculation.
- `2D-raw.tgz`: The initial json file containing the information of 2D materials and the corresponding cif files.
- `bulk.tgz`: The POSCAR files of metals and alloys for VASP calculation.
- `cluster.tgz`: The POSCAR files of clusters for VASP calculation.
- `cluster-raw.tgz`: The initial xyz files of clusters.
- `molecule.tgz`: The total POSCAR files of molecules for VASP calculation.
- `molecule-raw.tgz`: The initial xyz files of molecules.
- `slab.tgz`: The POSCAR files of slabs for VASP calculation.

## 4. License

This project is licensed under the GNU General Public License v3.0. See the [LICENSE](LICENSE) file for details.

## 5. Acknowledgements

* ASE: [https://wiki.fysik.dtu.dk/ase/](https://wiki.fysik.dtu.dk/ase/)
* VASP: [https://www.vasp.at/](https://www.vasp.at/)
* DeePMD-kit: [https://github.com/deepmodeling/deepmd-kit](https://github.com/deepmodeling/deepmd-kit)
* dpdata: [https://github.com/deepmodeling/dpdata](https://github.com/deepmodeling/dpdata)
* fairchem: [https://fair-chem.github.io/index.html](https://fair-chem.github.io/index.html)
* VASPKIT: [https://vaspkit.com/](https://vaspkit.com/)
* Pymatgen: [https://pymatgen.org/](https://pymatgen.org/)
* SLURM: [https://slurm.schedmd.com/](https://slurm.schedmd.com/)

## 6. Citation

Please cite the works below if this repository is helpful.

* [CLAM](https://chemrxiv.org/engage/chemrxiv/article-details/66c56756a4e53c48760ee61c)
```
Wu Z, Zhou L, Hou P, Liu Y, Guo T, Liu J-C. Catalytic Large Atomic Model (CLAM): A Machine-Learning-Based Interatomic Potential Universal Model. ChemRxiv. 2024; doi:10.26434/chemrxiv-2024-2xzct 
```
