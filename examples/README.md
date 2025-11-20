# Examples

This directory contains VASP calculation workflow examples for different types of materials, demonstrating how to use this project for automated catalytic material computations.

## Directory Structure

Each subdirectory contains the following components:
- [flow.py](examples/cluster/flow.py): Main workflow execution script
- [input](examples/2D/input): Input parameter configuration file
- [record.txt](vaspworkflow/record.txt): Running record file (in some directories)
- `utils/`: Utility functions directory

## Example Types Description

### 1. 2D Materials (2D)

This example demonstrates how to perform VASP calculations on two-dimensional materials. It includes structural optimization and molecular dynamics simulations.

### 2. Bulk Materials (bulk)

This example shows how to perform VASP calculations on bulk materials (such as metals, alloys, etc.). It includes the complete process from retrieving structures from the structure database, setting calculation parameters, submitting jobs, to collecting results.

### 3. Clusters (cluster)

This example focuses on the calculation process for atomic cluster materials. It is suitable for structural optimization and property calculations of small-sized atomic cluster systems.

### 4. Molecules (molecule)

This example is specifically designed for molecular system calculations. Compared to other types, this example may include special processing steps, such as using dpdata to handle molecular data.

### 5. Surfaces (slab)

This example demonstrates the calculation process for surface systems (such as catalyst surfaces). It typically involves surface model construction, adsorbate placement, and surface reaction studies.

## Usage

The [flow.py](examples/cluster/flow.py) script in each example directory is the main execution program. Calculation parameters can be adjusted by modifying the [input](examples/2D/input) configuration file. Run it as follows:

```bash
python flow.py
```

Before running, please ensure:

VASP-related environment is properly configured
Modify parameters in the input configuration file as needed
Sufficient computational resources are available
Notes
Each example is independent and can be run separately
Choose the appropriate example as a template based on your needs
Examples can be appropriately modified according to the characteristics of your specific research object