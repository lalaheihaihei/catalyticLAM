from fairchem.core import OCPCalculator
from ase.io import read, write
from ase.optimize import BFGS
import sys

# Path to the checkpoint file
checkpoint_path = './finetune4/checkpoint.pt'

# Initialize the OCPCalculator with the checkpoint file
calculator = OCPCalculator(checkpoint_path=checkpoint_path)

# Path to the POSCAR file
poscar_path = './opt4/CONTCAR'

# Read the POSCAR file
slab = read(poscar_path, format='vasp')

# Set tags: Ru atoms to 1, others to 0
tags = []
for atom in slab:
    if atom.symbol == 'Ru':
        tags.append(1)
    else:
        tags.append(0)
slab.set_tags(tags)
print(slab.get_tags())
# Set the calculator
slab.set_calculator(calculator)

# Redirect output to file
with open('aseopt.out', 'w') as f:
    sys.stdout = f
    sys.stderr = f

    # Initialize the BFGS optimizer
    optimizer = BFGS(slab)

    # Run the optimization
    optimizer.run(fmax=1.5e-1)

    # Reset output redirection
    sys.stdout = sys.__stdout__
    sys.stderr = sys.__stderr__

# Save the final structure as CONTCAR
contcar_path = 'CONTCAR'
write(contcar_path, slab)

# Print the final optimized positions
print("Optimized positions:")
print(slab.get_positions())
