from ase.io import write
from ase.cluster import wulff_construction as wc

metals = {
    'Li': [0.383, 0.289, 0.750],
    'Na': [0.290, 0.197, 0.546],
    'K':  [0.249, 0.167, 0.462],
    'Rb': [0.229, 0.150, 0.417],
    'Cs': [0.228, 0.142, 0.390],
    'Ba': [0.616, 0.464, 1.199],
    'V':  [1.725, 1.312, 3.494],
    'Cr': [2.020, 1.258, 3.626],
    'Fe': [1.265, 0.978, 2.694],
    'Nb': [1.987, 1.320, 3.668],
    'Mo': [2.410, 1.534, 4.068],
    'Ta': [2.174, 1.531, 4.201],
    'W':  [2.955, 1.806, 4.916]
}

surfaces = [(1, 0, 0), (1, 1, 0), (1, 1, 1)]
sizes = [8, 20, 50, 80, 150, 180]

for metal, esurf in metals.items():
    for size in sizes:
        atoms = wc(metal, surfaces, esurf, 
                   size, 'bcc', rounding='closest')
        num_atoms = len(atoms)
        filename = f'{metal}{num_atoms}.xyz'
        write(filename, atoms)
