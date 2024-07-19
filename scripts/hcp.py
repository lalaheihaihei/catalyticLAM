from ase.io import write
from ase.cluster import wulff_construction as wc

metals = {
    'Be': [1.796, 1.904, 2.473],
    'Mg': [0.510, 0.597, 0.716],
    'Sc': [1.272, 1.200, 1.261],
    'Ti': [1.975, 2.031, 1.887],
    'Co': [2.108, 2.257, 2.462],
    'Zn': [0.334, 0.528, 0.933],
    'Y':  [0.996, 0.961, 1.019],
    'Zr': [1.599, 1.660, 1.717],
    'Ru': [2.598, 2.906, 3.257],
    'Cd': [0.159, 0.379, 0.483],
    'Hf': [1.721, 1.871, 1.836],
    'Re': [2.566, 3.287, 3.075],
    'Os': [2.950, 3.454, 4.110],
    'Tl': [0.266, 0.332, 0.287],
    'La': [0.695, 0.738, 0.828]
}

surfaces = [(0, 0, 0, 1), (1, 0, -1, 0), (1, 1, -2, 0)]
sizes = [5, 10, 20, 35, 50, 65, 80, 100, 120, 150, 175, 200]

for metal, esurf in metals.items():
    for size in sizes:
        atoms = wc(metal, surfaces, esurf, 
                   size, 'hcp', rounding='closest')
        num_atoms = len(atoms)
        filename = f'{metal}{num_atoms}.xyz'
        write(filename, atoms)
