from ase.io import write
from ase.cluster import wulff_construction as wc

metals = {
    'Cu': [0.906, 1.323, 0.707],
    'Ag': [0.653, 0.953, 0.553],
    'Au': [0.895, 1.321, 0.611],
    'Ni': [0.969, 1.337, 0.695],
    'Rh': [1.310, 1.919, 1.002],
    'Pd': [1.152, 1.559, 0.824],
    'Ir': [1.772, 2.428, 1.225],
    'Pt': [1.378, 2.009, 1.004],
    'Al': [0.689, 0.919, 0.531],
    'Ca': [0.535, 0.811, 0.484],
    'Sr': [0.484, 0.725, 0.440],
    'Pb': [0.307, 0.513, 0.226],
    'Ce': [1.133, 1.129, 1.018]
}

surfaces = [(1, 0, 0), (1, 1, 0), (1, 1, 1)]
sizes = [8, 20, 50, 80, 200]

for metal, esurf in metals.items():
    for size in sizes:
        atoms = wc(metal, surfaces, esurf, 
                   size, 'fcc', rounding='closest')
        num_atoms = len(atoms)
        filename = f'{metal}{num_atoms}.xyz'
        write(filename, atoms)
