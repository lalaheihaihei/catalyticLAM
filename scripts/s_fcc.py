from ase.io import write
from ase.cluster import FaceCenteredCubic as fcc

surfaces = [(1, 0, 0), (1, 1, 1), (1, -1, 1)]
layers_list = [
    [2, 1, -1],
    [3, 2, -1],
    [4, 3, -1],
    [5, 4, -1]
]

metals = ['Ni', 'Cu', 'Rh',
          'Pd', 'Au', 'Ag',
          'Ir', 'Pt', 'Al',
          'Ca', 'Sr', 'Pb', 'Ce'
]

for layers in layers_list:
    for metal in metals:
        atoms = fcc(metal, surfaces, layers)
        num_atoms = len(atoms)
        filename = f's_{metal}{num_atoms}.xyz'
        write(filename, atoms)
