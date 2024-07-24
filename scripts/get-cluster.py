import argparse
from ase.io import write
from ase.cluster import FaceCenteredCubic as fcc
from ase.cluster import wulff_construction as wc

s_metals = ['Ni', 'Cu', 'Rh', 'Pd', 'Au', 'Ag', 'Ir', 'Pt', 'Al', 'Ca', 'Sr', 'Pb', 'Ce']
s_layers_list = [[2, 1, -1], [3, 2, -1], [4, 3, -1], [5, 4, -1]]
s_surfaces = [(1, 0, 0), (1, 1, 1), (1, -1, 1)]

fcc_metals = {
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

bcc_metals = {
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

hcp_metals = {
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

fcc_surfaces = [(1, 0, 0), (1, 1, 0), (1, 1, 1)]
fcc_sizes = [8, 20, 50, 80, 200]
bcc_surfaces = [(1, 0, 0), (1, 1, 0), (1, 1, 1)]
bcc_sizes = [8, 20, 50, 80, 150, 180]
hcp_surfaces = [(0, 0, 0, 1), (1, 0, -1, 0), (1, 1, -2, 0)]
hcp_sizes = [5, 10, 20, 35, 50, 65, 80, 100, 120, 150, 175, 200]

def generate_s_fcc_structures(metals):
    for layers in s_layers_list:
        for metal in metals:
            atoms = fcc(metal, s_surfaces, layers)
            num_atoms = len(atoms)
            filename = f's_{metal}{num_atoms}.xyz'
            write(filename, atoms)

def generate_structures(metals, surfaces, sizes, lattice):
    for metal, esurf in metals.items():
        for size in sizes:
            atoms = wc(metal, surfaces, esurf, size, lattice, rounding='closest')
            num_atoms = len(atoms)
            filename = f'{metal}{num_atoms}.xyz'
            write(filename, atoms)

def main():
    parser = argparse.ArgumentParser(description='Generate cluster structures.')
    parser.add_argument('--type', choices=['s_fcc', 'fcc', 'bcc', 'hcp', 'all'], default='all', help='Type of structure to generate')
    parser.add_argument('--element', choices=['all'] + list(set(fcc_metals) | set(bcc_metals) | set(hcp_metals)), default='all', help='Element or metal to use')
    args = parser.parse_args()

    if args.element == 'all':
        elements = set(s_metals) | set(fcc_metals) | set(bcc_metals) | set(hcp_metals)
    else:
        elements = [args.element]

    if args.type == 's_fcc':
        generate_s_fcc_structures(elements)
    
    elif args.type == 'fcc':
        if args.element == 'all':
            generate_structures(fcc_metals, fcc_surfaces, fcc_sizes, 'fcc')
        elif args.element in fcc_metals:
            generate_structures({args.element: fcc_metals[args.element]}, fcc_surfaces, fcc_sizes, 'fcc')
    
    elif args.type == 'bcc':
        if args.element == 'all':
            generate_structures(bcc_metals, bcc_surfaces, bcc_sizes, 'bcc')
        elif args.element in bcc_metals:
            generate_structures({args.element: bcc_metals[args.element]}, bcc_surfaces, bcc_sizes, 'bcc')
    
    elif args.type == 'hcp':
        if args.element == 'all':
            generate_structures(hcp_metals, hcp_surfaces, hcp_sizes, 'hcp')
        elif args.element in hcp_metals:
            generate_structures({args.element: hcp_metals[args.element]}, hcp_surfaces, hcp_sizes, 'hcp')
    
    elif args.type == 'all':
        if args.element == 'all':
            generate_s_fcc_structures(s_metals)
            generate_structures(fcc_metals, fcc_surfaces, fcc_sizes, 'fcc')
            generate_structures(bcc_metals, bcc_surfaces, bcc_sizes, 'bcc')
            generate_structures(hcp_metals, hcp_surfaces, hcp_sizes, 'hcp')
        else:
            if args.element in s_metals:
                generate_s_fcc_structures([args.element])
            if args.element in fcc_metals:
                generate_structures({args.element: fcc_metals[args.element]}, fcc_surfaces, fcc_sizes, 'fcc')
            if args.element in bcc_metals:
                generate_structures({args.element: bcc_metals[args.element]}, bcc_surfaces, bcc_sizes, 'bcc')
            if args.element in hcp_metals:
                generate_structures({args.element: hcp_metals[args.element]}, hcp_surfaces, hcp_sizes, 'hcp')

if __name__ == '__main__':
    main()
