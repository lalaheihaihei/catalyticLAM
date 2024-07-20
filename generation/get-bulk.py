import argparse
import itertools
import json
from mp_api.client import MPRester
from pymatgen.io.vasp.inputs import Poscar
from plotter import BulkPlotter
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from mp_api.client import MPRester
from element_list import metal_element_list,alloy_element_list
from material_db import type1,type4  # Replace 'some_library' with the actual library name

class BulkGenerator:
    def __init__(self, api_key):
        """
        Initialize the BulkGenerator with an API key for accessing the Materials Project database.
        
        :param api_key: API key for the Materials Project REST API.
        """
        self.mpr = MPRester(api_key)
        self.plotter = BulkPlotter(rows=7, cols=4)  # Adjust subplot layout and size as needed

    def fetch_structures(self, material_ids, min_lw):
        """
        Fetches and scales structures for a list of given material IDs.
        Ensures that a, b, and c lattice parameters are at least min_lw Å.
        
        :param material_ids: List of Materials Project IDs for the bulk materials.
        :param min_lw: Minimum length/width for the lattice vectors.
        :return: Dictionary of material ID to scaled structure.
        """
        structures = {}
        for material_id in material_ids:
            print(f"Fetching and scaling structure for {material_id['name']} with ID {material_id['mp_id']}...")
            structure = self.mpr.get_structure_by_material_id(material_id['mp_id'])
            magmom = structure.site_properties.get('magmom', 'N/A')  # Default to 'N/A' if magmom is not available
            structure = SpacegroupAnalyzer(structure).get_conventional_standard_structure()
            
            # Calculate scaling factors for each lattice parameter
            scaling_factors = []
            for a in structure.lattice.abc:
                if a < min_lw:
                    scaling_factors.append(int(min_lw / a)+1)
                else:
                    scaling_factors.append(1)
            
            # Only make supercell if at least one scaling factor is greater than 1
            if any(sf > 1 for sf in scaling_factors):
                structure.make_supercell(scaling_factors)

            structures[material_id['mp_id']] = {'structure': structure, 'magmom': magmom}   
        return structures

    def generate_bulk(self, structures):
        """
        Saves the given structures as POSCAR files with a 'mag' suffix in the filename if any magnetic moments are >= 0.1.
    
        :param structures: Dictionary of material ID to dict containing 'structure' and 'magmom'.
        """
        for material_id, data in structures.items():
            structure = data['structure']
            magmom = data['magmom']
    
            # Check if any magnetic moment values are greater than or equal to 0.1
            if isinstance(magmom, list):
                has_significant_magmom = any(float(m) >= 0.1 for m in magmom if isinstance(m, (float, int)))
            else:
                has_significant_magmom = False
    
            # Generate a filename with or without the 'mag' suffix
            magmom_suffix = "_mag" if has_significant_magmom else ""
            filename = f"POSCAR_bulk_{structure.composition.reduced_formula}_{material_id}{magmom_suffix}"
    
            with open(filename, 'w') as file:
                file.write(Poscar(structure).get_str())
            print(f"POSCAR file for {material_id} saved as {filename}")
            

    def plot_bulk(self, structures):
        """
        Plots the given structures using BulkPlotter, batching the plots into groups of 28 structures per image.
    
        :param structures: Dictionary of material ID to dict containing 'structure' and 'magmom'.
        """
        structures_per_plot = 28  # Number of structures per plot
        structure_ids = list(structures.keys())
        total_structures = len(structure_ids)
        batches = (total_structures + structures_per_plot - 1) // structures_per_plot  # Calculate total number of batches
    
        for batch in range(batches):
            start_index = batch * structures_per_plot
            end_index = start_index + structures_per_plot
            current_structure_ids = structure_ids[start_index:end_index]
    
            # Initialize a new plotter for each batch
            self.plotter = BulkPlotter(rows=7, cols=4)  # Adjust subplot layout and size as needed
    
            for material_id in current_structure_ids:
                if material_id in structures:
                    data = structures[material_id]
                    structure = data['structure']  # Extract the structure from the dict
                    self.plotter.plot_bulk(structure, title=f"Bulk: {material_id}")
    
            # Save the figure for the current batch
            filename = f"bulk_structures{batch + 1}.png"
            self.plotter.save_figure(filename=filename)
            print(f"Saved {filename} with structures {start_index + 1} to {min(end_index, total_structures)}.")
    
def parse_command_line_arguments():
    """
    Parses command line arguments and returns the parsed arguments.
    
    Returns:
        argparse.Namespace: An object containing the parsed command line arguments.
    """
    parser = argparse.ArgumentParser(description="Generate and save bulk structure POSCAR files, and optionally plot them.")
    parser.add_argument("--plot", action="store_true", help="Enable plotting of bulk structure based on bulkgenerate, default: False")
    parser.add_argument("--api-key", type=str, default="Cl1zzLwiyfrzcWoBz88veLs6NRrJYHT4", help="Materials Project API key")
    parser.add_argument("--min-lw", type=float, default=10.0, help="Minimum slab model a and b vector, default: 10.0")
    parser.add_argument("--bulktype", type=str, default="metal", help="input metal or alloy or oxide, default: metal")
    parser.add_argument("--elementNumber", type=int, default=2, help="Number of element types of alloy, default: 2")
    parser.add_argument("--task", type=str, default="generate", help="Task is bulk search or generate, default: generation")
    parser.add_argument("--abovehull", type=float, default=0.3, help="set the energy scope (0,energy_above_hull) for structure search, default: 0.3")
    parser.add_argument("--ificsd", action="store_true", help="icsd true means only experimental strcutures, theoretical means include all possible structures include the material does not match any known ICSD entries")
    
    return parser.parse_args()

def get_most_stable_metal_mp_id(api_key, elements, ificsd=True, abovehull=0.3):
    """
    Fetch the most stable mp_id for given elements from the Materials Project.

    :param api_key: Your Materials Project API key as a string.
    :param elements: List of element symbols, e.g., ["Cu", "Fe"] for copper and iron.
    :return: A dictionary with the most stable mp_id for the given elements or None if not found.
    """
    stable_ids = [] 
    with MPRester(api_key) as mpr:
        for element in elements:
            print('search',element)
            try:
                if ificsd == False:
                    results = mpr.materials.summary.search(elements=[element], num_elements=1,
                        energy_above_hull=(0,abovehull),fields=["material_id", "energy_above_hull"])
                else:
                    results = mpr.materials.summary.search(elements=[element], num_elements=1,
                        energy_above_hull=(0,abovehull),fields=["material_id", "energy_above_hull"], theoretical=False)
                if results:
                   for result in results:
                        stable_ids.append({'name': element, 'mp_id': result.material_id})
            except Exception as e:
                print(f"An error occurred for element {element}: {e}")
    return stable_ids

def get_most_stable_alloy_mp_id(api_key, elementNumber, elements, ificsd=True, abovehull=0.3):
    stable_ids = []
    element_pairs = list(itertools.combinations(alloy_element_list, elementNumber))
    with MPRester(api_key) as mpr:
        for element_pair in element_pairs:
            print('search',element_pair)
            try:
                if ificsd == False:
                    results = mpr.materials.summary.search(elements=list(element_pair), num_elements=elementNumber,
                              energy_above_hull=(0,abovehull), fields=["material_id", "energy_above_hull"]) 
                else:
                    results = mpr.materials.summary.search(elements=list(element_pair), num_elements=elementNumber,
                              energy_above_hull=(0,abovehull), fields=["material_id", "energy_above_hull"], theoretical=False)
                if results:
                    for result in results:
                        structure = mpr.get_structure_by_material_id(result.material_id)
                        formula = structure.composition.reduced_formula
                        stable_ids.append({'name': formula, 'mp_id': result.material_id})
            except Exception as e:
                print(f"An error occurred for element pair {element_pair}: {e}")
    return stable_ids

if __name__ == "__main__":
    args = parse_command_line_arguments()
    if args.task == "search":
        if args.bulktype == "metal": 
            elements = metal_element_list
            most_stable_mp_ids = get_most_stable_metal_mp_id(args.api_key, elements, args.ificsd)
            with open("metal_mp_ids.json", "w") as f:
                for entry in most_stable_mp_ids:
                    f.write(json.dumps(entry) + ",\n")
            for entry in most_stable_mp_ids:
                print(f"The most stable mp_id for {entry['name']} is: {entry['mp_id']}")
        elif args.bulktype == "alloy":
            elements = alloy_element_list
            most_stable_mp_ids = get_most_stable_alloy_mp_id(args.api_key, args.elementNumber, elements, args.ificsd)
            with open("alloy_mp_ids.json", "w") as f:
                for entry in most_stable_mp_ids:
                    f.write(json.dumps(entry) + ",\n")
            for entry in most_stable_mp_ids:
                print(f"The most stable mp_id for {entry['name']} is: {entry['mp_id']}")
        elif args.bulktype == "oxide":
            elements = oxide_element_list
        else:
            print("Error bulktype must be metal, alloy, or oxide")
        
    elif args.task == "generate": 
        bulk_generator = BulkGenerator(api_key=args.api_key)
        # Fetch structures once and reuse
        structures = bulk_generator.fetch_structures(type4,args.min_lw)
        
        # Generate POSCAR and plot using the fetched structures
        bulk_generator.generate_bulk(structures)
        if args.plot:
            bulk_generator.plot_bulk(structures)
