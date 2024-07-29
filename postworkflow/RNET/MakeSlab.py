import importlib.util
import os
import matplotlib.pyplot as plt
import argparse
from pymatgen.analysis.adsorption import AdsorbateSiteFinder, reorient_z
from pymatgen.core.surface import generate_all_slabs
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.core.structure import Molecule
from pymatgen.io.vasp.inputs import Poscar
from mp_api.client import MPRester
from ase import Atoms
from ase.visualize.plot import plot_atoms
from pymatgen.io.ase import AseAtomsAdaptor
from molecule_db import molecule_database
from material_db import type1
import numpy as np

class SlabPlotter:
    def __init__(self, rows=8, cols=4):
        """
        Initialize the plotter with a specific grid size.
        """
        self.rows = rows
        self.cols = cols
        self.fig = None
        self.axs = None
        self.current_ax_index = 0

    def get_ase_atoms(self, pmg_structure):
        """Convert a pymatgen Structure or Slab to an ASE Atoms object."""
        adaptor = AseAtomsAdaptor()
        return adaptor.get_atoms(pmg_structure)

    def plot_slabs_with_side_view_ase(self, ads_structs, material_id, miller_index):
        """
        Plot the structures on the current figure, creating a new one if necessary,
        and add Miller index above each subplot.

        :param ads_structs: A list of pymatgen Structure objects for plotting.
        :param material_id: The material ID, used for saving the figure.
        :param miller_index: The Miller index of the slab to be displayed on the plot.
        """
        if self.fig is None:
            self.fig, self.axs = plt.subplots(self.rows, self.cols, figsize=(self.cols * 3, self.rows * 3))
            self.axs = self.axs.flatten()  # Flatten the axes array for easier access
        
        for ads_struct in ads_structs:
            ase_atoms = self.get_ase_atoms(ads_struct)
            if self.current_ax_index < len(self.axs):
                ax = self.axs[self.current_ax_index]
                plot_atoms(ase_atoms, ax, rotation='90x')
                ax.title.set_text(f"Miller index: {miller_index}")  # 添加 Miller index 信息
                self.current_ax_index += 1
            if self.current_ax_index < len(self.axs):
                ax = self.axs[self.current_ax_index]
                plot_atoms(ase_atoms, ax, rotation='0x,0y,0z')
                ax.title.set_text(f"Miller index: {miller_index}")  # 添加 Miller index 信息
                self.current_ax_index += 1

    def save_figure(self, material_id,molecule_type):
        """Save the current figure to a file and reset the plotter for the next use."""
        if self.fig:
            plt.tight_layout()
            filename = f"slabs-{material_id['name']}-{material_id['mp_id']}-{molecule_type}.png"
            self.fig.savefig(filename, dpi=300)
            plt.close(self.fig)
            self.fig = None
            self.axs = None
            self.current_ax_index = 0


class SlabGenerator:
    def __init__(self, api_key, enable_plotting):
        self.mpr = MPRester(api_key)  # Materials Project REST API client
        self.enable_plotting = enable_plotting  # 控制是否绘图
        print("Plot:",self.enable_plotting)
        if self.enable_plotting:
            self.plotter = SlabPlotter()  # 如果启用绘图，创建SlabPlotter的实例

    def generate_slabs(self, material_ids, molecule_type, up_down, max_index, min_slab_size, min_vacuum_size, min_lw, distance, element,max_slabs):
        """
        Generate slabs from given material IDs and adsorb molecules on them.
        
        :param material_ids: List of material IDs to generate slabs from.
        :param molecule_type: Type of molecule to adsorb on the slabs.
        :param max_index: Maximum Miller index to consider for slab generation.
        :param min_slab_size: Minimum size of the slab.
        :param min_vacuum_size: Minimum size of the vacuum layer.
        :param min_lw: Minimum slab model a and b vector
        """
        # Iterate through all provided material IDs
        for mp_id in material_ids:
            if element == mp_id["name"]:
                struct = self.mpr.get_structure_by_material_id(mp_id["mp_id"])  # Get structure from Materials Project
                struct = SpacegroupAnalyzer(struct,symprec=0.03).get_conventional_standard_structure()  # Standardize structure
                slabs = generate_all_slabs(struct, max_index=max_index, min_slab_size=min_slab_size,
                                       min_vacuum_size=min_vacuum_size, center_slab=True,max_normal_search=2)  # Generate slabs
                
                #delete repeat miller index
                unique_slabs = {}
                for slab in slabs:
                    miller_index = tuple(slab.miller_index)  # Convert list to tuple for hashability
                    if miller_index not in unique_slabs:
                        unique_slabs[miller_index] = slab  # Store slab with unique Miller index
            
                # Limit the number of slabs if max_slabs is set
                unique_slabs = list(unique_slabs.values())[:max_slabs] if max_slabs is not None else list(unique_slabs.values())
            else:
                continue

            for slab in unique_slabs:
                #slab = slab.reorient_z()
                # Optionally adjust slab size
                slab_a, slab_b = slab.lattice.a, slab.lattice.b
                if min_lw/slab_a > 1:
                    slab.make_supercell([[int(min_lw/slab_a)+1, 0, 0], [0, 1, 0], [0, 0, 1]])
                if min_lw/slab_b > 1:
                    slab.make_supercell([[1, 0, 0], [0, int(min_lw/slab_b)+1, 0], [0, 0, 1]])
                self.process_slabs(slab, struct, mp_id, molecule_type, up_down, distance, mp_id)
            if self.enable_plotting:
                self.plotter.save_figure(mp_id,molecule_type)  # Save the figure


    def process_slabs(self, slab, struct, mp_id, molecule_type, up_down, distance, material_id):
        """
        Process each slab: adsorb molecules and handle plotting and structure saving.
        
        :param slab: The slab to process.
        :param struct: The original structure from which the slab was generated.
        :param mp_id: Material ID of the slab.
        :param molecule_type: Type of molecule to adsorb.
        :param up_down: Indicator of where to adsorb molecules: "U" for up, "D" for down, "UD" for both, "UUD" for two up and one down.
        :param distance: the distance between adsorbate and slab 
        """
        asf = AdsorbateSiteFinder(slab)
        # Find the top and bottom atoms of slab
        min_z_atom = min(slab, key=lambda atom: atom.z)
        max_z_atom = max(slab, key=lambda atom: atom.z)
        # Check molecule database for the specified molecule type
        if molecule_type in molecule_database:
            species, coords = molecule_database[molecule_type]
            adsorbate = Molecule(species, coords)  # Create Molecule object
        else:
            print(f"Molecule type {molecule_type} not found in database. Using default CO molecule.")
            adsorbate = Molecule(["C", "O"], [(0, 0, 0), (0, 0, -1.3)])
        
        # Calculate the The span of the molecule in the Z direction, to adjast the adsorption distance
        z_coords = [coord[2] for coord in adsorbate.cart_coords]
        z_span = max(z_coords) - min(z_coords)
 
        # Adsorb molecule based on the specified surface direction
        if up_down == "D":  # Adsorb on the bottom surface
            min_z_atom.z -= distance + z_span 
            adsorbate = self.flip_molecule_vertically(adsorbate)  # Flip molecule for adsorption on the bottom
            ads_structs = [reorient_z(asf.add_adsorbate(adsorbate, min_z_atom.coords))]  # Add adsorbate
            if self.enable_plotting:
                self.plotter.plot_slabs_with_side_view_ase(ads_structs,material_id)
        elif up_down == "U":  # Adsorb on the top surface
            max_z_atom.z += distance
            ads_structs = [reorient_z(asf.add_adsorbate(adsorbate, max_z_atom.coords))]  # Add adsorbate
            if self.enable_plotting:
                self.plotter.plot_slabs_with_side_view_ase(ads_structs,material_id)

        elif up_down == "UD":  # Adsorb on both surfaces
            # For both surfaces, first flip and adsorb on bottom, then adsorb on top
            min_z_atom.z -= distance + z_span
            max_z_atom.z += distance


            #Flip a direction of adsorbate and put one on the bottom surface
            adsorbate_D = self.flip_molecule_vertically(adsorbate)
            ads_structs_D = asf.add_adsorbate(adsorbate_D, min_z_atom.coords)

            #Put one adsorbate on top surface
            asf_D = AdsorbateSiteFinder(ads_structs_D)
            ads_structs = [reorient_z(asf_D.add_adsorbate(adsorbate, max_z_atom.coords))]
            
            if self.enable_plotting:
                self.plotter.plot_slabs_with_side_view_ase(ads_structs,material_id)

        elif up_down == "UUD":  # Adsorb on both surfaces
            # For both surfaces, first flip and adsorb on bottom, then adsorb on 2 top
            min_z_atom.z -= distance + z_span
            max_z_atom.z += distance


            #Flip a direction of adsorbate and put one on the bottom surface
            adsorbate_D = self.flip_molecule_vertically(adsorbate)
            ads_structs_D = asf.add_adsorbate(adsorbate_D, min_z_atom.coords)

            #Put one adsorbate on top surface
            asf_D = AdsorbateSiteFinder(ads_structs_D)
            ads_structs_DU = asf_D.add_adsorbate(adsorbate, max_z_atom.coords)
            asf_DU = AdsorbateSiteFinder(ads_structs_DU)
            print(slab.lattice.matrix[1]) 
            
            # Put another adsorbate on top surface
            max_newx = max_z_atom.x + slab.lattice.matrix[1][0]/2
            max_newy = max_z_atom.y + slab.lattice.matrix[1][1]/2
            max_newz = max_z_atom.z + slab.lattice.matrix[1][2]/2
            max_newx += slab.lattice.matrix[0][0]/2
            max_newy += slab.lattice.matrix[0][1]/2
            max_newz += slab.lattice.matrix[0][2]/2
            max_new = (max_newx, max_newy, max_newz)
            ads_structs = [reorient_z(asf_DU.add_adsorbate(adsorbate, max_new))]
            
            miller_index = slab.miller_index

            if self.enable_plotting:
                #self.plotter.plot_slabs(slab, [min_z_atom.coords, max_z_atom.coords, max_new])
                self.plotter.plot_slabs_with_side_view_ase(ads_structs, material_id, miller_index=miller_index)

        self.save_structures(ads_structs, mp_id, struct, slab.miller_index, molecule_type)  # Save the structure

    def flip_molecule_vertically(self, molecule, angle_degrees=180):
        """
        Flip a molecule vertically by rotating it around the x-axis.
        
        :param molecule: Molecule to flip.
        :param angle_degrees: Angle of rotation in degrees, default is 180 for a full flip.
        :return: A new Molecule object that is flipped.
        """
        species = molecule.species
        coords = molecule.cart_coords
        angle_radians = np.radians(angle_degrees)
        rotation_matrix = np.array([
            [1, 0, 0],
            [0, np.cos(angle_radians), -np.sin(angle_radians)],
            [0, np.sin(angle_radians), np.cos(angle_radians)]
        ])
    
        rotated_coords = np.dot(coords, rotation_matrix.T)
        flipped_molecule = Molecule(species, rotated_coords)
        print(flipped_molecule) 
        return flipped_molecule
    
    def save_structures(self, ads_structs, mp_id, struct, miller_index, molecule_type):
        """
        Save the structures with adsorbates to POSCAR files.
        
        :param ads_structs: Structures with adsorbates.
        :param mp_id: Material ID of the structure.
        :param struct: Original structure.
        :param miller_index: Miller index of the slab.
        """
        mi_string = "".join(map(str, miller_index))
        base_filename = f"{struct.composition.reduced_formula}{mi_string}"
        for i, ads_struct in enumerate(ads_structs):
            filename = f"{base_filename}_{molecule_type}.vasp"
            with open(filename, 'w') as file:
                file.write(str(Poscar(reorient_z(ads_struct))))

    def run(self, material_ids, molecule_type, up_down, max_index, min_slab_size, min_vacuum_size, min_lw, distance, element, max_slabs):
        """
        Entry point to generate slabs, adsorb molecules, and save structures.
        
        :param material_ids: List of material IDs for slab generation.
        :param molecule_type: Type of molecule to adsorb.
        """
        for mp_id in material_ids:  
            if element == mp_id["name"]:  # find the mp_id corresponding to element
                if molecule_type == "all":
                    # Fetch the list of adsorbates for the element
                    molecule_types = mp_id["ads"]
                    break
                else:
                    molecule_types = [molecule_type]
        print('##############',molecule_types)
        for molecule in molecule_types:
            self.generate_slabs(material_ids, molecule, up_down, max_index, min_slab_size, min_vacuum_size, min_lw, distance, element, max_slabs)

def rename_vasp(directory, id2chem_file):
    """
    Rename all .vasp files in the specified directory by replacing the number 
    in the filename with the corresponding molecule formula from the id2chem_file.

    :param directory: Directory containing .vasp files.
    :param id2chem_file: Path to the file containing id to molecule formula mapping.
    """
    def load_id2chem_mapping(file_path):
        """
        Load the id to molecule formula mapping from the specified file.
        
        :param file_path: Path to the id2chem file.
        :return: Dictionary mapping ids to molecule formulas.
        """
        id2chem = {}
        with open(file_path, 'r') as file:
            for line in file:
                parts = line.strip().split(': ')
                if len(parts) == 2:
                    id2chem[int(parts[0])] = parts[1]
        return id2chem
    
    # Load id2chem mapping
    id2chem = load_id2chem_mapping(id2chem_file)
    
    # Iterate over files in the directory
    for filename in os.listdir(directory):
        if filename.endswith('.vasp'):
            base_name, ext = os.path.splitext(filename)
            # Extract the number part
            if '_' in base_name:
                number_str = base_name.split('_')[-1]
                try:
                    number = int(number_str)
                    # Get the corresponding molecule formula
                    molecule = id2chem.get(number, '')
                    new_filename = f"{base_name}-{molecule}{ext}"
                    # Rename the file
                    old_file_path = os.path.join(directory, filename)
                    new_file_path = os.path.join(directory, new_filename)
                    os.rename(old_file_path, new_file_path)
                    #print(f"Renamed '{filename}' to '{new_filename}'")
                except ValueError:
                    print(f"Skipping file '{filename}' due to invalid number format.")
            else:
                print(f"Skipping file '{filename}' due to invalid format.")

def parse_command_line_arguments():
    """
    Parses command line arguments and returns the parsed arguments.

    Returns:
        argparse.Namespace: An object containing the parsed command line arguments.
    """
    parser = argparse.ArgumentParser(description="Process some integers.")

    # Add all required arguments
    parser.add_argument("--plot", action="store_true", help="Enable plotting of top view of slabs and its adsorption sites, default: False")
    parser.add_argument("--api-key", type=str, default="OSpXCDIm3zKRpUNhh27uTXcmqEkSXbky", help="Materials Project API key")
    # parser.add_argument("--api-key", type=str, required=True, help="Materials Project API key")
    parser.add_argument("--max-index", type=int, default=1, help="Maximum Miller index to consider for slab generation, default: 1")
    parser.add_argument("--min-slab-size", type=float, default=8.0, help="Minimum size of the slab, default: 8.0")
    parser.add_argument("--min-vacuum-size", type=float, default=15.0, help="Minimum size of the vacuum layer, default: 15.0")
    parser.add_argument("--min-lw", type=float, default=10.0, help="Minimum slab model a and b vector, default: 10.0")
    parser.add_argument("--distance", type=float, default=2.0, help="Distance between adsorbate and slab, default: 2.0")
    parser.add_argument("--element", type=str, default="Au", help="Chemical formula of the materials to process, ex: Pt3Ni or Pt")
    parser.add_argument("--max-slabs", type=int, default=None, help="Maximum number of slabs to process per material")

    return parser.parse_args()
    
if __name__ == "__main__":

    args = parse_command_line_arguments()
    api_key = args.api_key
    molecule_type = args.molecule_type
    enable_plotting = args.plot
    up_down = args.up_down
    slab_generator = SlabGenerator(api_key, enable_plotting)
    slab_generator.run(material_ids=type1, molecule_type="all", up_down='U',
                       max_index=args.max_index, min_slab_size=args.min_slab_size,
                       min_vacuum_size=args.min_vacuum_size, min_lw=args.min_lw, distance=args.distance,
                       element=args.element,max_slabs=args.max_slabs)    
    rename_vasp('.', 'id2chem.txt')