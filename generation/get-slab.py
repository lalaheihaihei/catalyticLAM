from itertools import combinations
import numpy as np
from pymatgen.analysis.adsorption import AdsorbateSiteFinder, plot_slab, reorient_z, get_rot
from pymatgen.core.surface import generate_all_slabs
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.core.structure import Molecule
from pymatgen.io.vasp.inputs import Poscar
from matplotlib import pyplot as plt
from mp_api.client import MPRester
from plotter import SlabPlotter
from molecule_db import molecule_database
import material_db
from material_db import *
import argparse

class SlabGenerator:
    def __init__(self, api_key, enable_plotting):
        self.mpr = MPRester(api_key)  # Materials Project REST API client
        self.enable_plotting = enable_plotting  # 控制是否绘图
        print("Plot:",self.enable_plotting)
        if self.enable_plotting:
            self.plotter = SlabPlotter()  # 如果启用绘图，创建SlabPlotter的实例

    def generate_slabs(self, material_ids, molecule_type, up_down, max_index, min_slab_size, min_vacuum_size, min_lw, distance, element,max_slabs,molecule_types):
        """
        Generate slabs from given material IDs and adsorb molecules on them.
        
        :param material_ids: List of material IDs to generate slabs from.
        :param molecule_type: Type of molecule to adsorb on the slabs.
        :param max_index: Maximum Miller index to consider for slab generation.
        :param min_slab_size: Minimum size of the slab.
        :param min_vacuum_size: Minimum size of the vacuum layer.
        :param min_lw: Minimum slab model a and b vector.
        :param molecule_types:for UUUUDDDD.
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
                self.process_slabs(slab, struct, mp_id, molecule_type, up_down, distance, mp_id, molecule_types)
            if self.enable_plotting:
                self.plotter.save_figure(mp_id,molecule_type)  # Save the figure


    def process_slabs(self, slab, struct, mp_id, molecule_type, up_down, distance, material_id, molecule_types):
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
            self.save_structures(ads_structs, mp_id, struct, slab.miller_index, molecule_type)  # Save the structure
       
        
        elif up_down == "U":  # Adsorb on the top surface
            max_z_atom.z += distance
            ads_structs = [reorient_z(asf.add_adsorbate(adsorbate, max_z_atom.coords))]  # Add adsorbate
            if self.enable_plotting:
                self.plotter.plot_slabs_with_side_view_ase(ads_structs,material_id)

            self.save_structures(ads_structs, mp_id, struct, slab.miller_index, molecule_type)  # Save the structure


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

            self.save_structures(ads_structs, mp_id, struct, slab.miller_index, molecule_type)  # Save the structure


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

        elif up_down == "UUUUDDDD":  # Adsorb four molecules, two on each surface
            # Ensure that there are enough molecule types to form combinations
            if len(molecule_types) < 4:
                print("Not enough molecule types provided.")
                return
            print("min and max z: ", min_z_atom.z,max_z_atom.z)
            min_z_atom.z -= distance
            max_z_atom.z += distance
 
            # Generate all possible combinations of four molecule types
            for comb in combinations(molecule_types, 4):
                # Create molecules for the current combination
                molecules = []
                for mol_type in comb:
                    if mol_type in molecule_database:
                        species, coords = molecule_database[mol_type]
                        molecule = Molecule(species, coords)
                    else:
                        print(f"Molecule type {mol_type} not found in database. Using default CO molecule.")
                        molecule = Molecule(["C", "O"], [(0, 0, 0), (0, 0, -1.3)])
                    molecules.append(molecule)
                    
                #Flip a direction of adsorbate and put one on the bottom surface
                adsorbate_D0 = self.flip_molecule_vertically(molecules[0])
                z_coords0 = [coord[2] for coord in molecules[0].cart_coords]
                z_span0 = max(z_coords0) - min(z_coords0)

                adsorbate_D1 = self.flip_molecule_vertically(molecules[1])
                z_coords1 = [coord[2] for coord in molecules[1].cart_coords]
                z_span1 = max(z_coords1) - min(z_coords1)

                adsorbate_D2 = self.flip_molecule_vertically(molecules[2])
                z_coords2 = [coord[2] for coord in molecules[2].cart_coords]
                z_span2 = max(z_coords2) - min(z_coords2)
                
                adsorbate_D3 = self.flip_molecule_vertically(molecules[3])
                z_coords3 = [coord[2] for coord in molecules[3].cart_coords]
                z_span3 = max(z_coords3) - min(z_coords3)
                
                min_z_atom.z -= z_span0
                ads_structs_D = asf.add_adsorbate(adsorbate_D0, min_z_atom.coords)
                asf_D = AdsorbateSiteFinder(ads_structs_D)
    
                min_newx1 = min_z_atom.x + slab.lattice.matrix[1][0]/2
                min_newy1 = min_z_atom.y + slab.lattice.matrix[1][1]/2
                min_newz1 = min_z_atom.z + slab.lattice.matrix[1][2]/2 + z_span0 - z_span1
                min_newx2 = min_newx1 + slab.lattice.matrix[0][0]/2
                min_newy2 = min_newy1 + slab.lattice.matrix[0][1]/2
                min_newz2 = min_newz1 + slab.lattice.matrix[0][2]/2 + z_span1 - z_span2
                min_newx3 = min_z_atom.x + slab.lattice.matrix[0][0]/2
                min_newy3 = min_z_atom.y + slab.lattice.matrix[0][1]/2
                min_newz3 = min_z_atom.z + slab.lattice.matrix[0][2]/2 + z_span0 - z_span3
                
                min_new1 = (min_newx1, min_newy1, min_newz1)
                min_new2 = (min_newx2, min_newy2, min_newz2)
                min_new3 = (min_newx3, min_newy3, min_newz3)
                ads_structs_DD = asf_D.add_adsorbate(adsorbate_D1, min_new1)
                asf_DD = AdsorbateSiteFinder(ads_structs_DD)
                ads_structs_DDD = asf_DD.add_adsorbate(adsorbate_D2, min_new2)
                asf_DDD = AdsorbateSiteFinder(ads_structs_DDD)
                ads_structs_DDDD = asf_DDD.add_adsorbate(adsorbate_D3, min_new3)
                asf_DDDD = AdsorbateSiteFinder(ads_structs_DDDD)
    
    
                #Put one adsorbate on top surface
                asf_DDDD = AdsorbateSiteFinder(ads_structs_DDDD)
                ads_structs_UDDDD = reorient_z((asf_DDDD.add_adsorbate(molecules[0], max_z_atom.coords)))
                asf_UDDDD = AdsorbateSiteFinder(ads_structs_UDDDD)
                print(slab.lattice.matrix[1])
    
                # Put another adsorbate on top surface
                max_newx1 = max_z_atom.x + slab.lattice.matrix[1][0]/2
                max_newy1 = max_z_atom.y + slab.lattice.matrix[1][1]/2
                max_newz1 = max_z_atom.z + slab.lattice.matrix[1][2]/2
                max_newx2 = max_newx1 + slab.lattice.matrix[0][0]/2
                max_newy2 = max_newy1 + slab.lattice.matrix[0][1]/2
                max_newz2 = max_newz1 + slab.lattice.matrix[0][2]/2
                max_newx3 = max_z_atom.x + slab.lattice.matrix[0][0]/2
                max_newy3 = max_z_atom.y + slab.lattice.matrix[0][1]/2
                max_newz3 = max_z_atom.z + slab.lattice.matrix[0][2]/2
    
                max_new1 = (max_newx1, max_newy1, max_newz1)
                max_new2 = (max_newx2, max_newy2, max_newz2)
                max_new3 = (max_newx3, max_newy3, max_newz3)
                ads_structs_UUDDDD = reorient_z(asf_UDDDD.add_adsorbate(molecules[1], max_new1))
                asf_UUDDDD = AdsorbateSiteFinder(ads_structs_UUDDDD)
                ads_structs_UUUDDDD = reorient_z(asf_UUDDDD.add_adsorbate(molecules[2], max_new2))
                asf_UUUDDDD = AdsorbateSiteFinder(ads_structs_UUUDDDD)
                ads_structs = [reorient_z(asf_UUUDDDD.add_adsorbate(molecules[3], max_new3))]
    
                miller_index = slab.miller_index
    
                if self.enable_plotting:
                    #self.plotter.plot_slabs(slab, [min_z_atom.coords, max_z_atom.coords, max_new])
                    self.plotter.plot_slabs_with_side_view_ase(ads_structs, material_id, miller_index)
                self.save_structures(ads_structs, mp_id, struct, slab.miller_index, comb[0]+'_'+comb[1]+'_'+comb[2]+'_'+comb[3])  # Save the structure
            self.plotter.save_figure(mp_id,str(slab.miller_index)+molecule_type)

    def flip_molecule_vertically(self, molecule, angle_degrees=180):
        """
        Flip a molecule vertically by rotating it around the x-axis and adjust the maximum z-coordinate to be 0.
    
        :param molecule: Molecule to flip.
        :param angle_degrees: Angle of rotation in degrees, default is 180 for a full flip.
        :return: A new Molecule object that is flipped and adjusted.
        """
        species = molecule.species
        coords = molecule.cart_coords
        angle_radians = np.radians(angle_degrees)
        rotation_matrix = np.array([
            [1, 0, 0],
            [0, np.cos(angle_radians), -np.sin(angle_radians)],
            [0, np.sin(angle_radians), np.cos(angle_radians)]
        ])
    
        # Rotate coordinates
        rotated_coords = np.dot(coords, rotation_matrix.T)
    
        # Adjust the maximum z-coordinate to be 0
        max_z = np.max(rotated_coords[:, 2])
        adjusted_coords = rotated_coords.copy()
        adjusted_coords[:, 2] -= max_z
    
        flipped_molecule = Molecule(species, adjusted_coords)
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
            filename = f"POSCAR-{base_filename}-{molecule_type}"
            with open(filename, 'w') as file:
                file.write(str(Poscar(reorient_z(ads_struct))))

    def run(self, material_ids, molecule_type, up_down, max_index, min_slab_size, min_vacuum_size, min_lw, distance, element, max_slabs):
        """
        Entry point to generate slabs, adsorb molecules, and save structures.
        
        :param material_ids: List of material IDs for slab generation.
        :param molecule_type: Type of molecule to adsorb.
        """
        print(material_ids)
        for mp_id in material_ids:  
            if element == mp_id["name"]:  # find the mp_id corresponding to element
                if molecule_type == "all":
                    # Fetch the list of adsorbates for the element
                    molecule_types = mp_id["ads"]
                    for molecule in molecule_types:
                        self.generate_slabs(material_ids, molecule, up_down, max_index, min_slab_size, min_vacuum_size, min_lw, distance, element, max_slabs, molecule_types)
                        if up_down == "UUUUDDDD":
                            break
                else:
                    molecule_types = [molecule_type]
                    print('##############',molecule_types)
                    for molecule in molecule_types:
                        self.generate_slabs(material_ids, molecule, up_down, max_index, min_slab_size, min_vacuum_size, min_lw, distance, element, max_slabs, molecule_types)
                        if up_down == "UUUUDDDD":
                            break

def parse_command_line_arguments():
    """
    Parses command line arguments and returns the parsed arguments.

    Returns:
        argparse.Namespace: An object containing the parsed command line arguments.
    """
    parser = argparse.ArgumentParser(description="Process some integers.")

    # Add all required arguments
    parser.add_argument("--plot", action="store_true", help="Enable plotting of top view of slabs and its adsorption sites, default: False")
    parser.add_argument("--api-key", type=str, default="rod4hw9iGEVqyrKAsfCHzZ7FhbYgYPgd", help="Materials Project API key")
    # parser.add_argument("--api-key", type=str, required=True, help="Materials Project API key")
    parser.add_argument("--molecule-type", type=str, default="NH3", help="Type of molecule to adsorb, default: NH3")
    parser.add_argument("--up-down", type=str, choices=['U', 'D', 'UD', 'UUD', 'UUUUDDDD'], default="UUD",
                        help="Indicator of where to adsorb molecules: 'U' for up, 'D' for down, 'UD' for both, 'UUD' for two up and one down, 'UUUUDDDD' for type3 to generate different adsorbates combination on top and bottom surface. Default: UUD")
    parser.add_argument("--max-index", type=int, default=2, help="Maximum Miller index to consider for slab generation, default: 2")
    parser.add_argument("--min-slab-size", type=float, default=8.0, help="Minimum size of the slab, default: 8.0")
    parser.add_argument("--min-vacuum-size", type=float, default=15.0, help="Minimum size of the vacuum layer, default: 15.0")
    parser.add_argument("--min-lw", type=float, default=10.0, help="Minimum slab model a and b vector, default: 10.0")
    parser.add_argument("--distance", type=float, default=2.0, help="Distance between adsorbate and slab, default: 2.0")
    parser.add_argument("--element", type=str, default="Au", help="Chemical formula of the materials to process, ex: Pt3Ni or Pt")
    parser.add_argument("--max-slabs", type=int, default=None, help="Maximum number of slabs to process per material")
    parser.add_argument("--type", type=str, default="type1", help="select which type of materials group in materials_db.py, default: type1")
    
    return parser.parse_args()


if __name__ == "__main__":
    # Example usage
    args = parse_command_line_arguments()
    api_key = args.api_key

    molecule_type = args.molecule_type
    enable_plotting = args.plot
    up_down = args.up_down

    material_ids = getattr(material_db, args.type)

    slab_generator = SlabGenerator(api_key, enable_plotting)
    slab_generator.run(material_ids=material_ids, molecule_type=args.molecule_type, up_down=args.up_down,
                       max_index=args.max_index, min_slab_size=args.min_slab_size,
                       min_vacuum_size=args.min_vacuum_size, min_lw=args.min_lw, distance=args.distance,
                       element=args.element,max_slabs=args.max_slabs)
