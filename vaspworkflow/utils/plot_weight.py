import os
import glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from pymatgen.core.periodic_table import Element
from pymatgen.io.vasp import Poscar
import dpdata as dp
from tqdm import tqdm
# from typing import Dict

class PeriodicTableWeightsVisualizer:
    """ the class to visualize the periodic table weights

    Args:
        poscar_dir (str): The directory containing POSCAR files.
        fmt (str): The format of the data files. only support npy_mix, npy_single, POSCAR three formats
        npy format is the output of dpdata
    """
    def __init__(self, poscar_dir="./",fmt='npy_mix'): 
        self.poscar_dir = poscar_dir
        self.fmt=fmt
        if self.fmt=='npy_mix':
            self.element_weights = self.calculate_element_weights_npy_mul()
        elif self.fmt=='npy_single':
            self.element_weights = self.calculate_element_weights_npy()
        elif self.fmt=='POSCAR':
            self.element_weights = self.calculate_element_weights()
    
    def calculate_element_weights(self):
        """ Calculate the element weights from POSCAR files."""

        # Initialize a dictionary to count element frequencies
        element_frequencies = {}

        # Iterate over POSCAR files to calculate frequencies
        for poscar_file in glob.glob(os.path.join(self.poscar_dir, "POSCAR*")):
            poscar = Poscar.from_file(poscar_file)
            structure = poscar.structure

            for element in structure.composition.elements:
                symbol = element.symbol
                element_frequencies[symbol] = element_frequencies.get(symbol, 0) + 1 # structure.composition[element]

        # Convert frequencies to proportions
        total_count = sum(element_frequencies.values())
        return {symbol: freq / total_count for symbol, freq in element_frequencies.items()}

    def calculate_element_weights_npy(self):
        """ Calculate the element weights from npy files deduced by dpdata with the most common data format."""
        element_frequencies = {}

        # Iterate over POSCAR files to calculate frequencies
        for poscar_file in glob.glob(pathname=os.path.join(self.poscar_dir,'**/type.raw'),recursive=True):

            ntype=open(os.path.join(os.path.dirname(poscar_file),'type_map.raw'))
            ntype1=[i.strip('\n') for i in ntype.readlines()]
            map=open(os.path.join(os.path.dirname(poscar_file),'type.raw'))
            map1=[i.strip('\n') for i in map.readlines()]

            structure=set([ntype1[int(i)] for i in map1])

            for element in structure:
                # symbol = element.symbol
                element_frequencies[element] = element_frequencies.get(element, 0) + 1 # structure.composition[element]

        # Convert frequencies to proportions
        total_count = sum(element_frequencies.values())
        return {symbol: freq / total_count for symbol, freq in element_frequencies.items()}

    def calculate_element_weights_npy_mul(self):
        """ Calculate the element weights from npy files. the mix_format deduced by dpdata, please refer to the doc of dpdata"""
        element_frequencies = {}

        # Iterate over POSCAR files to calculate frequencies
        for poscar_file in tqdm(glob.glob(pathname=os.path.join(self.poscar_dir,'**/type_map.raw'),recursive=True)):
            # ther may raise error for the older dataset deduced by dpdata ,the direction should be set.000/real_atom_types.npy
            type=np.load(os.path.join(os.path.dirname(poscar_file),'set.000000/real_atom_types.npy'))
            ntype=open(os.path.join(os.path.dirname(poscar_file),'type_map.raw'))
            e=[i.strip('\n') for i in ntype.readlines()]
            all_stru=[]
            for j in type:

                structure = [e[i] for i in set(j)]
                structure=sorted(structure)
                if structure not in all_stru:
                    all_stru.append(structure)
                    for element in structure:
                        element_frequencies[element] = element_frequencies.get(element, 0) + 1 # structure.composition[element]
                else:
                    pass

                for element in structure:
                # symbol = element.symbol
                    element_frequencies[element] = element_frequencies.get(element, 0) + 1 # structure.composition[element]

        # Convert frequencies to proportions
        total_count = sum(element_frequencies.values())
        return {symbol: freq / total_count for symbol, freq in element_frequencies.items()}
    def calculate_atoms_element_weights(self):
        """ Calculate the element weights from POSCAR files."""
        element_occurrences = {}

         # iterate over POSCAR files to calculate frequencies
        for poscar_file in glob.glob(os.path.join(self.poscar_dir, "POSCAR-*")):
            poscar = Poscar.from_file(poscar_file)
            structure = poscar.structure
            
            # use set to avoid duplicate elements
            elements_in_file = {element.symbol for element in structure.composition.elements}
            
           # updata the element occurrences
            for element in elements_in_file:
                if element in element_occurrences:
                    element_occurrences[element] += 1
                else:
                    element_occurrences[element] = 1

            # tatal number of files
        total_files = len(glob.glob(os.path.join(self.poscar_dir, "POSCAR-*")))

        # transform to proportions
        element_weights = {symbol: occurrences / total_files for symbol, occurrences in element_occurrences.items()}

        return element_weights

    def get_text_color(self,background_color):
        """ Get the text color based on the background color."""
        color = mcolors.to_rgb(background_color)

        if sum(color[:3]) / 3 < 0.5:
            return 'white'
        else:
            return 'black'


    def plot_2d(self,name):
        """ Plot the element weights as a 2D heatmap."""
        fig, ax = plt.subplots(figsize=(18, 8))  # set the figure size
        cmap = plt.get_cmap('YlGnBu')  # use the YlGnBu colormap
        norm = Normalize(vmin=-0.005, vmax=0.05)
        sm = ScalarMappable(norm=norm, cmap=cmap)
    
        max_row = max(el.row for el in Element if el.Z <= 86)
    
        for el in Element:
            if el.Z > 86 or el.Z in range(57, 72):  # exclude lanthanides and actinides
                continue
            weight = self.element_weights.get(el.symbol, 0)
            color = sm.to_rgba(weight)
            x, y = el.group - 1, max_row - el.row
            text_color = self.get_text_color(color)  

            ax.add_patch(plt.Rectangle((x + 0.05, y + 0.05), 0.9, 0.9, color=color, edgecolor=text_color, linewidth=0.5))
            # text element name and position
            ax.text(x + 0.5, y + 0.5, el.symbol, ha='center', va='center', color=text_color, fontsize=15)  
            ax.text(x + 0.1, y + 0.9, f'{el.Z}', ha='left', va='top', color=text_color, fontsize=10)  
            ax.text(x + 0.5, y + 0.05, f'{weight:.3f}', ha='center', va='bottom', color=text_color, fontsize=8)  
    
        # special treatment for lanthanides
        lanthanides = [el for el in Element if 57 <= el.Z <= 71]
        for i, el in enumerate(lanthanides):
            weight = self.element_weights.get(el.symbol, 0)
            color = sm.to_rgba(weight)
            x = i + 2
            y = -1  
            text_color = self.get_text_color(color) 
            ax.add_patch(plt.Rectangle((x + 0.05, y + 0.05), 0.9, 0.9, color=color, edgecolor=text_color, linewidth=0.5))
            ax.text(x + 0.5, y + 0.5, el.symbol, ha='center', va='center', color=text_color, fontsize=12)  
            ax.text(x + 0.1, y + 0.9, f'{el.Z}', ha='left', va='top', color=text_color, fontsize=8)  
            ax.text(x + 0.5, y + 0.05, f'{weight:.3f}', ha='center', va='bottom', color=text_color, fontsize=8) 
    
        ax.set_xlim(0, max(18, len(lanthanides)))
        ax.set_ylim(-2, max_row)
        ax.set_aspect('equal')
        ax.axis('off')
    
        cbar = plt.colorbar(sm, ax=ax, orientation='vertical', fraction=0.02, pad=0.04)
        cbar.set_label('Weight')
    
        plt.savefig(name)
    
    
    
    def plot_3d(self):
        """ Plot the element weights as a 3D bar chart."""
        fig = plt.figure(figsize=(20, 10))
        ax = fig.add_subplot(111, projection='3d')
        cmap = plt.get_cmap('YlGnBu')  # Use predefined YlGnBu colormap

        max_row = max(el.row for el in Element if el.Z <= 86)
        norm = Normalize(vmin=min(self.element_weights.values()), vmax=max(self.element_weights.values()))
        sm = ScalarMappable(norm=norm, cmap=cmap)

        dx = 0.8
        dy = 0.8
        spacing = 0.2
        text_z_offset = max(self.element_weights.values()) + 1

        for el in sorted(Element, key=lambda el: el.Z):
            if el.Z > 86:
                continue
            weight = self.element_weights.get(el.symbol, 0)
            color = sm.to_rgba(weight)
            x = (el.group - 1) * (dx + spacing) if el.group else 0
            y = (max_row - el.row) * (dy + spacing)
            z_height = weight
            ax.bar3d(x, y, 0, dx, dy, z_height, color=color, shade=True)
            ax.text(x + dx/2, y + dy/2, text_z_offset, f'{el.symbol}', color='black', ha='center', va='center', fontsize=10)

        ax.set_zlim(0, text_z_offset)
        ax.set_xlabel('Group')
        ax.set_ylabel('Period')
        ax.set_zlabel('Weight')
        ax.set_box_aspect([3,1,1])
        ax.view_init(elev=80, azim=270)
        plt.savefig("weight3d.png")

# Example usage



