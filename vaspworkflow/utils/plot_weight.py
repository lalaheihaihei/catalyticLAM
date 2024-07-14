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
    def __init__(self, poscar_dir="./",fmt='npy_mix'): 
        self.poscar_dir = poscar_dir
        self.fmt=fmt
        if self.fmt=='npy_mix':
            self.element_weights = self.calculate_element_weights_npy_mul()
        elif self.fmt=='npy_single':
            self.element_weights = self.calculate_element_weights_npy()
    
    def calculate_element_weights(self):
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
        element_frequencies = {}

        # Iterate over POSCAR files to calculate frequencies
        for poscar_file in glob.glob(pathname=os.path.join(self.poscar_dir,'**/type.raw'),recursive=True):
            # poscar = dp.LabeledSystem(os.path.dirname(poscar_file),fmt='deepmd/npy')
            # structure = poscar['atom_names']
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
        element_frequencies = {}

        # Iterate over POSCAR files to calculate frequencies
        for poscar_file in tqdm(glob.glob(pathname=os.path.join(self.poscar_dir,'**/type_map.raw'),recursive=True)):
            # poscar = dp.LabeledSystem(os.path.dirname(poscar_file),fmt='deepmd/npy')
            # poscar=dp.MultiSystems().load_systems_from_file(os.path.dirname(poscar_file), fmt="deepmd/npy/mixed")
            type=np.load(os.path.join(os.path.dirname(poscar_file),'set.000000/real_atom_types.npy'))
            ntype=open(os.path.join(os.path.dirname(poscar_file),'type_map.raw'))
            e=[i.strip('\n') for i in ntype.readlines()]
            all_stru=[]
            for j in type:
                # structure = set(j)
                structure = [e[i] for i in set(j)]
                structure=sorted(structure)
                if structure not in all_stru:
                    all_stru.append(structure)
                    for element in structure:
                    # symbol = element.symbol
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
        # 初始化字典来记录每种元素出现过的POSCAR文件数
        element_occurrences = {}

        # 遍历POSCAR文件
        for poscar_file in glob.glob(os.path.join(self.poscar_dir, "POSCAR-*")):
            poscar = Poscar.from_file(poscar_file)
            structure = poscar.structure
            
            # 使用集合来记录当前POSCAR文件中出现过的元素
            elements_in_file = {element.symbol for element in structure.composition.elements}
            
            # 更新每种元素的出现次数
            for element in elements_in_file:
                if element in element_occurrences:
                    element_occurrences[element] += 1
                else:
                    element_occurrences[element] = 1

        # 总的POSCAR文件数量
        total_files = len(glob.glob(os.path.join(self.poscar_dir, "POSCAR-*")))

        # 将出现次数转换为比例
        element_weights = {symbol: occurrences / total_files for symbol, occurrences in element_occurrences.items()}

        return element_weights

    def get_text_color(self,background_color):
        """
        根据背景颜色的亮度决定文本颜色。对于较暗的背景使用白色字体，对于较亮的背景使用黑色字体。
        """
        color = mcolors.to_rgb(background_color)
        # 使用亮度的简单估算来决定文本颜色
        # print(sum(color[:3]) / 3)
        if sum(color[:3]) / 3 < 0.5:
            return 'white'
        else:
            return 'black'


    def plot_2d(self,name):
        fig, ax = plt.subplots(figsize=(18, 8))  # 调整尺寸以适应镧系元素
        cmap = plt.get_cmap('YlGnBu')  # 使用YlGnBu颜色映射
        norm = Normalize(vmin=-0.005, vmax=0.05)
        sm = ScalarMappable(norm=norm, cmap=cmap)
    
        max_row = max(el.row for el in Element if el.Z <= 86)
    
        for el in Element:
            if el.Z > 86 or el.Z in range(57, 72):  # 排除镧系元素以外的元素
                continue
            weight = self.element_weights.get(el.symbol, 0)
            color = sm.to_rgba(weight)
            x, y = el.group - 1, max_row - el.row
            text_color = self.get_text_color(color)  # 动态选择文本颜色

            ax.add_patch(plt.Rectangle((x + 0.05, y + 0.05), 0.9, 0.9, color=color, edgecolor=text_color, linewidth=0.5))
            # 元素符号和原子序号
            ax.text(x + 0.5, y + 0.5, el.symbol, ha='center', va='center', color=text_color, fontsize=15)  # 增大字体
            ax.text(x + 0.1, y + 0.9, f'{el.Z}', ha='left', va='top', color=text_color, fontsize=10)  # 原子序号
            ax.text(x + 0.5, y + 0.05, f'{weight:.3f}', ha='center', va='bottom', color=text_color, fontsize=8)  # 权重
    
        # 特别处理镧系元素
        lanthanides = [el for el in Element if 57 <= el.Z <= 71]
        for i, el in enumerate(lanthanides):
            weight = self.element_weights.get(el.symbol, 0)
            color = sm.to_rgba(weight)
            x = i + 2
            y = -1  # 镧系元素所在行
            text_color = self.get_text_color(color)  # 动态选择文本颜色
            ax.add_patch(plt.Rectangle((x + 0.05, y + 0.05), 0.9, 0.9, color=color, edgecolor=text_color, linewidth=0.5))
            ax.text(x + 0.5, y + 0.5, el.symbol, ha='center', va='center', color=text_color, fontsize=12)  # 增大字体
            ax.text(x + 0.1, y + 0.9, f'{el.Z}', ha='left', va='top', color=text_color, fontsize=8)  # 原子序号
            ax.text(x + 0.5, y + 0.05, f'{weight:.3f}', ha='center', va='bottom', color=text_color, fontsize=8)  # 权重
    
        ax.set_xlim(0, max(18, len(lanthanides)))
        ax.set_ylim(-2, max_row)
        ax.set_aspect('equal')
        ax.axis('off')
    
        cbar = plt.colorbar(sm, ax=ax, orientation='vertical', fraction=0.02, pad=0.04)
        cbar.set_label('Weight')
    
        plt.savefig(name)
    
    
    
    def plot_3d(self):
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



