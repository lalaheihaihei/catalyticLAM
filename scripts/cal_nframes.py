# calculate numbers of frames of dataset
import os
from glob import glob
import numpy as np
import argparse

# 创建 ArgumentParser 对象
parser = argparse.ArgumentParser(description='calculate numbers of frames of dataset')

parser.add_argument('--dataset', type=str,
                    help='the path of the dataset')
args = parser.parse_args()

def find_folders_with_outcar(search_path='./', label='OUTCAR'):
    """
    Finds all the folders containing a file with the specified label in the given search path.
    Args:
        search_path (str, optional): The path to search for the folders. Defaults to './'.
        label (str, optional): The label of the file to search for. Defaults to 'OUTCAR'.
    Returns:
        list: A list of unique folder paths containing the file with the specified label.
    """
    seen_folders = []  
    for outcar_path in glob(os.path.join(search_path, '**', label), recursive=True):  
        folder_path = os.path.dirname(outcar_path)  
        if folder_path not in seen_folders:  
            seen_folders.append(folder_path)
    return seen_folders

<<<<<<< HEAD

d = find_folders_with_outcar(args.dataset, 'energy.npy') # At this point, d contains all folders with 'energy.npy'
=======
d = find_folders_with_outcar('./Your-Path/', 'energy.npy') # At this point, d contains all folders with 'energy.npy'
>>>>>>> upstream/master
num = 0
for i in d:
    num += len(np.load(f'{i}/energy.npy'))
print(num)
