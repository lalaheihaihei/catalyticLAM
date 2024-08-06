# 统计数据集中数据量总数
import os
from glob import glob
import numpy as np

def find_folders_with_outcar(search_path='./', label='OUTCAR'):
    seen_folders = []  
    # 遍历指定搜索路径下的所有OUTCAR文件  
    for outcar_path in glob(os.path.join(search_path, '**', label), recursive=True):  
    # 获取OUTCAR文件所在的文件夹路径  
        folder_path = os.path.dirname(outcar_path)  
    # 如果文件夹路径不是我们已经打印过的（避免重复打印）  
        if folder_path not in seen_folders:  
        # 将文件夹路径添加到已打印的文件夹集合中  
            seen_folders.append(folder_path)
    return seen_folders

d = find_folders_with_outcar('/home/ljcgroup/wzh/data/qm/data-0.9/data-qm-train0.9', 'energy.npy') # 这里我们要找的文件是energy.npy
num = 0
for i in d:
    num += len(np.load(f'{i}/energy.npy'))
print(num)
