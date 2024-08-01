import torch   # 此脚本用于删除模型中多余的参数，简化模型的finetune过程
checkpoint = torch.load('/home/ljcgroup/wzh/fairchem/fairchem/fairchem/ceshi/batch-4/checkpoints/2024-07-24-13-13-36/checkpoint.pt', map_location='cpu') # 这里的写原始checkpoint的路径
checkpoint_re = {i:j for i,j in checkpoint_re.items() if i not in ['epoch', 'step', 'optimizer', 'scheduler', 'ema', 'best_val_metric', 'primary_metric'] }
torch.save(checkpoint_re, '/home/ljcgroup/wzh/fairchem/fairchem/fairchem/ceshi/batch-4/checkpoints/2024-07-24-13-13-36/model-new.pt')  # 保持的新的checkpoint的路径
