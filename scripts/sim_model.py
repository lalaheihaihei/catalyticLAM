import torch
checkpoint = torch.load('./Your-Path/checkpoint.pt', map_location='cpu')
checkpoint_re = {i:j for i,j in checkpoint_re.items() if i not in ['epoch', 'step', 'optimizer', 'scheduler', 'ema', 'best_val_metric', 'primary_metric'] }
torch.save(checkpoint_re, './Your-New-Path/model-new.pt')
