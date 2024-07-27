"""For backwards compatibility"""

from __future__ import annotations
import sys
#sys.setrecursionlimit(100000)
import wandb
wandb_api_key = '697cee75f4ed39a1ffea28ab11f806b6ba88563f'
wandb.login(key=wandb_api_key)
#697cee75f4ed39a1ffea28ab11f806b6ba88563f

if __name__ == "__main__":
    from fairchem.core._cli import main

    main()
